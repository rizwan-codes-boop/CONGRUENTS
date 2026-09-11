/* Only galaxy-loop solver/emission work is exposed from this optional library. */
#include "include/solver.h"
#include "solver_runtime.h"
#include <stdatomic.h>
#include <omp.h>
#define hcubature_v cg_cubature
#include "../spectra_funcs.h"
#include "../CR_spectra/CR_funcs.h"
#define gen_funcs_h
#define free2D cg_free2D
#define malloc cg_allocate
#define free cg_release
#include "steady_state_native.h"
#undef malloc
#undef free
#undef free2D
#undef hcubature_v

double q_p_inject=2.2, q_e_inject=2.2;
double T_CR_lims__GeV[2]={1e-3,1e8};
double E_CRe_lims__GeV[2]={1e-3+5.10998950e-4,1e8+5.10998950e-4};
static atomic_flag cg_solver_busy=ATOMIC_FLAG_INIT;
unsigned cg_solver_abi(void) { return 2; }

static int cg_axis(size_t n,const double *v) {
    if(!v || n<2 || n>4096) return 0;
    for(size_t i=0;i<n;i++)
        if(!isfinite(v[i]) || v[i]<=0 || (i && v[i]<=v[i-1])) return 0;
    return 1;
}
static int cg_values(size_t n,const double *v,int positive) {
    if(!v) return 0;
    for(size_t i=0;i<n;i++) if(!isfinite(v[i]) || (positive?v[i]<=0:v[i]<0)) return 0;
    return 1;
}
static int cg_table_valid(const cg_table_input *t,int one) {
    return t && cg_axis(t->nx,t->x) && (one?t->ny==1:cg_axis(t->ny,t->y)) &&
           t->ny<=4096 && t->nx*t->ny<=4194304 && cg_values(t->nx*t->ny,t->values,0);
}
static gsl_spline_object_2D cg_bind_table(const cg_table_input *t) {
    return gsl_so2D(t->nx,t->ny,t->x,t->y,t->values);
}

/* Only the secondary-injection quadrature remains in this worker.
 * Calorimetry, diffusion, normalizations and primary/proton arrays are Python-owned. */
static void cg_secondary_galaxy(size_t ne,const double *t,double nh,double Cp,
    const double *fcal,double *out) {
    gsl_spline_object_1D fc=gsl_so1D(ne,t,fcal);
    for(size_t j=0;j<ne;j++) out[j]=q_e(t[j],nh,Cp,1e8,fc);
    if(!cg_values(ne,out,0)) cg_fail(CG_NUMERIC);
}

int cg_transport_batch(int threads,size_t n,size_t ne,const double *t,const double *nh,
    const double *cp,const double *fcal,double *out,int *statuses) {
    if(threads<1 || !n || n>100000 || !cg_axis(ne,t) ||
       !cg_values(n,nh,1) || !cg_values(n,cp,1) ||
       !cg_values(n*ne,fcal,0) || !out || !statuses) return CG_INVALID;
    if(atomic_flag_test_and_set(&cg_solver_busy)) return CG_INVALID;
    gsl_error_handler_t *previous=gsl_set_error_handler_off();
    const size_t team=(size_t)threads<n?(size_t)threads:n;
    cg_worker *workers=calloc(team,sizeof(cg_worker));
    if(!workers) { gsl_set_error_handler(previous);atomic_flag_clear(&cg_solver_busy);return CG_ALLOC; }
    #pragma omp parallel for num_threads(team) schedule(guided)
    for(size_t i=0;i<n;i++) {
        cg_worker *w=&workers[omp_get_thread_num()];
        w->count=0;w->status=0;w->phase="setup";
        cg_active=w;
        if(!setjmp(w->escape)) cg_secondary_galaxy(ne,t,nh[i],cp[i],fcal+ne*i,out+ne*i);
        statuses[i]=w->status;
        cg_cleanup(w);cg_active=NULL;
    }
    free(workers);gsl_set_error_handler(previous);atomic_flag_clear(&cg_solver_busy);
    return CG_OK;
}

static void cg_solve_galaxy(size_t ne,size_t np,size_t ns,const double *e,const double *ph,
    const double *p,const double *diff,const double *inj,const double *kinetic,const double *fcal,const cg_table_input *ic,
    const cg_table_input *gamma,const cg_table_input *bs,const cg_table_input *sy,
    double *electrons,double *emission) {
    gsl_spline_object_2D emit=cg_bind_table(ic), transition=cg_bind_table(gamma), brems=cg_bind_table(bs);
    gsl_spline_object_1D sync=gsl_so1D(sy->nx,sy->x,sy->values);
    gsl_spline_object_1D dd=gsl_so1D(ne,e,diff), dh=gsl_so1D(ne,e,diff+ne);
    gsl_spline_object_1D q1=gsl_so1D(ne,e,inj), q2=gsl_so1D(ne,e,inj+ne);
    gsl_spline_object_1D ss[4];
    gsl_spline_object_1D fc=gsl_so1D(ne,kinetic,fcal);
    const double Cp=p[7];
    double bounds[2]={E_CRe_lims__GeV[0],E_CRe_lims__GeV[1]};
    cg_active->phase="disc steady state";
    CRe_steadystate_solve(1,bounds,(int)ns,p[1],p[2],p[0],1,&transition,brems,dd,q1,q2,&ss[0],&ss[1]);
    double *escape=cg_allocate(2*ne*sizeof(double));
    for(size_t j=0;j<ne;j++) {
        escape[j]=gsl_so1D_eval(ss[0],e[j])/tau_diff__s(e[j],p[0],dd);
        escape[ne+j]=gsl_so1D_eval(ss[1],e[j])/tau_diff__s(e[j],p[0],dd);
    }
    gsl_spline_object_1D qh1=gsl_so1D(ne,e,escape), qh2=gsl_so1D(ne,e,escape+ne);
    cg_active->phase="halo steady state";
    CRe_steadystate_solve(2,bounds,(int)ns,p[1]/1000.,p[3],50*p[0],1,&transition,brems,dh,qh1,qh2,&ss[2],&ss[3]);
    for(size_t k=0;k<4;k++) for(size_t j=0;j<ne;j++) electrons[k*ne+j]=gsl_so1D_eval(ss[k],e[j]);
    cg_active->phase="emission";
    for(size_t j=0;j<np;j++) {
        double energy=ph[j];
        emission[0*np+j]=eps_IC_3(energy,emit,ss[0]);
        emission[1*np+j]=eps_IC_3(energy,emit,ss[1]);
        emission[2*np+j]=eps_BS_3(energy,p[1],brems,ss[0]);
        emission[3*np+j]=eps_BS_3(energy,p[1],brems,ss[1]);
        emission[4*np+j]=eps_SY_4(energy,p[2],sync,ss[0]);
        emission[5*np+j]=eps_SY_4(energy,p[2],sync,ss[1]);
        emission[6*np+j]=eps_IC_3(energy,emit,ss[2]);
        emission[7*np+j]=eps_IC_3(energy,emit,ss[3]);
        emission[8*np+j]=eps_SY_4(energy,p[3],sync,ss[2]);
        emission[9*np+j]=eps_SY_4(energy,p[3],sync,ss[3]);
        /* Python supplies free-free and applies disc attenuation after return. */
        emission[10*np+j]=0.;
        emission[11*np+j]=0.;
        emission[12*np+j]=eps_pi(energy,p[1],Cp,1e8,fc);
        emission[13*np+j]=eps_pi_fcal1(energy,p[1],Cp,1e8,fc);
        emission[14*np+j]=q_nu(energy,p[1],Cp,1e8,fc);
    }
    if(!cg_values(4*ne,electrons,0) || !cg_values(15*np,emission,0)) cg_fail(CG_NUMERIC);
}

int cg_solver_batch(int threads,size_t n,size_t ne,size_t np,size_t ns,
    const double *e,const double *ph,const double *props,const double *diff,const double *inj,
    const double *kinetic,const double *fcal,
    const cg_table_input *ic,const cg_table_input *gamma,const cg_table_input *bs,const cg_table_input *sy,
    double *electrons,double *emission,int *statuses) {
    if(threads<1 || !n || n>100000 || ns<4 || ns>500 || !cg_axis(ne,e) || !cg_axis(np,ph) ||
       !cg_values(8*n,props,1) || !cg_values(2*n*ne,diff,1) || !cg_values(2*n*ne,inj,0) ||
       !cg_axis(ne,kinetic) || !cg_values(n*ne,fcal,0) ||
       !ic || !gamma || !electrons || !emission || !statuses ||
       !cg_table_valid(bs,0) || !cg_table_valid(sy,1)) return CG_INVALID;
    /* Same domain as the production model, permitting logspace endpoint roundoff. */
    for(size_t k=0;k<2;k++)
        if(fabs(e[k?(ne-1):0]/E_CRe_lims__GeV[k]-1)>1e-12) return CG_INVALID;
    for(size_t i=0;i<n;i++)
        if(!cg_table_valid(ic+i,0) || !cg_table_valid(gamma+i,0)) return CG_INVALID;
    /* GSL's handler is process-global: serialize this backend's batches.
     * Python also serializes calls. Never install the legacy abort handler. */
    if(atomic_flag_test_and_set(&cg_solver_busy)) return CG_INVALID;
    gsl_error_handler_t *previous=gsl_set_error_handler_off();
    const size_t team=(size_t)threads<n?(size_t)threads:n;
    cg_worker *workers=calloc(team,sizeof(cg_worker));
    if(!workers) {
        gsl_set_error_handler(previous);atomic_flag_clear(&cg_solver_busy);return CG_ALLOC;
    }
    #pragma omp parallel for num_threads(team) schedule(guided)
    for(size_t i=0;i<n;i++) {
        cg_worker *w=&workers[omp_get_thread_num()];
        w->count=0;w->status=0;w->phase="setup";
        cg_active=w;
        if(!setjmp(w->escape))
            cg_solve_galaxy(ne,np,ns,e,ph,props+8*i,diff+2*ne*i,inj+2*ne*i,
                           kinetic,fcal+ne*i,ic+i,gamma+i,bs,sy,electrons+4*ne*i,emission+15*np*i);
        statuses[i]=w->status;
        cg_cleanup(w);
        cg_active=NULL;
    }
    free(workers);
    gsl_set_error_handler(previous);
    atomic_flag_clear(&cg_solver_busy);
    return CG_OK;
}
