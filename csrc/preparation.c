/* Boundary adapter only: physics helpers are included unchanged.
 * Driver expressions copied below retain their native constants/order.
 */
#include "include/congruents.h"
#include "internal.h"
#include <string.h>
#include <stdint.h>
#include "../spectra_funcs.h"
#include "../data_calc.h"
#include "../spec_integrate.h"

double q_p_inject=2.2, q_e_inject=2.2;
double T_CR_lims__GeV[2]={1.e-3,1.e8};
double E_CRe_lims__GeV[2]={1.e-3+5.10998950e-4,1.e8+5.10998950e-4};

static int valid_row(const double *r) {
    return r && isfinite(r[0]) && r[0]>=0 && r[0]<=20 &&
        isfinite(r[1]) && r[1]>0 && isfinite(r[2]) && r[2]>0 &&
        isfinite(r[3]) && r[3]>0;
}
static int derive(const double *r, double *p) {
    if (!valid_row(r)) return CG_INVALID;
    double z=r[0], mass=r[1], radius=r[2], sfr=r[3];
    double area=M_PI*pow(radius*1e3,2);
    double stars=mass/(2.*area), sfs=sfr/(2.*area);
    double gas=Sigma_gas_Shi_iKS__Msolpcm2(sfs,stars);
    double sigma=sigma_gas_Yu__kmsm1(sfr);
    double dust=Tdust__K(z,sfr,mass);
    double h=pow(sigma,2)/(M_PI*4.302e-3*
        (gas+sigma/sigma_star_Bezanson__kmsm1(mass,radius)*stars));
    double nh=gas/(1.4*m_H__kg*2.*h)*Msol__kg/pow(pc__cm,3);
    double n=nh*1.4/1.17;
    double vAi=1000.*((sigma/sqrt(2.))/10.)/(sqrt(1.e-4/1.e-4)*2.);
    double B=sqrt(4.*M_PI*1.e-4*n*1.17*m_H__kg*1e3)*vAi*1e5;
    double halo=log10(sfr/mass)>-10. ? B/3. : B/1.5;
    double values[10]={h,nh,B,sigma,area,gas,sfs,stars,dust,halo};
    for (int i=0;i<10;i++)
        if (!isfinite(values[i]) || values[i]<=0) return CG_NUMERIC;
    memcpy(p,values,sizeof(values));
    return CG_OK;
}
int cg_galaxies(const cg_context *ctx, size_t n,const double *rows,double *out) {
    if (!ctx || n>100000 || (n && (!rows || !out))) return CG_INVALID;
    for(size_t i=0;i<n;i++) if(!valid_row(rows+4*i)) return CG_INVALID;
    int bad=0;
    #pragma omp parallel for num_threads(ctx->threads) reduction(|:bad)
    for(size_t i=0;i<n;i++) bad |= derive(rows+4*i,out+10*i);
    return bad ? CG_NUMERIC : CG_OK;
}
static int parameters(const double *row,double *p,double *props) {
    int status=derive(row,props);
    if(status) return status;
    p[0]=T_0_CMB__K*(1.+row[0]); p[1]=props[8];
    p[2]=row[1]; p[3]=row[3]; p[4]=row[2]; p[5]=props[0]; p[6]=0;
    return CG_OK;
}
int cg_preparation_bounds(const double *rows,size_t count,double c[8],double t[3]) {
    if(!rows || !count || count>100000 || !c || !t) return CG_INVALID;
    double maxb=0,minb=INFINITY,maxz=0,mint=INFINITY,maxt=0,p[10];
    for(size_t i=0;i<count;i++) {
        int status=derive(rows+4*i,p);if(status) return status;
        maxb=fmax(maxb,p[2]);minb=fmin(minb,p[9]);maxz=fmax(maxz,rows[4*i]);
        mint=fmin(mint,p[8]);maxt=fmax(maxt,p[8]);
    }
    c[0]=1e-16;c[1]=1e8;c[2]=1e-3+m_e__GeV;c[3]=1e8+m_e__GeV;
    c[4]=E_BB_peak__GeV(T_0_CMB__K)*1e-4;c[5]=1e-7;
    c[6]=(2.*pow(M_PI,2)*pow(m_e__g,2)*pow(c__cmsm1,3))/
        (3.*e__esu*maxb*h__ergs)*(c[0]*m_e__GeV)/pow(c[3],2)*0.1;
    c[7]=(2.*pow(M_PI,2)*pow(m_e__g,2)*pow(c__cmsm1,3))/
        (3.*e__esu*minb/10.*h__ergs)*(c[1]*(1+maxz)*m_e__GeV)/pow(c[2],2)*10.;
    t[0]=(maxz+1.)*T_0_CMB__K;t[1]=mint;t[2]=maxt;
    for(int i=0;i<8;i++) if(!isfinite(c[i]) || c[i]<=0) return CG_NUMERIC;
    return CG_OK;
}
typedef double (*phot_func)(double *,double);
static phot_func field_functions[7]={
    dndEphot_CMB__cmm3GeVm1,dndEphot_FIR__cmm3GeVm1,
    dndEphot_3000__cmm3GeVm1,dndEphot_4000__cmm3GeVm1,
    dndEphot_7500__cmm3GeVm1,dndEphot_UV__cmm3GeVm1,
    dndEphot_total__cmm3GeVm1};
int cg_radiation(const double row[4],size_t n,const double *energies,
                 const double bounds[2],double *fields,double urad[8]) {
    if(!bounds || !urad || (n && (!energies || !fields)) ||
       n>1000000 || !isfinite(bounds[0]) || !isfinite(bounds[1]) ||
       bounds[0]<=0 || bounds[1]<=bounds[0]) return CG_INVALID;
    for(size_t i=0;i<n;i++)
        if(!isfinite(energies[i]) || energies[i]<=0) return CG_INVALID;
    double p[7],props[10]; int status=parameters(row,p,props);
    if(status) return status;
    for(size_t i=0;i<n;i++) for(int j=0;j<7;j++) {
        fields[i*7+j]=field_functions[j](p,energies[i]);
        if(!isfinite(fields[i*7+j])) return CG_NUMERIC;
    }
    double limits[2]={bounds[0],bounds[1]};
    urad[0]=ISRF_integrate__GeVcmm3(field_functions[6],p,limits)*1e9;
    urad[1]=pow(props[2],2)/(8.*M_PI)/GeV__erg*1e9;
    for(int j=0;j<6;j++)
        urad[j+2]=ISRF_integrate__GeVcmm3(field_functions[j],p,limits)*1e9;
    for(int j=0;j<8;j++) if(!isfinite(urad[j])) return CG_NUMERIC;
    return CG_OK;
}
struct cg_table {size_t nx,ny; double *x,*y,*z; int owns_arrays;};
void cg_table_destroy(cg_table *t) {
    if(t) {
        if(t->owns_arrays) {free(t->x);free(t->y);free(t->z);}
        free(t);
    }
}
static int shape_valid(size_t nx,size_t ny) {
    return nx>=2 && nx<=4096 && ny>=1 && ny<=4096 && nx*ny<=4194304;
}
static int axis_valid(size_t n,const double *x) {
    if(!x) return 0;
    for(size_t i=0;i<n;i++)
        if(!isfinite(x[i]) || x[i]<=0 || (i && x[i]<=x[i-1])) return 0;
    return 1;
}
int cg_table_import(size_t nx,size_t ny,const double *x,const double *y,
                    const double *z,cg_table **out) {
    if(!out) return CG_INVALID;
    *out=NULL;
    if(!shape_valid(nx,ny) || !axis_valid(nx,x) ||
       !axis_valid(ny,y) || !z) return CG_INVALID;
    for(size_t i=0;i<nx*ny;i++)
        if(!isfinite(z[i]) || z[i]<0) return CG_NUMERIC;
    cg_table *t=calloc(1,sizeof(*t));
    if(!t) return CG_ALLOC;
    t->nx=nx;t->ny=ny;t->owns_arrays=1;
    t->x=malloc(nx*sizeof(double));t->y=malloc(ny*sizeof(double));
    t->z=malloc(nx*ny*sizeof(double));
    if(!t->x || !t->y || !t->z) {cg_table_destroy(t);return CG_ALLOC;}
    memcpy(t->x,x,nx*sizeof(double));memcpy(t->y,y,ny*sizeof(double));
    memcpy(t->z,z,nx*ny*sizeof(double));*out=t;return CG_OK;
}
int cg_table_borrow(size_t nx,size_t ny,double *x,double *y,double *z,cg_table **out) {
    if(!out) return CG_INVALID;
    *out=NULL;
    if(!shape_valid(nx,ny) || !axis_valid(nx,x) ||
       !axis_valid(ny,y) || !z) return CG_INVALID;
    for(size_t i=0;i<nx*ny;i++)
        if(!isfinite(z[i]) || z[i]<0) return CG_NUMERIC;
    cg_table *t=calloc(1,sizeof(*t));
    if(!t) return CG_ALLOC;
    t->nx=nx;t->ny=ny;t->x=x;t->y=y;t->z=z;
    /* calloc leaves owns_arrays=0: destroy never frees caller buffers. */
    *out=t;
    return CG_OK;
}
int cg_table_shape(const cg_table *t,size_t *nx,size_t *ny) {
    if(!t || !nx || !ny) return CG_INVALID;
    *nx=t->nx;*ny=t->ny;return CG_OK;
}
int cg_table_copy(const cg_table *t,double *x,double *y,double *z) {
    if(!t || !x || !y || !z) return CG_INVALID;
    memcpy(x,t->x,t->nx*sizeof(double));memcpy(y,t->y,t->ny*sizeof(double));
    memcpy(z,t->z,t->nx*t->ny*sizeof(double));return CG_OK;
}
static size_t interval(size_t n,const double *a,double v) {
    size_t lo=0,hi=n-1;
    while(hi-lo>1) {size_t mid=(hi+lo)/2;if(a[mid]<=v) lo=mid;else hi=mid;}
    return lo;
}
int cg_table_eval(const cg_table *t,size_t n,const double *x,const double *y,double *out) {
    if(!t || (n && (!x || !out || (t->ny>1 && !y)))) return CG_INVALID;
    for(size_t k=0;k<n;k++)
        if(!isfinite(x[k]) || (t->ny>1 && !isfinite(y[k]))) return CG_INVALID;
    for(size_t k=0;k<n;k++) {
        if(x[k]<t->x[0] || x[k]>t->x[t->nx-1] ||
           (t->ny>1 && (y[k]<t->y[0] || y[k]>t->y[t->ny-1]))) {out[k]=0;continue;}
        size_t i=interval(t->nx,t->x,x[k]);
        double u=(x[k]-t->x[i])/(t->x[i+1]-t->x[i]);
        if(t->ny==1) {out[k]=(1-u)*t->z[i]+u*t->z[i+1];continue;}
        size_t j=interval(t->ny,t->y,y[k]);
        double v=(y[k]-t->y[j])/(t->y[j+1]-t->y[j]);
        out[k]=(1-v)*((1-u)*t->z[j*t->nx+i]+u*t->z[j*t->nx+i+1])+
            v*((1-u)*t->z[(j+1)*t->nx+i]+u*t->z[(j+1)*t->nx+i+1]);
    }
    return CG_OK;
}
int cg_table_create(int kind,int field,double temp,size_t nx,size_t ny,
                    const double c[8],cg_table **out) {
    if(!out) return CG_INVALID;
    *out=NULL;
    if(!c || !shape_valid(nx,ny) || ny<2 || kind<0 || kind>3 ||
       field<0 || field>5 || !isfinite(temp) || temp<0) return CG_INVALID;
    for(int i=0;i<8;i++) if(!isfinite(c[i]) || c[i]<=0) return CG_INVALID;
    for(int i=0;i<8;i+=2) if(c[i+1]<=c[i]) return CG_INVALID;
    /* Bounded tested numerical domain; avoids legacy GSL overflow paths.
     * Broader physical ranges require a separate numerical-domain review. */
    if(c[0]<1e-16 || c[1]>1e8 || c[2]<1e-3+m_e__GeV ||
       c[3]>1e8+m_e__GeV || c[4]<1e-20 || c[5]>1e-7 ||
       c[6]<1e-100 || c[7]>1e100 || (field>=4 && (temp<2 || temp>1000)))
        return CG_INVALID;
    double eg[2]={c[0],c[1]},ee[2]={c[2],c[3]},ep[2]={c[4],c[5]},sx[2]={c[6],c[7]};
    size_t dims[2]={nx,ny};
    if(kind==3) {
        data_object_1D a=init_do_1D_sync(sx,nx);
        double y=1;
        int status=cg_table_import(nx,1,a.x_data,&y,a.z_data,out);
        data_object_1D_free(a);return status;
    }
    data_object_2D a;
    if(kind==2) a=init_do_2D_BS(eg,ee,dims);
    else {
        double temps[6]={3000,4000,7500,0,temp,temp};
        phot_func funcs[6]={dndEphot_BB__cmm3GeVm1,dndEphot_BB__cmm3GeVm1,
            dndEphot_BB__cmm3GeVm1,dndEphot_UVMattis__cmm3GeVm1,
            dndEphot_BB__cmm3GeVm1,dndEphot_modBB__cmm3GeVm1};
        a=kind==0 ? init_do_2D_IC(funcs[field],&temps[field],eg,ee,ep,dims) :
                    init_do_2D_IC_Gamma(funcs[field],&temps[field],eg,ee,ep,dims);
    }
    int status=cg_table_import(nx,ny,a.x_data,a.y_data,a.z_data,out);
    data_object_2D_free(a);return status;
}
int cg_temperature_grid(int field,double low,double high,double step,
                        size_t capacity,double *values,size_t *count) {
    if(!count || (field!=4 && field!=5) || !isfinite(low) || !isfinite(high) ||
       !isfinite(step) || step<=0 || high<low || low<0 || high>1000) return CG_INVALID;
    size_t n=0; double lim[2];
    double span=field==4 ? high-T_0_CMB__K : high-low;
    if(!isfinite(span/step) || span/step>4095) return CG_INVALID;
    if(field==4) {
        if(high<T_0_CMB__K) return CG_INVALID;
        CMB_num(low,high,step,&n);CMB_lim(low,high,lim);
    } else {FIR_num(low,high,step,&n);FIR_lim(low,high,lim);}
    if(n<2 || n>4096 || lim[0]<=0 || lim[1]<=lim[0]) return CG_INVALID;
    *count=n;
    if(!values) return capacity==0 ? CG_OK : CG_INVALID;
    if(capacity<n) return CG_INVALID;
    values[0]=lim[0];values[n-1]=lim[1];
    double delta=(lim[1]-lim[0])/(n-1);
    for(size_t i=1;i<n-1;i++) values[i]=lim[0]+i*delta;
    return CG_OK;
}
int cg_combine_ic(const double row[4],cg_table *const a[8],
                  double fc,double ff,cg_table **out) {
    if(!out) return CG_INVALID;
    *out=NULL;
    if(!a || !isfinite(fc) || !isfinite(ff) || fc<0 || fc>=1 || ff<0 || ff>=1)
        return CG_INVALID;
    for(int k=0;k<8;k++) {
        if(!a[k] || a[k]->ny<2) return CG_INVALID;
        if(k && (a[k]->nx!=a[0]->nx || a[k]->ny!=a[0]->ny ||
           memcmp(a[k]->x,a[0]->x,a[0]->nx*sizeof(double)) ||
           memcmp(a[k]->y,a[0]->y,a[0]->ny*sizeof(double)))) return CG_INVALID;
    }
    double p[7],props[10];int status=parameters(row,p,props);if(status) return status;
    double dilution[5]={
        C_dil(u_rad_BB__GeVcmm3(3000),L3000K__Lsol(p[2]),p[4],p[5]),
        C_dil(u_rad_BB__GeVcmm3(4000),L4000K__Lsol(p[2]),p[4],p[5]),
        C_dil(u_rad_BB__GeVcmm3(7500),L7500K__Lsol(p[3]),p[4],p[5]),
        C_dil(u_rad_UVMattis__GeVcmm3(),LUV__Lsol(p[3]),p[4],p[5]),
        C_dil(u_rad_modBB__GeVcmm3(p[1]),LFIR__Lsol(p[3]),p[4],p[5])};
    size_t n=a[0]->nx*a[0]->ny;double *z=malloc(n*sizeof(double));
    if(!z) return CG_ALLOC;
    for(size_t i=0;i<n;i++) {
        double cmb=a[4]->z[i]*fc+a[5]->z[i]*(1.-fc);
        double fir=a[6]->z[i]*ff+a[7]->z[i]*(1.-ff);
        z[i]=cmb+dilution[4]*fir+dilution[0]*a[0]->z[i]+
             dilution[1]*a[1]->z[i]+dilution[2]*a[2]->z[i]+dilution[3]*a[3]->z[i];
    }
    status=cg_table_import(a[0]->nx,a[0]->ny,a[0]->x,a[0]->y,z,out);
    free(z);return status;
}
