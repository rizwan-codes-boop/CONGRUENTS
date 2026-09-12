/* Standalone checked-runtime audit: no Python interpreter or OpenMP needed. */
#define CG_TESTING
#include "../csrc/solver_runtime.h"
#include <assert.h>

static int invalid_integrand(unsigned dim,size_t n,const double *x,void *p,unsigned fd,double *v) {
    (void)dim;(void)x;(void)p;
    for(size_t i=0;i<n*fd;i++) v[i]=NAN;
    return 0;
}
static int exercise(size_t fail_at, int mode) {
    /* Heap storage remains defined across setjmp/longjmp. */
    cg_worker *w=calloc(1,sizeof(*w));
    assert(w);
    w->phase="runtime audit";w->fail_allocation=fail_at;cg_active=w;
    if(!setjmp(w->escape)) {
        double axis[2]={1.,2.}, values[2]={2.,4.}, plane[4]={2.,3.,3.,4.};
        gsl_spline_object_1D a=gsl_so1D(2,axis,values);
        gsl_spline_object_2D b=gsl_so2D(2,2,axis,axis,plane);
        assert(gsl_so1D_eval(a,1.5)==3.);
        assert(gsl_so2D_eval(b,1.5,1.5)==3.);
        assert(isnan(gsl_so1D_eval(a,.5)));
        assert(gsl_so2D_eval(b,.5,1.5)==0.);
        double matrix[4]={2.,0.,0.,4.}, rhs[2]={4.,8.}, solution[2];
        if(mode==1) matrix[3]=0.; /* Singular LU must not abort the process. */
        cg_linear_solve(2,matrix,rhs,solution);
        assert(solution[0]==2. && solution[1]==2.);
        void *extra=cg_allocate(128);
        cg_release(extra);
        gsl_so1D_free(a);gsl_so2D_free(b);
        if(mode==2) {
            double low=0.,high=1.,v,error;
            assert(cg_cubature(1,invalid_integrand,NULL,1,&low,&high,100,0.,1e-8,
                               ERROR_INDIVIDUAL,&v,&error)==1);
            assert(isnan(v));
        }
    }
    int status=w->status;
    cg_cleanup(w);
    assert(w->count==0);
    cg_cleanup(w); /* Cleanup is idempotent. */
    free(w);cg_active=NULL;
    return status;
}
int main(void) {
    gsl_error_handler_t *previous=gsl_set_error_handler_off();
    for(int repeat=0;repeat<100;repeat++) {
        assert(exercise(0,0)==0);
        /* 2 spline1D allocations, 3 spline2D, 2 LU, 1 raw buffer. */
        for(size_t fail=1;fail<=8;fail++) assert(exercise(fail,0)==2);
        assert(exercise(0,1)==3);
        assert(exercise(0,2)==3);
        assert(exercise(0,0)==0); /* Failure never poisons the next call. */
    }
    gsl_set_error_handler(previous);
    puts("Runtime audit passed: 1200 success/failure/recovery exercises.");
    return 0;
}
