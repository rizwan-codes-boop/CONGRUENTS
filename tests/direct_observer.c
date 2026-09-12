/* Independent distance function and GSL linear-resampling reference.
 * This fixture tests the output transform, not the source solver or EBL model. */
#include <stdio.h>
#include <stdlib.h>
#include <gsl/gsl_spline.h>
#include <gsl/gsl_interp2d.h>
#include "../cosmo_funcs.h"
int main(int argc,char **argv) {
    if(argc>=4) {
        FILE *f=fopen(argv[1],"r");
        if(!f) return 2;
        size_t nx,ny;
        fscanf(f,"%*[^\n]\n");fscanf(f,"%zu\n",&ny);
        fscanf(f,"%*[^\n]\n");fscanf(f,"%zu\n",&nx);
        fscanf(f,"%*[^\n]\n");
        double x[nx],y[ny],a[nx*ny];
        for(size_t i=0;i<ny;i++) fscanf(f,"%lf",y+i);
        fscanf(f,"\n%*[^\n]\n");
        for(size_t i=0;i<nx;i++) fscanf(f,"%lf",x+i);
        fscanf(f,"\n%*[^\n]\n");
        for(size_t i=0;i<nx*ny;i++) fscanf(f,"%lf",a+i);
        fclose(f);
        gsl_interp2d *interp=gsl_interp2d_alloc(gsl_interp2d_bilinear,nx,ny);
        gsl_interp_accel *xa=gsl_interp_accel_alloc(),*ya=gsl_interp_accel_alloc();
        gsl_interp2d_init(interp,x,y,a,nx,ny);
        double z=atof(argv[2]);
        for(int j=3;j<argc;j++) printf("%.17g ",fmax(0.,gsl_interp2d_eval_extrap(interp,x,y,a,
            fmin(atof(argv[j])/(1+z)*1e9,1e15),z,xa,ya)));
        printf("\n");gsl_interp2d_free(interp);gsl_interp_accel_free(xa);gsl_interp_accel_free(ya);
        return 0;
    }
    if(argc!=2) return 2;
    double z=atof(argv[1]),e[4]={1.,2.,4.,8.},shifted[4],rate[4]={3.,7.,11.,19.};
    double factor=pow(1+z,2)/(4*M_PI*pow(d_l_MPc(z)*Mpc__cm,2));
    printf("%.17g %.17g",d_l_MPc(z),factor);
    for(int i=0;i<4;i++) shifted[i]=e[i]/(1+z);
    gsl_interp_accel *acc=gsl_interp_accel_alloc();
    gsl_spline *spline=gsl_spline_alloc(gsl_interp_linear,4);
    gsl_spline_init(spline,shifted,rate,4);
    for(int j=0;j<4;j++) {
        double value=e[j]<shifted[0] || e[j]>shifted[3] ? 0. : gsl_spline_eval(spline,e[j],acc);
        printf(" %.17g",value*exp(-.1*j)*exp(-.2*j)*factor*pow(e[j],2));
    }
    printf("\n");gsl_spline_free(spline);gsl_interp_accel_free(acc);
    return 0;
}
