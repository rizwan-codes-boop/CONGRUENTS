/* Independent legacy reference: no cg_* API or adapter source is used. */
#include <string.h>
#include <stdlib.h>
#include "../spectra_funcs.h"
#include "../data_calc.h"
double q_p_inject=2.2,q_e_inject=2.2;
double T_CR_lims__GeV[2]={1.e-3,1.e8};
double E_CRe_lims__GeV[2]={1.e-3+5.10998950e-4,1.e8+5.10998950e-4};

int main(int argc,char **argv) {
    if(argc!=14 && argc!=21) return 2;
    int kind=atoi(argv[1]),field=atoi(argv[2]);
    double temp=atof(argv[3]);
    size_t dims[2]={(size_t)atoi(argv[4]),(size_t)atoi(argv[5])};
    double c[8];
    for(int i=0;i<8;i++) c[i]=atof(argv[i+6]);
    double eg[2]={c[0],c[1]},ee[2]={c[2],c[3]},ep[2]={c[4],c[5]},sx[2]={c[6],c[7]};
    typedef data_object_2D (*builder)(double (*)(double *,double),double *,double *,double *,double *,size_t *);
    builder make=kind==0 ? init_do_2D_IC : init_do_2D_IC_Gamma;
    typedef double (*photon)(double *,double);
    photon funcs[6]={dndEphot_BB__cmm3GeVm1,dndEphot_BB__cmm3GeVm1,
        dndEphot_BB__cmm3GeVm1,dndEphot_UVMattis__cmm3GeVm1,
        dndEphot_BB__cmm3GeVm1,dndEphot_modBB__cmm3GeVm1};
    double temps[6]={3000,4000,7500,0,temp,temp};
    if(kind==3) {
        data_object_1D a=init_do_1D_sync(sx,dims[0]);
        printf("%zu 1\n",a.nx);
        for(size_t i=0;i<a.nx;i++) printf("%.17g ",a.x_data[i]);
        printf("1 ");
        for(size_t i=0;i<a.nx;i++) printf("%.17g ",a.z_data[i]);
        data_object_1D_free(a);return 0;
    }
    if(field==6) {
        IC_object obj={0};
        for(int f=0;f<4;f++) obj.do_2D_IC[kind][f]=make(funcs[f],&temps[f],eg,ee,ep,dims);
        obj.do_3D_IC[kind][0]=do3D_IC(make,funcs[4],CMB_num,CMB_lim,0,atof(argv[18]),.5,eg,ee,ep,dims);
        obj.do_3D_IC[kind][1]=do3D_IC(make,funcs[5],FIR_num,FIR_lim,atof(argv[19]),atof(argv[20]),5,eg,ee,ep,dims);
        double z=atof(argv[14]),mass=atof(argv[15]),radius=atof(argv[16]),sfr=atof(argv[17]);
        /* h is unused in C_dil, hence any positive height gives the same IC. */
        double p[7]={T_0_CMB__K*(1+z),Tdust__K(z,sfr,mass),mass,sfr,radius,1,0};
        gsl_spline_object_2D result=construct_IC_gso2D(kind,obj,p);
        data_object_3D a=obj.do_3D_IC[kind][0];
        printf("%zu %zu\n",a.nx,a.ny);
        for(size_t i=0;i<a.nx;i++) printf("%.17g ",a.x_data[i]);
        for(size_t j=0;j<a.ny;j++) printf("%.17g ",a.y_data[j]);
        for(size_t j=0;j<a.ny;j++) for(size_t i=0;i<a.nx;i++)
            printf("%.17g ",gsl_so2D_eval(result,a.x_data[i],a.y_data[j]));
        gsl_so2D_free(result);
        for(int f=0;f<4;f++) data_object_2D_free(obj.do_2D_IC[kind][f]);
        for(int f=0;f<2;f++) data_object_3D_free(obj.do_3D_IC[kind][f]);
        return 0;
    }
    data_object_2D a=kind==2 ? init_do_2D_BS(eg,ee,dims) : make(funcs[field],&temps[field],eg,ee,ep,dims);
    printf("%zu %zu\n",a.nx,a.ny);
    for(size_t i=0;i<a.nx;i++) printf("%.17g ",a.x_data[i]);
    for(size_t j=0;j<a.ny;j++) printf("%.17g ",a.y_data[j]);
    for(size_t i=0;i<a.nx*a.ny;i++) printf("%.17g ",a.z_data[i]);
    data_object_2D_free(a);
    return 0;
}
