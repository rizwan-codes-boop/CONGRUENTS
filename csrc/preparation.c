/* Native galaxy-loop helpers only; no serial table generation in this library. */
#include "include/congruents.h"
#include "internal.h"
#include <string.h>
#include <math.h>
#include "../physical_constants.h"
static double Sigma_gas_Shi_iKS__Msolpcm2(double sf,double stars) {
    return pow(10.,10.28)*sf*pow(stars,-.48);
}
static double sigma_gas_Yu__kmsm1(double sfr) {
    return pow(10.,.2*log10(sfr)+1.6);
}
static double Tdust__K(double z,double sfr,double mass) {
    return 98.*pow(1.+z,-.065)+6.9*log10(sfr/mass);
}
static double sigma_star_Bezanson__kmsm1(double mass,double radius) {
    double kn=73.32/(10.465+pow(1.-.94,2))+.954;
    return sqrt(G__pcMsolm1km2sm2*mass/(.557*kn*radius*1e3));
}
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
    const int threads=ctx->threads;
    #pragma omp parallel for num_threads(threads) reduction(|:bad)
    for(size_t i=0;i<n;i++) bad |= derive(rows+4*i,out+10*i);
    return bad ? CG_NUMERIC : CG_OK;
}
