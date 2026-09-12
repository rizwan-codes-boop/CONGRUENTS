#ifndef CG_SOLVER_H
#define CG_SOLVER_H
#include "congruents.h"
typedef struct {
    size_t nx, ny;
    const double *x, *y, *values;
} cg_table_input;
/* Synchronous borrowed float64 buffers, row-major, caller-owned.
 * props[n][8]: h pc, nH cm^-3, B G, Bhalo G, SigmaSFR Msun/yr/pc², Re kpc, SFR Msun/yr, proton normalization Cp.
 * diffusion[n][2][ne] cm²/s: disc, halo.
 * injection[n][2][ne] GeV^-1 s^-1: primary, secondary disc injection.
 * electrons[n][4][ne] GeV^-1: primary/secondary disc, primary/secondary halo.
 * emission[n][15][np] GeV^-1 s^-1: raw, before Python disc attenuation.
 * Order: IC1d,IC2d,BS1d,BS2d,SY1d,SY2d,IC1h,IC2h,SY1h,SY2h,FF,tau_FF,pi,pi_fcal1,nu.
 * statuses[n]: discard a galaxy's output unless its status is CG_OK.
 * Optional radio[n][4]: raw SY primary/secondary disc/halo at 1.49 GHz.
 * Independent preparation API remains ABI 2; this optional solver ABI is 3.
 */
CG_API unsigned cg_solver_abi(void);
/* Secondary injection only: density[n], Cp[n], fcal[n][ne], output[n][ne].
 * Python prepares all other transport quantities before this call.
 * Solver emission slots FF and tau_FF are reserved zeros; Python fills them
 * and applies disc free-free attenuation after the native batch returns.
 */
CG_API int cg_transport_batch(int threads,size_t n,size_t ne,const double *kinetic,
    const double *density,const double *proton_normalization,const double *fcal,
    double *secondary_injection,int *statuses);
CG_API int cg_solver_batch(int threads,size_t n,size_t ne,size_t np,size_t ns,
    const double *electron_energy,const double *photon_energy,const double *props,
    const double *diffusion,const double *injection,const double *kinetic,const double *fcal,const cg_table_input *ic,
    const cg_table_input *gamma,const cg_table_input *bs,const cg_table_input *sy,
    double *electrons,double *emission,double *radio,int *statuses);
#endif
