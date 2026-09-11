#ifndef CG_SOLVER_H
#define CG_SOLVER_H
#include "congruents.h"
typedef struct {
    size_t nx, ny;
    const double *x, *y, *values;
} cg_table_input;
/* Synchronous borrowed float64 buffers, row-major, caller-owned.
 * props[n][7]: h pc, nH cm^-3, B G, Bhalo G, SigmaSFR Msun/yr/pc², Re kpc, SFR Msun/yr.
 * diffusion[n][2][ne] cm²/s: disc, halo.
 * injection[n][2][ne] GeV^-1 s^-1: primary, secondary disc injection.
 * electrons[n][4][ne] GeV^-1: primary/secondary disc, primary/secondary halo.
 * emission[n][15][np] GeV^-1 s^-1, except dimensionless tau_FF.
 * Order: IC1d,IC2d,BS1d,BS2d,SY1d,SY2d,IC1h,IC2h,SY1h,SY2h,FF,tau_FF,pi,pi_fcal1,nu.
 * statuses[n]: discard a galaxy's output unless its status is CG_OK.
 * Independent preparation API remains ABI 2; this optional solver ABI is 1.
 */
CG_API unsigned cg_solver_abi(void);
/* transport[n][7][ne]: fcal, Dp, De_disc, De_halo, Q1_disc, Q2_disc, proton_SS.
 * rows[n][4] and properties[n][10] follow cg_galaxies. Catalogue order is
 * significant: the reference's index-10 calorimetry experiment is preserved.
 */
CG_API int cg_transport_batch(int threads,size_t n,size_t ne,const double *kinetic,
    const double *electron,const double *rows,const double *properties,
    double *transport,int *statuses);
CG_API int cg_solver_batch(int threads,size_t n,size_t ne,size_t np,size_t ns,
    const double *electron_energy,const double *photon_energy,const double *props,
    const double *diffusion,const double *injection,const double *kinetic,const double *fcal,const cg_table_input *ic,
    const cg_table_input *gamma,const cg_table_input *bs,const cg_table_input *sy,
    double *electrons,double *emission,int *statuses);
#endif
