#ifndef CONGRUENTS_API_H
#define CONGRUENTS_API_H
#include <stddef.h>
#if defined(__GNUC__)
#define CG_API __attribute__((visibility("default")))
#else
#define CG_API
#endif
#ifdef __cplusplus
extern "C" {
#endif
/* ABI 1. All energies TOTAL GeV; densities cm^-3; outputs negative GeV/s.
 * Caller owns buffers, which must be distinct and have count doubles.
 * Context owns only configuration, never caller arrays. Calls borrow context.
 * No calls may race context destruction; NULL destruction is safe.
 * Invalid/destroyed non-NULL pointers are caller errors (not detectable).
 */
typedef struct cg_context cg_context;
enum cg_status { CG_OK=0, CG_INVALID=1, CG_ALLOC=2, CG_NUMERIC=3 };
CG_API const char *cg_version(void);
CG_API unsigned cg_abi_version(void);
CG_API const char *cg_status_message(int status); /* borrowed static string */
CG_API int cg_openmp_enabled(void);
CG_API int cg_context_create(int threads, cg_context **out);
CG_API void cg_context_destroy(cg_context *context);
/* Invalid arguments leave output unchanged. Numerical overflow can leave
 * partially calculated output; discard the entire output on any error.
 * count=0 permits NULL buffers. Nonzero count requires valid buffers.
 */
CG_API int cg_ionisation(const cg_context *context, size_t count,
                        const double *energy, double density, double *loss);
/* Week 2 additions, ABI 1 remains backwards compatible.
 * Galaxy input rows: z, Mstar [Msun], Re [kpc], SFR [Msun/yr].
 * Output rows (10): h [pc], nH [cm^-3], B [G], sigma_g [km/s],
 * area [pc²], Sigma_g [Msun/pc²], Sigma_SFR [Msun/yr/pc²],
 * Sigma_star [Msun/pc²], Tdust [K], Bhalo [G].
 */
CG_API int cg_galaxies(const cg_context *, size_t count,
                      const double *rows, double *properties);
CG_API int cg_preparation_bounds(const double *rows, size_t count,
                                double config[8], double temperatures[3]);
/* energies GeV; fields output row-major [energy][7], cm^-3 GeV^-1:
 * CMB,FIR,3000,4000,7500,UV,total. urad[8] matches Urad_Ub.txt in eV/cm³.
 */
CG_API int cg_radiation(const double row[4], size_t count,
                       const double *energies, const double photon_bounds[2],
                       double *fields, double urad[8]);
typedef struct cg_table cg_table;
/* kind: 0 IC emission, 1 IC Gamma, 2 BS, 3 SY.
 * field: 0=3000,1=4000,2=7500,3=UV,4=CMB,5=FIR (ignored for BS/SY).
 * config[8]: photon-output min/max, total-electron min/max,
 * target-photon min/max (GeV), synchrotron x min/max (dimensionless).
 * n sizes >=2; SY stores ny=1. Gamma x is DeltaE, not output energy.
 */
CG_API int cg_table_create(int kind, int field, double temperature,
                          size_t nx, size_t ny, const double config[8], cg_table **out);
CG_API int cg_table_import(size_t nx, size_t ny, const double *x,
                          const double *y, const double *z, cg_table **out);
/* Borrow without copying. Caller owns all three native contiguous double
 * buffers and keeps them alive/unchanged until handle destruction. C never
 * modifies or frees them. Buffer lengths are nx, ny, nx*ny, respectively.
 * Unlike import/create, destroy frees only the small descriptor.
 */
CG_API int cg_table_borrow(size_t nx, size_t ny, double *x,
                          double *y, double *z, cg_table **out);
CG_API void cg_table_destroy(cg_table *);
CG_API int cg_table_shape(const cg_table *, size_t *nx, size_t *ny);
CG_API int cg_table_copy(const cg_table *, double *x, double *y, double *z);
/* Linear/bilinear interpolation in native coordinates; zero outside domain. */
CG_API int cg_table_eval(const cg_table *, size_t n, const double *x,
                        const double *y, double *out);
/* Existing C field-grid count/bounds policy; fails safely for fewer than
 * two planes instead of entering legacy unsigned-loop underflow. */
CG_API int cg_temperature_grid(int field, double low, double high, double step,
                               size_t capacity, double *values, size_t *count);
/* Eight planes: 3000,4000,7500,UV,CMB lower/upper,FIR lower/upper.
 * Preserve legacy lower*frac+upper*(1-frac); output independently owned.
 */
CG_API int cg_combine_ic(const double row[4], cg_table *const planes[8],
                        double cmb_fraction, double fir_fraction, cg_table **out);
#ifdef __cplusplus
}
#endif
#endif
