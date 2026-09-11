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
/* ABI 2: native galaxy loops only. Serial table/diagnostic APIs were removed.
 * Caller owns row-major float64 buffers, alive and unchanged through the call.
 * All conversion/allocation happens before entry to the parallel region.
 */
typedef struct cg_context cg_context;
enum cg_status { CG_OK=0, CG_INVALID=1, CG_ALLOC=2, CG_NUMERIC=3 };
CG_API const char *cg_version(void);
CG_API unsigned cg_abi_version(void);
CG_API const char *cg_status_message(int status);
CG_API int cg_openmp_enabled(void);
CG_API int cg_context_create(int threads,cg_context **out);
CG_API void cg_context_destroy(cg_context *context);
/* rows[count][4]: z, Mstar Msun, Re kpc, SFR Msun/yr.
 * output[count][10]: h pc, nH cm^-3, B G, sigma_g km/s, area pc²,
 * Sigma_g Msun/pc², Sigma_SFR Msun/yr/pc², Sigma_star Msun/pc², Tdust K,
 * Bhalo G. Distinct valid buffers; discard output on any nonzero status.
 */
CG_API int cg_galaxies(const cg_context *,size_t count,const double *rows,double *output);
#ifdef __cplusplus
}
#endif
#endif
