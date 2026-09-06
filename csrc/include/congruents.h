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
#ifdef __cplusplus
}
#endif
#endif

