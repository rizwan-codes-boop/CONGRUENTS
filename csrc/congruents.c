#include "include/congruents.h"
#include <stdlib.h>
#include <math.h>
#include <stdint.h>
#include "../CR_spectra/ionisation.h"

#include "internal.h"
const char *cg_version(void) { return "0.1.0"; }
unsigned cg_abi_version(void) { return 1; }
int cg_openmp_enabled(void) {
#ifdef _OPENMP
    return 1;
#else
    return 0;
#endif
}
const char *cg_status_message(int status) {
    switch (status) {
    case CG_OK: return "success";
    case CG_INVALID: return "invalid argument: check context, buffers, total energy and density";
    case CG_ALLOC: return "context allocation failed";
    case CG_NUMERIC: return "non-finite numerical result";
    default: return "unknown status";
    }
}
int cg_context_create(int threads, cg_context **out) {
    if (!out) return CG_INVALID;
    *out = NULL;
    if (threads < 1) return CG_INVALID;
    cg_context *context = malloc(sizeof(*context));
    if (!context) return CG_ALLOC;
    context->threads = threads;
    *out = context;
    return CG_OK;
}
void cg_context_destroy(cg_context *context) { free(context); }
int cg_ionisation(const cg_context *context, size_t count,
                  const double *energy, double density, double *loss) {
    if (!context || !isfinite(density) || density < 0 ||
        count > (size_t) PTRDIFF_MAX ||
        (count && (!energy || !loss || energy == loss))) return CG_INVALID;
    for (size_t i=0; i<count; ++i)
        if (!isfinite(energy[i]) || energy[i] <= m_e__GeV) return CG_INVALID;
    int bad=0;
    #pragma omp parallel for num_threads(context->threads) reduction(|:bad)
    for (size_t i=0; i<count; ++i) {
        loss[i] = dEdtm1_ion__GeVsm1(energy[i], density);
        bad |= !isfinite(loss[i]);
    }
    return bad ? CG_NUMERIC : CG_OK;
}
