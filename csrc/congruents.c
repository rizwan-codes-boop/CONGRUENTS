#include "include/congruents.h"
#include <stdlib.h>
#include <math.h>
#include <stdint.h>


#include "internal.h"
const char *cg_version(void) { return "0.2.0"; }
unsigned cg_abi_version(void) { return 2; }
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
    case CG_INVALID: return "invalid argument: check context, buffers and galaxy inputs";
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
