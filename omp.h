#ifndef CONGRUENTS_OMP_COMPAT_H
#define CONGRUENTS_OMP_COMPAT_H

/*
 * Serial OpenMP compatibility layer for compilers without libomp.
 *
 * The scientific expressions are unchanged; OpenMP pragmas are ignored by
 * the compiler and these small query functions report a single worker.  A
 * true parallel build can instead define CONGRUENTS_USE_SYSTEM_OPENMP and
 * supply an OpenMP-capable compiler/runtime.
 */
#ifdef CONGRUENTS_USE_SYSTEM_OPENMP
#  include_next <omp.h>
#else
static inline int omp_get_max_threads(void) { return 1; }
static inline int omp_get_num_threads(void) { return 1; }
static inline int omp_get_thread_num(void) { return 0; }
static inline void omp_set_num_threads(int n) { (void)n; }
#endif

#endif
