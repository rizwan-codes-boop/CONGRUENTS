#ifndef CONGRUENTS_MATH_FUNCS_H
#define CONGRUENTS_MATH_FUNCS_H

#include <math.h>
#include <stddef.h>

static inline int logspace_array(size_t n, double low, double high, double *out)
{
    size_t i;
    if (!out || n < 2 || low <= 0.0 || high <= low) return 1;
    const double lo = log(low);
    const double step = (log(high) - lo) / (double)(n - 1);
    for (i = 0; i < n; ++i) out[i] = exp(lo + step * (double)i);
    return 0;
}

static inline int linspace_array(size_t n, double low, double high, double *out)
{
    size_t i;
    if (!out || n < 2 || high <= low) return 1;
    const double step = (high - low) / (double)(n - 1);
    for (i = 0; i < n; ++i) out[i] = low + step * (double)i;
    return 0;
}

static inline double minval(size_t n, const double *values)
{
    size_t i;
    double result = values[0];
    for (i = 1; i < n; ++i) if (values[i] < result) result = values[i];
    return result;
}

static inline double maxval(size_t n, const double *values)
{
    size_t i;
    double result = values[0];
    for (i = 1; i < n; ++i) if (values[i] > result) result = values[i];
    return result;
}

static inline void double_integer_array(int low, int high, double *out)
{
    int i;
    for (i = low; i <= high; ++i) out[i - low] = (double)i;
}

#endif
