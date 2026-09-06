#ifndef CONGRUENTS_GSL_DECS_H
#define CONGRUENTS_GSL_DECS_H

/* Reconstructed ownership wrappers around GSL's 1-D and 2-D splines. */

#include <gsl/gsl_spline.h>
#include <gsl/gsl_spline2d.h>
#include <stddef.h>
#include <stdio.h>
#include <stdlib.h>
#include <execinfo.h>
#include <unistd.h>
#include <dlfcn.h>

typedef struct {
    gsl_spline *spline;
    gsl_interp_accel *acc;
    double x_lim[2];
} gsl_spline_object_1D;

typedef struct {
    gsl_spline2d *spline;
    gsl_interp_accel *xacc;
    gsl_interp_accel *yacc;
    double x_lim[2];
    double y_lim[2];
} gsl_spline_object_2D;

static inline gsl_spline_object_1D gsl_so1D(
    size_t n, const double *x, const double *z)
{
    gsl_spline_object_1D object;
    object.acc = gsl_interp_accel_alloc();
    object.spline = gsl_spline_alloc(gsl_interp_linear, n);
    gsl_spline_init(object.spline, x, z, n);
    object.x_lim[0] = x[0];
    object.x_lim[1] = x[n - 1];
    return object;
}

static inline double congruents_gsl_so1D_eval(
    gsl_spline_object_1D object, double x, const char *caller_file,
    int caller_line)
{
    if (x < object.x_lim[0] || x > object.x_lim[1])
    {
        fprintf(stderr,
            "CONGRuENTS interpolation range error at %s:%d: "
            "x=%.17g, range=[%.17g, %.17g]\n",
            caller_file, caller_line, x, object.x_lim[0], object.x_lim[1]);
        void *frames[64];
        int count = backtrace(frames, 64);
        fprintf(stderr, "Backtrace frame count: %d, caller=%p\n",
            count, __builtin_return_address(0));
        backtrace_symbols_fd(frames, count, STDERR_FILENO);
        abort();
    }
    return gsl_spline_eval(object.spline, x, object.acc);
}

#define gsl_so1D_eval(object, x) \
    congruents_gsl_so1D_eval((object), (x), __FILE__, __LINE__)

static inline void gsl_so1D_free(gsl_spline_object_1D object)
{
    gsl_spline_free(object.spline);
    gsl_interp_accel_free(object.acc);
}

static inline gsl_spline_object_2D gsl_so2D(
    size_t nx, size_t ny, const double *x, const double *y, const double *z)
{
    gsl_spline_object_2D object;
    object.xacc = gsl_interp_accel_alloc();
    object.yacc = gsl_interp_accel_alloc();
    object.spline = gsl_spline2d_alloc(gsl_interp2d_bilinear, nx, ny);
    gsl_spline2d_init(object.spline, x, y, z, nx, ny);
    object.x_lim[0] = x[0];
    object.x_lim[1] = x[nx - 1];
    object.y_lim[0] = y[0];
    object.y_lim[1] = y[ny - 1];
    return object;
}

static inline gsl_spline_object_2D gsl_so2D_temp(
    size_t nx, size_t ny, const double *x, const double *y,
    const double x_lim[2], const double y_lim[2], const double *z)
{
    gsl_spline_object_2D object = gsl_so2D(nx, ny, x, y, z);
    object.x_lim[0] = x_lim[0];
    object.x_lim[1] = x_lim[1];
    object.y_lim[0] = y_lim[0];
    object.y_lim[1] = y_lim[1];
    return object;
}

static inline double gsl_so2D_eval(
    gsl_spline_object_2D object, double x, double y)
{
    if (x < object.x_lim[0] || x > object.x_lim[1]
        || y < object.y_lim[0] || y > object.y_lim[1]) return 0.0;
    return gsl_spline2d_eval(object.spline, x, y, object.xacc, object.yacc);
}

static inline void gsl_so2D_free(gsl_spline_object_2D object)
{
    gsl_spline2d_free(object.spline);
    gsl_interp_accel_free(object.xacc);
    gsl_interp_accel_free(object.yacc);
}

#endif
