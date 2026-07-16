/*
* Author: Rasheed Ajala
* Contact: abioyeajala@gmail.com
*/

#ifndef BLEND_WINDOW_H
#define BLEND_WINDOW_H

#include "blend_error.h"
#include "blend_polygon.h"

#include <stdio.h>

/* Window functions */
typedef enum blend_window_function {
    WFUNC_BOXCAR = 0,
    WFUNC_COSINE,
    WFUNC_TRAPEZOID,
    WFUNC_HAMMING,
    WFUNC_BLACKMAN,
    WFUNC_BLACKMANHARRIS,
    WFUNC_WELCH,
    WFUNC_PARZEN,
    WFUNC_GAUSSIAN,
    WFUNC_SMOOTHSTEP,
    WFUNC_SMOOTHERSTEP,
    WFUNC_EXPONENTIAL,
    WFUNC_SINE,
    WFUNC_BOHMAN,
    WFUNC_NUTTALL,
    WFUNC_KAISER,
    WFUNC_CAUCHY,
    WFUNC_QUADRATIC,
    WFUNC_CUBIC,
    WFUNC_POISSON,
    WFUNC_BARTLETT,
    WFUNC_BARTLETTHANN,
    WFUNC_EXACTBLACKMAN,
    WFUNC_BLACKMANNUTTALL,
    WFUNC_FLATTOP,
    WFUNC_LANCZOS,
    WFUNC_RIESZ,
    WFUNC_RIEMANN,
    WFUNC_FEJER,
    WFUNC_CONNES,
    WFUNC_HANNINGPOISSON,
    WFUNC_KAISERBESSEL,
    WFUNC_PLANCKTAPER,
    WFUNC_QUARTIC,
    WFUNC_QUINTIC,
    WFUNC_SEPTIC,
    WFUNC_NONIC,
    WFUNC_LOGISTIC,
    WFUNC_TANH,
    WFUNC_ERF,
    WFUNC_ARCTAN,
    WFUNC_GOMPERTZ,
    WFUNC_SOFTSIGN,
    WFUNC_AGNESI,
    WFUNC_INVERSEQUADRATIC,
    WFUNC_INVERSEMULTIQUADRIC,
    WFUNC_POWERLAW,
    WFUNC_ROOT,
    WFUNC_CIRCULAR,
    WFUNC_SECH,
    WFUNC_SECH2,
    WFUNC_STUDENT,
    WFUNC_LAPLACE,
    WFUNC_INVALID
} blend_window_function;

/* Supported assembly type */
typedef enum { ROW_STEP = 0,
               COLUMN_STEP } assembly_type;

/* Data structure for boundary and vector assembly */
typedef struct window {
    int row_size;
    double **vertices;
    assembly_type amode;
    int nx;
    int ny;
    int nz;
    int *nnx1;
    int *nnx2;
    int *nny1;
    int *nny2;
    int minnx;
    int minny;
    double ratio_x1;
    double ratio_x2;
    double ratio_y1;
    double ratio_y2;
    double ratio_z1;
    double ratio_z2;
    blend_window_function x_function;
    blend_window_function y_function;
    blend_window_function z_function;
    double contribution;
} window;

/* Window function name */
const char *blend_window_function_name(blend_window_function function);

/* Window function lookup */
int blend_window_function_from_name(const char *name, blend_window_function *function);

/* Print supported window function names */
void blend_print_window_function_names(FILE *fp);

/* Assign polygon support to a window */
int blend_window_set_polygon(window *data, const polygon *poly);

/* Clear polygon support assigned to a window */
void blend_window_clear_polygon(window *data);

/* Window/Taper function */
double window_function(int x, int x1, int x2, int n, int nmin,
                       double r1, double r2, blend_window_function func);

double boxcar(void);
double cosine(int x, int n, int nmax, int nmin, double r1, double r2);
double trapezoid(int x, int n, int nmax, int nmin, double r1, double r2);
double hamming(int x, int n, int nmax, int nmin, double r1, double r2);
double blackman(int x, int n, int nmax, int nmin, double r1, double r2);
double blackmanharris(int x, int n, int nmax, int nmin, double r1, double r2);
double welch(int x, int n, int nmax, int nmin, double r1, double r2);
double parzen(int x, int n, int nmax, int nmin, double r1, double r2);
double gaussian(int x, int n, int nmax, int nmin, double r1, double r2);
double smoothstep(int x, int n, int nmax, int nmin, double r1, double r2);
double smootherstep(int x, int n, int nmax, int nmin, double r1, double r2);
double exponential(int x, int n, int nmax, int nmin, double r1, double r2);
double sine(int x, int n, int nmax, int nmin, double r1, double r2);
double bohman(int x, int n, int nmax, int nmin, double r1, double r2);
double nuttall(int x, int n, int nmax, int nmin, double r1, double r2);
double kaiser(int x, int n, int nmax, int nmin, double r1, double r2);
double cauchy(int x, int n, int nmax, int nmin, double r1, double r2);
double quadratic(int x, int n, int nmax, int nmin, double r1, double r2);
double cubic(int x, int n, int nmax, int nmin, double r1, double r2);
double poisson(int x, int n, int nmax, int nmin, double r1, double r2);
double bartlett(int x, int n, int nmax, int nmin, double r1, double r2);
double bartletthann(int x, int n, int nmax, int nmin, double r1, double r2);
double exactblackman(int x, int n, int nmax, int nmin, double r1, double r2);
double blackmannuttall(int x, int n, int nmax, int nmin, double r1, double r2);
double flattop(int x, int n, int nmax, int nmin, double r1, double r2);
double lanczos(int x, int n, int nmax, int nmin, double r1, double r2);
double riesz(int x, int n, int nmax, int nmin, double r1, double r2);
double riemann(int x, int n, int nmax, int nmin, double r1, double r2);
double fejer(int x, int n, int nmax, int nmin, double r1, double r2);
double connes(int x, int n, int nmax, int nmin, double r1, double r2);
double hanningpoisson(int x, int n, int nmax, int nmin, double r1, double r2);
double kaiserbessel(int x, int n, int nmax, int nmin, double r1, double r2);
double plancktaper(int x, int n, int nmax, int nmin, double r1, double r2);
double quartic(int x, int n, int nmax, int nmin, double r1, double r2);
double quintic(int x, int n, int nmax, int nmin, double r1, double r2);
double septic(int x, int n, int nmax, int nmin, double r1, double r2);
double nonic(int x, int n, int nmax, int nmin, double r1, double r2);
double logistic(int x, int n, int nmax, int nmin, double r1, double r2);
double tanhwindow(int x, int n, int nmax, int nmin, double r1, double r2);
double erfwindow(int x, int n, int nmax, int nmin, double r1, double r2);
double arctanwindow(int x, int n, int nmax, int nmin, double r1, double r2);
double gompertz(int x, int n, int nmax, int nmin, double r1, double r2);
double softsign(int x, int n, int nmax, int nmin, double r1, double r2);
double agnesi(int x, int n, int nmax, int nmin, double r1, double r2);
double inversequadratic(int x, int n, int nmax, int nmin, double r1, double r2);
double inversemultiquadric(int x, int n, int nmax, int nmin, double r1, double r2);
double powerlaw(int x, int n, int nmax, int nmin, double r1, double r2);
double root(int x, int n, int nmax, int nmin, double r1, double r2);
double circular(int x, int n, int nmax, int nmin, double r1, double r2);
double sechwindow(int x, int n, int nmax, int nmin, double r1, double r2);
double sech2window(int x, int n, int nmax, int nmin, double r1, double r2);
double student(int x, int n, int nmax, int nmin, double r1, double r2);
double laplace(int x, int n, int nmax, int nmin, double r1, double r2);

#endif /* BLEND_WINDOW_H */
