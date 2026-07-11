/*
* Author: Rasheed Ajala
* Contact: abioyeajala@gmail.com
*
*/

/*
* Window function operates on the domain.
* Taper functions operate on the support.
* The support is a subset of the domain.
*/

#include "blend.h"

static const blend_window_function blend_window_functions[] = {
    WFUNC_BOXCAR,
    WFUNC_COSINE,
    WFUNC_TRAPEZOID
};

const char *blend_window_function_name(blend_window_function function)
{
    switch (function) {
        case WFUNC_BOXCAR:
            return "boxcar";
        case WFUNC_COSINE:
            return "cosine";
        case WFUNC_TRAPEZOID:
            return "trapezoid";
        case WFUNC_INVALID:
        default:
            return "invalid";
    }
}

int blend_window_function_from_name(const char *name, blend_window_function *function)
{
    if (name == NULL || function == NULL) {
        return FAIL;
    }

    if (strcmp(name, "boxcar") == 0) {
        *function = WFUNC_BOXCAR;
    }
    else if (strcmp(name, "cosine") == 0) {
        *function = WFUNC_COSINE;
    }
    else if (strcmp(name, "trapezoid") == 0) {
        *function = WFUNC_TRAPEZOID;
    }
    else {
        *function = WFUNC_INVALID;
        return FAIL;
    }

    return SUCCESS;
}

void blend_print_window_function_names(FILE *fp)
{
    size_t i;

    for (i = 0; i < sizeof(blend_window_functions) / sizeof(blend_window_functions[0]); i++) {
        fprintf(fp, "%s\n", blend_window_function_name(blend_window_functions[i]));
    }
}

/* boxcar function */
double boxcar(void)
{
    return 1.0;
}

/* cosine function */
double cosine(int x, int n, int nmax, int nmin, double r1, double r2)
{
    double ni, nf;
    double value, t;

    ni = r1 * (double)nmax;
    nf = r2 * (double)nmax;

    if (ni > floor(nmin / 2)) {
        ni = r1 * (double)n;
    }
    if (nf > floor(nmin / 2)) {
        nf = r2 * (double)n;
    }

    if (x > ni && x <= n - nf) {
        value = 1.0;
    }
    else if (x <= ni) {
        t = ((2.0 * (double)x) - 1.0) / (2.0 * ni);
        value = 0.5 * (1.0 - cos(M_PI * t));
    }
    else {
        t = (2.0 * ((double)x - (double)n + 2.0 * nf) - 1.0) / (2.0 * nf);
        value = 0.5 * (1.0 - cos(M_PI * t));
    }

    return value;
}

/* trapezoid function */
double trapezoid(int x, int n, int nmax, int nmin, double r1, double r2)
{
    double ni, nf;
    double value, t;

    ni = r1 * (double)nmax;
    nf = r2 * (double)nmax;

    if (ni > floor(nmin / 2)) {
        ni = r1 * (double)n;
    }
    if (nf > floor(nmin / 2)) {
        nf = r2 * (double)n;
    }

    if (x > ni && x <= n - nf) {
        value = 1.0;
    }
    else if (x <= ni) {
        t = ((2.0 * (double)x) - 1.0) / (2.0 * ni);
        value = t;
    }
    else {
        t = (2.0 * ((double)x - (double)n + 2.0 * nf) - 1.0) / (2.0 * nf);
        value = 2.0 - t;
    }

    return value;
}

/* Smoothing function */
double window_function(int x, int x1, int x2, int n, int nmin,
                       double r1, double r2, blend_window_function func)
{
    int n1, n2;
    double taper;

    if (r1 < 0.0 || r2 < 0.0 || r1 >= 0.5 || r2 >= 0.5) {
        BLEND_Report(BLEND_MSG_ERROR, "window: taper ratios must be >= 0 and < 0.5\n");
        return WFUNC_ERROR;
    }

    if (x < 1 || x > n) {
        BLEND_Report(BLEND_MSG_ERROR, "window: query point is outside model domain\n");
        return WFUNC_ERROR;
    }

    n1 = x1 - 1;
    n2 = x2 - x1 + 1;

    if (x < x1) {
        return 0.0;
    }

    if (x > x2) {
        return 0.0;
    }

    x = x - n1;

    switch (func) {
        case WFUNC_BOXCAR:
            taper = boxcar();
            break;
        case WFUNC_COSINE:
            taper = cosine(x, n2, n, nmin, r1, r2);
            break;
        case WFUNC_TRAPEZOID:
            taper = trapezoid(x, n2, n, nmin, r1, r2);
            break;
        case WFUNC_INVALID:
        default:
            BLEND_Report(BLEND_MSG_ERROR, "window: requested window function does not exist\n");
            return WFUNC_ERROR;
    }

    return taper;
}
