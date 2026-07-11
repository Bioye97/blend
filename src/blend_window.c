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

#include <float.h>

static const char *blend_window_function_names[] = {
    "boxcar",
    "cosine",
    "hann",
    "tukey",
    "trapezoid",
    "linear",
    "hamming",
    "blackman",
    "blackmanharris",
    "welch",
    "parzen",
    "gaussian",
    "smoothstep",
    "smootherstep",
    "exponential",
    "sine",
    "bohman",
    "nuttall",
    "kaiser",
    "cauchy",
    "quadratic",
    "cubic",
    "poisson",
    "bartlett",
    "barthann",
    "bartletthann",
    "exactblackman",
    "blackmannuttall",
    "flattop",
    "lanczos",
    "riesz",
    "riemann",
    "fejer",
    "connes",
    "hanningpoisson",
    "kaiserbessel",
    "plancktaper",
    "quartic",
    "quintic",
    "septic",
    "nonic",
    "logistic",
    "tanh",
    "erf",
    "arctan",
    "gompertz",
    "softsign",
    "agnesi",
    "inversequadratic",
    "inversemultiquadric",
    "powerlaw",
    "root",
    "circular",
    "sech",
    "sech2",
    "student",
    "laplace"
};

static double blend_clamp_unit(double value)
{
    if (value < 0.0) {
        return 0.0;
    }
    if (value > 1.0) {
        return 1.0;
    }
    return value;
}

const char *blend_window_function_name(blend_window_function function)
{
    switch (function) {
        case WFUNC_BOXCAR:
            return "boxcar";
        case WFUNC_COSINE:
            return "cosine";
        case WFUNC_TRAPEZOID:
            return "trapezoid";
        case WFUNC_HAMMING:
            return "hamming";
        case WFUNC_BLACKMAN:
            return "blackman";
        case WFUNC_BLACKMANHARRIS:
            return "blackmanharris";
        case WFUNC_WELCH:
            return "welch";
        case WFUNC_PARZEN:
            return "parzen";
        case WFUNC_GAUSSIAN:
            return "gaussian";
        case WFUNC_SMOOTHSTEP:
            return "smoothstep";
        case WFUNC_SMOOTHERSTEP:
            return "smootherstep";
        case WFUNC_EXPONENTIAL:
            return "exponential";
        case WFUNC_SINE:
            return "sine";
        case WFUNC_BOHMAN:
            return "bohman";
        case WFUNC_NUTTALL:
            return "nuttall";
        case WFUNC_KAISER:
            return "kaiser";
        case WFUNC_CAUCHY:
            return "cauchy";
        case WFUNC_QUADRATIC:
            return "quadratic";
        case WFUNC_CUBIC:
            return "cubic";
        case WFUNC_POISSON:
            return "poisson";
        case WFUNC_BARTLETT:
            return "bartlett";
        case WFUNC_BARTLETTHANN:
            return "bartletthann";
        case WFUNC_EXACTBLACKMAN:
            return "exactblackman";
        case WFUNC_BLACKMANNUTTALL:
            return "blackmannuttall";
        case WFUNC_FLATTOP:
            return "flattop";
        case WFUNC_LANCZOS:
            return "lanczos";
        case WFUNC_RIESZ:
            return "riesz";
        case WFUNC_RIEMANN:
            return "riemann";
        case WFUNC_FEJER:
            return "fejer";
        case WFUNC_CONNES:
            return "connes";
        case WFUNC_HANNINGPOISSON:
            return "hanningpoisson";
        case WFUNC_KAISERBESSEL:
            return "kaiserbessel";
        case WFUNC_PLANCKTAPER:
            return "plancktaper";
        case WFUNC_QUARTIC:
            return "quartic";
        case WFUNC_QUINTIC:
            return "quintic";
        case WFUNC_SEPTIC:
            return "septic";
        case WFUNC_NONIC:
            return "nonic";
        case WFUNC_LOGISTIC:
            return "logistic";
        case WFUNC_TANH:
            return "tanh";
        case WFUNC_ERF:
            return "erf";
        case WFUNC_ARCTAN:
            return "arctan";
        case WFUNC_GOMPERTZ:
            return "gompertz";
        case WFUNC_SOFTSIGN:
            return "softsign";
        case WFUNC_AGNESI:
            return "agnesi";
        case WFUNC_INVERSEQUADRATIC:
            return "inversequadratic";
        case WFUNC_INVERSEMULTIQUADRIC:
            return "inversemultiquadric";
        case WFUNC_POWERLAW:
            return "powerlaw";
        case WFUNC_ROOT:
            return "root";
        case WFUNC_CIRCULAR:
            return "circular";
        case WFUNC_SECH:
            return "sech";
        case WFUNC_SECH2:
            return "sech2";
        case WFUNC_STUDENT:
            return "student";
        case WFUNC_LAPLACE:
            return "laplace";
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
    else if (strcmp(name, "cosine") == 0 || strcmp(name, "hann") == 0 ||
             strcmp(name, "tukey") == 0) {
        *function = WFUNC_COSINE;
    }
    else if (strcmp(name, "trapezoid") == 0 || strcmp(name, "linear") == 0) {
        *function = WFUNC_TRAPEZOID;
    }
    else if (strcmp(name, "hamming") == 0) {
        *function = WFUNC_HAMMING;
    }
    else if (strcmp(name, "blackman") == 0) {
        *function = WFUNC_BLACKMAN;
    }
    else if (strcmp(name, "blackmanharris") == 0) {
        *function = WFUNC_BLACKMANHARRIS;
    }
    else if (strcmp(name, "welch") == 0) {
        *function = WFUNC_WELCH;
    }
    else if (strcmp(name, "parzen") == 0) {
        *function = WFUNC_PARZEN;
    }
    else if (strcmp(name, "gaussian") == 0) {
        *function = WFUNC_GAUSSIAN;
    }
    else if (strcmp(name, "smoothstep") == 0) {
        *function = WFUNC_SMOOTHSTEP;
    }
    else if (strcmp(name, "smootherstep") == 0) {
        *function = WFUNC_SMOOTHERSTEP;
    }
    else if (strcmp(name, "exponential") == 0) {
        *function = WFUNC_EXPONENTIAL;
    }
    else if (strcmp(name, "sine") == 0) {
        *function = WFUNC_SINE;
    }
    else if (strcmp(name, "bohman") == 0) {
        *function = WFUNC_BOHMAN;
    }
    else if (strcmp(name, "nuttall") == 0) {
        *function = WFUNC_NUTTALL;
    }
    else if (strcmp(name, "kaiser") == 0) {
        *function = WFUNC_KAISER;
    }
    else if (strcmp(name, "cauchy") == 0) {
        *function = WFUNC_CAUCHY;
    }
    else if (strcmp(name, "quadratic") == 0) {
        *function = WFUNC_QUADRATIC;
    }
    else if (strcmp(name, "cubic") == 0) {
        *function = WFUNC_CUBIC;
    }
    else if (strcmp(name, "poisson") == 0) {
        *function = WFUNC_POISSON;
    }
    else if (strcmp(name, "bartlett") == 0) {
        *function = WFUNC_BARTLETT;
    }
    else if (strcmp(name, "barthann") == 0 || strcmp(name, "bartletthann") == 0) {
        *function = WFUNC_BARTLETTHANN;
    }
    else if (strcmp(name, "exactblackman") == 0) {
        *function = WFUNC_EXACTBLACKMAN;
    }
    else if (strcmp(name, "blackmannuttall") == 0) {
        *function = WFUNC_BLACKMANNUTTALL;
    }
    else if (strcmp(name, "flattop") == 0) {
        *function = WFUNC_FLATTOP;
    }
    else if (strcmp(name, "lanczos") == 0) {
        *function = WFUNC_LANCZOS;
    }
    else if (strcmp(name, "riesz") == 0) {
        *function = WFUNC_RIESZ;
    }
    else if (strcmp(name, "riemann") == 0) {
        *function = WFUNC_RIEMANN;
    }
    else if (strcmp(name, "fejer") == 0) {
        *function = WFUNC_FEJER;
    }
    else if (strcmp(name, "connes") == 0) {
        *function = WFUNC_CONNES;
    }
    else if (strcmp(name, "hanningpoisson") == 0 || strcmp(name, "hannpoisson") == 0) {
        *function = WFUNC_HANNINGPOISSON;
    }
    else if (strcmp(name, "kaiserbessel") == 0) {
        *function = WFUNC_KAISERBESSEL;
    }
    else if (strcmp(name, "plancktaper") == 0 || strcmp(name, "planck") == 0) {
        *function = WFUNC_PLANCKTAPER;
    }
    else if (strcmp(name, "quartic") == 0) {
        *function = WFUNC_QUARTIC;
    }
    else if (strcmp(name, "quintic") == 0) {
        *function = WFUNC_QUINTIC;
    }
    else if (strcmp(name, "septic") == 0) {
        *function = WFUNC_SEPTIC;
    }
    else if (strcmp(name, "nonic") == 0) {
        *function = WFUNC_NONIC;
    }
    else if (strcmp(name, "logistic") == 0) {
        *function = WFUNC_LOGISTIC;
    }
    else if (strcmp(name, "tanh") == 0) {
        *function = WFUNC_TANH;
    }
    else if (strcmp(name, "erf") == 0) {
        *function = WFUNC_ERF;
    }
    else if (strcmp(name, "arctan") == 0) {
        *function = WFUNC_ARCTAN;
    }
    else if (strcmp(name, "gompertz") == 0) {
        *function = WFUNC_GOMPERTZ;
    }
    else if (strcmp(name, "softsign") == 0) {
        *function = WFUNC_SOFTSIGN;
    }
    else if (strcmp(name, "agnesi") == 0) {
        *function = WFUNC_AGNESI;
    }
    else if (strcmp(name, "inversequadratic") == 0) {
        *function = WFUNC_INVERSEQUADRATIC;
    }
    else if (strcmp(name, "inversemultiquadric") == 0) {
        *function = WFUNC_INVERSEMULTIQUADRIC;
    }
    else if (strcmp(name, "powerlaw") == 0) {
        *function = WFUNC_POWERLAW;
    }
    else if (strcmp(name, "root") == 0) {
        *function = WFUNC_ROOT;
    }
    else if (strcmp(name, "circular") == 0) {
        *function = WFUNC_CIRCULAR;
    }
    else if (strcmp(name, "sech") == 0) {
        *function = WFUNC_SECH;
    }
    else if (strcmp(name, "sech2") == 0) {
        *function = WFUNC_SECH2;
    }
    else if (strcmp(name, "student") == 0) {
        *function = WFUNC_STUDENT;
    }
    else if (strcmp(name, "laplace") == 0) {
        *function = WFUNC_LAPLACE;
    }
    else if (strcmp(name, "normal") == 0) {
        *function = WFUNC_GAUSSIAN;
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

    for (i = 0; i < sizeof(blend_window_function_names) / sizeof(blend_window_function_names[0]); i++) {
        fprintf(fp, "%s\n", blend_window_function_names[i]);
    }
}

/* boxcar function */
double boxcar(void)
{
    return 1.0;
}

static void blend_taper_lengths(int n, int nmax, int nmin, double r1, double r2,
                                double *ni, double *nf)
{
    *ni = r1 * (double)nmax;
    *nf = r2 * (double)nmax;

    if (*ni > floor(nmin / 2)) {
        *ni = r1 * (double)n;
    }
    if (*nf > floor(nmin / 2)) {
        *nf = r2 * (double)n;
    }
}

static double blend_taper_position(int x, int n, int nmax, int nmin,
                                   double r1, double r2, int *in_taper)
{
    double ni, nf;
    double j;

    blend_taper_lengths(n, nmax, nmin, r1, r2, &ni, &nf);

    if (x > ni && x <= n - nf) {
        *in_taper = 0;
        return 1.0;
    }

    *in_taper = 1;
    if (x <= ni) {
        if (ni <= 0.0) return 1.0;
        return ((2.0 * (double)x) - 1.0) / (2.0 * ni);
    }

    if (nf <= 0.0) return 1.0;
    j = (double)n - (double)x + 1.0;
    return ((2.0 * j) - 1.0) / (2.0 * nf);
}

static double blend_window_from_edge(int x, int n, int nmax, int nmin,
                                     double r1, double r2,
                                     double (*edge)(double))
{
    int in_taper = 0;
    double t = blend_taper_position(x, n, nmax, nmin, r1, r2, &in_taper);

    if (!in_taper) {
        return 1.0;
    }

    return blend_clamp_unit(edge(blend_clamp_unit(t)));
}

static double cosine_edge(double t)
{
    return 0.5 * (1.0 - cos(M_PI * t));
}

static double trapezoid_edge(double t)
{
    return t;
}

static double hamming_edge(double t)
{
    return 0.54 - 0.46 * cos(M_PI * t);
}

static double blackman_edge(double t)
{
    return 0.42 - 0.50 * cos(M_PI * t) + 0.08 * cos(2.0 * M_PI * t);
}

static double blackmanharris_edge(double t)
{
    return 0.35875 - 0.48829 * cos(M_PI * t) +
           0.14128 * cos(2.0 * M_PI * t) -
           0.01168 * cos(3.0 * M_PI * t);
}

static double welch_edge(double t)
{
    return 1.0 - (1.0 - t) * (1.0 - t);
}

static double parzen_edge(double t)
{
    double u;

    if (t <= 0.5) {
        return 2.0 * t * t * t;
    }

    u = 1.0 - t;
    return 1.0 - 6.0 * u * u + 6.0 * u * u * u;
}

static double gaussian_edge(double t)
{
    const double alpha = 3.0;
    double edge = exp(-0.5 * alpha * alpha);
    double value = exp(-0.5 * alpha * alpha * (1.0 - t) * (1.0 - t));

    return (value - edge) / (1.0 - edge);
}

static double smoothstep_edge(double t)
{
    return t * t * (3.0 - 2.0 * t);
}

static double smootherstep_edge(double t)
{
    return t * t * t * (t * (6.0 * t - 15.0) + 10.0);
}

static double exponential_edge(double t)
{
    const double alpha = 5.0;

    return (1.0 - exp(-alpha * t)) / (1.0 - exp(-alpha));
}

static double sine_edge(double t)
{
    return sin(0.5 * M_PI * t);
}

static double bohman_edge(double t)
{
    double u = 1.0 - t;

    return (1.0 - u) * cos(M_PI * u) + sin(M_PI * u) / M_PI;
}

static double nuttall_edge(double t)
{
    return 0.355768 - 0.487396 * cos(M_PI * t) +
           0.144232 * cos(2.0 * M_PI * t) -
           0.012604 * cos(3.0 * M_PI * t);
}

static double bessel_i0(double x)
{
    double ax = fabs(x);
    double y;

    if (ax < 3.75) {
        y = x / 3.75;
        y *= y;
        return 1.0 + y * (3.5156229 + y * (3.0899424 + y * (1.2067492 +
               y * (0.2659732 + y * (0.0360768 + y * 0.0045813)))));
    }

    y = 3.75 / ax;
    return (exp(ax) / sqrt(ax)) *
           (0.39894228 + y * (0.01328592 + y * (0.00225319 +
           y * (-0.00157565 + y * (0.00916281 + y * (-0.02057706 +
           y * (0.02635537 + y * (-0.01647633 + y * 0.00392377))))))));
}

static double kaiser_edge(double t)
{
    const double beta = 6.0;
    double edge = 1.0 / bessel_i0(beta);
    double u = 1.0 - t;
    double value = bessel_i0(beta * sqrt(1.0 - u * u)) / bessel_i0(beta);

    return (value - edge) / (1.0 - edge);
}

static double cauchy_edge(double t)
{
    const double alpha = 5.0;
    double edge = 1.0 / (1.0 + alpha * alpha);
    double u = 1.0 - t;
    double value = 1.0 / (1.0 + alpha * alpha * u * u);

    return (value - edge) / (1.0 - edge);
}

static double quadratic_edge(double t)
{
    return t * t;
}

static double cubic_edge(double t)
{
    return t * t * t;
}

static double poisson_edge(double t)
{
    const double alpha = 5.0;
    double edge = exp(-alpha);
    double value = exp(-alpha * (1.0 - t));

    return (value - edge) / (1.0 - edge);
}

static double bartlett_edge(double t)
{
    return t;
}

static double bartletthann_edge(double t)
{
    return 0.38 + 0.24 * t - 0.38 * cos(M_PI * t);
}

static double exactblackman_edge(double t)
{
    const double a0 = 7938.0 / 18608.0;
    const double a1 = 9240.0 / 18608.0;
    const double a2 = 1430.0 / 18608.0;

    return a0 - a1 * cos(M_PI * t) + a2 * cos(2.0 * M_PI * t);
}

static double blackmannuttall_edge(double t)
{
    return 0.3635819 - 0.4891775 * cos(M_PI * t) +
           0.1365995 * cos(2.0 * M_PI * t) -
           0.0106411 * cos(3.0 * M_PI * t);
}

static double flattop_edge(double t)
{
    return 0.21557895 - 0.41663158 * cos(M_PI * t) +
           0.277263158 * cos(2.0 * M_PI * t) -
           0.083578947 * cos(3.0 * M_PI * t) +
           0.006947368 * cos(4.0 * M_PI * t);
}

static double lanczos_edge(double t)
{
    double u = 1.0 - t;

    if (fabs(u) < DBL_EPSILON) return 1.0;
    return sin(M_PI * u) / (M_PI * u);
}

static double riesz_edge(double t)
{
    double u = 1.0 - t;

    return 1.0 - u * u;
}

static double riemann_edge(double t)
{
    return lanczos_edge(t);
}

static double fejer_edge(double t)
{
    return t;
}

static double connes_edge(double t)
{
    double r = riesz_edge(t);

    return r * r;
}

static double hanningpoisson_edge(double t)
{
    const double alpha = 2.0;

    return cosine_edge(t) * exp(-alpha * (1.0 - t));
}

static double kaiserbessel_edge(double t)
{
    return sqrt(blend_clamp_unit(kaiser_edge(t)));
}

static double plancktaper_edge(double t)
{
    if (t <= 0.0) return 0.0;
    if (t >= 1.0) return 1.0;

    return 1.0 / (exp((1.0 / t) - (1.0 / (1.0 - t))) + 1.0);
}

static double quartic_edge(double t)
{
    return t * t * t * t;
}

static double quintic_edge(double t)
{
    return t * t * t * t * t;
}

static double septic_edge(double t)
{
    double t2 = t * t;
    double t4 = t2 * t2;

    return t4 * t2 * t;
}

static double nonic_edge(double t)
{
    double t2 = t * t;
    double t4 = t2 * t2;

    return t4 * t4 * t;
}

static double logistic_edge(double t)
{
    const double alpha = 12.0;
    double edge = 1.0 / (1.0 + exp(0.5 * alpha));
    double center = 1.0 / (1.0 + exp(-0.5 * alpha));
    double value = 1.0 / (1.0 + exp(-alpha * (t - 0.5)));

    return (value - edge) / (center - edge);
}

static double tanh_edge(double t)
{
    const double alpha = 3.0;

    return (tanh(alpha * (t - 0.5)) + tanh(0.5 * alpha)) /
           (2.0 * tanh(0.5 * alpha));
}

static double erf_edge(double t)
{
    const double alpha = 3.0;

    return (erf(alpha * (t - 0.5)) + erf(0.5 * alpha)) /
           (2.0 * erf(0.5 * alpha));
}

static double arctan_edge(double t)
{
    const double alpha = 8.0;

    return (atan(alpha * (t - 0.5)) + atan(0.5 * alpha)) /
           (2.0 * atan(0.5 * alpha));
}

static double gompertz_edge(double t)
{
    const double b = 5.0;
    const double c = 6.0;
    double edge = exp(-b);
    double center = exp(-b * exp(-c));
    double value = exp(-b * exp(-c * t));

    return (value - edge) / (center - edge);
}

static double softsign_edge(double t)
{
    const double alpha = 2.0;

    return t / (t + alpha * (1.0 - t));
}

static double agnesi_edge(double t)
{
    const double alpha = 4.0;
    double edge = 1.0 / (1.0 + alpha * alpha);
    double u = 1.0 - t;
    double value = 1.0 / (1.0 + alpha * alpha * u * u);

    return (value - edge) / (1.0 - edge);
}

static double inversequadratic_edge(double t)
{
    return agnesi_edge(t);
}

static double inversemultiquadric_edge(double t)
{
    const double alpha = 4.0;
    double edge = 1.0 / sqrt(1.0 + alpha * alpha);
    double u = 1.0 - t;
    double value = 1.0 / sqrt(1.0 + alpha * alpha * u * u);

    return (value - edge) / (1.0 - edge);
}

static double powerlaw_edge(double t)
{
    return pow(t, 1.5);
}

static double root_edge(double t)
{
    return sqrt(t);
}

static double circular_edge(double t)
{
    double u = 1.0 - t;

    return sqrt(1.0 - u * u);
}

static double sech_value(double x)
{
    return 1.0 / cosh(x);
}

static double sech_edge(double t)
{
    const double alpha = 4.0;
    double edge = sech_value(alpha);
    double u = 1.0 - t;
    double value = sech_value(alpha * u);

    return (value - edge) / (1.0 - edge);
}

static double sech2_edge(double t)
{
    const double alpha = 4.0;
    double edge = sech_value(alpha) * sech_value(alpha);
    double u = 1.0 - t;
    double value = sech_value(alpha * u) * sech_value(alpha * u);

    return (value - edge) / (1.0 - edge);
}

static double student_edge(double t)
{
    const double nu = 3.0;
    double u = 1.0 - t;
    double edge = pow(1.0 + 1.0 / nu, -0.5 * (nu + 1.0));
    double value = pow(1.0 + u * u / nu, -0.5 * (nu + 1.0));

    return (value - edge) / (1.0 - edge);
}

static double laplace_edge(double t)
{
    return poisson_edge(t);
}

/* cosine function */
double cosine(int x, int n, int nmax, int nmin, double r1, double r2)
{
    return blend_window_from_edge(x, n, nmax, nmin, r1, r2, cosine_edge);
}

/* trapezoid function */
double trapezoid(int x, int n, int nmax, int nmin, double r1, double r2)
{
    return blend_window_from_edge(x, n, nmax, nmin, r1, r2, trapezoid_edge);
}

/* hamming function */
double hamming(int x, int n, int nmax, int nmin, double r1, double r2)
{
    return blend_window_from_edge(x, n, nmax, nmin, r1, r2, hamming_edge);
}

/* blackman function */
double blackman(int x, int n, int nmax, int nmin, double r1, double r2)
{
    return blend_window_from_edge(x, n, nmax, nmin, r1, r2, blackman_edge);
}

/* blackman-harris function */
double blackmanharris(int x, int n, int nmax, int nmin, double r1, double r2)
{
    return blend_window_from_edge(x, n, nmax, nmin, r1, r2, blackmanharris_edge);
}

/* welch function */
double welch(int x, int n, int nmax, int nmin, double r1, double r2)
{
    return blend_window_from_edge(x, n, nmax, nmin, r1, r2, welch_edge);
}

/* parzen function */
double parzen(int x, int n, int nmax, int nmin, double r1, double r2)
{
    return blend_window_from_edge(x, n, nmax, nmin, r1, r2, parzen_edge);
}

/* gaussian function */
double gaussian(int x, int n, int nmax, int nmin, double r1, double r2)
{
    return blend_window_from_edge(x, n, nmax, nmin, r1, r2, gaussian_edge);
}

/* smoothstep function */
double smoothstep(int x, int n, int nmax, int nmin, double r1, double r2)
{
    return blend_window_from_edge(x, n, nmax, nmin, r1, r2, smoothstep_edge);
}

/* smootherstep function */
double smootherstep(int x, int n, int nmax, int nmin, double r1, double r2)
{
    return blend_window_from_edge(x, n, nmax, nmin, r1, r2, smootherstep_edge);
}

/* exponential function */
double exponential(int x, int n, int nmax, int nmin, double r1, double r2)
{
    return blend_window_from_edge(x, n, nmax, nmin, r1, r2, exponential_edge);
}

/* sine function */
double sine(int x, int n, int nmax, int nmin, double r1, double r2)
{
    return blend_window_from_edge(x, n, nmax, nmin, r1, r2, sine_edge);
}

/* bohman function */
double bohman(int x, int n, int nmax, int nmin, double r1, double r2)
{
    return blend_window_from_edge(x, n, nmax, nmin, r1, r2, bohman_edge);
}

/* nuttall function */
double nuttall(int x, int n, int nmax, int nmin, double r1, double r2)
{
    return blend_window_from_edge(x, n, nmax, nmin, r1, r2, nuttall_edge);
}

/* kaiser function */
double kaiser(int x, int n, int nmax, int nmin, double r1, double r2)
{
    return blend_window_from_edge(x, n, nmax, nmin, r1, r2, kaiser_edge);
}

/* cauchy function */
double cauchy(int x, int n, int nmax, int nmin, double r1, double r2)
{
    return blend_window_from_edge(x, n, nmax, nmin, r1, r2, cauchy_edge);
}

/* quadratic function */
double quadratic(int x, int n, int nmax, int nmin, double r1, double r2)
{
    return blend_window_from_edge(x, n, nmax, nmin, r1, r2, quadratic_edge);
}

/* cubic function */
double cubic(int x, int n, int nmax, int nmin, double r1, double r2)
{
    return blend_window_from_edge(x, n, nmax, nmin, r1, r2, cubic_edge);
}

/* poisson function */
double poisson(int x, int n, int nmax, int nmin, double r1, double r2)
{
    return blend_window_from_edge(x, n, nmax, nmin, r1, r2, poisson_edge);
}

/* bartlett function */
double bartlett(int x, int n, int nmax, int nmin, double r1, double r2)
{
    return blend_window_from_edge(x, n, nmax, nmin, r1, r2, bartlett_edge);
}

/* bartlett-hann function */
double bartletthann(int x, int n, int nmax, int nmin, double r1, double r2)
{
    return blend_window_from_edge(x, n, nmax, nmin, r1, r2, bartletthann_edge);
}

/* exact blackman function */
double exactblackman(int x, int n, int nmax, int nmin, double r1, double r2)
{
    return blend_window_from_edge(x, n, nmax, nmin, r1, r2, exactblackman_edge);
}

/* blackman-nuttall function */
double blackmannuttall(int x, int n, int nmax, int nmin, double r1, double r2)
{
    return blend_window_from_edge(x, n, nmax, nmin, r1, r2, blackmannuttall_edge);
}

/* flat-top function */
double flattop(int x, int n, int nmax, int nmin, double r1, double r2)
{
    return blend_window_from_edge(x, n, nmax, nmin, r1, r2, flattop_edge);
}

/* lanczos function */
double lanczos(int x, int n, int nmax, int nmin, double r1, double r2)
{
    return blend_window_from_edge(x, n, nmax, nmin, r1, r2, lanczos_edge);
}

/* riesz function */
double riesz(int x, int n, int nmax, int nmin, double r1, double r2)
{
    return blend_window_from_edge(x, n, nmax, nmin, r1, r2, riesz_edge);
}

/* riemann function */
double riemann(int x, int n, int nmax, int nmin, double r1, double r2)
{
    return blend_window_from_edge(x, n, nmax, nmin, r1, r2, riemann_edge);
}

/* fejer function */
double fejer(int x, int n, int nmax, int nmin, double r1, double r2)
{
    return blend_window_from_edge(x, n, nmax, nmin, r1, r2, fejer_edge);
}

/* connes function */
double connes(int x, int n, int nmax, int nmin, double r1, double r2)
{
    return blend_window_from_edge(x, n, nmax, nmin, r1, r2, connes_edge);
}

/* hanning-poisson function */
double hanningpoisson(int x, int n, int nmax, int nmin, double r1, double r2)
{
    return blend_window_from_edge(x, n, nmax, nmin, r1, r2, hanningpoisson_edge);
}

/* kaiser-bessel function */
double kaiserbessel(int x, int n, int nmax, int nmin, double r1, double r2)
{
    return blend_window_from_edge(x, n, nmax, nmin, r1, r2, kaiserbessel_edge);
}

/* planck-taper function */
double plancktaper(int x, int n, int nmax, int nmin, double r1, double r2)
{
    return blend_window_from_edge(x, n, nmax, nmin, r1, r2, plancktaper_edge);
}

/* quartic function */
double quartic(int x, int n, int nmax, int nmin, double r1, double r2)
{
    return blend_window_from_edge(x, n, nmax, nmin, r1, r2, quartic_edge);
}

/* quintic function */
double quintic(int x, int n, int nmax, int nmin, double r1, double r2)
{
    return blend_window_from_edge(x, n, nmax, nmin, r1, r2, quintic_edge);
}

/* septic function */
double septic(int x, int n, int nmax, int nmin, double r1, double r2)
{
    return blend_window_from_edge(x, n, nmax, nmin, r1, r2, septic_edge);
}

/* nonic function */
double nonic(int x, int n, int nmax, int nmin, double r1, double r2)
{
    return blend_window_from_edge(x, n, nmax, nmin, r1, r2, nonic_edge);
}

/* logistic function */
double logistic(int x, int n, int nmax, int nmin, double r1, double r2)
{
    return blend_window_from_edge(x, n, nmax, nmin, r1, r2, logistic_edge);
}

/* tanh function */
double tanhwindow(int x, int n, int nmax, int nmin, double r1, double r2)
{
    return blend_window_from_edge(x, n, nmax, nmin, r1, r2, tanh_edge);
}

/* erf function */
double erfwindow(int x, int n, int nmax, int nmin, double r1, double r2)
{
    return blend_window_from_edge(x, n, nmax, nmin, r1, r2, erf_edge);
}

/* arctan function */
double arctanwindow(int x, int n, int nmax, int nmin, double r1, double r2)
{
    return blend_window_from_edge(x, n, nmax, nmin, r1, r2, arctan_edge);
}

/* gompertz function */
double gompertz(int x, int n, int nmax, int nmin, double r1, double r2)
{
    return blend_window_from_edge(x, n, nmax, nmin, r1, r2, gompertz_edge);
}

/* softsign function */
double softsign(int x, int n, int nmax, int nmin, double r1, double r2)
{
    return blend_window_from_edge(x, n, nmax, nmin, r1, r2, softsign_edge);
}

/* agnesi function */
double agnesi(int x, int n, int nmax, int nmin, double r1, double r2)
{
    return blend_window_from_edge(x, n, nmax, nmin, r1, r2, agnesi_edge);
}

/* inverse quadratic function */
double inversequadratic(int x, int n, int nmax, int nmin, double r1, double r2)
{
    return blend_window_from_edge(x, n, nmax, nmin, r1, r2, inversequadratic_edge);
}

/* inverse multiquadric function */
double inversemultiquadric(int x, int n, int nmax, int nmin, double r1, double r2)
{
    return blend_window_from_edge(x, n, nmax, nmin, r1, r2, inversemultiquadric_edge);
}

/* power-law function */
double powerlaw(int x, int n, int nmax, int nmin, double r1, double r2)
{
    return blend_window_from_edge(x, n, nmax, nmin, r1, r2, powerlaw_edge);
}

/* root function */
double root(int x, int n, int nmax, int nmin, double r1, double r2)
{
    return blend_window_from_edge(x, n, nmax, nmin, r1, r2, root_edge);
}

/* circular function */
double circular(int x, int n, int nmax, int nmin, double r1, double r2)
{
    return blend_window_from_edge(x, n, nmax, nmin, r1, r2, circular_edge);
}

/* hyperbolic secant function */
double sechwindow(int x, int n, int nmax, int nmin, double r1, double r2)
{
    return blend_window_from_edge(x, n, nmax, nmin, r1, r2, sech_edge);
}

/* squared hyperbolic secant function */
double sech2window(int x, int n, int nmax, int nmin, double r1, double r2)
{
    return blend_window_from_edge(x, n, nmax, nmin, r1, r2, sech2_edge);
}

/* student function */
double student(int x, int n, int nmax, int nmin, double r1, double r2)
{
    return blend_window_from_edge(x, n, nmax, nmin, r1, r2, student_edge);
}

/* laplace function */
double laplace(int x, int n, int nmax, int nmin, double r1, double r2)
{
    return blend_window_from_edge(x, n, nmax, nmin, r1, r2, laplace_edge);
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
        case WFUNC_HAMMING:
            taper = hamming(x, n2, n, nmin, r1, r2);
            break;
        case WFUNC_BLACKMAN:
            taper = blackman(x, n2, n, nmin, r1, r2);
            break;
        case WFUNC_BLACKMANHARRIS:
            taper = blackmanharris(x, n2, n, nmin, r1, r2);
            break;
        case WFUNC_WELCH:
            taper = welch(x, n2, n, nmin, r1, r2);
            break;
        case WFUNC_PARZEN:
            taper = parzen(x, n2, n, nmin, r1, r2);
            break;
        case WFUNC_GAUSSIAN:
            taper = gaussian(x, n2, n, nmin, r1, r2);
            break;
        case WFUNC_SMOOTHSTEP:
            taper = smoothstep(x, n2, n, nmin, r1, r2);
            break;
        case WFUNC_SMOOTHERSTEP:
            taper = smootherstep(x, n2, n, nmin, r1, r2);
            break;
        case WFUNC_EXPONENTIAL:
            taper = exponential(x, n2, n, nmin, r1, r2);
            break;
        case WFUNC_SINE:
            taper = sine(x, n2, n, nmin, r1, r2);
            break;
        case WFUNC_BOHMAN:
            taper = bohman(x, n2, n, nmin, r1, r2);
            break;
        case WFUNC_NUTTALL:
            taper = nuttall(x, n2, n, nmin, r1, r2);
            break;
        case WFUNC_KAISER:
            taper = kaiser(x, n2, n, nmin, r1, r2);
            break;
        case WFUNC_CAUCHY:
            taper = cauchy(x, n2, n, nmin, r1, r2);
            break;
        case WFUNC_QUADRATIC:
            taper = quadratic(x, n2, n, nmin, r1, r2);
            break;
        case WFUNC_CUBIC:
            taper = cubic(x, n2, n, nmin, r1, r2);
            break;
        case WFUNC_POISSON:
            taper = poisson(x, n2, n, nmin, r1, r2);
            break;
        case WFUNC_BARTLETT:
            taper = bartlett(x, n2, n, nmin, r1, r2);
            break;
        case WFUNC_BARTLETTHANN:
            taper = bartletthann(x, n2, n, nmin, r1, r2);
            break;
        case WFUNC_EXACTBLACKMAN:
            taper = exactblackman(x, n2, n, nmin, r1, r2);
            break;
        case WFUNC_BLACKMANNUTTALL:
            taper = blackmannuttall(x, n2, n, nmin, r1, r2);
            break;
        case WFUNC_FLATTOP:
            taper = flattop(x, n2, n, nmin, r1, r2);
            break;
        case WFUNC_LANCZOS:
            taper = lanczos(x, n2, n, nmin, r1, r2);
            break;
        case WFUNC_RIESZ:
            taper = riesz(x, n2, n, nmin, r1, r2);
            break;
        case WFUNC_RIEMANN:
            taper = riemann(x, n2, n, nmin, r1, r2);
            break;
        case WFUNC_FEJER:
            taper = fejer(x, n2, n, nmin, r1, r2);
            break;
        case WFUNC_CONNES:
            taper = connes(x, n2, n, nmin, r1, r2);
            break;
        case WFUNC_HANNINGPOISSON:
            taper = hanningpoisson(x, n2, n, nmin, r1, r2);
            break;
        case WFUNC_KAISERBESSEL:
            taper = kaiserbessel(x, n2, n, nmin, r1, r2);
            break;
        case WFUNC_PLANCKTAPER:
            taper = plancktaper(x, n2, n, nmin, r1, r2);
            break;
        case WFUNC_QUARTIC:
            taper = quartic(x, n2, n, nmin, r1, r2);
            break;
        case WFUNC_QUINTIC:
            taper = quintic(x, n2, n, nmin, r1, r2);
            break;
        case WFUNC_SEPTIC:
            taper = septic(x, n2, n, nmin, r1, r2);
            break;
        case WFUNC_NONIC:
            taper = nonic(x, n2, n, nmin, r1, r2);
            break;
        case WFUNC_LOGISTIC:
            taper = logistic(x, n2, n, nmin, r1, r2);
            break;
        case WFUNC_TANH:
            taper = tanhwindow(x, n2, n, nmin, r1, r2);
            break;
        case WFUNC_ERF:
            taper = erfwindow(x, n2, n, nmin, r1, r2);
            break;
        case WFUNC_ARCTAN:
            taper = arctanwindow(x, n2, n, nmin, r1, r2);
            break;
        case WFUNC_GOMPERTZ:
            taper = gompertz(x, n2, n, nmin, r1, r2);
            break;
        case WFUNC_SOFTSIGN:
            taper = softsign(x, n2, n, nmin, r1, r2);
            break;
        case WFUNC_AGNESI:
            taper = agnesi(x, n2, n, nmin, r1, r2);
            break;
        case WFUNC_INVERSEQUADRATIC:
            taper = inversequadratic(x, n2, n, nmin, r1, r2);
            break;
        case WFUNC_INVERSEMULTIQUADRIC:
            taper = inversemultiquadric(x, n2, n, nmin, r1, r2);
            break;
        case WFUNC_POWERLAW:
            taper = powerlaw(x, n2, n, nmin, r1, r2);
            break;
        case WFUNC_ROOT:
            taper = root(x, n2, n, nmin, r1, r2);
            break;
        case WFUNC_CIRCULAR:
            taper = circular(x, n2, n, nmin, r1, r2);
            break;
        case WFUNC_SECH:
            taper = sechwindow(x, n2, n, nmin, r1, r2);
            break;
        case WFUNC_SECH2:
            taper = sech2window(x, n2, n, nmin, r1, r2);
            break;
        case WFUNC_STUDENT:
            taper = student(x, n2, n, nmin, r1, r2);
            break;
        case WFUNC_LAPLACE:
            taper = laplace(x, n2, n, nmin, r1, r2);
            break;
        case WFUNC_INVALID:
        default:
            BLEND_Report(BLEND_MSG_ERROR, "window: requested window function does not exist\n");
            return WFUNC_ERROR;
    }

    return blend_clamp_unit(taper);
}
