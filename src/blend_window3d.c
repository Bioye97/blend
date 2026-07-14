/*
* Author: Rasheed Ajala
* Contact: abioyeajala@gmail.com
*
*/

#include "blend.h"

#include <ctype.h>
#include <errno.h>
#include <limits.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#ifndef PATH_MAX
#define PATH_MAX 4096
#endif

typedef struct window3d_options {
    double xmin;
    double xmax;
    double ymin;
    double ymax;
    double zmin;
    double zmax;
    double dx;
    double dy;
    double dz;
    double ratio_x1;
    double ratio_x2;
    double ratio_y1;
    double ratio_y2;
    double ratio_z1;
    double ratio_z2;
    int nx;
    int ny;
    int nz;
    int has_region;
    int has_increment;
    int has_function;
    int has_taper;
    int has_monotone_method;
    int has_write_monotone;
    blend_window_function x_function;
    blend_window_function y_function;
    blend_window_function z_function;
    const char *blendfile;
    char clobber;
    char monotone_method;
} window3d_options;

typedef struct window3d_support {
    double xmin;
    double xmax;
    double ymin;
    double ymax;
    double zlo;
    double zhi;
    size_t n_vertices;
    vertex *vertices;
    window data;
} window3d_support;

typedef struct window3d_support_list {
    window3d_support *items;
    size_t count;
    size_t capacity;
} window3d_support_list;

static void window3d_usage(FILE *fp)
{
    fprintf(fp, "blend window3d - Generate 3-D blending weights\n\n");
    fprintf(fp, "usage: blend window3d -R<xmin>/<xmax>/<ymin>/<ymax>/<zmin>/<zmax> -I<dx>[/<dy>[/<dz>]]\n");
    fprintf(fp, "       [-F<xfunction>[/<yfunction>[/<zfunction>]]] [-T<rx1>/<rx2>/<ry1>/<ry2>/<rz1>/<rz2>]\n");
    fprintf(fp, "       [-B<blendfile>] [-C<f|l|o|u|a|g|p>] [-M<method>] [-N] [-V[q|e|w|t|i|c|d]]\n\n");

    fprintf(fp, "Generate 3-D blending weights and write four columns: x y z weight. If x/y/z\n");
    fprintf(fp, "query coordinates are provided on standard input, weights are returned at those\n");
    fprintf(fp, "locations using trilinear interpolation from neighboring grid-point weights. If\n");
    fprintf(fp, "no query points are provided, weights are written for the complete -R/-I grid.\n\n");

    fprintf(fp, "REQUIRED ARGUMENTS:\n\n");
    fprintf(fp, "  -R<xmin>/<xmax>/<ymin>/<ymax>/<zmin>/<zmax>, --region=<...>\n");
    fprintf(fp, "     Specify the full 3-D output domain in user coordinates. Each dimension\n");
    fprintf(fp, "     must satisfy min < max. Upper bounds are adjusted upward when needed so\n");
    fprintf(fp, "     the domain lands exactly on the requested increment.\n\n");

    fprintf(fp, "  -I<dx>[/<dy>[/<dz>]], --increment=<dx>[/<dy>[/<dz>]]\n");
    fprintf(fp, "     Specify grid increments. If dy is omitted then dy = dx. If dz is omitted\n");
    fprintf(fp, "     then dz = dx. All increments must be positive.\n\n");

    fprintf(fp, "OPTIONAL ARGUMENTS:\n\n");
    fprintf(fp, "  -F<xfunction>[/<yfunction>[/<zfunction>]], --function=<...>\n");
    fprintf(fp, "     Select window functions used along x, y, and z for one support over the\n");
    fprintf(fp, "     full -R domain [Default is cosine/cosine/cosine]. Omitted functions copy\n");
    fprintf(fp, "     the x function. Run blend --show-windows to list available functions.\n");
    fprintf(fp, "     This option is ignored when -B is given.\n\n");

    fprintf(fp, "  -T<rx1>/<rx2>/<ry1>/<ry2>/<rz1>/<rz2>, --taper_ratio=<...>\n");
    fprintf(fp, "     Set beginning and ending taper ratios for x, y, and z. Ratios must be\n");
    fprintf(fp, "     >= 0 and < 0.5 [Default is 0.2/0.2/0.2/0.2/0.2/0.2]. This option is\n");
    fprintf(fp, "     ignored for boxcar dimensions and when -B is given.\n\n");

    fprintf(fp, "  -B<blendfile>, --blendfile=<blendfile>\n");
    fprintf(fp, "     Read one or more 3-D supports from <blendfile>. Each non-empty,\n");
    fprintf(fp, "     non-comment row must contain exactly five fields:\n");
    fprintf(fp, "       <polygonfile> <zlo> <zhi> <xfunction>/<yfunction>/<zfunction> <rx1>/<rx2>/<ry1>/<ry2>/<rz1>/<rz2>\n");
    fprintf(fp, "     where <polygonfile> is a two-column x/y vertex file contained inside -R,\n");
    fprintf(fp, "     and zlo/zhi define the vertical support. Grid points outside all\n");
    fprintf(fp, "     blendfile supports receive weight 0.\n\n");

    fprintf(fp, "  -C<f|l|o|u|a|g|p>, --clobber=<mode>\n");
    fprintf(fp, "     Select how overlapping support weights from -B are combined [Default is p].\n");
    fprintf(fp, "       f  Use the first support listed in the blendfile.\n");
    fprintf(fp, "       l  Use the lowest weight among overlapping supports.\n");
    fprintf(fp, "       o  Use the last support listed in the blendfile.\n");
    fprintf(fp, "       u  Use the highest weight among overlapping supports.\n");
    fprintf(fp, "       a  Use the arithmetic average of overlapping weights.\n");
    fprintf(fp, "       g  Use the geometric average of overlapping weights.\n");
    fprintf(fp, "       p  Use the product of overlapping weights.\n\n");

    fprintf(fp, "  -M<method>, --monotone=<method>\n");
    fprintf(fp, "     Here only the 'strictly increasing' monotone condition is allowed.\n");
    fprintf(fp, "     Select how non-strict xy-monotone blendfile polygons are converted before\n");
    fprintf(fp, "     window assembly. If this option is omitted, polygons that are not strict\n");
    fprintf(fp, "     enough for function localization are rejected. -MB is experimental and can\n");
    fprintf(fp, "     produce irregular polygons that make function localization awkward.\n");
    fprintf(fp, "       E  Use the strict xy-monotone envelope (i.e., convex hull) method.\n");
    fprintf(fp, "       B  Use the strict highest-IoU piecewise envelope from both traversal directions,\n");
    fprintf(fp, "          sampled on the -R/-I grid.\n\n");

    fprintf(fp, "  -N, --write-monotone\n");
    fprintf(fp, "     When -M modifies a polygon, write the modified polygon to a file named\n");
    fprintf(fp, "     by inserting _monotone before the file extension, for example\n");
    fprintf(fp, "     south_america_monotone.txt. If -N is omitted, modified polygons are used for\n");
    fprintf(fp, "     weights but are not written.\n\n");

    fprintf(fp, "  -V[level], --verbose=<level>\n");
    fprintf(fp, "     Select verbosity level [w]. Choose among q, e, w, t, i, c, and d:\n");
    fprintf(fp, "       q  Quiet; suppress all diagnostic messages.\n");
    fprintf(fp, "       e  Error messages only.\n");
    fprintf(fp, "       w  Warnings and errors [Default].\n");
    fprintf(fp, "       t  Timings, warnings, and errors.\n");
    fprintf(fp, "       i  Informational messages, timings, warnings, and errors.\n");
    fprintf(fp, "       c  Compatibility messages and all lower verbosity messages.\n");
    fprintf(fp, "       d  Debug messages and all lower verbosity messages.\n\n");

    fprintf(fp, "  -?, --help\n");
    fprintf(fp, "     Print this usage message and exit.\n\n");

    fprintf(fp, "EXAMPLES:\n\n");
    fprintf(fp, "  blend window3d -R0/10/0/5/0/2 -I1/1/0.5 -Fcosine/cosine/cosine\n");
    fprintf(fp, "     Generate one 3-D cosine-tapered window over the full domain.\n\n");

    fprintf(fp, "  printf '2.5 1.5 1\\n' | blend window3d -R0/10/0/5/0/2 -I1 -Fboxcar/boxcar/boxcar\n");
    fprintf(fp, "     Query the interpolated weight at x = 2.5, y = 1.5, z = 1.\n\n");

    fprintf(fp, "  blend window3d -R0/10/0/10/0/5 -I0.25/0.25/0.5 -Bsupports.txt -Ca\n");
    fprintf(fp, "     Read multiple polygon slabs from supports.txt and average overlaps.\n");
}

static int window3d_parse_region(const char *text, window3d_options *options)
{
    char *end = NULL;
    double values[6];
    int i;

    for (i = 0; i < 6; i++) {
        errno = 0;
        values[i] = strtod(text, &end);
        if (text == end || errno == ERANGE || !isfinite(values[i])) {
            BLEND_Report(BLEND_MSG_ERROR, "window3d: invalid region\n");
            return FAIL;
        }
        if (i < 5) {
            if (*end != '/') {
                BLEND_Report(BLEND_MSG_ERROR, "window3d: region must be xmin/xmax/ymin/ymax/zmin/zmax\n");
                return FAIL;
            }
            text = end + 1;
        }
        else if (*end != '\0') {
            BLEND_Report(BLEND_MSG_ERROR, "window3d: region must be xmin/xmax/ymin/ymax/zmin/zmax\n");
            return FAIL;
        }
    }

    if (values[0] >= values[1] || values[2] >= values[3] || values[4] >= values[5]) {
        BLEND_Report(BLEND_MSG_ERROR, "window3d: region must satisfy xmin < xmax, ymin < ymax, and zmin < zmax\n");
        return FAIL;
    }

    options->xmin = values[0];
    options->xmax = values[1];
    options->ymin = values[2];
    options->ymax = values[3];
    options->zmin = values[4];
    options->zmax = values[5];
    options->has_region = 1;
    return SUCCESS;
}

static int window3d_parse_increment(const char *text, window3d_options *options)
{
    char *end = NULL;
    double values[3];
    int count = 1;

    errno = 0;
    values[0] = strtod(text, &end);
    if (text == end || errno == ERANGE || !isfinite(values[0])) {
        BLEND_Report(BLEND_MSG_ERROR, "window3d: invalid increment: %s\n", text);
        return FAIL;
    }

    while (*end == '/' && count < 3) {
        text = end + 1;
        errno = 0;
        values[count] = strtod(text, &end);
        if (text == end || errno == ERANGE || !isfinite(values[count])) {
            BLEND_Report(BLEND_MSG_ERROR, "window3d: invalid increment: %s\n", text);
            return FAIL;
        }
        count++;
    }
    if (*end != '\0') {
        BLEND_Report(BLEND_MSG_ERROR, "window3d: increment must be dx[/dy[/dz]]\n");
        return FAIL;
    }

    options->dx = values[0];
    options->dy = count > 1 ? values[1] : values[0];
    options->dz = count > 2 ? values[2] : values[0];
    options->has_increment = 1;
    return SUCCESS;
}

static int window3d_parse_taper_ratio_values(const char *text, const char *context,
                                             double *rx1, double *rx2,
                                             double *ry1, double *ry2,
                                             double *rz1, double *rz2)
{
    char *end = NULL;
    double values[6];
    int i;

    for (i = 0; i < 6; i++) {
        errno = 0;
        values[i] = strtod(text, &end);
        if (text == end || errno == ERANGE || !isfinite(values[i])) {
            BLEND_Report(BLEND_MSG_ERROR, "window3d: invalid taper ratio in %s\n", context);
            return FAIL;
        }
        if (i < 5) {
            if (*end != '/') {
                BLEND_Report(BLEND_MSG_ERROR, "window3d: taper ratio in %s must be rx1/rx2/ry1/ry2/rz1/rz2\n", context);
                return FAIL;
            }
            text = end + 1;
        }
        else if (*end != '\0') {
            BLEND_Report(BLEND_MSG_ERROR, "window3d: taper ratio in %s must be rx1/rx2/ry1/ry2/rz1/rz2\n", context);
            return FAIL;
        }
    }

    for (i = 0; i < 6; i++) {
        if (values[i] < 0.0 || values[i] >= 0.5) {
            BLEND_Report(BLEND_MSG_ERROR, "window3d: taper ratios in %s must be >= 0 and < 0.5\n", context);
            return FAIL;
        }
    }

    *rx1 = values[0];
    *rx2 = values[1];
    *ry1 = values[2];
    *ry2 = values[3];
    *rz1 = values[4];
    *rz2 = values[5];
    return SUCCESS;
}

static int window3d_parse_taper_ratio(const char *text, window3d_options *options)
{
    return window3d_parse_taper_ratio_values(text, "option -T",
                                             &options->ratio_x1, &options->ratio_x2,
                                             &options->ratio_y1, &options->ratio_y2,
                                             &options->ratio_z1, &options->ratio_z2);
}

static int window3d_parse_function_triple(const char *text,
                                          blend_window_function *x_function,
                                          blend_window_function *y_function,
                                          blend_window_function *z_function)
{
    char x_name[64], y_name[64], z_name[64];
    const char *first_slash = strchr(text, '/');
    const char *second_slash = first_slash == NULL ? NULL : strchr(first_slash + 1, '/');
    size_t x_len, y_len, z_len;

    if (first_slash == NULL) {
        x_len = strlen(text);
        y_len = x_len;
        z_len = x_len;
        if (x_len == 0 || x_len >= sizeof(x_name)) {
            BLEND_Report(BLEND_MSG_ERROR, "window3d: invalid window function: %s\n", text);
            return FAIL;
        }
        memcpy(x_name, text, x_len + 1);
        memcpy(y_name, text, y_len + 1);
        memcpy(z_name, text, z_len + 1);
    }
    else if (second_slash == NULL) {
        x_len = (size_t)(first_slash - text);
        y_len = strlen(first_slash + 1);
        z_len = x_len;
        if (x_len == 0 || y_len == 0 || x_len >= sizeof(x_name) || y_len >= sizeof(y_name)) {
            BLEND_Report(BLEND_MSG_ERROR, "window3d: function must be xfunction[/yfunction[/zfunction]]\n");
            return FAIL;
        }
        memcpy(x_name, text, x_len);
        x_name[x_len] = '\0';
        memcpy(y_name, first_slash + 1, y_len + 1);
        memcpy(z_name, x_name, x_len + 1);
    }
    else {
        x_len = (size_t)(first_slash - text);
        y_len = (size_t)(second_slash - first_slash - 1);
        z_len = strlen(second_slash + 1);
        if (x_len == 0 || y_len == 0 || z_len == 0 ||
            x_len >= sizeof(x_name) || y_len >= sizeof(y_name) || z_len >= sizeof(z_name)) {
            BLEND_Report(BLEND_MSG_ERROR, "window3d: function must be xfunction/yfunction/zfunction\n");
            return FAIL;
        }
        memcpy(x_name, text, x_len);
        x_name[x_len] = '\0';
        memcpy(y_name, first_slash + 1, y_len);
        y_name[y_len] = '\0';
        memcpy(z_name, second_slash + 1, z_len + 1);
    }

    if (blend_window_function_from_name(x_name, x_function) != SUCCESS ||
        blend_window_function_from_name(y_name, y_function) != SUCCESS ||
        blend_window_function_from_name(z_name, z_function) != SUCCESS) {
        BLEND_Report(BLEND_MSG_ERROR, "window3d: unknown window function in %s\n", text);
        return FAIL;
    }

    return SUCCESS;
}

static int window3d_parse_clobber(const char *value, window3d_options *options)
{
    if (value == NULL || value[0] == '\0' || value[1] != '\0') {
        BLEND_Report(BLEND_MSG_ERROR, "window3d: clobber mode must be one of f, l, o, u, a, g, p\n");
        return FAIL;
    }

    if (strchr("flouagp", value[0]) == NULL) {
        BLEND_Report(BLEND_MSG_ERROR, "window3d: unknown clobber mode: %s\n", value);
        return FAIL;
    }

    options->clobber = value[0];
    return SUCCESS;
}

static int window3d_parse_monotone_method(const char *value, window3d_options *options)
{
    if (value == NULL || value[0] == '\0' || value[1] != '\0') {
        BLEND_Report(BLEND_MSG_ERROR, "window3d: monotone method must be one of E, B\n");
        return FAIL;
    }

    if (value[0] != 'E' && value[0] != 'B') {
        BLEND_Report(BLEND_MSG_ERROR, "window3d: unknown monotone method: %s\n", value);
        return FAIL;
    }

    options->monotone_method = value[0];
    options->has_monotone_method = 1;
    return SUCCESS;
}

static const char *window3d_option_value(int argc, char **argv, int *i, const char *option)
{
    if (*i + 1 >= argc) {
        BLEND_Report(BLEND_MSG_ERROR, "window3d: option %s requires an argument\n", option);
        return NULL;
    }

    *i += 1;
    return argv[*i];
}

static int window3d_adjust_axis(const char *axis, double start, double *end, double increment, int *n)
{
    double old_end = *end;
    double span = *end - start;
    double intervals_exact = span / increment;
    double intervals_nearest = floor(intervals_exact + 0.5);
    double tolerance = 1.0e-9 * fmax(1.0, fabs(intervals_exact));
    double intervals;

    if (fabs(intervals_exact - intervals_nearest) <= tolerance) {
        intervals = intervals_nearest;
        *end = start + intervals * increment;
    }
    else {
        intervals = ceil(intervals_exact - tolerance);
        *end = start + intervals * increment;
        BLEND_Report(
            BLEND_MSG_WARNING,
            "window3d: adjusted %s region end from %.12g to %.12g to match increment %.12g\n",
            axis, old_end, *end, increment
        );
    }

    if (intervals < 1.0 || intervals > (double)(INT_MAX - 1)) {
        BLEND_Report(BLEND_MSG_ERROR, "window3d: region and increment produce an invalid %s grid size\n", axis);
        return FAIL;
    }

    *n = (int)intervals + 1;
    return SUCCESS;
}

static int window3d_parse_options(int argc, char **argv, window3d_options *options)
{
    int i;

    if (argc == 0) {
        window3d_usage(stdout);
        return 2;
    }

    memset(options, 0, sizeof(*options));
    options->ratio_x1 = 0.2;
    options->ratio_x2 = 0.2;
    options->ratio_y1 = 0.2;
    options->ratio_y2 = 0.2;
    options->ratio_z1 = 0.2;
    options->ratio_z2 = 0.2;
    options->x_function = WFUNC_COSINE;
    options->y_function = WFUNC_COSINE;
    options->z_function = WFUNC_COSINE;
    options->clobber = 'p';

    for (i = 0; i < argc; i++) {
        const char *arg = argv[i];
        const char *value = NULL;

        if (strcmp(arg, "-?") == 0 || strcmp(arg, "--help") == 0) {
            window3d_usage(stdout);
            return 2;
        }
        else if (strcmp(arg, "-R") == 0 || strcmp(arg, "--region") == 0) {
            value = window3d_option_value(argc, argv, &i, arg);
            if (value == NULL || window3d_parse_region(value, options) != SUCCESS) return FAIL;
        }
        else if (strncmp(arg, "-R", 2) == 0 && arg[2] != '\0') {
            if (window3d_parse_region(arg + 2, options) != SUCCESS) return FAIL;
        }
        else if (strncmp(arg, "--region=", 9) == 0) {
            if (window3d_parse_region(arg + 9, options) != SUCCESS) return FAIL;
        }
        else if (strcmp(arg, "-I") == 0 || strcmp(arg, "--increment") == 0) {
            value = window3d_option_value(argc, argv, &i, arg);
            if (value == NULL || window3d_parse_increment(value, options) != SUCCESS) return FAIL;
        }
        else if (strncmp(arg, "-I", 2) == 0 && arg[2] != '\0') {
            if (window3d_parse_increment(arg + 2, options) != SUCCESS) return FAIL;
        }
        else if (strncmp(arg, "--increment=", 12) == 0) {
            if (window3d_parse_increment(arg + 12, options) != SUCCESS) return FAIL;
        }
        else if (strcmp(arg, "-F") == 0 || strcmp(arg, "--function") == 0) {
            value = window3d_option_value(argc, argv, &i, arg);
            if (value == NULL || window3d_parse_function_triple(value, &options->x_function,
                                                                &options->y_function,
                                                                &options->z_function) != SUCCESS) return FAIL;
            options->has_function = 1;
        }
        else if (strncmp(arg, "-F", 2) == 0 && arg[2] != '\0') {
            if (window3d_parse_function_triple(arg + 2, &options->x_function,
                                               &options->y_function,
                                               &options->z_function) != SUCCESS) return FAIL;
            options->has_function = 1;
        }
        else if (strncmp(arg, "--function=", 11) == 0) {
            if (window3d_parse_function_triple(arg + 11, &options->x_function,
                                               &options->y_function,
                                               &options->z_function) != SUCCESS) return FAIL;
            options->has_function = 1;
        }
        else if (strcmp(arg, "-T") == 0 || strcmp(arg, "--taper_ratio") == 0) {
            value = window3d_option_value(argc, argv, &i, arg);
            if (value == NULL || window3d_parse_taper_ratio(value, options) != SUCCESS) return FAIL;
            options->has_taper = 1;
        }
        else if (strncmp(arg, "-T", 2) == 0 && arg[2] != '\0') {
            if (window3d_parse_taper_ratio(arg + 2, options) != SUCCESS) return FAIL;
            options->has_taper = 1;
        }
        else if (strncmp(arg, "--taper_ratio=", 14) == 0) {
            if (window3d_parse_taper_ratio(arg + 14, options) != SUCCESS) return FAIL;
            options->has_taper = 1;
        }
        else if (strcmp(arg, "-B") == 0 || strcmp(arg, "--blendfile") == 0) {
            value = window3d_option_value(argc, argv, &i, arg);
            if (value == NULL) return FAIL;
            options->blendfile = value;
        }
        else if (strncmp(arg, "-B", 2) == 0 && arg[2] != '\0') {
            options->blendfile = arg + 2;
        }
        else if (strncmp(arg, "--blendfile=", 12) == 0) {
            options->blendfile = arg + 12;
        }
        else if (strcmp(arg, "-C") == 0 || strcmp(arg, "--clobber") == 0) {
            value = window3d_option_value(argc, argv, &i, arg);
            if (value == NULL || window3d_parse_clobber(value, options) != SUCCESS) return FAIL;
        }
        else if (strncmp(arg, "-C", 2) == 0 && arg[2] != '\0') {
            if (window3d_parse_clobber(arg + 2, options) != SUCCESS) return FAIL;
        }
        else if (strncmp(arg, "--clobber=", 10) == 0) {
            if (window3d_parse_clobber(arg + 10, options) != SUCCESS) return FAIL;
        }
        else if (strcmp(arg, "-M") == 0 || strcmp(arg, "--monotone") == 0) {
            value = window3d_option_value(argc, argv, &i, arg);
            if (value == NULL || window3d_parse_monotone_method(value, options) != SUCCESS) return FAIL;
        }
        else if (strncmp(arg, "-M", 2) == 0 && arg[2] != '\0') {
            if (window3d_parse_monotone_method(arg + 2, options) != SUCCESS) return FAIL;
        }
        else if (strncmp(arg, "--monotone=", 11) == 0) {
            if (window3d_parse_monotone_method(arg + 11, options) != SUCCESS) return FAIL;
        }
        else if (strcmp(arg, "-N") == 0 || strcmp(arg, "--write-monotone") == 0) {
            options->has_write_monotone = 1;
        }
        else {
            BLEND_Report(BLEND_MSG_ERROR, "window3d: unknown option: %s\n", arg);
            return FAIL;
        }
    }

    if (!options->has_region) {
        BLEND_Report(BLEND_MSG_ERROR, "window3d: -R/--region is required\n");
        return FAIL;
    }
    if (!options->has_increment) {
        BLEND_Report(BLEND_MSG_ERROR, "window3d: -I/--increment is required\n");
        return FAIL;
    }
    if (options->dx <= 0.0 || options->dy <= 0.0 || options->dz <= 0.0) {
        BLEND_Report(BLEND_MSG_ERROR, "window3d: increments must be positive\n");
        return FAIL;
    }

    if (options->blendfile != NULL) {
        if (options->has_function) {
            BLEND_Report(BLEND_MSG_WARNING, "window3d: -F/--function is ignored when -B/--blendfile is given\n");
        }
        if (options->has_taper) {
            BLEND_Report(BLEND_MSG_WARNING, "window3d: -T/--taper_ratio is ignored when -B/--blendfile is given\n");
        }
    }
    else if (options->has_taper &&
             (options->x_function == WFUNC_BOXCAR ||
              options->y_function == WFUNC_BOXCAR ||
              options->z_function == WFUNC_BOXCAR)) {
        BLEND_Report(BLEND_MSG_WARNING, "window3d: taper ratios are ignored for boxcar window dimensions\n");
    }
    if (options->blendfile == NULL && options->has_monotone_method) {
        BLEND_Report(BLEND_MSG_WARNING, "window3d: -M/--monotone is ignored unless -B/--blendfile is given\n");
    }
    if (options->has_write_monotone && (options->blendfile == NULL || !options->has_monotone_method)) {
        BLEND_Report(BLEND_MSG_WARNING, "window3d: -N/--write-monotone is ignored unless -B/--blendfile and -M/--monotone are given\n");
    }

    if (window3d_adjust_axis("x", options->xmin, &options->xmax, options->dx, &options->nx) != SUCCESS ||
        window3d_adjust_axis("y", options->ymin, &options->ymax, options->dy, &options->ny) != SUCCESS ||
        window3d_adjust_axis("z", options->zmin, &options->zmax, options->dz, &options->nz) != SUCCESS) {
        return FAIL;
    }

    return SUCCESS;
}

static char *window3d_trim_line(char *line)
{
    char *end;

    while (isspace((unsigned char)*line)) {
        line++;
    }
    if (*line == '\0') {
        return line;
    }
    end = line + strlen(line) - 1;
    while (end > line && isspace((unsigned char)*end)) {
        *end = '\0';
        end--;
    }
    return line;
}

static void window3d_support_clear(window3d_support *support)
{
    if (support == NULL) {
        return;
    }
    free(support->vertices);
    support->vertices = NULL;
    support->n_vertices = 0;
    blend_window_boundary_clear(&support->data);
}

static void window3d_support_list_free(window3d_support_list *supports)
{
    size_t i;

    for (i = 0; i < supports->count; i++) {
        window3d_support_clear(&supports->items[i]);
    }
    free(supports->items);
    supports->items = NULL;
    supports->count = 0;
    supports->capacity = 0;
}

static int window3d_support_list_append(window3d_support_list *supports, window3d_support *support)
{
    window3d_support *items;
    size_t capacity;

    if (supports->count == supports->capacity) {
        capacity = supports->capacity == 0 ? 8 : supports->capacity * 2;
        items = (window3d_support *)realloc(supports->items, capacity * sizeof(*supports->items));
        if (items == NULL) {
            BLEND_Report(BLEND_MSG_ERROR, "window3d: could not allocate polygon support storage\n");
            return FAIL;
        }
        supports->items = items;
        supports->capacity = capacity;
    }

    supports->items[supports->count] = *support;
    memset(support, 0, sizeof(*support));
    supports->count++;
    return SUCCESS;
}

static int window3d_support_grid_size(double min, double max, double increment, int *n)
{
    double span = max - min;
    double intervals = floor(span / increment + 0.5);

    if (intervals < 1.0 || intervals > (double)(INT_MAX - 1)) {
        BLEND_Report(BLEND_MSG_ERROR, "window3d: support produces an invalid grid size\n");
        return FAIL;
    }

    *n = (int)intervals + 1;
    return SUCCESS;
}

static int window3d_monotone_output_path(const char *polygon_file, char *path, size_t path_size)
{
    int n;
    const char *slash;
    const char *dot;
    size_t stem_len;

    slash = strrchr(polygon_file, '/');
    dot = strrchr(polygon_file, '.');
    if (dot != NULL && (slash == NULL || dot > slash)) {
        stem_len = (size_t)(dot - polygon_file);
        n = snprintf(path, path_size, "%.*s_monotone%s", (int)stem_len, polygon_file, dot);
    }
    else {
        n = snprintf(path, path_size, "%s_monotone", polygon_file);
    }
    if (n < 0 || (size_t)n >= path_size) {
        BLEND_Report(BLEND_MSG_ERROR, "window3d: monotone polygon output path is too long for %s\n", polygon_file);
        return FAIL;
    }

    return SUCCESS;
}

static int window3d_write_monotone_polygon(const char *polygon_file, const polygon *poly)
{
    char path[PATH_MAX];

    if (window3d_monotone_output_path(polygon_file, path, sizeof(path)) != SUCCESS) {
        return FAIL;
    }
    if (blend_polygon_write(path, poly) != SUCCESS) {
        return FAIL;
    }

    BLEND_Report(BLEND_MSG_WARNING, "window3d: wrote modified xy-monotone polygon to %s\n", path);
    return SUCCESS;
}

static int window3d_append_snapped_vertex(vertex **vertices, size_t *count, size_t *capacity,
                                          double x, double y)
{
    vertex *items;
    size_t new_capacity;

    if (*count > 0 &&
        fabs((*vertices)[*count - 1].x - x) <= 1.0e-12 &&
        fabs((*vertices)[*count - 1].y - y) <= 1.0e-12) {
        return SUCCESS;
    }

    if (*count == *capacity) {
        new_capacity = *capacity == 0 ? 16 : *capacity * 2;
        items = (vertex *)realloc(*vertices, new_capacity * sizeof(**vertices));
        if (items == NULL) {
            BLEND_Report(BLEND_MSG_ERROR, "window3d: could not allocate local polygon vertices\n");
            return FAIL;
        }
        *vertices = items;
        *capacity = new_capacity;
    }

    (*vertices)[*count].x = x;
    (*vertices)[*count].y = y;
    *count += 1;
    return SUCCESS;
}

static double window3d_snap_local_coordinate(double value, int max_index)
{
    double snapped = floor(value + 0.5);

    if (fabs(value) < 1.0e-10) {
        snapped = 0.0;
    }
    else if (fabs(value - (double)max_index) < 1.0e-10) {
        snapped = (double)max_index;
    }

    if (snapped < 0.0) {
        return 0.0;
    }
    if (snapped > (double)max_index) {
        return (double)max_index;
    }
    return snapped;
}

static int window3d_polygon_geometry(const window3d_options *options, const polygon *poly,
                                     const char *polygon_file,
                                     double *xmin, double *xmax, double *ymin, double *ymax,
                                     int *nx, int *ny)
{
    size_t i;

    if (poly == NULL || poly->vertices == NULL || poly->n_vertices == 0 ||
        xmin == NULL || xmax == NULL || ymin == NULL || ymax == NULL ||
        nx == NULL || ny == NULL) {
        return FAIL;
    }

    *xmin = poly->vertices[0].x;
    *xmax = poly->vertices[0].x;
    *ymin = poly->vertices[0].y;
    *ymax = poly->vertices[0].y;

    for (i = 0; i < poly->n_vertices; i++) {
        double x = poly->vertices[i].x;
        double y = poly->vertices[i].y;

        if (x < options->xmin || x > options->xmax ||
            y < options->ymin || y > options->ymax) {
            BLEND_Report(BLEND_MSG_ERROR, "window3d: polygon %s has a vertex outside -R\n", polygon_file);
            return FAIL;
        }
        if (x < *xmin) *xmin = x;
        if (x > *xmax) *xmax = x;
        if (y < *ymin) *ymin = y;
        if (y > *ymax) *ymax = y;
    }

    if (*xmin >= *xmax || *ymin >= *ymax) {
        BLEND_Report(BLEND_MSG_ERROR, "window3d: polygon %s has an invalid bounding box\n", polygon_file);
        return FAIL;
    }

    return window3d_support_grid_size(*xmin, *xmax, options->dx, nx) == SUCCESS &&
           window3d_support_grid_size(*ymin, *ymax, options->dy, ny) == SUCCESS
               ? SUCCESS
               : FAIL;
}

static int window3d_polygon_to_local_grid(const polygon *real_poly,
                                          double xmin, double xmax,
                                          double ymin, double ymax,
                                          int nx, int ny, polygon *local_poly)
{
    vertex *vertices = NULL;
    size_t count = 0;
    size_t capacity = 0;
    polygon tmp = {0};
    size_t i;

    if (real_poly == NULL || local_poly == NULL || nx < 2 || ny < 2 ||
        xmin >= xmax || ymin >= ymax) {
        return FAIL;
    }

    for (i = 0; i < real_poly->n_vertices; i++) {
        double local_x = (real_poly->vertices[i].x - xmin) * (double)(nx - 1) / (xmax - xmin);
        double local_y = (real_poly->vertices[i].y - ymin) * (double)(ny - 1) / (ymax - ymin);

        local_x = window3d_snap_local_coordinate(local_x, nx - 1);
        local_y = window3d_snap_local_coordinate(local_y, ny - 1);
        if (window3d_append_snapped_vertex(&vertices, &count, &capacity, local_x, local_y) != SUCCESS) {
            free(vertices);
            return FAIL;
        }
    }

    if (count > 3 &&
        fabs(vertices[0].x - vertices[count - 1].x) <= 1.0e-12 &&
        fabs(vertices[0].y - vertices[count - 1].y) <= 1.0e-12) {
        count--;
    }

    tmp.n_vertices = count;
    tmp.vertices = vertices;
    if (blend_polygon_validate(&tmp) != SUCCESS) {
        free(vertices);
        return FAIL;
    }

    blend_polygon_free(local_poly);
    local_poly->n_vertices = tmp.n_vertices;
    local_poly->vertices = tmp.vertices;
    return SUCCESS;
}

static int window3d_local_grid_to_real_polygon(const polygon *local_poly,
                                               double xmin, double xmax,
                                               double ymin, double ymax,
                                               int nx, int ny, polygon *real_poly)
{
    polygon tmp = {0};
    size_t i;

    if (local_poly == NULL || real_poly == NULL || nx < 2 || ny < 2 ||
        xmin >= xmax || ymin >= ymax) {
        return FAIL;
    }
    if (blend_polygon_alloc(&tmp, local_poly->n_vertices) != SUCCESS) {
        return FAIL;
    }

    for (i = 0; i < local_poly->n_vertices; i++) {
        double x = xmin + local_poly->vertices[i].x * (xmax - xmin) / (double)(nx - 1);
        double y = ymin + local_poly->vertices[i].y * (ymax - ymin) / (double)(ny - 1);

        if (blend_polygon_set_vertex(&tmp, i, x, y) != SUCCESS) {
            blend_polygon_free(&tmp);
            return FAIL;
        }
    }

    if (blend_polygon_validate(&tmp) != SUCCESS) {
        blend_polygon_free(&tmp);
        return FAIL;
    }

    blend_polygon_free(real_poly);
    real_poly->n_vertices = tmp.n_vertices;
    real_poly->vertices = tmp.vertices;
    return SUCCESS;
}

static int window3d_build_support_from_local(const char *polygon_file,
                                             const polygon *real_poly,
                                             const polygon *local_poly,
                                             double xmin, double xmax,
                                             double ymin, double ymax,
                                             int nx, int ny,
                                             double zlo, double zhi, int nz,
                                             blend_window_function x_function,
                                             blend_window_function y_function,
                                             blend_window_function z_function,
                                             double rx1, double rx2,
                                             double ry1, double ry2,
                                             double rz1, double rz2,
                                             window3d_support *support)
{
    permuted_vertex pv;

    if (real_poly == NULL || local_poly == NULL || support == NULL ||
        real_poly->vertices == NULL || local_poly->vertices == NULL) {
        return FAIL;
    }

    memset(support, 0, sizeof(*support));
    memset(&pv, 0, sizeof(pv));

    support->vertices = (vertex *)calloc(real_poly->n_vertices, sizeof(*support->vertices));
    if (support->vertices == NULL) {
        BLEND_Report(BLEND_MSG_ERROR, "window3d: could not allocate polygon support vertices\n");
        return FAIL;
    }
    memcpy(support->vertices, real_poly->vertices, real_poly->n_vertices * sizeof(*support->vertices));
    support->n_vertices = real_poly->n_vertices;
    support->xmin = xmin;
    support->xmax = xmax;
    support->ymin = ymin;
    support->ymax = ymax;
    support->zlo = zlo;
    support->zhi = zhi;
    support->data.nx = nx;
    support->data.ny = ny;
    support->data.nz = nz;
    support->data.ratio_x1 = rx1;
    support->data.ratio_x2 = rx2;
    support->data.ratio_y1 = ry1;
    support->data.ratio_y2 = ry2;
    support->data.ratio_z1 = rz1;
    support->data.ratio_z2 = rz2;
    support->data.x_function = x_function;
    support->data.y_function = y_function;
    support->data.z_function = z_function;

    if (blend_window_set_polygon(&support->data, local_poly) != SUCCESS) {
        BLEND_Report(BLEND_MSG_ERROR, "window3d: could not assign local polygon support\n");
        window3d_support_clear(support);
        return FAIL;
    }

    if (boundary_assembly(&support->data, &pv) != SUCCESS) {
        BLEND_Report(BLEND_MSG_ERROR, "window3d: could not assemble boundary for polygon %s\n", polygon_file);
        blend_permuted_vertex_free(&pv);
        window3d_support_clear(support);
        return FAIL;
    }
    blend_permuted_vertex_free(&pv);
    return SUCCESS;
}

static int window3d_make_polygon(vertex *real_vertices, size_t n_vertices, window3d_support *support,
                                 double *xmin, double *xmax, double *ymin, double *ymax)
{
    size_t i;
    polygon local_poly;

    memset(&local_poly, 0, sizeof(local_poly));
    if (blend_polygon_alloc(&local_poly, n_vertices) != SUCCESS) {
        BLEND_Report(BLEND_MSG_ERROR, "window3d: could not allocate local polygon\n");
        return FAIL;
    }

    for (i = 0; i < n_vertices; i++) {
        double x = real_vertices[i].x;
        double y = real_vertices[i].y;
        double local_x = (x - *xmin) * (double)(support->data.nx - 1) / (*xmax - *xmin);
        double local_y = (y - *ymin) * (double)(support->data.ny - 1) / (*ymax - *ymin);

        if (fabs(local_x) < 1.0e-10) local_x = 0.0;
        if (fabs(local_y) < 1.0e-10) local_y = 0.0;
        if (fabs(local_x - (double)(support->data.nx - 1)) < 1.0e-10) local_x = (double)(support->data.nx - 1);
        if (fabs(local_y - (double)(support->data.ny - 1)) < 1.0e-10) local_y = (double)(support->data.ny - 1);

        if (blend_polygon_set_vertex(&local_poly, i, local_x, local_y) != SUCCESS) {
            blend_polygon_free(&local_poly);
            return FAIL;
        }
    }

    if (blend_window_set_polygon(&support->data, &local_poly) != SUCCESS) {
        BLEND_Report(BLEND_MSG_ERROR, "window3d: could not assign polygon support\n");
        blend_polygon_free(&local_poly);
        return FAIL;
    }

    blend_polygon_free(&local_poly);
    return SUCCESS;
}

static int window3d_build_support(const window3d_options *options, const char *polygon_file,
                                  const vertex *vertices, size_t n_vertices,
                                  double zlo, double zhi,
                                  blend_window_function x_function,
                                  blend_window_function y_function,
                                  blend_window_function z_function,
                                  double rx1, double rx2,
                                  double ry1, double ry2,
                                  double rz1, double rz2,
                                  window3d_support *support)
{
    size_t i;
    permuted_vertex pv;

    memset(support, 0, sizeof(*support));
    memset(&pv, 0, sizeof(pv));

    support->vertices = (vertex *)calloc(n_vertices, sizeof(*support->vertices));
    if (support->vertices == NULL) {
        BLEND_Report(BLEND_MSG_ERROR, "window3d: could not allocate polygon support vertices\n");
        return FAIL;
    }
    memcpy(support->vertices, vertices, n_vertices * sizeof(*support->vertices));
    support->n_vertices = n_vertices;

    support->xmin = vertices[0].x;
    support->xmax = vertices[0].x;
    support->ymin = vertices[0].y;
    support->ymax = vertices[0].y;
    support->zlo = zlo;
    support->zhi = zhi;

    if (zlo < options->zmin || zhi > options->zmax || zlo >= zhi) {
        BLEND_Report(BLEND_MSG_ERROR, "window3d: support %s has invalid zlo/zhi outside -R\n", polygon_file);
        window3d_support_clear(support);
        return FAIL;
    }

    for (i = 0; i < n_vertices; i++) {
        if (vertices[i].x < options->xmin || vertices[i].x > options->xmax ||
            vertices[i].y < options->ymin || vertices[i].y > options->ymax) {
            BLEND_Report(BLEND_MSG_ERROR, "window3d: polygon %s has a vertex outside -R\n", polygon_file);
            window3d_support_clear(support);
            return FAIL;
        }
        if (vertices[i].x < support->xmin) support->xmin = vertices[i].x;
        if (vertices[i].x > support->xmax) support->xmax = vertices[i].x;
        if (vertices[i].y < support->ymin) support->ymin = vertices[i].y;
        if (vertices[i].y > support->ymax) support->ymax = vertices[i].y;
    }

    if (support->xmin >= support->xmax || support->ymin >= support->ymax) {
        BLEND_Report(BLEND_MSG_ERROR, "window3d: polygon %s has an invalid bounding box\n", polygon_file);
        window3d_support_clear(support);
        return FAIL;
    }

    if (window3d_support_grid_size(support->xmin, support->xmax, options->dx, &support->data.nx) != SUCCESS ||
        window3d_support_grid_size(support->ymin, support->ymax, options->dy, &support->data.ny) != SUCCESS ||
        window3d_support_grid_size(support->zlo, support->zhi, options->dz, &support->data.nz) != SUCCESS) {
        window3d_support_clear(support);
        return FAIL;
    }

    support->data.ratio_x1 = rx1;
    support->data.ratio_x2 = rx2;
    support->data.ratio_y1 = ry1;
    support->data.ratio_y2 = ry2;
    support->data.ratio_z1 = rz1;
    support->data.ratio_z2 = rz2;
    support->data.x_function = x_function;
    support->data.y_function = y_function;
    support->data.z_function = z_function;

    if (window3d_make_polygon(support->vertices, support->n_vertices, support,
                              &support->xmin, &support->xmax, &support->ymin, &support->ymax) != SUCCESS) {
        window3d_support_clear(support);
        return FAIL;
    }

    if (boundary_assembly(&support->data, &pv) != SUCCESS) {
        BLEND_Report(BLEND_MSG_ERROR, "window3d: could not assemble boundary for polygon %s\n", polygon_file);
        blend_permuted_vertex_free(&pv);
        window3d_support_clear(support);
        return FAIL;
    }
    blend_permuted_vertex_free(&pv);

    return SUCCESS;
}

static int window3d_parse_blendfile_line(const char *line, long line_number,
                                         const window3d_options *options,
                                         window3d_support *support)
{
    char polygon_file[PATH_MAX];
    char function_triple[128];
    char taper_ratio[160];
    char extra[2];
    polygon input_polygon = {0};
    polygon support_polygon = {0};
    polygon local_input_polygon = {0};
    polygon local_support_polygon = {0};
    blend_window_function x_function, y_function, z_function;
    double xmin, xmax, ymin, ymax;
    double zlo, zhi;
    double rx1, rx2, ry1, ry2, rz1, rz2;
    int nx, ny, nz;
    int is_xy_monotone = 0;
    int fields;
    int status;

    fields = sscanf(line, "%1023s %lf %lf %127s %159s %1s",
                    polygon_file, &zlo, &zhi, function_triple, taper_ratio, extra);
    if (fields != 5 || !isfinite(zlo) || !isfinite(zhi)) {
        BLEND_Report(BLEND_MSG_ERROR,
                     "window3d: blendfile line %ld must have five fields: polygon zlo zhi function_triple taper_ratio\n",
                     line_number);
        return FAIL;
    }

    if (window3d_parse_function_triple(function_triple, &x_function, &y_function, &z_function) != SUCCESS) {
        BLEND_Report(BLEND_MSG_ERROR, "window3d: invalid window function triple on blendfile line %ld\n", line_number);
        return FAIL;
    }
    if (window3d_parse_taper_ratio_values(taper_ratio, "blendfile",
                                          &rx1, &rx2, &ry1, &ry2, &rz1, &rz2) != SUCCESS) {
        BLEND_Report(BLEND_MSG_ERROR, "window3d: invalid taper ratios on blendfile line %ld\n", line_number);
        return FAIL;
    }
    if ((x_function == WFUNC_BOXCAR || y_function == WFUNC_BOXCAR || z_function == WFUNC_BOXCAR) &&
        (rx1 != 0.0 || rx2 != 0.0 || ry1 != 0.0 || ry2 != 0.0 || rz1 != 0.0 || rz2 != 0.0)) {
        BLEND_Report(BLEND_MSG_WARNING,
                     "window3d: taper ratios on blendfile line %ld are ignored for boxcar window dimensions\n",
                     line_number);
    }

    if (blend_polygon_read(polygon_file, &input_polygon) != SUCCESS) {
        return FAIL;
    }

    if (zlo < options->zmin || zhi > options->zmax || zlo >= zhi) {
        BLEND_Report(BLEND_MSG_ERROR, "window3d: support %s has invalid zlo/zhi outside -R\n", polygon_file);
        blend_polygon_free(&input_polygon);
        return FAIL;
    }

    if (window3d_polygon_geometry(options, &input_polygon, polygon_file,
                                  &xmin, &xmax, &ymin, &ymax, &nx, &ny) != SUCCESS ||
        window3d_support_grid_size(zlo, zhi, options->dz, &nz) != SUCCESS ||
        window3d_polygon_to_local_grid(&input_polygon, xmin, xmax, ymin, ymax, nx, ny,
                                       &local_input_polygon) != SUCCESS) {
        blend_polygon_free(&input_polygon);
        return FAIL;
    }

    if (blend_polygon_is_xy_monotone_strict(&local_input_polygon, &is_xy_monotone) != SUCCESS) {
        blend_polygon_free(&input_polygon);
        blend_polygon_free(&local_input_polygon);
        return FAIL;
    }
    if (!is_xy_monotone) {
        if (options->monotone_method != 'E' && options->monotone_method != 'B') {
            BLEND_Report(BLEND_MSG_ERROR,
                         "window3d: grid-snapped polygon %s is not strict enough for boundary assembly; use -ME or -MB if you want BLEND to refine it\n",
                         polygon_file);
            blend_polygon_free(&input_polygon);
            blend_polygon_free(&local_input_polygon);
            return FAIL;
        }
        if (options->monotone_method == 'E') {
            BLEND_Report(BLEND_MSG_WARNING,
                         "window3d: grid-snapped polygon %s is not strict enough for boundary assembly; using its strict xy-monotone envelope from -ME\n",
                         polygon_file);
        }
        else if (options->monotone_method == 'B') {
            BLEND_Report(BLEND_MSG_WARNING,
                         "window3d: grid-snapped polygon %s is not strict enough for boundary assembly; using its strict best IoU piecewise xy-monotone envelope from -MB\n",
                         polygon_file);
        }
    }

    if (is_xy_monotone) {
        status = blend_polygon_copy(&local_input_polygon, &local_support_polygon);
    }
    else {
        status = options->monotone_method == 'B'
                     ? blend_polygon_xy_monotone_best_piecewise_envelope_strict(
                           &local_input_polygon, &local_support_polygon,
                           0.0, (double)(nx - 1), 0.0, (double)(ny - 1), nx, ny)
                     : blend_polygon_xy_monotone_envelope_strict(&local_input_polygon,
                                                                 &local_support_polygon);
    }
    if (status != SUCCESS) {
        blend_polygon_free(&input_polygon);
        blend_polygon_free(&local_input_polygon);
        blend_polygon_free(&local_support_polygon);
        return FAIL;
    }

    if (window3d_local_grid_to_real_polygon(&local_support_polygon, xmin, xmax, ymin, ymax,
                                            nx, ny, &support_polygon) != SUCCESS) {
        blend_polygon_free(&input_polygon);
        blend_polygon_free(&local_input_polygon);
        blend_polygon_free(&local_support_polygon);
        return FAIL;
    }

    status = window3d_build_support_from_local(polygon_file, &support_polygon,
                                               &local_support_polygon,
                                               xmin, xmax, ymin, ymax, nx, ny,
                                               zlo, zhi, nz,
                                               x_function, y_function, z_function,
                                               rx1, rx2, ry1, ry2, rz1, rz2, support);
    if (status != SUCCESS) {
        blend_polygon_free(&input_polygon);
        blend_polygon_free(&support_polygon);
        blend_polygon_free(&local_input_polygon);
        blend_polygon_free(&local_support_polygon);
        return FAIL;
    }

    if (!is_xy_monotone) {
        BLEND_Report(BLEND_MSG_WARNING,
                     "%s: original number of vertices = %zu, final number of vertices = %zu.\n",
                     polygon_file, input_polygon.n_vertices, support_polygon.n_vertices);
        if (options->has_write_monotone) {
            if (window3d_write_monotone_polygon(polygon_file, &support_polygon) != SUCCESS) {
                blend_polygon_free(&input_polygon);
                blend_polygon_free(&support_polygon);
                blend_polygon_free(&local_input_polygon);
                blend_polygon_free(&local_support_polygon);
                return FAIL;
            }
        }
        else {
            BLEND_Report(BLEND_MSG_WARNING,
                         "window3d: modified polygon for %s will not be written; use -N to write it\n",
                         polygon_file);
        }
    }

    blend_polygon_free(&input_polygon);
    blend_polygon_free(&support_polygon);
    blend_polygon_free(&local_input_polygon);
    blend_polygon_free(&local_support_polygon);
    return status;
}

static int window3d_read_blendfile(const window3d_options *options, window3d_support_list *supports)
{
    FILE *fp;
    char line[2048];
    long line_number = 0;

    memset(supports, 0, sizeof(*supports));

    fp = fopen(options->blendfile, "r");
    if (fp == NULL) {
        BLEND_Report(BLEND_MSG_ERROR, "%s: %s\n", options->blendfile, strerror(errno));
        return FAIL;
    }

    while (fgets(line, sizeof(line), fp) != NULL) {
        char *text;
        char *comment;
        window3d_support support;

        line_number++;
        comment = strchr(line, '#');
        if (comment != NULL) {
            *comment = '\0';
        }
        text = window3d_trim_line(line);
        if (*text == '\0') {
            continue;
        }

        if (window3d_parse_blendfile_line(text, line_number, options, &support) != SUCCESS) {
            fclose(fp);
            window3d_support_list_free(supports);
            return FAIL;
        }
        if (window3d_support_list_append(supports, &support) != SUCCESS) {
            window3d_support_clear(&support);
            fclose(fp);
            window3d_support_list_free(supports);
            return FAIL;
        }
    }

    if (ferror(fp)) {
        BLEND_Report(BLEND_MSG_ERROR, "%s: %s\n", options->blendfile, strerror(errno));
        fclose(fp);
        window3d_support_list_free(supports);
        return FAIL;
    }
    fclose(fp);

    if (supports->count == 0) {
        BLEND_Report(BLEND_MSG_ERROR, "window3d: blendfile contains no supports\n");
        return FAIL;
    }

    return SUCCESS;
}

static int window3d_make_full_support(const window3d_options *options, window3d_support *support)
{
    vertex vertices[4];

    vertices[0].x = options->xmin;
    vertices[0].y = options->ymin;
    vertices[1].x = options->xmax;
    vertices[1].y = options->ymin;
    vertices[2].x = options->xmax;
    vertices[2].y = options->ymax;
    vertices[3].x = options->xmin;
    vertices[3].y = options->ymax;

    return window3d_build_support(options, "full -R domain", vertices, 4,
                                  options->zmin, options->zmax,
                                  options->x_function, options->y_function, options->z_function,
                                  options->ratio_x1, options->ratio_x2,
                                  options->ratio_y1, options->ratio_y2,
                                  options->ratio_z1, options->ratio_z2,
                                  support);
}

static double window3d_coordinate_x(const window3d_options *options, int index)
{
    return options->xmin + (double)index * options->dx;
}

static double window3d_coordinate_y(const window3d_options *options, int index)
{
    return options->ymin + (double)index * options->dy;
}

static double window3d_coordinate_z(const window3d_options *options, int index)
{
    return options->zmin + (double)index * options->dz;
}

static int window3d_point_on_segment(double x, double y, const vertex *a, const vertex *b)
{
    double dx = b->x - a->x;
    double dy = b->y - a->y;
    double cross = (x - a->x) * dy - (y - a->y) * dx;
    double scale = fmax(1.0, fmax(fabs(dx), fabs(dy)));
    double tol = 1.0e-10 * scale;

    if (fabs(cross) > tol) {
        return 0;
    }
    if (x < fmin(a->x, b->x) - tol || x > fmax(a->x, b->x) + tol ||
        y < fmin(a->y, b->y) - tol || y > fmax(a->y, b->y) + tol) {
        return 0;
    }
    return 1;
}

static int window3d_point_in_polygon(double x, double y, const window3d_support *support)
{
    int inside = 0;
    size_t i, j;

    for (i = 0, j = support->n_vertices - 1; i < support->n_vertices; j = i++) {
        const vertex *vi = &support->vertices[i];
        const vertex *vj = &support->vertices[j];

        if (window3d_point_on_segment(x, y, vj, vi)) {
            return 1;
        }
        if (((vi->y > y) != (vj->y > y)) &&
            (x < (vj->x - vi->x) * (y - vi->y) / (vj->y - vi->y) + vi->x)) {
            inside = !inside;
        }
    }

    return inside;
}

static int window3d_support_grid_weight(int ix, int iy, int iz, window3d_support *support, double *weight)
{
    if (ix < 0 || iy < 0 || iz < 0 ||
        ix >= support->data.nx || iy >= support->data.ny || iz >= support->data.nz) {
        *weight = 0.0;
        return SUCCESS;
    }

    if (embedding_contribution3d(ix, iy, iz, &support->data) != SUCCESS) {
        return FAIL;
    }

    *weight = support->data.contribution;
    return SUCCESS;
}

static int window3d_support_weight(double x, double y, double z, window3d_support *support, double *weight)
{
    double sx, sy, sz;
    double w000, w100, w010, w110, w001, w101, w011, w111;
    int ix0, iy0, iz0, ix1, iy1, iz1;

    if (x < support->xmin || x > support->xmax ||
        y < support->ymin || y > support->ymax ||
        z < support->zlo || z > support->zhi ||
        !window3d_point_in_polygon(x, y, support)) {
        *weight = 0.0;
        return SUCCESS;
    }

    sx = (x - support->xmin) * (double)(support->data.nx - 1) / (support->xmax - support->xmin);
    sy = (y - support->ymin) * (double)(support->data.ny - 1) / (support->ymax - support->ymin);
    sz = (z - support->zlo) * (double)(support->data.nz - 1) / (support->zhi - support->zlo);

    ix0 = (int)floor(sx);
    iy0 = (int)floor(sy);
    iz0 = (int)floor(sz);

    if (ix0 < 0) ix0 = 0;
    if (iy0 < 0) iy0 = 0;
    if (iz0 < 0) iz0 = 0;
    if (ix0 >= support->data.nx - 1) ix0 = support->data.nx - 2;
    if (iy0 >= support->data.ny - 1) iy0 = support->data.ny - 2;
    if (iz0 >= support->data.nz - 1) iz0 = support->data.nz - 2;

    ix1 = ix0 < support->data.nx - 1 ? ix0 + 1 : ix0;
    iy1 = iy0 < support->data.ny - 1 ? iy0 + 1 : iy0;
    iz1 = iz0 < support->data.nz - 1 ? iz0 + 1 : iz0;

    if (window3d_support_grid_weight(ix0, iy0, iz0, support, &w000) != SUCCESS ||
        window3d_support_grid_weight(ix1, iy0, iz0, support, &w100) != SUCCESS ||
        window3d_support_grid_weight(ix0, iy1, iz0, support, &w010) != SUCCESS ||
        window3d_support_grid_weight(ix1, iy1, iz0, support, &w110) != SUCCESS ||
        window3d_support_grid_weight(ix0, iy0, iz1, support, &w001) != SUCCESS ||
        window3d_support_grid_weight(ix1, iy0, iz1, support, &w101) != SUCCESS ||
        window3d_support_grid_weight(ix0, iy1, iz1, support, &w011) != SUCCESS ||
        window3d_support_grid_weight(ix1, iy1, iz1, support, &w111) != SUCCESS) {
        return FAIL;
    }

    return interpolate_trilinear((double)ix0, (double)ix1,
                                 (double)iy0, (double)iy1,
                                 (double)iz0, (double)iz1,
                                 w000, w100, w010, w110,
                                 w001, w101, w011, w111,
                                 sx, sy, sz, weight);
}

static int window3d_combine_weights(const double *weights, size_t count, char clobber, double *weight)
{
    size_t i;

    if (count == 0) {
        *weight = 0.0;
        return SUCCESS;
    }

    switch (clobber) {
        case 'f':
            *weight = weights[0];
            break;
        case 'o':
            *weight = weights[count - 1];
            break;
        case 'l':
            *weight = weights[0];
            for (i = 1; i < count; i++) if (weights[i] < *weight) *weight = weights[i];
            break;
        case 'u':
            *weight = weights[0];
            for (i = 1; i < count; i++) if (weights[i] > *weight) *weight = weights[i];
            break;
        case 'a':
            *weight = 0.0;
            for (i = 0; i < count; i++) *weight += weights[i];
            *weight /= (double)count;
            break;
        case 'g':
            *weight = 1.0;
            for (i = 0; i < count; i++) *weight *= weights[i];
            *weight = pow(*weight, 1.0 / (double)count);
            break;
        case 'p':
            *weight = 1.0;
            for (i = 0; i < count; i++) *weight *= weights[i];
            break;
        default:
            BLEND_Report(BLEND_MSG_ERROR, "window3d: unknown clobber mode: %c\n", clobber);
            return FAIL;
    }

    return SUCCESS;
}

static int window3d_blendfile_weight(double x, double y, double z,
                                     const window3d_options *options,
                                     window3d_support_list *supports,
                                     double *weight)
{
    double *weights = NULL;
    size_t i, count = 0;
    int status;

    weights = (double *)calloc(supports->count, sizeof(*weights));
    if (weights == NULL) {
        BLEND_Report(BLEND_MSG_ERROR, "window3d: could not allocate overlap weights\n");
        return FAIL;
    }

    for (i = 0; i < supports->count; i++) {
        double support_weight;

        if (window3d_support_weight(x, y, z, &supports->items[i], &support_weight) != SUCCESS) {
            free(weights);
            return FAIL;
        }
        if (support_weight > 0.0 ||
            (z >= supports->items[i].zlo && z <= supports->items[i].zhi &&
             window3d_point_in_polygon(x, y, &supports->items[i]))) {
            weights[count] = support_weight;
            count++;
        }
    }

    status = window3d_combine_weights(weights, count, options->clobber, weight);
    free(weights);
    return status;
}

static int window3d_get_weight(double x, double y, double z, int ix, int iy, int iz,
                               const window3d_options *options,
                               window3d_support *full_support,
                               window3d_support_list *supports,
                               double *weight)
{
    if (supports != NULL) {
        return window3d_blendfile_weight(x, y, z, options, supports, weight);
    }

    if (ix >= 0 && iy >= 0 && iz >= 0) {
        return window3d_support_grid_weight(ix, iy, iz, full_support, weight);
    }

    return window3d_support_weight(x, y, z, full_support, weight);
}

static int window3d_write_weight(FILE *fp, double x, double y, double z, double weight)
{
    return fprintf(fp, "%.12g %.12g %.12g %.12f\n", x, y, z, weight) < 0 ? FAIL : SUCCESS;
}

static int window3d_write_grid(FILE *fp, const window3d_options *options,
                               window3d_support *full_support,
                               window3d_support_list *supports)
{
    int ix, iy, iz;

    for (iz = 0; iz < options->nz; iz++) {
        double z = window3d_coordinate_z(options, iz);
        for (iy = 0; iy < options->ny; iy++) {
            double y = window3d_coordinate_y(options, iy);
            for (ix = 0; ix < options->nx; ix++) {
                double x = window3d_coordinate_x(options, ix);
                double weight;

                if (window3d_get_weight(x, y, z, ix, iy, iz, options,
                                        full_support, supports, &weight) != SUCCESS) {
                    return FAIL;
                }
                if (window3d_write_weight(fp, x, y, z, weight) != SUCCESS) {
                    return FAIL;
                }
            }
        }
    }

    return SUCCESS;
}

static int window3d_write_queries(FILE *fp, const window3d_options *options,
                                  window3d_support *full_support,
                                  window3d_support_list *supports,
                                  int *n_queries)
{
    double x, y, z;
    int status;

    *n_queries = 0;
    while ((status = scanf("%lf %lf %lf", &x, &y, &z)) == 3) {
        double weight;

        if (!isfinite(x) || !isfinite(y) || !isfinite(z)) {
            BLEND_Report(BLEND_MSG_ERROR, "window3d: standard input contains a non-finite query\n");
            return FAIL;
        }
        if (window3d_get_weight(x, y, z, -1, -1, -1, options,
                                full_support, supports, &weight) != SUCCESS) {
            return FAIL;
        }
        if (window3d_write_weight(fp, x, y, z, weight) != SUCCESS) {
            return FAIL;
        }
        *n_queries += 1;
    }

    if (status == 1 || status == 2 || status == 0) {
        BLEND_Report(BLEND_MSG_ERROR, "window3d: standard input must contain x y z query triples\n");
        return FAIL;
    }

    return SUCCESS;
}

int blend_window3d_module(int argc, char **argv)
{
    window3d_options options;
    window3d_support full_support;
    window3d_support_list supports = {0};
    window3d_support_list *support_list = NULL;
    int parse_status;
    int n_queries = 0;
    int status = SUCCESS;

    memset(&full_support, 0, sizeof(full_support));

    parse_status = window3d_parse_options(argc, argv, &options);
    if (parse_status == 2) {
        return SUCCESS;
    }
    if (parse_status != SUCCESS) {
        if (blend_get_verbosity() >= BLEND_MSG_ERROR) {
            window3d_usage(stderr);
        }
        return FAIL;
    }

    if (options.blendfile != NULL) {
        double start = blend_elapsed_seconds();

        if (window3d_read_blendfile(&options, &supports) != SUCCESS) {
            return FAIL;
        }
        support_list = &supports;
        BLEND_Report(BLEND_MSG_TIMING,
                     "window3d: prepared %zu blendfile supports in %.3f s\n",
                     supports.count, blend_elapsed_seconds() - start);
    }
    else if (window3d_make_full_support(&options, &full_support) != SUCCESS) {
        return FAIL;
    }

    if (!isatty(STDIN_FILENO)) {
        double start = blend_elapsed_seconds();

        status = window3d_write_queries(stdout, &options, &full_support, support_list, &n_queries);
        if (status == SUCCESS && n_queries > 0) {
            BLEND_Report(BLEND_MSG_TIMING,
                         "window3d: evaluated %d query points in %.3f s\n",
                         n_queries, blend_elapsed_seconds() - start);
        }
    }

    if (status == SUCCESS && n_queries == 0) {
        double start = blend_elapsed_seconds();
        size_t node_count = (size_t)options.nx * (size_t)options.ny * (size_t)options.nz;

        BLEND_Report(BLEND_MSG_TIMING,
                     "window3d: writing %zu grid nodes (%d x %d x %d)\n",
                     node_count, options.nx, options.ny, options.nz);
        status = window3d_write_grid(stdout, &options, &full_support, support_list);
        if (status == SUCCESS) {
            BLEND_Report(BLEND_MSG_TIMING,
                         "window3d: wrote %zu grid nodes in %.3f s\n",
                         node_count, blend_elapsed_seconds() - start);
        }
    }

    window3d_support_clear(&full_support);
    window3d_support_list_free(&supports);
    return status;
}
