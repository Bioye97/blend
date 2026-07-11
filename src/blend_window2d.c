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

typedef struct window2d_options {
    double xmin;
    double xmax;
    double ymin;
    double ymax;
    double dx;
    double dy;
    double ratio_x1;
    double ratio_x2;
    double ratio_y1;
    double ratio_y2;
    int nx;
    int ny;
    int has_region;
    int has_increment;
    int has_function;
    int has_taper;
    int has_monotone_method;
    blend_window_function x_function;
    blend_window_function y_function;
    const char *blendfile;
    char clobber;
    char monotone_method;
} window2d_options;

typedef struct window2d_support {
    double xmin;
    double xmax;
    double ymin;
    double ymax;
    size_t n_vertices;
    vertex *vertices;
    window data;
} window2d_support;

typedef struct window2d_support_list {
    window2d_support *items;
    size_t count;
    size_t capacity;
} window2d_support_list;

static void window2d_usage(FILE *fp)
{
    fprintf(fp, "blend window2d - Generate 2-D blending weights\n\n");
    fprintf(fp, "usage: blend window2d -R<xmin>/<xmax>/<ymin>/<ymax> -I<dx>[/<dy>]\n");
    fprintf(fp, "       [-F<xfunction>[/<yfunction>]] [-T<rx1>/<rx2>/<ry1>/<ry2>]\n");
    fprintf(fp, "       [-B<blendfile>] [-C<f|l|o|u|a|g|p>] [-M<method>] [-V[q|e|w|t|i|c|d]]\n\n");

    fprintf(fp, "Generate 2-D blending weights and write three columns: x y weight. If x/y\n");
    fprintf(fp, "query coordinates are provided on standard input, weights are returned at those\n");
    fprintf(fp, "locations using bilinear interpolation from neighboring grid-point weights. If\n");
    fprintf(fp, "no query points are provided, weights are written for the complete -R/-I grid.\n\n");

    fprintf(fp, "REQUIRED ARGUMENTS:\n\n");
    fprintf(fp, "  -R<xmin>/<xmax>/<ymin>/<ymax>, --region=<xmin>/<xmax>/<ymin>/<ymax>\n");
    fprintf(fp, "     Specify the full 2-D output domain in user coordinates. The domain must\n");
    fprintf(fp, "     satisfy xmin < xmax and ymin < ymax. Internally, each dimension is mapped\n");
    fprintf(fp, "     to zero-based grid indices, but output always uses user coordinates. If an\n");
    fprintf(fp, "     upper bound does not fall exactly on its increment, BLEND adjusts it upward\n");
    fprintf(fp, "     to the next increment and reports a warning.\n\n");

    fprintf(fp, "  -I<dx>[/<dy>], --increment=<dx>[/<dy>]\n");
    fprintf(fp, "     Specify grid increments in user coordinates. Both increments must be\n");
    fprintf(fp, "     positive. If dy is omitted then dy = dx.\n\n");

    fprintf(fp, "OPTIONAL ARGUMENTS:\n\n");
    fprintf(fp, "  -F<xfunction>[/<yfunction>], --function=<xfunction>[/<yfunction>]\n");
    fprintf(fp, "     Select the window functions used along x and y for a single support over\n");
    fprintf(fp, "     the full -R domain [Default is cosine/cosine]. If yfunction is omitted,\n");
    fprintf(fp, "     yfunction = xfunction. Run blend --show-windows to list available window\n");
    fprintf(fp, "     functions. This option is ignored when -B is given.\n\n");

    fprintf(fp, "  -T<rx1>/<rx2>/<ry1>/<ry2>, --taper_ratio=<rx1>/<rx2>/<ry1>/<ry2>\n");
    fprintf(fp, "     Set beginning and ending taper ratios for x and y. Ratios must be >= 0\n");
    fprintf(fp, "     and < 0.5 [Default is 0.2/0.2/0.2/0.2]. This option is ignored for boxcar\n");
    fprintf(fp, "     dimensions and when -B is given.\n\n");

    fprintf(fp, "  -B<blendfile>, --blendfile=<blendfile>\n");
    fprintf(fp, "     Read one or more xy-monotone polygon supports from <blendfile>. Each\n");
    fprintf(fp, "     non-empty, non-comment row must contain exactly three fields:\n");
    fprintf(fp, "       <polygonfile> <xfunction>/<yfunction> <rx1>/<rx2>/<ry1>/<ry2>\n");
    fprintf(fp, "     where <polygonfile> is a two-column x/y vertex file contained inside -R.\n");
    fprintf(fp, "     Lines may contain comments introduced by '#'. Grid points outside all\n");
    fprintf(fp, "     blendfile polygons receive weight 0.\n\n");

    fprintf(fp, "  -C<f|l|o|u|a|g|p>, --clobber=<mode>\n");
    fprintf(fp, "     Select how overlapping polygon weights from -B are combined [Default is p].\n");
    fprintf(fp, "       f  Use the first polygon listed in the blendfile.\n");
    fprintf(fp, "       l  Use the lowest weight among overlapping polygons.\n");
    fprintf(fp, "       o  Use the last polygon listed in the blendfile.\n");
    fprintf(fp, "       u  Use the highest weight among overlapping polygons.\n");
    fprintf(fp, "       a  Use the arithmetic average of overlapping weights.\n");
    fprintf(fp, "       g  Use the geometric average of overlapping weights.\n");
    fprintf(fp, "       p  Use the product of overlapping weights.\n\n");

    fprintf(fp, "  -M<method>, --monotone=<method>\n");
    fprintf(fp, "     Select how non-xy-monotone blendfile polygons are converted before window\n");
    fprintf(fp, "     assembly. If this option is omitted, non-xy-monotone polygons are rejected.\n");
    fprintf(fp, "       e  Use the xy-monotone envelope method.\n");
    fprintf(fp, "       b  Use the highest-IoU piecewise envelope from both traversal directions,\n");
    fprintf(fp, "          sampled on the -R/-I grid.\n\n");

    fprintf(fp, "  -V[level], --verbose=<level>\n");
    fprintf(fp, "     Select verbosity level [w]. Choose among q, e, w, t, i, c, and d.\n\n");

    fprintf(fp, "  -?, --help\n");
    fprintf(fp, "     Print this usage message and exit.\n\n");

    fprintf(fp, "EXAMPLES:\n\n");
    fprintf(fp, "  blend window2d -R0/10/0/5 -I0.5 -Fcosine/cosine -T0.2/0.2/0.2/0.2\n");
    fprintf(fp, "     Generate one 2-D cosine-tapered window over the full domain.\n\n");
    fprintf(fp, "  printf '2.5 1.5\\n' | blend window2d -R0/10/0/5 -I1 -Fcosine/cosine\n");
    fprintf(fp, "     Query the interpolated weight at x = 2.5, y = 1.5.\n\n");
    fprintf(fp, "  blend window2d -R0/10/0/10 -I0.25 -Bsupports.txt -Ca\n");
    fprintf(fp, "     Read multiple polygon supports from supports.txt and average overlaps.\n");
}

static int window2d_parse_region(const char *text, window2d_options *options)
{
    char *end = NULL;
    double values[4];
    int i;

    for (i = 0; i < 4; i++) {
        errno = 0;
        values[i] = strtod(text, &end);
        if (text == end || errno == ERANGE || !isfinite(values[i])) {
            BLEND_Report(BLEND_MSG_ERROR, "window2d: invalid region\n");
            return FAIL;
        }
        if (i < 3) {
            if (*end != '/') {
                BLEND_Report(BLEND_MSG_ERROR, "window2d: region must be xmin/xmax/ymin/ymax\n");
                return FAIL;
            }
            text = end + 1;
        }
        else if (*end != '\0') {
            BLEND_Report(BLEND_MSG_ERROR, "window2d: region must be xmin/xmax/ymin/ymax\n");
            return FAIL;
        }
    }

    if (values[0] >= values[1] || values[2] >= values[3]) {
        BLEND_Report(BLEND_MSG_ERROR, "window2d: region must satisfy xmin < xmax and ymin < ymax\n");
        return FAIL;
    }

    options->xmin = values[0];
    options->xmax = values[1];
    options->ymin = values[2];
    options->ymax = values[3];
    options->has_region = 1;
    return SUCCESS;
}

static int window2d_parse_increment(const char *text, window2d_options *options)
{
    char *end = NULL;
    double dx, dy;

    errno = 0;
    dx = strtod(text, &end);
    if (text == end || errno == ERANGE || !isfinite(dx)) {
        BLEND_Report(BLEND_MSG_ERROR, "window2d: invalid increment: %s\n", text);
        return FAIL;
    }

    if (*end == '\0') {
        dy = dx;
    }
    else if (*end == '/') {
        text = end + 1;
        errno = 0;
        dy = strtod(text, &end);
        if (text == end || *end != '\0' || errno == ERANGE || !isfinite(dy)) {
            BLEND_Report(BLEND_MSG_ERROR, "window2d: invalid increment: %s\n", text);
            return FAIL;
        }
    }
    else {
        BLEND_Report(BLEND_MSG_ERROR, "window2d: invalid increment: %s\n", text);
        return FAIL;
    }

    options->dx = dx;
    options->dy = dy;
    options->has_increment = 1;
    return SUCCESS;
}

static int window2d_parse_taper_ratio_values(const char *text, const char *context,
                                             double *rx1, double *rx2, double *ry1, double *ry2)
{
    char *end = NULL;
    double values[4];
    int i;

    for (i = 0; i < 4; i++) {
        errno = 0;
        values[i] = strtod(text, &end);
        if (text == end || errno == ERANGE || !isfinite(values[i])) {
            BLEND_Report(BLEND_MSG_ERROR, "window2d: invalid taper ratio in %s\n", context);
            return FAIL;
        }
        if (i < 3) {
            if (*end != '/') {
                BLEND_Report(BLEND_MSG_ERROR, "window2d: taper ratio in %s must be rx1/rx2/ry1/ry2\n", context);
                return FAIL;
            }
            text = end + 1;
        }
        else if (*end != '\0') {
            BLEND_Report(BLEND_MSG_ERROR, "window2d: taper ratio in %s must be rx1/rx2/ry1/ry2\n", context);
            return FAIL;
        }
    }

    for (i = 0; i < 4; i++) {
        if (values[i] < 0.0 || values[i] >= 0.5) {
            BLEND_Report(BLEND_MSG_ERROR, "window2d: taper ratios in %s must be >= 0 and < 0.5\n", context);
            return FAIL;
        }
    }

    *rx1 = values[0];
    *rx2 = values[1];
    *ry1 = values[2];
    *ry2 = values[3];
    return SUCCESS;
}

static int window2d_parse_taper_ratio(const char *text, window2d_options *options)
{
    return window2d_parse_taper_ratio_values(text, "option -T",
                                             &options->ratio_x1, &options->ratio_x2,
                                             &options->ratio_y1, &options->ratio_y2);
}

static int window2d_parse_function_pair(const char *text, blend_window_function *x_function,
                                        blend_window_function *y_function)
{
    char x_name[64];
    char y_name[64];
    const char *slash;
    size_t x_len, y_len;

    slash = strchr(text, '/');
    if (slash == NULL) {
        x_len = strlen(text);
        y_len = x_len;
        if (x_len == 0 || x_len >= sizeof(x_name)) {
            BLEND_Report(BLEND_MSG_ERROR, "window2d: invalid window function: %s\n", text);
            return FAIL;
        }
        memcpy(x_name, text, x_len + 1);
        memcpy(y_name, text, y_len + 1);
    }
    else {
        x_len = (size_t)(slash - text);
        y_len = strlen(slash + 1);
        if (x_len == 0 || y_len == 0 || x_len >= sizeof(x_name) || y_len >= sizeof(y_name)) {
            BLEND_Report(BLEND_MSG_ERROR, "window2d: invalid window function pair: %s\n", text);
            return FAIL;
        }
        memcpy(x_name, text, x_len);
        x_name[x_len] = '\0';
        memcpy(y_name, slash + 1, y_len + 1);
    }

    if (blend_window_function_from_name(x_name, x_function) != SUCCESS) {
        BLEND_Report(BLEND_MSG_ERROR, "window2d: unknown x window function: %s\n", x_name);
        return FAIL;
    }
    if (blend_window_function_from_name(y_name, y_function) != SUCCESS) {
        BLEND_Report(BLEND_MSG_ERROR, "window2d: unknown y window function: %s\n", y_name);
        return FAIL;
    }

    return SUCCESS;
}

static int window2d_parse_clobber(const char *value, window2d_options *options)
{
    if (value == NULL || value[0] == '\0' || value[1] != '\0') {
        BLEND_Report(BLEND_MSG_ERROR, "window2d: clobber mode must be one of f, l, o, u, a, g, p\n");
        return FAIL;
    }

    if (strchr("flouagp", value[0]) == NULL) {
        BLEND_Report(BLEND_MSG_ERROR, "window2d: unknown clobber mode: %s\n", value);
        return FAIL;
    }

    options->clobber = value[0];
    return SUCCESS;
}

static int window2d_parse_monotone_method(const char *value, window2d_options *options)
{
    if (value == NULL || value[0] == '\0' || value[1] != '\0') {
        BLEND_Report(BLEND_MSG_ERROR, "window2d: monotone method must be one of e, b\n");
        return FAIL;
    }

    if (value[0] != 'e' && value[0] != 'b') {
        BLEND_Report(BLEND_MSG_ERROR, "window2d: unknown monotone method: %s\n", value);
        return FAIL;
    }

    options->monotone_method = value[0];
    options->has_monotone_method = 1;
    return SUCCESS;
}

static const char *window2d_option_value(int argc, char **argv, int *i, const char *option)
{
    if (*i + 1 >= argc) {
        BLEND_Report(BLEND_MSG_ERROR, "window2d: option %s requires an argument\n", option);
        return NULL;
    }

    *i += 1;
    return argv[*i];
}

static int window2d_adjust_axis(const char *axis, double start, double *end, double increment, int *n)
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
            "window2d: adjusted %s region end from %.12g to %.12g to match increment %.12g\n",
            axis, old_end, *end, increment
        );
    }

    if (intervals < 1.0 || intervals > (double)(INT_MAX - 1)) {
        BLEND_Report(BLEND_MSG_ERROR, "window2d: region and increment produce an invalid %s grid size\n", axis);
        return FAIL;
    }

    *n = (int)intervals + 1;
    return SUCCESS;
}

static int window2d_parse_options(int argc, char **argv, window2d_options *options)
{
    int i;

    if (argc == 0) {
        window2d_usage(stdout);
        return 2;
    }

    memset(options, 0, sizeof(*options));
    options->ratio_x1 = 0.2;
    options->ratio_x2 = 0.2;
    options->ratio_y1 = 0.2;
    options->ratio_y2 = 0.2;
    options->x_function = WFUNC_COSINE;
    options->y_function = WFUNC_COSINE;
    options->clobber = 'p';

    for (i = 0; i < argc; i++) {
        const char *arg = argv[i];
        const char *value = NULL;

        if (strcmp(arg, "-?") == 0 || strcmp(arg, "--help") == 0) {
            window2d_usage(stdout);
            return 2;
        }
        else if (strcmp(arg, "-R") == 0 || strcmp(arg, "--region") == 0) {
            value = window2d_option_value(argc, argv, &i, arg);
            if (value == NULL || window2d_parse_region(value, options) != SUCCESS) return FAIL;
        }
        else if (strncmp(arg, "-R", 2) == 0 && arg[2] != '\0') {
            if (window2d_parse_region(arg + 2, options) != SUCCESS) return FAIL;
        }
        else if (strncmp(arg, "--region=", 9) == 0) {
            if (window2d_parse_region(arg + 9, options) != SUCCESS) return FAIL;
        }
        else if (strcmp(arg, "-I") == 0 || strcmp(arg, "--increment") == 0) {
            value = window2d_option_value(argc, argv, &i, arg);
            if (value == NULL || window2d_parse_increment(value, options) != SUCCESS) return FAIL;
        }
        else if (strncmp(arg, "-I", 2) == 0 && arg[2] != '\0') {
            if (window2d_parse_increment(arg + 2, options) != SUCCESS) return FAIL;
        }
        else if (strncmp(arg, "--increment=", 12) == 0) {
            if (window2d_parse_increment(arg + 12, options) != SUCCESS) return FAIL;
        }
        else if (strcmp(arg, "-F") == 0 || strcmp(arg, "--function") == 0) {
            value = window2d_option_value(argc, argv, &i, arg);
            if (value == NULL || window2d_parse_function_pair(value, &options->x_function, &options->y_function) != SUCCESS) return FAIL;
            options->has_function = 1;
        }
        else if (strncmp(arg, "-F", 2) == 0 && arg[2] != '\0') {
            if (window2d_parse_function_pair(arg + 2, &options->x_function, &options->y_function) != SUCCESS) return FAIL;
            options->has_function = 1;
        }
        else if (strncmp(arg, "--function=", 11) == 0) {
            if (window2d_parse_function_pair(arg + 11, &options->x_function, &options->y_function) != SUCCESS) return FAIL;
            options->has_function = 1;
        }
        else if (strcmp(arg, "-T") == 0 || strcmp(arg, "--taper_ratio") == 0) {
            value = window2d_option_value(argc, argv, &i, arg);
            if (value == NULL || window2d_parse_taper_ratio(value, options) != SUCCESS) return FAIL;
            options->has_taper = 1;
        }
        else if (strncmp(arg, "-T", 2) == 0 && arg[2] != '\0') {
            if (window2d_parse_taper_ratio(arg + 2, options) != SUCCESS) return FAIL;
            options->has_taper = 1;
        }
        else if (strncmp(arg, "--taper_ratio=", 14) == 0) {
            if (window2d_parse_taper_ratio(arg + 14, options) != SUCCESS) return FAIL;
            options->has_taper = 1;
        }
        else if (strcmp(arg, "-B") == 0 || strcmp(arg, "--blendfile") == 0) {
            value = window2d_option_value(argc, argv, &i, arg);
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
            value = window2d_option_value(argc, argv, &i, arg);
            if (value == NULL || window2d_parse_clobber(value, options) != SUCCESS) return FAIL;
        }
        else if (strncmp(arg, "-C", 2) == 0 && arg[2] != '\0') {
            if (window2d_parse_clobber(arg + 2, options) != SUCCESS) return FAIL;
        }
        else if (strncmp(arg, "--clobber=", 10) == 0) {
            if (window2d_parse_clobber(arg + 10, options) != SUCCESS) return FAIL;
        }
        else if (strcmp(arg, "-M") == 0 || strcmp(arg, "--monotone") == 0) {
            value = window2d_option_value(argc, argv, &i, arg);
            if (value == NULL || window2d_parse_monotone_method(value, options) != SUCCESS) return FAIL;
        }
        else if (strncmp(arg, "-M", 2) == 0 && arg[2] != '\0') {
            if (window2d_parse_monotone_method(arg + 2, options) != SUCCESS) return FAIL;
        }
        else if (strncmp(arg, "--monotone=", 11) == 0) {
            if (window2d_parse_monotone_method(arg + 11, options) != SUCCESS) return FAIL;
        }
        else {
            BLEND_Report(BLEND_MSG_ERROR, "window2d: unknown option: %s\n", arg);
            return FAIL;
        }
    }

    if (!options->has_region) {
        BLEND_Report(BLEND_MSG_ERROR, "window2d: -R/--region is required\n");
        return FAIL;
    }
    if (!options->has_increment) {
        BLEND_Report(BLEND_MSG_ERROR, "window2d: -I/--increment is required\n");
        return FAIL;
    }
    if (options->dx <= 0.0 || options->dy <= 0.0) {
        BLEND_Report(BLEND_MSG_ERROR, "window2d: increments must be positive\n");
        return FAIL;
    }

    if (options->blendfile != NULL) {
        if (options->has_function) {
            BLEND_Report(BLEND_MSG_WARNING, "window2d: -F/--function is ignored when -B/--blendfile is given\n");
        }
        if (options->has_taper) {
            BLEND_Report(BLEND_MSG_WARNING, "window2d: -T/--taper_ratio is ignored when -B/--blendfile is given\n");
        }
    }
    else if (options->has_taper && (options->x_function == WFUNC_BOXCAR || options->y_function == WFUNC_BOXCAR)) {
        BLEND_Report(BLEND_MSG_WARNING, "window2d: taper ratios are ignored for boxcar window dimensions\n");
    }
    if (options->blendfile == NULL && options->has_monotone_method) {
        BLEND_Report(BLEND_MSG_WARNING, "window2d: -M/--monotone is ignored unless -B/--blendfile is given\n");
    }

    if (window2d_adjust_axis("x", options->xmin, &options->xmax, options->dx, &options->nx) != SUCCESS) {
        return FAIL;
    }
    if (window2d_adjust_axis("y", options->ymin, &options->ymax, options->dy, &options->ny) != SUCCESS) {
        return FAIL;
    }

    return SUCCESS;
}

static char *window2d_trim_line(char *line)
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

static void window2d_support_clear(window2d_support *support)
{
    if (support == NULL) {
        return;
    }
    free(support->vertices);
    support->vertices = NULL;
    support->n_vertices = 0;
    blend_window_boundary_clear(&support->data);
}

static void window2d_support_list_free(window2d_support_list *supports)
{
    size_t i;

    for (i = 0; i < supports->count; i++) {
        window2d_support_clear(&supports->items[i]);
    }
    free(supports->items);
    supports->items = NULL;
    supports->count = 0;
    supports->capacity = 0;
}

static int window2d_support_list_append(window2d_support_list *supports, window2d_support *support)
{
    window2d_support *items;
    size_t capacity;

    if (supports->count == supports->capacity) {
        capacity = supports->capacity == 0 ? 8 : supports->capacity * 2;
        items = (window2d_support *)realloc(supports->items, capacity * sizeof(*supports->items));
        if (items == NULL) {
            BLEND_Report(BLEND_MSG_ERROR, "window2d: could not allocate polygon support storage\n");
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

static int window2d_support_grid_size(double min, double max, double increment, int *n)
{
    double span = max - min;
    double intervals = floor(span / increment + 0.5);

    if (intervals < 1.0 || intervals > (double)(INT_MAX - 1)) {
        BLEND_Report(BLEND_MSG_ERROR, "window2d: polygon support produces an invalid grid size\n");
        return FAIL;
    }

    *n = (int)intervals + 1;
    return SUCCESS;
}

static int window2d_make_polygon(vertex *real_vertices, size_t n_vertices, window2d_support *support,
                                 double *xmin, double *xmax, double *ymin, double *ymax)
{
    size_t i;
    polygon local_poly;

    memset(&local_poly, 0, sizeof(local_poly));
    if (blend_polygon_alloc(&local_poly, n_vertices) != SUCCESS) {
        BLEND_Report(BLEND_MSG_ERROR, "window2d: could not allocate local polygon\n");
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
        BLEND_Report(BLEND_MSG_ERROR, "window2d: could not assign polygon support\n");
        blend_polygon_free(&local_poly);
        return FAIL;
    }

    blend_polygon_free(&local_poly);
    return SUCCESS;
}

static int window2d_build_support(const window2d_options *options, const char *polygon_file,
                                  const vertex *vertices, size_t n_vertices,
                                  blend_window_function x_function, blend_window_function y_function,
                                  double rx1, double rx2, double ry1, double ry2,
                                  window2d_support *support)
{
    size_t i;
    permuted_vertex pv;

    memset(support, 0, sizeof(*support));
    memset(&pv, 0, sizeof(pv));

    support->vertices = (vertex *)calloc(n_vertices, sizeof(*support->vertices));
    if (support->vertices == NULL) {
        BLEND_Report(BLEND_MSG_ERROR, "window2d: could not allocate polygon support vertices\n");
        return FAIL;
    }
    memcpy(support->vertices, vertices, n_vertices * sizeof(*support->vertices));
    support->n_vertices = n_vertices;

    support->xmin = vertices[0].x;
    support->xmax = vertices[0].x;
    support->ymin = vertices[0].y;
    support->ymax = vertices[0].y;
    for (i = 0; i < n_vertices; i++) {
        if (vertices[i].x < options->xmin || vertices[i].x > options->xmax ||
            vertices[i].y < options->ymin || vertices[i].y > options->ymax) {
            BLEND_Report(BLEND_MSG_ERROR, "window2d: polygon %s has a vertex outside -R\n", polygon_file);
            window2d_support_clear(support);
            return FAIL;
        }
        if (vertices[i].x < support->xmin) support->xmin = vertices[i].x;
        if (vertices[i].x > support->xmax) support->xmax = vertices[i].x;
        if (vertices[i].y < support->ymin) support->ymin = vertices[i].y;
        if (vertices[i].y > support->ymax) support->ymax = vertices[i].y;
    }

    if (support->xmin >= support->xmax || support->ymin >= support->ymax) {
        BLEND_Report(BLEND_MSG_ERROR, "window2d: polygon %s has an invalid bounding box\n", polygon_file);
        window2d_support_clear(support);
        return FAIL;
    }

    if (window2d_support_grid_size(support->xmin, support->xmax, options->dx, &support->data.nx) != SUCCESS ||
        window2d_support_grid_size(support->ymin, support->ymax, options->dy, &support->data.ny) != SUCCESS) {
        window2d_support_clear(support);
        return FAIL;
    }

    support->data.ratio_x1 = rx1;
    support->data.ratio_x2 = rx2;
    support->data.ratio_y1 = ry1;
    support->data.ratio_y2 = ry2;
    support->data.x_function = x_function;
    support->data.y_function = y_function;

    if (window2d_make_polygon(support->vertices, support->n_vertices, support,
                              &support->xmin, &support->xmax, &support->ymin, &support->ymax) != SUCCESS) {
        window2d_support_clear(support);
        return FAIL;
    }

    if (boundary_assembly(&support->data, &pv) != SUCCESS) {
        BLEND_Report(BLEND_MSG_ERROR, "window2d: could not assemble boundary for polygon %s\n", polygon_file);
        blend_permuted_vertex_free(&pv);
        window2d_support_clear(support);
        return FAIL;
    }
    blend_permuted_vertex_free(&pv);

    return SUCCESS;
}

static int window2d_parse_blendfile_line(const char *line, long line_number,
                                         const window2d_options *options,
                                         window2d_support *support)
{
    char polygon_file[PATH_MAX];
    char function_pair[128];
    char taper_ratio[128];
    char extra[2];
    polygon input_polygon = {0};
    polygon support_polygon = {0};
    blend_window_function x_function, y_function;
    double rx1, rx2, ry1, ry2;
    int is_xy_monotone = 0;
    int fields;
    int status;

    fields = sscanf(line, "%1023s %127s %127s %1s", polygon_file, function_pair, taper_ratio, extra);
    if (fields != 3) {
        BLEND_Report(BLEND_MSG_ERROR, "window2d: blendfile line %ld must have three fields: polygon function_pair taper_ratio\n",
                     line_number);
        return FAIL;
    }

    if (window2d_parse_function_pair(function_pair, &x_function, &y_function) != SUCCESS) {
        BLEND_Report(BLEND_MSG_ERROR, "window2d: invalid window function pair on blendfile line %ld\n", line_number);
        return FAIL;
    }
    if (window2d_parse_taper_ratio_values(taper_ratio, "blendfile", &rx1, &rx2, &ry1, &ry2) != SUCCESS) {
        BLEND_Report(BLEND_MSG_ERROR, "window2d: invalid taper ratios on blendfile line %ld\n", line_number);
        return FAIL;
    }
    if ((x_function == WFUNC_BOXCAR || y_function == WFUNC_BOXCAR) &&
        (rx1 != 0.0 || rx2 != 0.0 || ry1 != 0.0 || ry2 != 0.0)) {
        BLEND_Report(BLEND_MSG_WARNING, "window2d: taper ratios on blendfile line %ld are ignored for boxcar window dimensions\n",
                     line_number);
    }

    if (blend_polygon_read(polygon_file, &input_polygon) != SUCCESS) {
        return FAIL;
    }

    if (blend_polygon_is_xy_monotone(&input_polygon, &is_xy_monotone) != SUCCESS) {
        blend_polygon_free(&input_polygon);
        return FAIL;
    }
    if (!is_xy_monotone) {
        if (options->monotone_method != 'e' && options->monotone_method != 'b') {
            BLEND_Report(BLEND_MSG_ERROR,
                         "window2d: polygon %s is not xy-monotone; check -M if you want BLEND to refine it with -Me or -Mb\n",
                         polygon_file);
            blend_polygon_free(&input_polygon);
            return FAIL;
        }
        if (options->monotone_method == 'e') {
            BLEND_Report(BLEND_MSG_WARNING,
                         "window2d: polygon %s is not xy-monotone; using its xy-monotone envelope from -Me\n",
                         polygon_file);
        }
        else if (options->monotone_method == 'b') {
            BLEND_Report(BLEND_MSG_WARNING,
                         "window2d: polygon %s is not xy-monotone; using its best IoU piecewise xy-monotone envelope from -Mb\n",
                         polygon_file);
        }
    }

    if ((options->monotone_method == 'b'
             ? blend_polygon_xy_monotone_best_piecewise_envelope(
                   &input_polygon, &support_polygon,
                   options->xmin, options->xmax, options->ymin, options->ymax,
                   options->nx, options->ny)
             : blend_polygon_xy_monotone_envelope(&input_polygon, &support_polygon)) != SUCCESS) {
        blend_polygon_free(&input_polygon);
        return FAIL;
    }

    status = window2d_build_support(options, polygon_file, support_polygon.vertices, support_polygon.n_vertices,
                                    x_function, y_function, rx1, rx2, ry1, ry2, support);
    blend_polygon_free(&input_polygon);
    blend_polygon_free(&support_polygon);
    return status;
}

static int window2d_read_blendfile(const window2d_options *options, window2d_support_list *supports)
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
        window2d_support support;

        line_number++;
        comment = strchr(line, '#');
        if (comment != NULL) {
            *comment = '\0';
        }
        text = window2d_trim_line(line);
        if (*text == '\0') {
            continue;
        }

        if (window2d_parse_blendfile_line(text, line_number, options, &support) != SUCCESS) {
            fclose(fp);
            window2d_support_list_free(supports);
            return FAIL;
        }
        if (window2d_support_list_append(supports, &support) != SUCCESS) {
            window2d_support_clear(&support);
            fclose(fp);
            window2d_support_list_free(supports);
            return FAIL;
        }
    }

    if (ferror(fp)) {
        BLEND_Report(BLEND_MSG_ERROR, "%s: %s\n", options->blendfile, strerror(errno));
        fclose(fp);
        window2d_support_list_free(supports);
        return FAIL;
    }
    fclose(fp);

    if (supports->count == 0) {
        BLEND_Report(BLEND_MSG_ERROR, "window2d: blendfile contains no polygon supports\n");
        return FAIL;
    }

    return SUCCESS;
}

static int window2d_make_full_support(const window2d_options *options, window2d_support *support)
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

    return window2d_build_support(options, "full -R domain", vertices, 4,
                                  options->x_function, options->y_function,
                                  options->ratio_x1, options->ratio_x2,
                                  options->ratio_y1, options->ratio_y2,
                                  support);
}

static double window2d_coordinate_x(const window2d_options *options, int index)
{
    return options->xmin + (double)index * options->dx;
}

static double window2d_coordinate_y(const window2d_options *options, int index)
{
    return options->ymin + (double)index * options->dy;
}

static int window2d_point_on_segment(double x, double y, const vertex *a, const vertex *b)
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

static int window2d_point_in_polygon(double x, double y, const window2d_support *support)
{
    int inside = 0;
    size_t i, j;

    for (i = 0, j = support->n_vertices - 1; i < support->n_vertices; j = i++) {
        const vertex *vi = &support->vertices[i];
        const vertex *vj = &support->vertices[j];

        if (window2d_point_on_segment(x, y, vj, vi)) {
            return 1;
        }
        if (((vi->y > y) != (vj->y > y)) &&
            (x < (vj->x - vi->x) * (y - vi->y) / (vj->y - vi->y) + vi->x)) {
            inside = !inside;
        }
    }

    return inside;
}

static int window2d_support_grid_weight(int ix, int iy, window2d_support *support, double *weight)
{
    if (ix < 0 || iy < 0 || ix >= support->data.nx || iy >= support->data.ny) {
        *weight = 0.0;
        return SUCCESS;
    }

    if (embedding_contribution2d(ix, iy, &support->data) != SUCCESS) {
        return FAIL;
    }

    *weight = support->data.contribution;
    return SUCCESS;
}

static int window2d_support_weight(double x, double y, window2d_support *support, double *weight)
{
    double sx, sy;
    double x0, x1, y0, y1;
    double w00, w10, w01, w11;
    int ix0, iy0, ix1, iy1;

    if (x < support->xmin || x > support->xmax || y < support->ymin || y > support->ymax ||
        !window2d_point_in_polygon(x, y, support)) {
        *weight = 0.0;
        return SUCCESS;
    }

    sx = (x - support->xmin) * (double)(support->data.nx - 1) / (support->xmax - support->xmin);
    sy = (y - support->ymin) * (double)(support->data.ny - 1) / (support->ymax - support->ymin);
    ix0 = (int)floor(sx);
    iy0 = (int)floor(sy);

    if (ix0 < 0) ix0 = 0;
    if (iy0 < 0) iy0 = 0;
    if (ix0 >= support->data.nx - 1) ix0 = support->data.nx - 1;
    if (iy0 >= support->data.ny - 1) iy0 = support->data.ny - 1;

    ix1 = ix0 < support->data.nx - 1 ? ix0 + 1 : ix0;
    iy1 = iy0 < support->data.ny - 1 ? iy0 + 1 : iy0;

    if (window2d_support_grid_weight(ix0, iy0, support, &w00) != SUCCESS ||
        window2d_support_grid_weight(ix1, iy0, support, &w10) != SUCCESS ||
        window2d_support_grid_weight(ix0, iy1, support, &w01) != SUCCESS ||
        window2d_support_grid_weight(ix1, iy1, support, &w11) != SUCCESS) {
        return FAIL;
    }

    if (ix0 == ix1 && iy0 == iy1) {
        *weight = w00;
        return SUCCESS;
    }
    if (ix0 == ix1) {
        return interpolate_linear((double)iy0, w00, (double)iy1, w01, sy, weight);
    }
    if (iy0 == iy1) {
        return interpolate_linear((double)ix0, w00, (double)ix1, w10, sx, weight);
    }

    x0 = (double)ix0;
    x1 = (double)ix1;
    y0 = (double)iy0;
    y1 = (double)iy1;
    return interpolate_bilinear(x0, x1, y0, y1, w00, w10, w01, w11, sx, sy, weight);
}

static int window2d_combine_weights(const double *weights, size_t count, char clobber, double *weight)
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
            BLEND_Report(BLEND_MSG_ERROR, "window2d: unknown clobber mode: %c\n", clobber);
            return FAIL;
    }

    return SUCCESS;
}

static int window2d_blendfile_weight(double x, double y, const window2d_options *options,
                                     window2d_support_list *supports, double *weight)
{
    double *weights = NULL;
    size_t i, count = 0;
    int status;

    weights = (double *)calloc(supports->count, sizeof(*weights));
    if (weights == NULL) {
        BLEND_Report(BLEND_MSG_ERROR, "window2d: could not allocate overlap weights\n");
        return FAIL;
    }

    for (i = 0; i < supports->count; i++) {
        double support_weight;

        if (window2d_support_weight(x, y, &supports->items[i], &support_weight) != SUCCESS) {
            free(weights);
            return FAIL;
        }
        if (support_weight > 0.0 || window2d_point_in_polygon(x, y, &supports->items[i])) {
            weights[count] = support_weight;
            count++;
        }
    }

    status = window2d_combine_weights(weights, count, options->clobber, weight);
    free(weights);
    return status;
}

static int window2d_get_weight(double x, double y, int ix, int iy,
                               const window2d_options *options,
                               window2d_support *full_support,
                               window2d_support_list *supports,
                               double *weight)
{
    if (supports != NULL) {
        return window2d_blendfile_weight(x, y, options, supports, weight);
    }

    if (ix >= 0 && iy >= 0) {
        return window2d_support_grid_weight(ix, iy, full_support, weight);
    }

    return window2d_support_weight(x, y, full_support, weight);
}

static int window2d_write_weight(FILE *fp, double x, double y, double weight)
{
    return fprintf(fp, "%.12g %.12g %.12f\n", x, y, weight) < 0 ? FAIL : SUCCESS;
}

static int window2d_write_grid(FILE *fp, const window2d_options *options,
                               window2d_support *full_support,
                               window2d_support_list *supports)
{
    int ix, iy;

    for (iy = 0; iy < options->ny; iy++) {
        double y = window2d_coordinate_y(options, iy);
        for (ix = 0; ix < options->nx; ix++) {
            double x = window2d_coordinate_x(options, ix);
            double weight;

            if (window2d_get_weight(x, y, ix, iy, options, full_support, supports, &weight) != SUCCESS) {
                return FAIL;
            }
            if (window2d_write_weight(fp, x, y, weight) != SUCCESS) {
                return FAIL;
            }
        }
    }

    return SUCCESS;
}

static int window2d_write_queries(FILE *fp, const window2d_options *options,
                                  window2d_support *full_support,
                                  window2d_support_list *supports,
                                  int *n_queries)
{
    double x, y;
    int status;

    *n_queries = 0;
    while ((status = scanf("%lf %lf", &x, &y)) == 2) {
        double weight;

        if (!isfinite(x) || !isfinite(y)) {
            BLEND_Report(BLEND_MSG_ERROR, "window2d: standard input contains a non-finite query\n");
            return FAIL;
        }
        if (window2d_get_weight(x, y, -1, -1, options, full_support, supports, &weight) != SUCCESS) {
            return FAIL;
        }
        if (window2d_write_weight(fp, x, y, weight) != SUCCESS) {
            return FAIL;
        }
        *n_queries += 1;
    }

    if (status == 1 || status == 0) {
        BLEND_Report(BLEND_MSG_ERROR, "window2d: standard input must contain x y query pairs\n");
        return FAIL;
    }

    return SUCCESS;
}

int blend_window2d_module(int argc, char **argv)
{
    window2d_options options;
    window2d_support full_support;
    window2d_support_list supports = {0};
    window2d_support_list *support_list = NULL;
    int parse_status;
    int n_queries = 0;
    int status = SUCCESS;

    memset(&full_support, 0, sizeof(full_support));

    parse_status = window2d_parse_options(argc, argv, &options);
    if (parse_status == 2) {
        return SUCCESS;
    }
    if (parse_status != SUCCESS) {
        if (blend_get_verbosity() >= BLEND_MSG_ERROR) {
            window2d_usage(stderr);
        }
        return FAIL;
    }

    if (options.blendfile != NULL) {
        if (window2d_read_blendfile(&options, &supports) != SUCCESS) {
            return FAIL;
        }
        support_list = &supports;
    }
    else if (window2d_make_full_support(&options, &full_support) != SUCCESS) {
        return FAIL;
    }

    if (!isatty(STDIN_FILENO)) {
        status = window2d_write_queries(stdout, &options, &full_support, support_list, &n_queries);
    }

    if (status == SUCCESS && n_queries == 0) {
        status = window2d_write_grid(stdout, &options, &full_support, support_list);
    }

    window2d_support_clear(&full_support);
    window2d_support_list_free(&supports);
    return status;
}
