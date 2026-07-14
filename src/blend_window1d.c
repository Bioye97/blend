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

typedef struct window1d_support {
    double xmin;
    double xmax;
    double ratio1;
    double ratio2;
    int nx;
    blend_window_function function;
} window1d_support;

typedef struct window1d_support_list {
    window1d_support *items;
    size_t count;
    size_t capacity;
} window1d_support_list;

typedef struct window1d_options {
    double xmin;
    double xmax;
    double dx;
    double ratio1;
    double ratio2;
    int nx;
    int has_region;
    int has_increment;
    int has_function;
    int has_taper;
    blend_window_function function;
    const char *blendfile;
    char clobber;
} window1d_options;

static void window1d_usage(FILE *fp)
{
    fprintf(fp, "blend window1d - Generate 1-D blending weights\n\n");
    fprintf(fp, "usage: blend window1d -R<xmin>/<xmax> -I<dx> [-F<function>] [-T<r1>[/<r2>]]\n");
    fprintf(fp, "       [-B<blendfile>] [-C<f|l|o|u|a|g|p>] [-V[q|e|w|t|i|c|d]]\n");
    fprintf(fp, "\n");

    fprintf(fp, "Generate 1-D blending weights and write two columns: x weight. If x coordinates\n");
    fprintf(fp, "are provided on standard input, weights are returned at those query locations.\n");
    fprintf(fp, "Query locations may be arbitrary real coordinates; off-grid weights are computed\n");
    fprintf(fp, "by linear interpolation from neighboring grid-point weights. If no query points\n");
    fprintf(fp, "are provided, weights are written for the complete -R/-I grid.\n\n");

    fprintf(fp, "REQUIRED ARGUMENTS:\n\n");
    fprintf(fp, "  -R<xmin>/<xmax>, --region=<xmin>/<xmax>\n");
    fprintf(fp, "     Specify the full 1-D output domain in user coordinates. The domain must\n");
    fprintf(fp, "     satisfy xmin < xmax. If xmax does not fall exactly on the -I increment,\n");
    fprintf(fp, "     BLEND adjusts xmax upward to the next increment and reports a warning\n");
    fprintf(fp, "     at the default verbosity level.\n\n");

    fprintf(fp, "  -I<dx>, --increment=<dx>\n");
    fprintf(fp, "     Specify the grid increment in user coordinates. The increment must be\n");
    fprintf(fp, "     positive. For example, -R10/14 -I0.1 creates real coordinates\n");
    fprintf(fp, "     10, 10.1, 10.2, ..., 14.\n\n");

    fprintf(fp, "OPTIONAL ARGUMENTS:\n\n");
    fprintf(fp, "  -F<function>, --function=<function>\n");
    fprintf(fp, "     Select the window function used for a single support over the full -R\n");
    fprintf(fp, "     domain [Default is cosine]. Run blend --show-windows to list available\n");
    fprintf(fp, "     window functions.\n");
    fprintf(fp, "     This option is ignored when -B is given, since each blendfile row supplies\n");
    fprintf(fp, "     its own function.\n\n");

    fprintf(fp, "  -T<r1>[/<r2>], --taper_ratio=<r1>[/<r2>]\n");
    fprintf(fp, "     Set beginning and ending taper ratios for a single support over the full\n");
    fprintf(fp, "     -R domain. Ratios must be >= 0 and < 0.5. If r2 is omitted then r2 = r1.\n");
    fprintf(fp, "     [Default is 0.2/0.2]. This option is ignored for boxcar windows and when\n");
    fprintf(fp, "     -B is given, since each blendfile row supplies its own taper ratios.\n\n");

    fprintf(fp, "  -B<blendfile>, --blendfile=<blendfile>\n");
    fprintf(fp, "     Read one or more interval supports from <blendfile>. Each non-empty,\n");
    fprintf(fp, "     non-comment row must contain exactly four fields:\n");
    fprintf(fp, "       <left> <right> <function> <r1>[/<r2>]\n");
    fprintf(fp, "     where <left>/<right> define an interval contained inside the -R domain,\n");
    fprintf(fp, "     <function> is any supported window function listed by blend --show-windows,\n");
    fprintf(fp, "     and <r1>/<r2> are taper ratios for that interval. Lines may contain\n");
    fprintf(fp, "     comments introduced by '#'. Taper ratios are ignored for boxcar windows.\n");
    fprintf(fp, "     Example row for -R0/10:\n");
    fprintf(fp, "       2 3 cosine 0.2/0.2\n");
    fprintf(fp, "     Grid points not covered by any blendfile interval within -R receive weight 0.\n\n");

    fprintf(fp, "  -C<f|l|o|u|a|g|p>, --clobber=<mode>\n");
    fprintf(fp, "     Select how overlapping interval weights from -B are combined [Default is p].\n");
    fprintf(fp, "       f  Use the first interval listed in the blendfile.\n");
    fprintf(fp, "       l  Use the lowest weight among overlapping intervals.\n");
    fprintf(fp, "       o  Use the last interval listed in the blendfile.\n");
    fprintf(fp, "       u  Use the highest weight among overlapping intervals.\n");
    fprintf(fp, "       a  Use the arithmetic average of overlapping weights.\n");
    fprintf(fp, "       g  Use the geometric average of overlapping weights.\n");
    fprintf(fp, "       p  Use the product of overlapping weights.\n\n");

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
    fprintf(fp, "  blend window1d -R0/10 -I1 -Fcosine -T0.2/0.2\n");
    fprintf(fp, "     Generate one cosine-tapered window over the full domain.\n\n");
    fprintf(fp, "  printf '2.5\\n' | blend window1d -R0/10 -I1 -Fcosine -T0.2\n");
    fprintf(fp, "     Query the interpolated weight at x = 2.5.\n\n");
    fprintf(fp, "  blend window1d -R0/10 -I0.5 -Bsupports.txt -Ca\n");
    fprintf(fp, "     Read multiple interval supports from supports.txt and average overlaps.\n");
    fprintf(fp, "     Use shell redirection to write the output to a file, e.g., > weights.txt.\n");
}

static int window1d_parse_double(const char *text, const char *name, double *value)
{
    char *end = NULL;
    double parsed;

    errno = 0;
    parsed = strtod(text, &end);
    if (text == end || *end != '\0' || errno == ERANGE || !isfinite(parsed)) {
        BLEND_Report(BLEND_MSG_ERROR, "window1d: invalid %s: %s\n", name, text);
        return FAIL;
    }

    *value = parsed;
    return SUCCESS;
}

static int window1d_parse_region(const char *text, window1d_options *options)
{
    char *end = NULL;
    double xmin, xmax;

    errno = 0;
    xmin = strtod(text, &end);
    if (text == end || *end != '/' || errno == ERANGE || !isfinite(xmin)) {
        BLEND_Report(BLEND_MSG_ERROR, "window1d: invalid region: %s\n", text);
        return FAIL;
    }

    text = end + 1;
    errno = 0;
    xmax = strtod(text, &end);
    if (text == end || *end != '\0' || errno == ERANGE || !isfinite(xmax)) {
        BLEND_Report(BLEND_MSG_ERROR, "window1d: invalid region: %s\n", text);
        return FAIL;
    }

    if (xmin >= xmax) {
        BLEND_Report(BLEND_MSG_ERROR, "window1d: region must satisfy xmin < xmax\n");
        return FAIL;
    }

    options->xmin = xmin;
    options->xmax = xmax;
    options->has_region = 1;
    return SUCCESS;
}

static int window1d_parse_taper_ratio_values(const char *text, const char *context, double *ratio1, double *ratio2)
{
    char *end = NULL;
    double r1, r2;

    errno = 0;
    r1 = strtod(text, &end);
    if (text == end || errno == ERANGE || !isfinite(r1)) {
        BLEND_Report(BLEND_MSG_ERROR, "window1d: invalid taper ratio in %s: %s\n", context, text);
        return FAIL;
    }

    if (*end == '\0') {
        r2 = r1;
    }
    else if (*end == '/') {
        text = end + 1;
        errno = 0;
        r2 = strtod(text, &end);
        if (text == end || *end != '\0' || errno == ERANGE || !isfinite(r2)) {
            BLEND_Report(BLEND_MSG_ERROR, "window1d: invalid taper ratio in %s: %s\n", context, text);
            return FAIL;
        }
    }
    else {
        BLEND_Report(BLEND_MSG_ERROR, "window1d: invalid taper ratio in %s: %s\n", context, text);
        return FAIL;
    }

    if (r1 < 0.0 || r2 < 0.0 || r1 >= 0.5 || r2 >= 0.5) {
        BLEND_Report(BLEND_MSG_ERROR, "window1d: taper ratios in %s must be >= 0 and < 0.5\n", context);
        return FAIL;
    }

    *ratio1 = r1;
    *ratio2 = r2;
    return SUCCESS;
}

static int window1d_parse_taper_ratio(const char *text, window1d_options *options)
{
    return window1d_parse_taper_ratio_values(text, "option -T", &options->ratio1, &options->ratio2);
}

static int window1d_parse_function(blend_window_function *function, const char *value)
{
    if (blend_window_function_from_name(value, function) != SUCCESS) {
        BLEND_Report(BLEND_MSG_ERROR, "window1d: unknown window function: %s\n", value);
        return FAIL;
    }

    return SUCCESS;
}

static int window1d_parse_clobber(const char *value, window1d_options *options)
{
    if (value == NULL || value[0] == '\0' || value[1] != '\0') {
        BLEND_Report(BLEND_MSG_ERROR, "window1d: clobber mode must be one of f, l, o, u, a, g, p\n");
        return FAIL;
    }

    if (strchr("flouagp", value[0]) == NULL) {
        BLEND_Report(BLEND_MSG_ERROR, "window1d: unknown clobber mode: %s\n", value);
        return FAIL;
    }

    options->clobber = value[0];
    return SUCCESS;
}

static const char *window1d_option_value(int argc, char **argv, int *i, const char *option)
{
    if (*i + 1 >= argc) {
        BLEND_Report(BLEND_MSG_ERROR, "window1d: option %s requires an argument\n", option);
        return NULL;
    }

    *i += 1;
    return argv[*i];
}

static int window1d_parse_options(int argc, char **argv, window1d_options *options)
{
    int i;

    if (argc == 0) {
        window1d_usage(stdout);
        return 2;
    }

    options->xmin = 0.0;
    options->xmax = 0.0;
    options->dx = 0.0;
    options->ratio1 = 0.2;
    options->ratio2 = 0.2;
    options->nx = 0;
    options->has_region = 0;
    options->has_increment = 0;
    options->has_function = 0;
    options->has_taper = 0;
    options->blendfile = NULL;
    options->clobber = 'p';
    options->function = WFUNC_COSINE;

    for (i = 0; i < argc; i++) {
        const char *arg = argv[i];
        const char *value = NULL;

        if (strcmp(arg, "-?") == 0 || strcmp(arg, "--help") == 0) {
            window1d_usage(stdout);
            return 2;
        }
        else if (strcmp(arg, "-R") == 0 || strcmp(arg, "--region") == 0) {
            value = window1d_option_value(argc, argv, &i, arg);
            if (value == NULL || window1d_parse_region(value, options) != SUCCESS) return FAIL;
        }
        else if (strncmp(arg, "-R", 2) == 0 && arg[2] != '\0') {
            if (window1d_parse_region(arg + 2, options) != SUCCESS) return FAIL;
        }
        else if (strncmp(arg, "--region=", 9) == 0) {
            if (window1d_parse_region(arg + 9, options) != SUCCESS) return FAIL;
        }
        else if (strcmp(arg, "-I") == 0 || strcmp(arg, "--increment") == 0) {
            value = window1d_option_value(argc, argv, &i, arg);
            if (value == NULL || window1d_parse_double(value, "increment", &options->dx) != SUCCESS) return FAIL;
            options->has_increment = 1;
        }
        else if (strncmp(arg, "-I", 2) == 0 && arg[2] != '\0') {
            if (window1d_parse_double(arg + 2, "increment", &options->dx) != SUCCESS) return FAIL;
            options->has_increment = 1;
        }
        else if (strncmp(arg, "--increment=", 12) == 0) {
            if (window1d_parse_double(arg + 12, "increment", &options->dx) != SUCCESS) return FAIL;
            options->has_increment = 1;
        }
        else if (strcmp(arg, "-F") == 0 || strcmp(arg, "--function") == 0) {
            value = window1d_option_value(argc, argv, &i, arg);
            if (value == NULL || window1d_parse_function(&options->function, value) != SUCCESS) return FAIL;
            options->has_function = 1;
        }
        else if (strncmp(arg, "-F", 2) == 0 && arg[2] != '\0') {
            if (window1d_parse_function(&options->function, arg + 2) != SUCCESS) return FAIL;
            options->has_function = 1;
        }
        else if (strncmp(arg, "--function=", 11) == 0) {
            if (window1d_parse_function(&options->function, arg + 11) != SUCCESS) return FAIL;
            options->has_function = 1;
        }
        else if (strcmp(arg, "-T") == 0 || strcmp(arg, "--taper_ratio") == 0) {
            value = window1d_option_value(argc, argv, &i, arg);
            if (value == NULL || window1d_parse_taper_ratio(value, options) != SUCCESS) return FAIL;
            options->has_taper = 1;
        }
        else if (strncmp(arg, "-T", 2) == 0 && arg[2] != '\0') {
            if (window1d_parse_taper_ratio(arg + 2, options) != SUCCESS) return FAIL;
            options->has_taper = 1;
        }
        else if (strncmp(arg, "--taper_ratio=", 14) == 0) {
            if (window1d_parse_taper_ratio(arg + 14, options) != SUCCESS) return FAIL;
            options->has_taper = 1;
        }
        else if (strcmp(arg, "-B") == 0 || strcmp(arg, "--blendfile") == 0) {
            value = window1d_option_value(argc, argv, &i, arg);
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
            value = window1d_option_value(argc, argv, &i, arg);
            if (value == NULL || window1d_parse_clobber(value, options) != SUCCESS) return FAIL;
        }
        else if (strncmp(arg, "-C", 2) == 0 && arg[2] != '\0') {
            if (window1d_parse_clobber(arg + 2, options) != SUCCESS) return FAIL;
        }
        else if (strncmp(arg, "--clobber=", 10) == 0) {
            if (window1d_parse_clobber(arg + 10, options) != SUCCESS) return FAIL;
        }
        else {
            BLEND_Report(BLEND_MSG_ERROR, "window1d: unknown option: %s\n", arg);
            return FAIL;
        }
    }

    if (!options->has_region) {
        BLEND_Report(BLEND_MSG_ERROR, "window1d: -R/--region is required\n");
        return FAIL;
    }

    if (!options->has_increment) {
        BLEND_Report(BLEND_MSG_ERROR, "window1d: -I/--increment is required\n");
        return FAIL;
    }

    if (options->dx <= 0.0) {
        BLEND_Report(BLEND_MSG_ERROR, "window1d: increment must be positive\n");
        return FAIL;
    }

    if (options->ratio1 < 0.0 || options->ratio2 < 0.0 ||
        options->ratio1 >= 0.5 || options->ratio2 >= 0.5) {
        BLEND_Report(BLEND_MSG_ERROR, "window1d: taper ratios must be >= 0 and < 0.5\n");
        return FAIL;
    }

    if (options->blendfile != NULL) {
        if (options->has_function) {
            BLEND_Report(BLEND_MSG_WARNING, "window1d: -F/--function is ignored when -B/--blendfile is given\n");
        }
        if (options->has_taper) {
            BLEND_Report(BLEND_MSG_WARNING, "window1d: -T/--taper_ratio is ignored when -B/--blendfile is given\n");
        }
    }
    else if (options->function == WFUNC_BOXCAR && options->has_taper) {
        BLEND_Report(BLEND_MSG_WARNING, "window1d: -T/--taper_ratio is ignored for boxcar windows\n");
    }

    {
        double old_xmax = options->xmax;
        double span = options->xmax - options->xmin;
        double intervals_exact = span / options->dx;
        double intervals_nearest = floor(intervals_exact + 0.5);
        double tolerance = 1.0e-9 * fmax(1.0, fabs(intervals_exact));
        double intervals;

        if (fabs(intervals_exact - intervals_nearest) <= tolerance) {
            intervals = intervals_nearest;
            options->xmax = options->xmin + intervals * options->dx;
        }
        else {
            intervals = ceil(intervals_exact - tolerance);
            options->xmax = options->xmin + intervals * options->dx;
            BLEND_Report(
                BLEND_MSG_WARNING,
                "window1d: adjusted region end from %.12g to %.12g to match increment %.12g\n",
                old_xmax, options->xmax, options->dx
            );
        }

        if (intervals < 1.0 || intervals > (double)(INT_MAX - 1)) {
            BLEND_Report(BLEND_MSG_ERROR, "window1d: region and increment produce an invalid grid size\n");
            return FAIL;
        }
        options->nx = (int)intervals + 1;
    }

    return SUCCESS;
}

static void window1d_support_list_free(window1d_support_list *supports)
{
    free(supports->items);
    supports->items = NULL;
    supports->count = 0;
    supports->capacity = 0;
}

static int window1d_support_list_append(window1d_support_list *supports, const window1d_support *support)
{
    window1d_support *items;
    size_t capacity;

    if (supports->count == supports->capacity) {
        capacity = supports->capacity == 0 ? 8 : supports->capacity * 2;
        items = (window1d_support *)realloc(supports->items, capacity * sizeof(*supports->items));
        if (items == NULL) {
            BLEND_Report(BLEND_MSG_ERROR, "window1d: could not allocate blendfile interval storage\n");
            return FAIL;
        }
        supports->items = items;
        supports->capacity = capacity;
    }

    supports->items[supports->count] = *support;
    supports->count++;
    return SUCCESS;
}

static char *window1d_trim_line(char *line)
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

static int window1d_support_grid_size(const window1d_options *options, const window1d_support *support, int *nx)
{
    double span = support->xmax - support->xmin;
    double intervals = floor(span / options->dx + 0.5);

    if (intervals < 1.0 || intervals > (double)(INT_MAX - 1)) {
        BLEND_Report(BLEND_MSG_ERROR, "window1d: blendfile interval %.12g/%.12g produces an invalid grid size\n",
                     support->xmin, support->xmax);
        return FAIL;
    }

    *nx = (int)intervals + 1;
    return SUCCESS;
}

static int window1d_parse_blendfile_line(const char *line, long line_number,
                                         const window1d_options *options,
                                         window1d_support *support)
{
    char function_name[64];
    char taper_ratio[64];
    char extra[2];
    int fields;

    fields = sscanf(line, "%lf %lf %63s %63s %1s",
                    &support->xmin, &support->xmax, function_name, taper_ratio, extra);
    if (fields != 4) {
        BLEND_Report(BLEND_MSG_ERROR, "window1d: blendfile line %ld must have four fields: xmin xmax function taper_ratio\n",
                     line_number);
        return FAIL;
    }

    if (!isfinite(support->xmin) || !isfinite(support->xmax) || support->xmin >= support->xmax) {
        BLEND_Report(BLEND_MSG_ERROR, "window1d: invalid interval on blendfile line %ld\n", line_number);
        return FAIL;
    }

    if (support->xmin < options->xmin || support->xmax > options->xmax) {
        BLEND_Report(BLEND_MSG_ERROR, "window1d: blendfile line %ld interval %.12g/%.12g is outside -R domain %.12g/%.12g\n",
                     line_number, support->xmin, support->xmax, options->xmin, options->xmax);
        return FAIL;
    }

    if (blend_window_function_from_name(function_name, &support->function) != SUCCESS) {
        BLEND_Report(BLEND_MSG_ERROR, "window1d: unknown window function on blendfile line %ld: %s\n",
                     line_number, function_name);
        return FAIL;
    }

    if (window1d_parse_taper_ratio_values(taper_ratio, "blendfile", &support->ratio1, &support->ratio2) != SUCCESS) {
        BLEND_Report(BLEND_MSG_ERROR, "window1d: invalid taper ratio on blendfile line %ld\n", line_number);
        return FAIL;
    }

    if (support->function == WFUNC_BOXCAR && (support->ratio1 != 0.0 || support->ratio2 != 0.0)) {
        BLEND_Report(BLEND_MSG_WARNING, "window1d: taper ratios on blendfile line %ld are ignored for boxcar windows\n",
                     line_number);
    }

    return window1d_support_grid_size(options, support, &support->nx);
}

static int window1d_read_blendfile(const window1d_options *options, window1d_support_list *supports)
{
    FILE *fp;
    char line[1024];
    long line_number = 0;

    supports->items = NULL;
    supports->count = 0;
    supports->capacity = 0;

    fp = fopen(options->blendfile, "r");
    if (fp == NULL) {
        BLEND_Report(BLEND_MSG_ERROR, "%s: %s\n", options->blendfile, strerror(errno));
        return FAIL;
    }

    while (fgets(line, sizeof(line), fp) != NULL) {
        char *text;
        char *comment;
        window1d_support support;

        line_number++;
        comment = strchr(line, '#');
        if (comment != NULL) {
            *comment = '\0';
        }

        text = window1d_trim_line(line);
        if (*text == '\0') {
            continue;
        }

        if (window1d_parse_blendfile_line(text, line_number, options, &support) != SUCCESS) {
            fclose(fp);
            window1d_support_list_free(supports);
            return FAIL;
        }

        if (window1d_support_list_append(supports, &support) != SUCCESS) {
            fclose(fp);
            window1d_support_list_free(supports);
            return FAIL;
        }
    }

    if (ferror(fp)) {
        BLEND_Report(BLEND_MSG_ERROR, "%s: %s\n", options->blendfile, strerror(errno));
        fclose(fp);
        window1d_support_list_free(supports);
        return FAIL;
    }

    if (fclose(fp) != 0) {
        BLEND_Report(BLEND_MSG_ERROR, "%s: %s\n", options->blendfile, strerror(errno));
        window1d_support_list_free(supports);
        return FAIL;
    }

    if (supports->count == 0) {
        BLEND_Report(BLEND_MSG_ERROR, "window1d: blendfile contains no intervals\n");
        return FAIL;
    }

    return SUCCESS;
}

static void window1d_init_window(window *data, const window1d_options *options)
{
    memset(data, 0, sizeof(*data));
    data->nx = options->nx;
    data->ratio_x1 = options->ratio1;
    data->ratio_x2 = options->ratio2;
    data->x_function = options->function;
}

static double window1d_coordinate(const window1d_options *options, int index)
{
    return options->xmin + (double)index * options->dx;
}

static int window1d_bracket_coordinate(double x, const window1d_options *options, int *left, int *right)
{
    double scaled;

    if (x < options->xmin || x > options->xmax) {
        *left = -1;
        *right = -1;
        return SUCCESS;
    }

    scaled = (x - options->xmin) / options->dx;
    *left = (int)floor(scaled);

    if (*left < 0) {
        *left = 0;
    }
    else if (*left >= options->nx - 1) {
        *left = options->nx - 1;
    }

    if (*left == options->nx - 1) {
        *right = *left;
    }
    else {
        *right = *left + 1;
    }

    return SUCCESS;
}

static int window1d_grid_weight(int index, window *data, double *weight)
{
    if (index < 0) {
        *weight = 0.0;
        return SUCCESS;
    }

    if (embedding_contribution1d(index, data) != SUCCESS) {
        return FAIL;
    }

    *weight = data->contribution;
    return SUCCESS;
}

static int window1d_interpolate_weight(double x, const window1d_options *options, window *data, double *weight)
{
    int left, right;
    double x_left, x_right, weight_left, weight_right;

    if (window1d_bracket_coordinate(x, options, &left, &right) != SUCCESS) {
        return FAIL;
    }

    if (left < 0) {
        *weight = 0.0;
        return SUCCESS;
    }

    if (window1d_grid_weight(left, data, &weight_left) != SUCCESS) {
        return FAIL;
    }

    if (left == right) {
        *weight = weight_left;
        return SUCCESS;
    }

    if (window1d_grid_weight(right, data, &weight_right) != SUCCESS) {
        return FAIL;
    }

    x_left = window1d_coordinate(options, left);
    x_right = window1d_coordinate(options, right);
    return interpolate_linear(x_left, weight_left, x_right, weight_right, x, weight);
}

static double window1d_support_coordinate(const window1d_support *support, int index)
{
    double scale;

    if (support->nx <= 1) {
        return support->xmin;
    }

    scale = (double)index / (double)(support->nx - 1);
    return support->xmin + scale * (support->xmax - support->xmin);
}

static int window1d_support_grid_weight(int index, const window1d_support *support, double *weight)
{
    window data;

    memset(&data, 0, sizeof(data));
    data.nx = support->nx;
    data.ratio_x1 = support->ratio1;
    data.ratio_x2 = support->ratio2;
    data.x_function = support->function;

    if (embedding_contribution1d(index, &data) != SUCCESS) {
        return FAIL;
    }

    *weight = data.contribution;
    return SUCCESS;
}

static int window1d_support_weight(double x, const window1d_support *support, double *weight)
{
    int left, right;
    double scaled, x_left, x_right, weight_left, weight_right;

    if (x < support->xmin || x > support->xmax) {
        *weight = 0.0;
        return SUCCESS;
    }

    scaled = (x - support->xmin) * (double)(support->nx - 1) / (support->xmax - support->xmin);
    left = (int)floor(scaled);
    if (left < 0) {
        left = 0;
    }
    else if (left >= support->nx - 1) {
        left = support->nx - 1;
    }

    if (left == support->nx - 1) {
        return window1d_support_grid_weight(left, support, weight);
    }

    right = left + 1;
    if (window1d_support_grid_weight(left, support, &weight_left) != SUCCESS) {
        return FAIL;
    }
    if (window1d_support_grid_weight(right, support, &weight_right) != SUCCESS) {
        return FAIL;
    }

    x_left = window1d_support_coordinate(support, left);
    x_right = window1d_support_coordinate(support, right);
    return interpolate_linear(x_left, weight_left, x_right, weight_right, x, weight);
}

static int window1d_combine_weights(const double *weights, size_t count, char clobber, double *weight)
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
            for (i = 1; i < count; i++) {
                if (weights[i] < *weight) {
                    *weight = weights[i];
                }
            }
            break;
        case 'u':
            *weight = weights[0];
            for (i = 1; i < count; i++) {
                if (weights[i] > *weight) {
                    *weight = weights[i];
                }
            }
            break;
        case 'a':
            *weight = 0.0;
            for (i = 0; i < count; i++) {
                *weight += weights[i];
            }
            *weight /= (double)count;
            break;
        case 'g':
            *weight = 1.0;
            for (i = 0; i < count; i++) {
                *weight *= weights[i];
            }
            *weight = pow(*weight, 1.0 / (double)count);
            break;
        case 'p':
            *weight = 1.0;
            for (i = 0; i < count; i++) {
                *weight *= weights[i];
            }
            break;
        default:
            BLEND_Report(BLEND_MSG_ERROR, "window1d: unknown clobber mode: %c\n", clobber);
            return FAIL;
    }

    return SUCCESS;
}

static int window1d_blendfile_weight(double x, const window1d_options *options,
                                     const window1d_support_list *supports, double *weight)
{
    double *weights = NULL;
    size_t i, count = 0;
    int status;

    weights = (double *)calloc(supports->count, sizeof(*weights));
    if (weights == NULL) {
        BLEND_Report(BLEND_MSG_ERROR, "window1d: could not allocate overlap weights\n");
        return FAIL;
    }

    for (i = 0; i < supports->count; i++) {
        double support_weight;

        if (x < supports->items[i].xmin || x > supports->items[i].xmax) {
            continue;
        }

        if (window1d_support_weight(x, &supports->items[i], &support_weight) != SUCCESS) {
            free(weights);
            return FAIL;
        }

        weights[count] = support_weight;
        count++;
    }

    status = window1d_combine_weights(weights, count, options->clobber, weight);
    free(weights);
    return status;
}

static int window1d_write_weight(FILE *fp, double x, double weight)
{
    return fprintf(fp, "%.12g %.12f\n", x, weight) < 0 ? FAIL : SUCCESS;
}

static int window1d_get_weight(double x, int index, const window1d_options *options, window *data,
                               const window1d_support_list *supports, double *weight)
{
    if (supports != NULL) {
        return window1d_blendfile_weight(x, options, supports, weight);
    }

    if (index >= 0) {
        return window1d_grid_weight(index, data, weight);
    }

    return window1d_interpolate_weight(x, options, data, weight);
}

static int window1d_write_grid(FILE *fp, const window1d_options *options, window *data,
                               const window1d_support_list *supports)
{
    int i;

    for (i = 0; i < options->nx; i++) {
        double x = window1d_coordinate(options, i);
        double weight;

        if (window1d_get_weight(x, i, options, data, supports, &weight) != SUCCESS) {
            return FAIL;
        }
        if (window1d_write_weight(fp, x, weight) != SUCCESS) {
            return FAIL;
        }
    }

    return SUCCESS;
}

static int window1d_write_queries(FILE *fp, const window1d_options *options, window *data,
                                  const window1d_support_list *supports, int *n_queries)
{
    double x;
    int status;

    *n_queries = 0;
    while ((status = scanf("%lf", &x)) == 1) {
        if (!isfinite(x)) {
            BLEND_Report(BLEND_MSG_ERROR, "window1d: standard input contains a non-finite query\n");
            return FAIL;
        }
        {
            double weight;

            if (window1d_get_weight(x, -1, options, data, supports, &weight) != SUCCESS) {
                return FAIL;
            }
            if (window1d_write_weight(fp, x, weight) != SUCCESS) {
                return FAIL;
            }
        }
        *n_queries += 1;
    }

    if (status == 0) {
        BLEND_Report(BLEND_MSG_ERROR, "window1d: standard input contains a non-numeric query\n");
        return FAIL;
    }

    return SUCCESS;
}

int blend_window1d_module(int argc, char **argv)
{
    window1d_options options;
    window1d_support_list supports = {0};
    window1d_support_list *support_list = NULL;
    window data;
    int parse_status;
    int n_queries = 0;
    int status;

    parse_status = window1d_parse_options(argc, argv, &options);
    if (parse_status == 2) {
        return SUCCESS;
    }
    if (parse_status != SUCCESS) {
        if (blend_get_verbosity() >= BLEND_MSG_ERROR) {
            window1d_usage(stderr);
        }
        return FAIL;
    }

    window1d_init_window(&data, &options);

    if (options.blendfile != NULL) {
        double start = blend_elapsed_seconds();

        if (window1d_read_blendfile(&options, &supports) != SUCCESS) {
            return FAIL;
        }
        support_list = &supports;
        BLEND_Report(BLEND_MSG_TIMING,
                     "window1d: read %zu blendfile intervals in %.3f s\n",
                     supports.count, blend_elapsed_seconds() - start);
    }

    if (!isatty(STDIN_FILENO)) {
        double start = blend_elapsed_seconds();

        status = window1d_write_queries(stdout, &options, &data, support_list, &n_queries);
        if (status == SUCCESS && n_queries > 0) {
            BLEND_Report(BLEND_MSG_TIMING,
                         "window1d: evaluated %d query points in %.3f s\n",
                         n_queries, blend_elapsed_seconds() - start);
        }
        if (status != SUCCESS) {
            window1d_support_list_free(&supports);
            return FAIL;
        }
    }

    if (n_queries == 0) {
        double start = blend_elapsed_seconds();

        BLEND_Report(BLEND_MSG_TIMING,
                     "window1d: writing %d grid nodes\n",
                     options.nx);
        status = window1d_write_grid(stdout, &options, &data, support_list);
        if (status == SUCCESS) {
            BLEND_Report(BLEND_MSG_TIMING,
                         "window1d: wrote %d grid nodes in %.3f s\n",
                         options.nx, blend_elapsed_seconds() - start);
        }
        if (status != SUCCESS) {
            window1d_support_list_free(&supports);
            return FAIL;
        }
    }

    window1d_support_list_free(&supports);
    return SUCCESS;
}
