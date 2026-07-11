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

typedef struct monotone_options {
    const char *polygonfile;
    int has_monotone_method;
    int has_grid;
    int nx;
    int ny;
    char monotone_method;
} monotone_options;

static void monotone_usage(FILE *fp)
{
    fprintf(fp, "blend monotone - Check or enforce xy-monotone polygon boundaries\n\n");
    fprintf(fp, "usage: blend monotone [<polygonfile>] [-M<method>] [-G<nx>[/<ny>]] [-V[q|e|w|t|i|c|d]]\n\n");

    fprintf(fp, "Read a two-column polygon vertex file and check whether the polygon is simple\n");
    fprintf(fp, "and xy-monotone. If <polygonfile> is omitted, vertices are read from standard\n");
    fprintf(fp, "input. Default check mode writes a short status message to standard output.\n\n");

    fprintf(fp, "OPTIONAL ARGUMENTS:\n\n");
    fprintf(fp, "  -M<method>, --monotone=<method>\n");
    fprintf(fp, "     Convert a non-xy-monotone polygon using the selected method and write the\n");
    fprintf(fp, "     resulting polygon vertices to standard output. If the input polygon is\n");
    fprintf(fp, "     already xy-monotone, it is written unchanged.\n");
    fprintf(fp, "       e  Use the xy-monotone envelope method.\n");
    fprintf(fp, "       b  Use the highest-IoU piecewise envelope from both traversal directions.\n");
    fprintf(fp, "          This is sampled on the -G grid [Default is 256/256].\n\n");

    fprintf(fp, "  -G<nx>[/<ny>], --grid=<nx>[/<ny>]\n");
    fprintf(fp, "     Set the grid used to estimate IoU for -Mb. nx and ny must be positive\n");
    fprintf(fp, "     integers. If ny is omitted then ny = nx [Default is 256/256].\n\n");

    fprintf(fp, "  -?, --help\n");
    fprintf(fp, "     Print this usage message and exit.\n\n");

    fprintf(fp, "EXAMPLES:\n\n");
    fprintf(fp, "  blend monotone polygon.txt\n");
    fprintf(fp, "     Check whether polygon.txt is simple and xy-monotone.\n\n");
    fprintf(fp, "  blend monotone polygon.txt -Mb -G512/512 > polygon_monotone.txt\n");
    fprintf(fp, "     Convert polygon.txt using the best-IoU piecewise method and write the result.\n\n");
    fprintf(fp, "  cat polygon.txt | blend monotone\n");
    fprintf(fp, "     Read polygon vertices from standard input and check monotonicity.\n");
}

static int monotone_parse_method(const char *value, monotone_options *options)
{
    if (value == NULL || value[0] == '\0' || value[1] != '\0') {
        BLEND_Report(BLEND_MSG_ERROR, "monotone: monotone method must be one of e, b\n");
        return FAIL;
    }

    if (value[0] != 'e' && value[0] != 'b') {
        BLEND_Report(BLEND_MSG_ERROR, "monotone: unknown monotone method: %s\n", value);
        return FAIL;
    }

    options->monotone_method = value[0];
    options->has_monotone_method = 1;
    return SUCCESS;
}

static int monotone_parse_grid(const char *value, monotone_options *options)
{
    char *end = NULL;
    long nx, ny;

    if (value == NULL || *value == '\0') {
        BLEND_Report(BLEND_MSG_ERROR, "monotone: grid must be nx[/ny]\n");
        return FAIL;
    }

    errno = 0;
    nx = strtol(value, &end, 10);
    if (value == end || errno == ERANGE || nx <= 0 || nx > INT_MAX) {
        BLEND_Report(BLEND_MSG_ERROR, "monotone: grid nx must be a positive integer\n");
        return FAIL;
    }

    if (*end == '\0') {
        ny = nx;
    }
    else if (*end == '/') {
        const char *text = end + 1;

        errno = 0;
        ny = strtol(text, &end, 10);
        if (text == end || errno == ERANGE || ny <= 0 || ny > INT_MAX || *end != '\0') {
            BLEND_Report(BLEND_MSG_ERROR, "monotone: grid ny must be a positive integer\n");
            return FAIL;
        }
    }
    else {
        BLEND_Report(BLEND_MSG_ERROR, "monotone: grid must be nx[/ny]\n");
        return FAIL;
    }

    options->nx = (int)nx;
    options->ny = (int)ny;
    options->has_grid = 1;
    return SUCCESS;
}

static const char *monotone_option_value(int argc, char **argv, int *i, const char *option)
{
    if (*i + 1 >= argc) {
        BLEND_Report(BLEND_MSG_ERROR, "monotone: option %s requires an argument\n", option);
        return NULL;
    }

    *i += 1;
    return argv[*i];
}

static int monotone_parse_options(int argc, char **argv, monotone_options *options)
{
    int i;

    if (argc == 0 && isatty(STDIN_FILENO)) {
        monotone_usage(stdout);
        return 2;
    }

    memset(options, 0, sizeof(*options));
    options->nx = 256;
    options->ny = 256;

    for (i = 0; i < argc; i++) {
        const char *arg = argv[i];
        const char *value = NULL;

        if (strcmp(arg, "-?") == 0 || strcmp(arg, "--help") == 0) {
            monotone_usage(stdout);
            return 2;
        }
        else if (strcmp(arg, "-M") == 0 || strcmp(arg, "--monotone") == 0) {
            value = monotone_option_value(argc, argv, &i, arg);
            if (value == NULL || monotone_parse_method(value, options) != SUCCESS) return FAIL;
        }
        else if (strncmp(arg, "-M", 2) == 0 && arg[2] != '\0') {
            if (monotone_parse_method(arg + 2, options) != SUCCESS) return FAIL;
        }
        else if (strncmp(arg, "--monotone=", 11) == 0) {
            if (monotone_parse_method(arg + 11, options) != SUCCESS) return FAIL;
        }
        else if (strcmp(arg, "-G") == 0 || strcmp(arg, "--grid") == 0) {
            value = monotone_option_value(argc, argv, &i, arg);
            if (value == NULL || monotone_parse_grid(value, options) != SUCCESS) return FAIL;
        }
        else if (strncmp(arg, "-G", 2) == 0 && arg[2] != '\0') {
            if (monotone_parse_grid(arg + 2, options) != SUCCESS) return FAIL;
        }
        else if (strncmp(arg, "--grid=", 7) == 0) {
            if (monotone_parse_grid(arg + 7, options) != SUCCESS) return FAIL;
        }
        else if (arg[0] == '-') {
            BLEND_Report(BLEND_MSG_ERROR, "monotone: unknown option: %s\n", arg);
            return FAIL;
        }
        else {
            if (options->polygonfile != NULL) {
                BLEND_Report(BLEND_MSG_ERROR, "monotone: only one polygon file may be specified\n");
                return FAIL;
            }
            options->polygonfile = arg;
        }
    }

    if (options->has_grid && options->monotone_method != 'b') {
        BLEND_Report(BLEND_MSG_WARNING, "monotone: -G/--grid is ignored unless -Mb is used\n");
    }

    return SUCCESS;
}

static char *monotone_trim_line(char *line)
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

static int monotone_append_vertex(vertex **vertices, size_t *count, size_t *capacity, double x, double y)
{
    vertex *items;
    size_t new_capacity;

    if (*count == *capacity) {
        new_capacity = *capacity == 0 ? 16 : *capacity * 2;
        items = (vertex *)realloc(*vertices, new_capacity * sizeof(**vertices));
        if (items == NULL) {
            BLEND_Report(BLEND_MSG_ERROR, "monotone: could not allocate polygon vertices\n");
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

static int monotone_read_stream(FILE *fp, const char *name, polygon *poly)
{
    char line[1024];
    long line_number = 0;
    vertex *vertices = NULL;
    size_t count = 0, capacity = 0;
    polygon tmp;

    while (fgets(line, sizeof(line), fp) != NULL) {
        char *text;
        char *comment;
        char extra[2];
        double x, y;
        int fields;

        line_number++;
        comment = strchr(line, '#');
        if (comment != NULL) {
            *comment = '\0';
        }

        text = monotone_trim_line(line);
        if (*text == '\0') {
            continue;
        }

        fields = sscanf(text, "%lf %lf %1s", &x, &y, extra);
        if (fields != 2 || !isfinite(x) || !isfinite(y)) {
            BLEND_Report(BLEND_MSG_ERROR, "monotone: %s line %ld must contain two finite columns: x y\n",
                         name, line_number);
            free(vertices);
            return FAIL;
        }

        if (monotone_append_vertex(&vertices, &count, &capacity, x, y) != SUCCESS) {
            free(vertices);
            return FAIL;
        }
    }

    if (ferror(fp)) {
        BLEND_Report(BLEND_MSG_ERROR, "%s: %s\n", name, strerror(errno));
        free(vertices);
        return FAIL;
    }

    tmp.n_vertices = count;
    tmp.vertices = vertices;
    if (blend_polygon_validate(&tmp) != SUCCESS) {
        free(vertices);
        return FAIL;
    }

    blend_polygon_free(poly);
    poly->n_vertices = tmp.n_vertices;
    poly->vertices = tmp.vertices;
    return SUCCESS;
}

static int monotone_read_polygon(const monotone_options *options, polygon *poly)
{
    if (options->polygonfile != NULL) {
        return blend_polygon_read(options->polygonfile, poly);
    }

    if (isatty(STDIN_FILENO)) {
        BLEND_Report(BLEND_MSG_ERROR, "monotone: no polygon file or standard input was provided\n");
        return FAIL;
    }

    return monotone_read_stream(stdin, "standard input", poly);
}

static int monotone_write_polygon(FILE *fp, const polygon *poly)
{
    size_t i;

    for (i = 0; i < poly->n_vertices; i++) {
        if (fprintf(fp, "%.17g %.17g\n", poly->vertices[i].x, poly->vertices[i].y) < 0) {
            return FAIL;
        }
    }

    return SUCCESS;
}

int blend_monotone_module(int argc, char **argv)
{
    monotone_options options;
    polygon input = {0};
    polygon output = {0};
    double xmin, xmax, ymin, ymax;
    int parse_status;
    int is_xy_monotone = 0;
    int status = SUCCESS;

    parse_status = monotone_parse_options(argc, argv, &options);
    if (parse_status == 2) {
        return SUCCESS;
    }
    if (parse_status != SUCCESS) {
        if (blend_get_verbosity() >= BLEND_MSG_ERROR) {
            monotone_usage(stderr);
        }
        return FAIL;
    }

    if (monotone_read_polygon(&options, &input) != SUCCESS) {
        return FAIL;
    }

    if (blend_polygon_is_xy_monotone(&input, &is_xy_monotone) != SUCCESS) {
        blend_polygon_free(&input);
        return FAIL;
    }

    if (!options.has_monotone_method) {
        printf("monotone: polygon is %sxy-monotone\n", is_xy_monotone ? "" : "not ");
        blend_polygon_free(&input);
        return SUCCESS;
    }

    if (options.monotone_method == 'e' || options.monotone_method == 'b') {
        if (is_xy_monotone) {
            printf("monotone: polygon is already xy-monotone\n");
            blend_polygon_free(&input);
            return SUCCESS;
        }

        if (options.monotone_method == 'e') {
            BLEND_Report(BLEND_MSG_WARNING, "monotone: polygon is not xy-monotone; using envelope method\n");
            status = blend_polygon_xy_monotone_envelope(&input, &output);
        }
        else if (options.monotone_method == 'b') {
            BLEND_Report(BLEND_MSG_WARNING,
                         "monotone: polygon is not xy-monotone; using best IoU piecewise envelope method\n");
            if (blend_polygon_bounds(&input, &xmin, &xmax, &ymin, &ymax) != SUCCESS) {
                blend_polygon_free(&input);
                return FAIL;
            }
            status = blend_polygon_xy_monotone_best_piecewise_envelope(&input, &output,
                                                                       xmin, xmax, ymin, ymax,
                                                                       options.nx, options.ny);
        }
        if (status != SUCCESS) {
            blend_polygon_free(&input);
            return FAIL;
        }
    }
    else {
        BLEND_Report(BLEND_MSG_ERROR, "monotone: unknown monotone method: %c\n", options.monotone_method);
        blend_polygon_free(&input);
        return FAIL;
    }

    if (monotone_write_polygon(stdout, &output) != SUCCESS) {
        BLEND_Report(BLEND_MSG_ERROR, "monotone: could not write polygon to standard output\n");
        status = FAIL;
    }

    blend_polygon_free(&input);
    blend_polygon_free(&output);
    return status;
}
