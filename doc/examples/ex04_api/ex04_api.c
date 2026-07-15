#include "blend.h"

#include <limits.h>

#define REGION_XMIN -85.0
#define REGION_XMAX -30.0
#define REGION_YMIN -60.0
#define REGION_YMAX 15.0
#define INCREMENT 0.1

static int grid_size(double min, double max, double increment, int *n)
{
    double intervals = floor((max - min) / increment + 0.5);

    if (intervals < 1.0 || intervals > (double)(INT_MAX - 1)) {
        return FAIL;
    }

    *n = (int)intervals + 1;
    return SUCCESS;
}

static int output_axis_size(double min, double max, double increment, int *n)
{
    double intervals_exact = (max - min) / increment;
    double intervals_nearest = floor(intervals_exact + 0.5);

    if (fabs(intervals_exact - intervals_nearest) > 1.0e-9) {
        return FAIL;
    }

    *n = (int)intervals_nearest + 1;
    return SUCCESS;
}

static double snap_local_coordinate(double value, int max_index)
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

static int polygon_bounds(const polygon *poly, double *xmin, double *xmax,
                          double *ymin, double *ymax)
{
    size_t i;

    if (poly == NULL || poly->vertices == NULL || poly->n_vertices == 0) {
        return FAIL;
    }

    *xmin = poly->vertices[0].x;
    *xmax = poly->vertices[0].x;
    *ymin = poly->vertices[0].y;
    *ymax = poly->vertices[0].y;

    for (i = 0; i < poly->n_vertices; i++) {
        double x = poly->vertices[i].x;
        double y = poly->vertices[i].y;

        if (x < REGION_XMIN || x > REGION_XMAX || y < REGION_YMIN || y > REGION_YMAX) {
            return FAIL;
        }
        if (x < *xmin) *xmin = x;
        if (x > *xmax) *xmax = x;
        if (y < *ymin) *ymin = y;
        if (y > *ymax) *ymax = y;
    }

    return SUCCESS;
}

static int map_polygon_to_local_grid(const polygon *src,
                                     double xmin, double xmax,
                                     double ymin, double ymax,
                                     int nx, int ny, polygon *dst)
{
    polygon tmp = {0};
    size_t i;

    if (blend_polygon_alloc(&tmp, src->n_vertices) != SUCCESS) {
        return FAIL;
    }

    for (i = 0; i < src->n_vertices; i++) {
        double x = (src->vertices[i].x - xmin) * (double)(nx - 1) / (xmax - xmin);
        double y = (src->vertices[i].y - ymin) * (double)(ny - 1) / (ymax - ymin);

        if (blend_polygon_set_vertex(&tmp, i,
                                     snap_local_coordinate(x, nx - 1),
                                     snap_local_coordinate(y, ny - 1)) != SUCCESS) {
            blend_polygon_free(&tmp);
            return FAIL;
        }
    }

    if (blend_polygon_validate(&tmp) != SUCCESS) {
        blend_polygon_free(&tmp);
        return FAIL;
    }

    blend_polygon_free(dst);
    dst->n_vertices = tmp.n_vertices;
    dst->vertices = tmp.vertices;
    return SUCCESS;
}

static int map_polygon_to_real_coordinates(const polygon *src,
                                           double xmin, double xmax,
                                           double ymin, double ymax,
                                           int nx, int ny, polygon *dst)
{
    polygon tmp = {0};
    size_t i;

    if (blend_polygon_alloc(&tmp, src->n_vertices) != SUCCESS) {
        return FAIL;
    }

    for (i = 0; i < src->n_vertices; i++) {
        double x = xmin + src->vertices[i].x * (xmax - xmin) / (double)(nx - 1);
        double y = ymin + src->vertices[i].y * (ymax - ymin) / (double)(ny - 1);

        if (blend_polygon_set_vertex(&tmp, i, x, y) != SUCCESS) {
            blend_polygon_free(&tmp);
            return FAIL;
        }
    }

    if (blend_polygon_validate(&tmp) != SUCCESS) {
        blend_polygon_free(&tmp);
        return FAIL;
    }

    blend_polygon_free(dst);
    dst->n_vertices = tmp.n_vertices;
    dst->vertices = tmp.vertices;
    return SUCCESS;
}

static int point_on_segment(double x, double y, const vertex *a, const vertex *b)
{
    double cross = (x - a->x) * (b->y - a->y) - (y - a->y) * (b->x - a->x);
    double dot;

    if (fabs(cross) > 1.0e-10) {
        return 0;
    }

    dot = (x - a->x) * (x - b->x) + (y - a->y) * (y - b->y);
    return dot <= 1.0e-10;
}

static int point_in_polygon(double x, double y, const polygon *poly)
{
    int inside = 0;
    size_t i, j;

    for (i = 0, j = poly->n_vertices - 1; i < poly->n_vertices; j = i++) {
        const vertex *vi = &poly->vertices[i];
        const vertex *vj = &poly->vertices[j];

        if (point_on_segment(x, y, vj, vi)) {
            return 1;
        }
        if (((vi->y > y) != (vj->y > y)) &&
            (x < (vj->x - vi->x) * (y - vi->y) / (vj->y - vi->y) + vi->x)) {
            inside = !inside;
        }
    }

    return inside;
}

static int support_grid_weight(window *data, int ix, int iy, double *weight)
{
    if (ix < 0 || iy < 0 || ix >= data->nx || iy >= data->ny) {
        *weight = 0.0;
        return SUCCESS;
    }
    if (embedding_contribution2d(ix, iy, data) != SUCCESS) {
        return FAIL;
    }

    *weight = data->contribution;
    return SUCCESS;
}

static int support_weight(double x, double y, window *data, const polygon *support,
                          double xmin, double xmax, double ymin, double ymax,
                          double *weight)
{
    double sx, sy;
    double w00, w10, w01, w11;
    int ix0, iy0, ix1, iy1;

    if (x < xmin || x > xmax || y < ymin || y > ymax ||
        !point_in_polygon(x, y, support)) {
        *weight = 0.0;
        return SUCCESS;
    }

    sx = (x - xmin) * (double)(data->nx - 1) / (xmax - xmin);
    sy = (y - ymin) * (double)(data->ny - 1) / (ymax - ymin);
    ix0 = (int)floor(sx);
    iy0 = (int)floor(sy);

    if (ix0 < 0) ix0 = 0;
    if (iy0 < 0) iy0 = 0;
    if (ix0 >= data->nx - 1) ix0 = data->nx - 1;
    if (iy0 >= data->ny - 1) iy0 = data->ny - 1;

    ix1 = ix0 < data->nx - 1 ? ix0 + 1 : ix0;
    iy1 = iy0 < data->ny - 1 ? iy0 + 1 : iy0;

    if (support_grid_weight(data, ix0, iy0, &w00) != SUCCESS ||
        support_grid_weight(data, ix1, iy0, &w10) != SUCCESS ||
        support_grid_weight(data, ix0, iy1, &w01) != SUCCESS ||
        support_grid_weight(data, ix1, iy1, &w11) != SUCCESS) {
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

    return interpolate_bilinear((double)ix0, (double)ix1, (double)iy0, (double)iy1,
                                w00, w10, w01, w11, sx, sy, weight);
}

int main(void)
{
    const char *polygon_path = "south_america.txt";
    const char *grid_path = "ex04_api_grid.txt";
    polygon input = {0};
    polygon local_input = {0};
    polygon local_support = {0};
    polygon support = {0};
    window data = {0};
    permuted_vertex boundary = {0};
    FILE *fp = NULL;
    double xmin, xmax, ymin, ymax;
    int nx, ny, out_nx, out_ny;
    int is_strict = 0;
    int ix, iy;
    int status = FAIL;

    if (blend_polygon_read(polygon_path, &input) != SUCCESS ||
        blend_polygon_validate(&input) != SUCCESS ||
        polygon_bounds(&input, &xmin, &xmax, &ymin, &ymax) != SUCCESS ||
        grid_size(xmin, xmax, INCREMENT, &nx) != SUCCESS ||
        grid_size(ymin, ymax, INCREMENT, &ny) != SUCCESS ||
        map_polygon_to_local_grid(&input, xmin, xmax, ymin, ymax, nx, ny,
                                  &local_input) != SUCCESS) {
        goto cleanup;
    }

    if (blend_polygon_is_xy_monotone_strict(&local_input, &is_strict) != SUCCESS) {
        goto cleanup;
    }

    if (is_strict) {
        if (blend_polygon_copy(&local_input, &local_support) != SUCCESS) {
            goto cleanup;
        }
    }
    else if (blend_polygon_xy_monotone_envelope_strict(&local_input,
                                                       &local_support) != SUCCESS) {
        goto cleanup;
    }

    if (blend_polygon_is_xy_monotone_strict(&local_support, &is_strict) != SUCCESS ||
        !is_strict ||
        map_polygon_to_real_coordinates(&local_support, xmin, xmax, ymin, ymax,
                                        nx, ny, &support) != SUCCESS ||
        blend_polygon_write("south_america_monotone.txt", &support) != SUCCESS) {
        goto cleanup;
    }

    data.nx = nx;
    data.ny = ny;
    data.ratio_x1 = 0.49;
    data.ratio_x2 = 0.49;
    data.ratio_y1 = 0.49;
    data.ratio_y2 = 0.49;
    data.x_function = WFUNC_COSINE;
    data.y_function = WFUNC_COSINE;

    if (blend_window_set_polygon(&data, &local_support) != SUCCESS ||
        boundary_assembly(&data, &boundary) != SUCCESS ||
        output_axis_size(REGION_XMIN, REGION_XMAX, INCREMENT, &out_nx) != SUCCESS ||
        output_axis_size(REGION_YMIN, REGION_YMAX, INCREMENT, &out_ny) != SUCCESS) {
        goto cleanup;
    }

    fp = fopen(grid_path, "w");
    if (fp == NULL) {
        perror(grid_path);
        goto cleanup;
    }

    for (iy = 0; iy < out_ny; iy++) {
        double y = REGION_YMIN + (double)iy * INCREMENT;

        for (ix = 0; ix < out_nx; ix++) {
            double x = REGION_XMIN + (double)ix * INCREMENT;
            double weight;

            if (support_weight(x, y, &data, &support, xmin, xmax, ymin, ymax,
                               &weight) != SUCCESS ||
                fprintf(fp, "%.12g %.12g %.12f\n", x, y, weight) < 0) {
                goto cleanup;
            }
        }
        if (fprintf(fp, "\n") < 0) {
            goto cleanup;
        }
    }

    if (fclose(fp) != 0) {
        fp = NULL;
        perror(grid_path);
        goto cleanup;
    }
    fp = NULL;
    status = SUCCESS;

cleanup:
    if (fp != NULL) {
        fclose(fp);
    }
    blend_permuted_vertex_free(&boundary);
    blend_window_boundary_clear(&data);
    blend_polygon_free(&support);
    blend_polygon_free(&local_support);
    blend_polygon_free(&local_input);
    blend_polygon_free(&input);
    return status;
}
