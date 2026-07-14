/*
* Author: Rasheed Ajala
* Contact: abioyeajala@gmail.com
*
*/

#ifndef BLEND_POLYGON_H
#define BLEND_POLYGON_H

#include <stddef.h>

typedef struct vertex {
    double x;
    double y;
} vertex;

typedef struct polygon {
    size_t n_vertices;
    vertex *vertices;
} polygon;

int blend_polygon_alloc(polygon *poly, size_t n_vertices);
void blend_polygon_free(polygon *poly);
int blend_polygon_set_vertex(polygon *poly, size_t index, double x, double y);
int blend_polygon_get_vertex(const polygon *poly, size_t index, vertex *out);
int blend_polygon_copy(const polygon *src, polygon *dst);
int blend_polygon_from_array(polygon *poly, double **vertices, size_t n_vertices);
int blend_polygon_read(const char *path, polygon *poly);
int blend_polygon_write(const char *path, const polygon *poly);
int blend_polygon_validate(const polygon *poly);
int blend_polygon_is_simple(const polygon *poly, int *is_simple);
int blend_polygon_bounds(const polygon *poly, double *xmin, double *xmax, double *ymin, double *ymax);
int blend_polygon_contains_point(const polygon *poly, double x, double y, int *inside);
int blend_polygon_is_closed(const polygon *poly, int *is_closed);
int blend_polygon_close(polygon *poly);
int blend_polygon_map_to_grid(const polygon *src, double xmin, double xmax, double ymin, double ymax,
                              int nx, int ny, polygon *dst);
int blend_polygon_is_xy_monotone(const polygon *poly, int *is_xy_monotone);
int blend_polygon_is_xy_monotone_strict(const polygon *poly, int *is_xy_monotone);
int blend_polygon_xy_monotone_envelope(const polygon *src, polygon *dst);
int blend_polygon_xy_monotone_envelope_strict(const polygon *src, polygon *dst);
int blend_polygon_xy_monotone_best_piecewise_envelope(const polygon *src, polygon *dst,
                                                      double xmin, double xmax,
                                                      double ymin, double ymax,
                                                      int nx, int ny);
int blend_polygon_xy_monotone_best_piecewise_envelope_strict(const polygon *src, polygon *dst,
                                                             double xmin, double xmax,
                                                             double ymin, double ymax,
                                                             int nx, int ny);

#endif /* BLEND_POLYGON_H */
