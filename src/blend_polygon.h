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

#endif /* BLEND_POLYGON_H */
