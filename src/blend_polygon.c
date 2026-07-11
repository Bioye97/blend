/*
* Author: Rasheed Ajala
* Contact: abioyeajala@gmail.com
*
*/

#include "blend.h"

#include <limits.h>

int blend_polygon_alloc(polygon *poly, size_t n_vertices)
{
    if (poly == NULL || n_vertices == 0) {
        return FAIL;
    }

    poly->vertices = (vertex *)calloc(n_vertices, sizeof(vertex));
    if (poly->vertices == NULL) {
        poly->n_vertices = 0;
        return FAIL;
    }

    poly->n_vertices = n_vertices;
    return SUCCESS;
}

void blend_polygon_free(polygon *poly)
{
    if (poly == NULL) {
        return;
    }

    free(poly->vertices);
    poly->vertices = NULL;
    poly->n_vertices = 0;
}

int blend_polygon_set_vertex(polygon *poly, size_t index, double x, double y)
{
    if (poly == NULL || poly->vertices == NULL || index >= poly->n_vertices) {
        return FAIL;
    }

    poly->vertices[index].x = x;
    poly->vertices[index].y = y;
    return SUCCESS;
}

int blend_polygon_get_vertex(const polygon *poly, size_t index, vertex *out)
{
    if (poly == NULL || poly->vertices == NULL || out == NULL || index >= poly->n_vertices) {
        return FAIL;
    }

    out->x = poly->vertices[index].x;
    out->y = poly->vertices[index].y;
    return SUCCESS;
}

int blend_polygon_copy(const polygon *src, polygon *dst)
{
    size_t i;

    if (src == NULL || dst == NULL || src->vertices == NULL || src->n_vertices == 0) {
        return FAIL;
    }

    if (blend_polygon_alloc(dst, src->n_vertices) != SUCCESS) {
        return FAIL;
    }

    for (i = 0; i < src->n_vertices; i++) {
        dst->vertices[i] = src->vertices[i];
    }

    return SUCCESS;
}

int blend_polygon_from_array(polygon *poly, double **vertices, size_t n_vertices)
{
    size_t i;

    if (poly == NULL || vertices == NULL || n_vertices == 0) {
        return FAIL;
    }

    if (blend_polygon_alloc(poly, n_vertices) != SUCCESS) {
        return FAIL;
    }

    for (i = 0; i < n_vertices; i++) {
        if (vertices[i] == NULL) {
            blend_polygon_free(poly);
            return FAIL;
        }
        poly->vertices[i].x = vertices[i][0];
        poly->vertices[i].y = vertices[i][1];
    }

    return SUCCESS;
}

int blend_window_set_polygon(window *data, const polygon *poly)
{
    size_t i;

    if (data == NULL || poly == NULL || poly->vertices == NULL || poly->n_vertices == 0) {
        return FAIL;
    }
    if (poly->n_vertices > (size_t)INT_MAX) {
        return FAIL;
    }

    if (data->vertices != NULL) {
        blend_window_clear_polygon(data);
    }

    data->vertices = (double **)calloc(poly->n_vertices, sizeof(double *));
    if (data->vertices == NULL) {
        data->row_size = 0;
        return FAIL;
    }

    for (i = 0; i < poly->n_vertices; i++) {
        data->vertices[i] = (double *)calloc(2, sizeof(double));
        if (data->vertices[i] == NULL) {
            blend_window_clear_polygon(data);
            return FAIL;
        }
        data->vertices[i][0] = poly->vertices[i].x;
        data->vertices[i][1] = poly->vertices[i].y;
    }

    data->row_size = (int)poly->n_vertices;
    return SUCCESS;
}

void blend_window_clear_polygon(window *data)
{
    int i;

    if (data == NULL || data->vertices == NULL) {
        return;
    }

    for (i = 0; i < data->row_size; i++) {
        free(data->vertices[i]);
    }
    free(data->vertices);
    data->vertices = NULL;
    data->row_size = 0;
}
