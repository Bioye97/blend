/*
* Author: Rasheed Ajala
* Contact: abioyeajala@gmail.com
*
*/

#include "blend.h"

#include <ctype.h>
#include <errno.h>
#include <limits.h>

static double blend_polygon_coordinate_tolerance(const polygon *poly)
{
    double scale = 1.0;
    size_t i;

    if (poly == NULL || poly->vertices == NULL) {
        return 1.0e-12;
    }

    for (i = 0; i < poly->n_vertices; i++) {
        scale = fmax(scale, fabs(poly->vertices[i].x));
        scale = fmax(scale, fabs(poly->vertices[i].y));
    }

    return 1.0e-12 * scale;
}

static int blend_polygon_vertices_equal(vertex a, vertex b, double tolerance)
{
    return fabs(a.x - b.x) <= tolerance && fabs(a.y - b.y) <= tolerance;
}

static double blend_polygon_signed_area(const polygon *poly)
{
    double area = 0.0;
    size_t i, j;

    for (i = 0, j = poly->n_vertices - 1; i < poly->n_vertices; j = i++) {
        area += poly->vertices[j].x * poly->vertices[i].y - poly->vertices[i].x * poly->vertices[j].y;
    }

    return 0.5 * area;
}

static int blend_polygon_append_vertex(vertex **vertices, size_t *count, size_t *capacity, double x, double y)
{
    vertex *items;
    size_t new_capacity;

    if (*count == *capacity) {
        new_capacity = *capacity == 0 ? 16 : *capacity * 2;
        items = (vertex *)realloc(*vertices, new_capacity * sizeof(**vertices));
        if (items == NULL) {
            BLEND_Report(BLEND_MSG_ERROR, "polygon: could not allocate vertices\n");
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

static char *blend_polygon_trim_line(char *line)
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

static int blend_polygon_point_on_segment(double x, double y, const vertex *a, const vertex *b)
{
    double dx = b->x - a->x;
    double dy = b->y - a->y;
    double cross = (x - a->x) * dy - (y - a->y) * dx;
    double scale = fmax(1.0, fmax(fabs(dx), fabs(dy)));
    double tolerance = 1.0e-10 * scale;

    if (fabs(cross) > tolerance) {
        return 0;
    }
    if (x < fmin(a->x, b->x) - tolerance || x > fmax(a->x, b->x) + tolerance ||
        y < fmin(a->y, b->y) - tolerance || y > fmax(a->y, b->y) + tolerance) {
        return 0;
    }

    return 1;
}

static double blend_polygon_orientation(vertex a, vertex b, vertex c)
{
    return (b.x - a.x) * (c.y - a.y) - (b.y - a.y) * (c.x - a.x);
}

static int blend_polygon_segments_intersect(vertex a, vertex b, vertex c, vertex d, double tolerance)
{
    double o1 = blend_polygon_orientation(a, b, c);
    double o2 = blend_polygon_orientation(a, b, d);
    double o3 = blend_polygon_orientation(c, d, a);
    double o4 = blend_polygon_orientation(c, d, b);

    if (((o1 > tolerance && o2 < -tolerance) || (o1 < -tolerance && o2 > tolerance)) &&
        ((o3 > tolerance && o4 < -tolerance) || (o3 < -tolerance && o4 > tolerance))) {
        return 1;
    }

    if (fabs(o1) <= tolerance && blend_polygon_point_on_segment(c.x, c.y, &a, &b)) return 1;
    if (fabs(o2) <= tolerance && blend_polygon_point_on_segment(d.x, d.y, &a, &b)) return 1;
    if (fabs(o3) <= tolerance && blend_polygon_point_on_segment(a.x, a.y, &c, &d)) return 1;
    if (fabs(o4) <= tolerance && blend_polygon_point_on_segment(b.x, b.y, &c, &d)) return 1;

    return 0;
}

static int blend_polygon_edges_are_adjacent(size_t i, size_t j, size_t n_edges)
{
    if (i == j) {
        return 1;
    }
    if (i + 1 == j || j + 1 == i) {
        return 1;
    }
    if ((i == 0 && j == n_edges - 1) || (j == 0 && i == n_edges - 1)) {
        return 1;
    }

    return 0;
}

static int blend_polygon_compare_vertex(const void *a, const void *b)
{
    const vertex *va = (const vertex *)a;
    const vertex *vb = (const vertex *)b;

    if (va->x < vb->x) return -1;
    if (va->x > vb->x) return 1;
    if (va->y < vb->y) return -1;
    if (va->y > vb->y) return 1;
    return 0;
}

static int blend_polygon_compare_x_desc(const void *a, const void *b)
{
    const vertex *va = (const vertex *)a;
    const vertex *vb = (const vertex *)b;

    if (va->x > vb->x) return -1;
    if (va->x < vb->x) return 1;
    if (va->y < vb->y) return -1;
    if (va->y > vb->y) return 1;
    return 0;
}

static int blend_polygon_compare_y_asc_x_asc(const void *a, const void *b)
{
    const vertex *va = (const vertex *)a;
    const vertex *vb = (const vertex *)b;

    if (va->y < vb->y) return -1;
    if (va->y > vb->y) return 1;
    if (va->x < vb->x) return -1;
    if (va->x > vb->x) return 1;
    return 0;
}

static int blend_polygon_compare_y_asc_x_desc(const void *a, const void *b)
{
    const vertex *va = (const vertex *)a;
    const vertex *vb = (const vertex *)b;

    if (va->y < vb->y) return -1;
    if (va->y > vb->y) return 1;
    if (va->x > vb->x) return -1;
    if (va->x < vb->x) return 1;
    return 0;
}

static double blend_polygon_axis_value(vertex point, int axis)
{
    return axis == 0 ? point.x : point.y;
}

static int blend_polygon_value_between_open(double value, double a, double b, double tolerance)
{
    double min_value = fmin(a, b);
    double max_value = fmax(a, b);

    return value > min_value + tolerance && value < max_value - tolerance;
}

static int blend_polygon_point_on_bbox_side(vertex point, double xmin, double xmax,
                                            double ymin, double ymax, double tolerance)
{
    return fabs(point.x - xmin) <= tolerance || fabs(point.x - xmax) <= tolerance ||
           fabs(point.y - ymin) <= tolerance || fabs(point.y - ymax) <= tolerance;
}

static int blend_polygon_append_clean(vertex **vertices, size_t *count, size_t *capacity,
                                      vertex point, double tolerance)
{
    vertex *items;
    size_t new_capacity;

    if (*count > 0 && blend_polygon_vertices_equal((*vertices)[*count - 1], point, tolerance)) {
        return SUCCESS;
    }

    if (*count == *capacity) {
        new_capacity = *capacity == 0 ? 16 : *capacity * 2;
        items = (vertex *)realloc(*vertices, new_capacity * sizeof(**vertices));
        if (items == NULL) {
            BLEND_Report(BLEND_MSG_ERROR, "polygon: could not allocate piecewise envelope vertices\n");
            return FAIL;
        }
        *vertices = items;
        *capacity = new_capacity;
    }

    (*vertices)[*count] = point;
    *count += 1;
    return SUCCESS;
}

static int blend_polygon_append_chain(vertex **vertices, size_t *count, size_t *capacity,
                                      const vertex *chain, size_t chain_count, double tolerance)
{
    size_t i;

    for (i = 0; i < chain_count; i++) {
        if (blend_polygon_append_clean(vertices, count, capacity, chain[i], tolerance) != SUCCESS) {
            return FAIL;
        }
    }

    return SUCCESS;
}

static int blend_polygon_chain_is_strictly_monotone(const vertex *chain, size_t count,
                                                    int x_direction, int y_direction, double tolerance)
{
    size_t i;

    if (count < 2) {
        return 1;
    }

    for (i = 1; i < count; i++) {
        double dx = chain[i].x - chain[i - 1].x;
        double dy = chain[i].y - chain[i - 1].y;

        if ((double)x_direction * dx <= tolerance || (double)y_direction * dy <= tolerance) {
            return 0;
        }
    }

    return 1;
}

static int blend_polygon_adjust_corner_chain(vertex *chain, size_t *count,
                                             int x_direction, int y_direction, double tolerance)
{
    vertex end;
    size_t i;
    size_t write_index = 1;

    if (chain == NULL || count == NULL || *count < 2) {
        return FAIL;
    }

    end = chain[*count - 1];
    for (i = 1; i + 1 < *count; i++) {
        double dx = chain[i].x - chain[write_index - 1].x;
        double dy = chain[i].y - chain[write_index - 1].y;

        if ((double)x_direction * dx > tolerance && (double)y_direction * dy > tolerance) {
            chain[write_index] = chain[i];
            write_index++;
        }
    }

    while (write_index > 1 &&
           ((double)x_direction * (end.x - chain[write_index - 1].x) <= tolerance ||
            (double)y_direction * (end.y - chain[write_index - 1].y) <= tolerance)) {
        write_index--;
    }
    chain[write_index] = end;
    write_index++;
    *count = write_index;

    return blend_polygon_chain_is_strictly_monotone(chain, *count, x_direction, y_direction, tolerance)
               ? SUCCESS
               : FAIL;
}

static int blend_polygon_build_convex_hull(const vertex *input, size_t input_count, double tolerance,
                                           vertex **hull_out, size_t *hull_count_out)
{
    vertex *points = NULL;
    vertex *hull = NULL;
    size_t n_points = 0;
    size_t k = 0;
    size_t i;

    if (input == NULL || input_count < 3 || hull_out == NULL || hull_count_out == NULL) {
        return FAIL;
    }

    points = (vertex *)calloc(input_count, sizeof(*points));
    hull = (vertex *)calloc(2 * input_count, sizeof(*hull));
    if (points == NULL || hull == NULL) {
        BLEND_Report(BLEND_MSG_ERROR, "polygon: could not allocate convex hull\n");
        free(points);
        free(hull);
        return FAIL;
    }

    for (i = 0; i < input_count; i++) {
        points[i] = input[i];
    }
    qsort(points, input_count, sizeof(*points), blend_polygon_compare_vertex);

    for (i = 0; i < input_count; i++) {
        if (n_points == 0 || !blend_polygon_vertices_equal(points[i], points[n_points - 1], tolerance)) {
            points[n_points] = points[i];
            n_points++;
        }
    }

    if (n_points < 3) {
        free(points);
        free(hull);
        return FAIL;
    }

    for (i = 0; i < n_points; i++) {
        while (k >= 2 && blend_polygon_orientation(hull[k - 2], hull[k - 1], points[i]) <= tolerance) {
            k--;
        }
        hull[k] = points[i];
        k++;
    }

    {
        size_t lower_size = k;
        size_t reverse_index;

        for (reverse_index = n_points - 1; reverse_index > 0; reverse_index--) {
            vertex point = points[reverse_index - 1];

            while (k > lower_size && blend_polygon_orientation(hull[k - 2], hull[k - 1], point) <= tolerance) {
                k--;
            }
            hull[k] = point;
            k++;
        }
    }

    if (k > 1) {
        k--;
    }
    if (k < 3) {
        free(points);
        free(hull);
        return FAIL;
    }

    *hull_out = hull;
    *hull_count_out = k;
    free(points);
    return SUCCESS;
}

static int blend_polygon_find_vertex(const vertex *vertices, size_t count, vertex point, double tolerance,
                                     size_t *index_out)
{
    size_t i;

    for (i = 0; i < count; i++) {
        if (blend_polygon_vertices_equal(vertices[i], point, tolerance)) {
            *index_out = i;
            return SUCCESS;
        }
    }

    return FAIL;
}

static int blend_polygon_hull_path(const vertex *hull, size_t hull_count, size_t start_index, size_t end_index,
                                   int step, vertex **path_out, size_t *path_count_out)
{
    vertex *path;
    size_t count = 1;
    size_t index = start_index;
    size_t i = 0;

    while (index != end_index) {
        count++;
        if (step > 0) {
            index = (index + 1) % hull_count;
        }
        else {
            index = index == 0 ? hull_count - 1 : index - 1;
        }
    }

    path = (vertex *)calloc(count, sizeof(*path));
    if (path == NULL) {
        BLEND_Report(BLEND_MSG_ERROR, "polygon: could not allocate convex hull path\n");
        return FAIL;
    }

    index = start_index;
    while (1) {
        path[i] = hull[index];
        i++;
        if (index == end_index) {
            break;
        }
        if (step > 0) {
            index = (index + 1) % hull_count;
        }
        else {
            index = index == 0 ? hull_count - 1 : index - 1;
        }
    }

    *path_out = path;
    *path_count_out = count;
    return SUCCESS;
}

static int blend_polygon_corner_chain(vertex start, vertex end, const vertex *points, size_t point_count,
                                      int x_direction, int y_direction, double tolerance,
                                      vertex **chain_out, size_t *chain_count_out)
{
    vertex *chain = NULL;
    vertex *hull = NULL;
    vertex *path_forward = NULL;
    vertex *path_backward = NULL;
    size_t chain_count = point_count + 2;
    size_t hull_count = 0;
    size_t forward_count = 0;
    size_t backward_count = 0;
    size_t start_index = 0;
    size_t end_index = 0;
    size_t i;
    int forward_ok = 0;
    int backward_ok = 0;

    if (chain_out == NULL || chain_count_out == NULL) {
        return FAIL;
    }

    if ((double)x_direction * (end.x - start.x) <= tolerance ||
        (double)y_direction * (end.y - start.y) <= tolerance) {
        *chain_out = NULL;
        *chain_count_out = 0;
        return SUCCESS;
    }

    chain = (vertex *)calloc(chain_count, sizeof(*chain));
    if (chain == NULL) {
        BLEND_Report(BLEND_MSG_ERROR, "polygon: could not allocate corner chain\n");
        return FAIL;
    }

    chain[0] = start;
    for (i = 0; i < point_count; i++) {
        chain[i + 1] = points[i];
    }
    chain[chain_count - 1] = end;

    qsort(chain, chain_count, sizeof(*chain),
          x_direction > 0 ? blend_polygon_compare_y_asc_x_asc : blend_polygon_compare_y_asc_x_desc);

    if (blend_polygon_chain_is_strictly_monotone(chain, chain_count, x_direction, y_direction, tolerance) ||
        blend_polygon_adjust_corner_chain(chain, &chain_count, x_direction, y_direction, tolerance) == SUCCESS) {
        *chain_out = chain;
        *chain_count_out = chain_count;
        return SUCCESS;
    }

    if (blend_polygon_build_convex_hull(chain, chain_count, tolerance, &hull, &hull_count) != SUCCESS ||
        blend_polygon_find_vertex(hull, hull_count, start, tolerance, &start_index) != SUCCESS ||
        blend_polygon_find_vertex(hull, hull_count, end, tolerance, &end_index) != SUCCESS) {
        chain[0] = start;
        chain[1] = end;
        *chain_out = chain;
        *chain_count_out = 2;
        free(hull);
        return SUCCESS;
    }

    if (blend_polygon_hull_path(hull, hull_count, start_index, end_index, 1,
                                &path_forward, &forward_count) == SUCCESS) {
        forward_ok = blend_polygon_chain_is_strictly_monotone(path_forward, forward_count,
                                                              x_direction, y_direction, tolerance);
    }
    if (blend_polygon_hull_path(hull, hull_count, start_index, end_index, -1,
                                &path_backward, &backward_count) == SUCCESS) {
        backward_ok = blend_polygon_chain_is_strictly_monotone(path_backward, backward_count,
                                                               x_direction, y_direction, tolerance);
    }

    if (forward_ok && (!backward_ok || forward_count >= backward_count)) {
        free(chain);
        free(path_backward);
        free(hull);
        *chain_out = path_forward;
        *chain_count_out = forward_count;
        return SUCCESS;
    }
    if (backward_ok) {
        free(chain);
        free(path_forward);
        free(hull);
        *chain_out = path_backward;
        *chain_count_out = backward_count;
        return SUCCESS;
    }

    free(path_forward);
    free(path_backward);
    free(hull);
    chain[0] = start;
    chain[1] = end;
    *chain_out = chain;
    *chain_count_out = 2;
    return SUCCESS;
}

static int blend_polygon_axis_path_is_nondecreasing(const polygon *poly, size_t start, size_t end,
                                                    int step, int axis, double tolerance)
{
    size_t n = poly->n_vertices;
    size_t index = start;
    double previous = blend_polygon_axis_value(poly->vertices[index], axis);

    while (index != end) {
        double current;

        if (step > 0) {
            index = (index + 1) % n;
        }
        else {
            index = index == 0 ? n - 1 : index - 1;
        }

        current = blend_polygon_axis_value(poly->vertices[index], axis);
        if (current + tolerance < previous) {
            return 0;
        }
        previous = current;
    }

    return 1;
}

static int blend_polygon_axis_is_monotone(const polygon *poly, int axis)
{
    double min_value, max_value;
    double tolerance = blend_polygon_coordinate_tolerance(poly);
    size_t i, min_index, max_index;

    min_value = blend_polygon_axis_value(poly->vertices[0], axis);
    max_value = min_value;
    for (i = 1; i < poly->n_vertices; i++) {
        double value = blend_polygon_axis_value(poly->vertices[i], axis);

        if (value < min_value) min_value = value;
        if (value > max_value) max_value = value;
    }

    if (fabs(max_value - min_value) <= tolerance) {
        return 0;
    }

    for (min_index = 0; min_index < poly->n_vertices; min_index++) {
        double min_candidate = blend_polygon_axis_value(poly->vertices[min_index], axis);

        if (fabs(min_candidate - min_value) > tolerance) {
            continue;
        }

        for (max_index = 0; max_index < poly->n_vertices; max_index++) {
            double max_candidate = blend_polygon_axis_value(poly->vertices[max_index], axis);

            if (fabs(max_candidate - max_value) > tolerance) {
                continue;
            }

            if (blend_polygon_axis_path_is_nondecreasing(poly, min_index, max_index, 1, axis, tolerance) &&
                blend_polygon_axis_path_is_nondecreasing(poly, min_index, max_index, -1, axis, tolerance)) {
                return 1;
            }
        }
    }

    return 0;
}

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
    if (!isfinite(x) || !isfinite(y)) {
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

static int blend_polygon_corner_chain_from_order(vertex start, vertex end, const vertex *points,
                                                 size_t point_count, int x_direction, int y_direction,
                                                 double tolerance, vertex **chain_out, size_t *chain_count_out)
{
    vertex *chain = NULL;
    size_t chain_count = point_count + 2;
    size_t i;

    if (chain_out == NULL || chain_count_out == NULL) {
        return FAIL;
    }

    if ((double)x_direction * (end.x - start.x) <= tolerance ||
        (double)y_direction * (end.y - start.y) <= tolerance) {
        *chain_out = NULL;
        *chain_count_out = 0;
        return SUCCESS;
    }

    chain = (vertex *)calloc(chain_count, sizeof(*chain));
    if (chain == NULL) {
        BLEND_Report(BLEND_MSG_ERROR, "polygon: could not allocate traversal corner chain\n");
        return FAIL;
    }

    chain[0] = start;
    for (i = 0; i < point_count; i++) {
        chain[i + 1] = points[i];
    }
    chain[chain_count - 1] = end;

    if (blend_polygon_chain_is_strictly_monotone(chain, chain_count, x_direction, y_direction, tolerance) ||
        blend_polygon_adjust_corner_chain(chain, &chain_count, x_direction, y_direction, tolerance) == SUCCESS) {
        *chain_out = chain;
        *chain_count_out = chain_count;
        return SUCCESS;
    }

    free(chain);
    return blend_polygon_corner_chain(start, end, points, point_count, x_direction, y_direction, tolerance,
                                      chain_out, chain_count_out);
}

static int blend_polygon_collect_traversal_sector(const polygon *poly, size_t start_index, size_t end_index,
                                                  vertex start, vertex end, double xmin, double xmax,
                                                  double ymin, double ymax, double tolerance,
                                                  vertex *points, size_t *point_count)
{
    size_t index;
    size_t guard = 0;

    if (poly == NULL || points == NULL || point_count == NULL || poly->n_vertices == 0) {
        return FAIL;
    }

    *point_count = 0;
    index = start_index;
    while (index != end_index) {
        vertex point;

        index = (index + 1) % poly->n_vertices;
        if (index == end_index) {
            break;
        }

        point = poly->vertices[index];
        if (!blend_polygon_point_on_bbox_side(point, xmin, xmax, ymin, ymax, tolerance) &&
            blend_polygon_value_between_open(point.x, start.x, end.x, tolerance) &&
            blend_polygon_value_between_open(point.y, start.y, end.y, tolerance)) {
            points[*point_count] = point;
            *point_count += 1;
        }

        guard++;
        if (guard > poly->n_vertices) {
            return FAIL;
        }
    }

    return SUCCESS;
}

typedef struct blend_polygon_path_candidate {
    int valid;
    size_t count;
    size_t *indices;
    vertex *points;
} blend_polygon_path_candidate;

static void blend_polygon_path_candidate_free(blend_polygon_path_candidate *candidate)
{
    if (candidate == NULL) {
        return;
    }

    free(candidate->indices);
    free(candidate->points);
    candidate->indices = NULL;
    candidate->points = NULL;
    candidate->count = 0;
    candidate->valid = 0;
}

static int blend_polygon_path_candidate_alloc(blend_polygon_path_candidate *candidate, size_t capacity)
{
    if (candidate == NULL) {
        return FAIL;
    }

    memset(candidate, 0, sizeof(*candidate));
    candidate->indices = (size_t *)calloc(capacity == 0 ? 1 : capacity, sizeof(*candidate->indices));
    candidate->points = (vertex *)calloc(capacity == 0 ? 1 : capacity, sizeof(*candidate->points));
    if (candidate->indices == NULL || candidate->points == NULL) {
        blend_polygon_path_candidate_free(candidate);
        BLEND_Report(BLEND_MSG_ERROR, "polygon: could not allocate traversal path candidate\n");
        return FAIL;
    }

    return SUCCESS;
}

static int blend_polygon_collect_strict_path_candidate(const polygon *poly, size_t n_input,
                                                       size_t start_index, size_t end_index, int step,
                                                       vertex start, vertex end, double xmin, double xmax,
                                                       double ymin, double ymax, int x_direction,
                                                       int y_direction, double tolerance,
                                                       blend_polygon_path_candidate *candidate)
{
    size_t index;
    size_t guard = 0;

    if (poly == NULL || candidate == NULL || poly->vertices == NULL || n_input == 0 ||
        (step != 1 && step != -1)) {
        return FAIL;
    }

    candidate->valid = 0;
    candidate->count = 0;

    if (start_index == end_index) {
        candidate->valid = 1;
        return SUCCESS;
    }

    index = start_index;
    while (index != end_index) {
        vertex point;

        if (step > 0) {
            index = (index + 1) % n_input;
        }
        else {
            index = (index == 0) ? n_input - 1 : index - 1;
        }
        if (index == end_index) {
            break;
        }

        point = poly->vertices[index];
        if (blend_polygon_point_on_bbox_side(point, xmin, xmax, ymin, ymax, tolerance)) {
            continue;
        }
        if (!blend_polygon_value_between_open(point.x, start.x, end.x, tolerance) ||
            !blend_polygon_value_between_open(point.y, start.y, end.y, tolerance)) {
            return SUCCESS;
        }

        candidate->indices[candidate->count] = index;
        candidate->points[candidate->count] = point;
        candidate->count++;

        guard++;
        if (guard > n_input) {
            return FAIL;
        }
    }

    if (candidate->count == 0) {
        if (blend_polygon_vertices_equal(start, end, tolerance) ||
            ((double)x_direction * (end.x - start.x) > tolerance &&
             (double)y_direction * (end.y - start.y) > tolerance)) {
            candidate->valid = 1;
        }
        return SUCCESS;
    }

    {
        vertex *chain = NULL;
        size_t i;
        int ok;

        chain = (vertex *)calloc(candidate->count + 2, sizeof(*chain));
        if (chain == NULL) {
            BLEND_Report(BLEND_MSG_ERROR, "polygon: could not allocate traversal path chain\n");
            return FAIL;
        }
        chain[0] = start;
        for (i = 0; i < candidate->count; i++) {
            chain[i + 1] = candidate->points[i];
        }
        chain[candidate->count + 1] = end;

        ok = blend_polygon_chain_is_strictly_monotone(chain, candidate->count + 2,
                                                      x_direction, y_direction, tolerance);
        free(chain);
        candidate->valid = ok;
    }

    return SUCCESS;
}

static int blend_polygon_nonbbox_vertex_count(const polygon *poly, size_t n_input,
                                              double xmin, double xmax,
                                              double ymin, double ymax,
                                              double tolerance, size_t *count)
{
    size_t i;

    if (poly == NULL || count == NULL) {
        return FAIL;
    }

    *count = 0;
    for (i = 0; i < n_input; i++) {
        if (!blend_polygon_point_on_bbox_side(poly->vertices[i], xmin, xmax, ymin, ymax, tolerance)) {
            *count += 1;
        }
    }

    return SUCCESS;
}

static int blend_polygon_choose_strict_path_candidates(const polygon *poly, size_t n_input,
                                                       double xmin, double xmax,
                                                       double ymin, double ymax,
                                                       double tolerance,
                                                       blend_polygon_path_candidate candidates[4][2],
                                                       int selected[4])
{
    int a, b, c, d;
    size_t nonbbox_count = 0;

    if (blend_polygon_nonbbox_vertex_count(poly, n_input, xmin, xmax, ymin, ymax,
                                           tolerance, &nonbbox_count) != SUCCESS) {
        return FAIL;
    }

    for (a = 0; a < 2; a++) {
        for (b = 0; b < 2; b++) {
            for (c = 0; c < 2; c++) {
                for (d = 0; d < 2; d++) {
                    int choices[4];
                    unsigned char *assigned = NULL;
                    size_t assigned_count = 0;
                    size_t i;
                    int sector;
                    int ok = 1;

                    choices[0] = a;
                    choices[1] = b;
                    choices[2] = c;
                    choices[3] = d;

                    for (sector = 0; sector < 4; sector++) {
                        if (!candidates[sector][choices[sector]].valid) {
                            ok = 0;
                            break;
                        }
                    }
                    if (!ok) {
                        continue;
                    }

                    assigned = (unsigned char *)calloc(n_input == 0 ? 1 : n_input, sizeof(*assigned));
                    if (assigned == NULL) {
                        BLEND_Report(BLEND_MSG_ERROR, "polygon: could not allocate traversal path assignments\n");
                        return FAIL;
                    }

                    for (sector = 0; sector < 4 && ok; sector++) {
                        blend_polygon_path_candidate *candidate = &candidates[sector][choices[sector]];

                        for (i = 0; i < candidate->count; i++) {
                            size_t index = candidate->indices[i];

                            if (assigned[index]) {
                                ok = 0;
                                break;
                            }
                            assigned[index] = 1;
                            assigned_count++;
                        }
                    }

                    if (ok && assigned_count == nonbbox_count) {
                        selected[0] = choices[0];
                        selected[1] = choices[1];
                        selected[2] = choices[2];
                        selected[3] = choices[3];
                        free(assigned);
                        return SUCCESS;
                    }

                    free(assigned);
                }
            }
        }
    }

    return FAIL;
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

int blend_polygon_read(const char *path, polygon *poly)
{
    FILE *fp;
    char line[1024];
    long line_number = 0;
    vertex *vertices = NULL;
    size_t count = 0, capacity = 0;
    polygon tmp;

    if (path == NULL || poly == NULL) {
        return FAIL;
    }

    fp = fopen(path, "r");
    if (fp == NULL) {
        BLEND_Report(BLEND_MSG_ERROR, "%s: %s\n", path, strerror(errno));
        return FAIL;
    }

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

        text = blend_polygon_trim_line(line);
        if (*text == '\0') {
            continue;
        }

        fields = sscanf(text, "%lf %lf %1s", &x, &y, extra);
        if (fields != 2 || !isfinite(x) || !isfinite(y)) {
            BLEND_Report(BLEND_MSG_ERROR, "polygon: %s line %ld must contain two finite columns: x y\n",
                         path, line_number);
            fclose(fp);
            free(vertices);
            return FAIL;
        }

        if (blend_polygon_append_vertex(&vertices, &count, &capacity, x, y) != SUCCESS) {
            fclose(fp);
            free(vertices);
            return FAIL;
        }
    }

    if (ferror(fp)) {
        BLEND_Report(BLEND_MSG_ERROR, "%s: %s\n", path, strerror(errno));
        fclose(fp);
        free(vertices);
        return FAIL;
    }
    fclose(fp);

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

int blend_polygon_write(const char *path, const polygon *poly)
{
    FILE *fp;
    size_t i;

    if (path == NULL || blend_polygon_validate(poly) != SUCCESS) {
        return FAIL;
    }

    fp = fopen(path, "w");
    if (fp == NULL) {
        BLEND_Report(BLEND_MSG_ERROR, "%s: %s\n", path, strerror(errno));
        return FAIL;
    }

    for (i = 0; i < poly->n_vertices; i++) {
        if (fprintf(fp, "%.17g %.17g\n", poly->vertices[i].x, poly->vertices[i].y) < 0) {
            BLEND_Report(BLEND_MSG_ERROR, "%s: %s\n", path, strerror(errno));
            fclose(fp);
            return FAIL;
        }
    }

    if (fclose(fp) != 0) {
        BLEND_Report(BLEND_MSG_ERROR, "%s: %s\n", path, strerror(errno));
        return FAIL;
    }

    return SUCCESS;
}

int blend_polygon_validate(const polygon *poly)
{
    double tolerance;
    double area;
    int is_simple = 0;
    size_t i, j;
    size_t unique_count = 0;

    if (poly == NULL || poly->vertices == NULL || poly->n_vertices < 3) {
        BLEND_Report(BLEND_MSG_ERROR, "polygon: at least three vertices are required\n");
        return FAIL;
    }

    tolerance = blend_polygon_coordinate_tolerance(poly);
    for (i = 0; i < poly->n_vertices; i++) {
        if (!isfinite(poly->vertices[i].x) || !isfinite(poly->vertices[i].y)) {
            BLEND_Report(BLEND_MSG_ERROR, "polygon: vertex coordinates must be finite\n");
            return FAIL;
        }
    }

    for (i = 0; i < poly->n_vertices; i++) {
        int seen = 0;

        for (j = 0; j < i; j++) {
            if (blend_polygon_vertices_equal(poly->vertices[i], poly->vertices[j], tolerance)) {
                seen = 1;
                break;
            }
        }
        if (!seen) {
            unique_count++;
        }
    }

    if (unique_count < 3) {
        BLEND_Report(BLEND_MSG_ERROR, "polygon: at least three unique vertices are required\n");
        return FAIL;
    }

    area = fabs(blend_polygon_signed_area(poly));
    if (area <= tolerance * tolerance) {
        BLEND_Report(BLEND_MSG_ERROR, "polygon: polygon area is zero\n");
        return FAIL;
    }

    if (blend_polygon_is_simple(poly, &is_simple) != SUCCESS || !is_simple) {
        BLEND_Report(BLEND_MSG_ERROR, "polygon: polygon edges must not self-intersect\n");
        return FAIL;
    }

    return SUCCESS;
}

int blend_polygon_is_simple(const polygon *poly, int *is_simple)
{
    double tolerance;
    size_t n_edges;
    size_t i, j;

    if (is_simple == NULL) {
        return FAIL;
    }
    *is_simple = 0;

    if (poly == NULL || poly->vertices == NULL || poly->n_vertices < 3) {
        return FAIL;
    }

    tolerance = blend_polygon_coordinate_tolerance(poly);
    n_edges = poly->n_vertices;
    if (poly->n_vertices > 3 &&
        blend_polygon_vertices_equal(poly->vertices[0], poly->vertices[poly->n_vertices - 1], tolerance)) {
        n_edges = poly->n_vertices - 1;
    }

    if (n_edges < 3) {
        return FAIL;
    }

    for (i = 0; i < n_edges; i++) {
        vertex a = poly->vertices[i];
        vertex b = poly->vertices[(i + 1) % n_edges];

        if (!isfinite(a.x) || !isfinite(a.y) || !isfinite(b.x) || !isfinite(b.y)) {
            return FAIL;
        }
        if (blend_polygon_vertices_equal(a, b, tolerance)) {
            return SUCCESS;
        }

        for (j = i + 1; j < n_edges; j++) {
            vertex c;
            vertex d;

            if (blend_polygon_edges_are_adjacent(i, j, n_edges)) {
                continue;
            }

            c = poly->vertices[j];
            d = poly->vertices[(j + 1) % n_edges];
            if (!isfinite(c.x) || !isfinite(c.y) || !isfinite(d.x) || !isfinite(d.y)) {
                return FAIL;
            }

            if (blend_polygon_segments_intersect(a, b, c, d, tolerance)) {
                return SUCCESS;
            }
        }
    }

    *is_simple = 1;
    return SUCCESS;
}

int blend_polygon_bounds(const polygon *poly, double *xmin, double *xmax, double *ymin, double *ymax)
{
    size_t i;

    if (poly == NULL || poly->vertices == NULL || poly->n_vertices == 0 ||
        xmin == NULL || xmax == NULL || ymin == NULL || ymax == NULL) {
        return FAIL;
    }

    *xmin = poly->vertices[0].x;
    *xmax = poly->vertices[0].x;
    *ymin = poly->vertices[0].y;
    *ymax = poly->vertices[0].y;

    for (i = 0; i < poly->n_vertices; i++) {
        if (!isfinite(poly->vertices[i].x) || !isfinite(poly->vertices[i].y)) {
            BLEND_Report(BLEND_MSG_ERROR, "polygon: vertex coordinates must be finite\n");
            return FAIL;
        }
        if (poly->vertices[i].x < *xmin) *xmin = poly->vertices[i].x;
        if (poly->vertices[i].x > *xmax) *xmax = poly->vertices[i].x;
        if (poly->vertices[i].y < *ymin) *ymin = poly->vertices[i].y;
        if (poly->vertices[i].y > *ymax) *ymax = poly->vertices[i].y;
    }

    return SUCCESS;
}

int blend_polygon_contains_point(const polygon *poly, double x, double y, int *inside)
{
    int is_inside = 0;
    size_t i, j;

    if (inside == NULL || !isfinite(x) || !isfinite(y)) {
        return FAIL;
    }
    if (blend_polygon_validate(poly) != SUCCESS) {
        return FAIL;
    }

    for (i = 0, j = poly->n_vertices - 1; i < poly->n_vertices; j = i++) {
        const vertex *vi = &poly->vertices[i];
        const vertex *vj = &poly->vertices[j];

        if (blend_polygon_point_on_segment(x, y, vj, vi)) {
            *inside = 1;
            return SUCCESS;
        }
        if (((vi->y > y) != (vj->y > y)) &&
            (x < (vj->x - vi->x) * (y - vi->y) / (vj->y - vi->y) + vi->x)) {
            is_inside = !is_inside;
        }
    }

    *inside = is_inside;
    return SUCCESS;
}

int blend_polygon_is_closed(const polygon *poly, int *is_closed)
{
    double tolerance;

    if (poly == NULL || poly->vertices == NULL || poly->n_vertices < 2 || is_closed == NULL) {
        return FAIL;
    }

    tolerance = blend_polygon_coordinate_tolerance(poly);
    *is_closed = blend_polygon_vertices_equal(poly->vertices[0], poly->vertices[poly->n_vertices - 1], tolerance);
    return SUCCESS;
}

int blend_polygon_close(polygon *poly)
{
    int is_closed = 0;
    vertex *vertices;

    if (blend_polygon_validate(poly) != SUCCESS) {
        return FAIL;
    }
    if (blend_polygon_is_closed(poly, &is_closed) != SUCCESS) {
        return FAIL;
    }
    if (is_closed) {
        return SUCCESS;
    }

    vertices = (vertex *)realloc(poly->vertices, (poly->n_vertices + 1) * sizeof(*poly->vertices));
    if (vertices == NULL) {
        BLEND_Report(BLEND_MSG_ERROR, "polygon: could not close polygon\n");
        return FAIL;
    }

    poly->vertices = vertices;
    poly->vertices[poly->n_vertices] = poly->vertices[0];
    poly->n_vertices++;
    return SUCCESS;
}

int blend_polygon_map_to_grid(const polygon *src, double xmin, double xmax, double ymin, double ymax,
                              int nx, int ny, polygon *dst)
{
    polygon tmp;
    size_t i;

    if (src == NULL || dst == NULL || src == dst || nx < 2 || ny < 2 ||
        !isfinite(xmin) || !isfinite(xmax) || !isfinite(ymin) || !isfinite(ymax) ||
        xmin >= xmax || ymin >= ymax) {
        return FAIL;
    }
    if (blend_polygon_validate(src) != SUCCESS) {
        return FAIL;
    }

    memset(&tmp, 0, sizeof(tmp));
    if (blend_polygon_alloc(&tmp, src->n_vertices) != SUCCESS) {
        return FAIL;
    }

    for (i = 0; i < src->n_vertices; i++) {
        double x = src->vertices[i].x;
        double y = src->vertices[i].y;
        double gx, gy;

        if (x < xmin || x > xmax || y < ymin || y > ymax) {
            BLEND_Report(BLEND_MSG_ERROR, "polygon: source vertex is outside the target domain\n");
            blend_polygon_free(&tmp);
            return FAIL;
        }

        gx = (x - xmin) * (double)(nx - 1) / (xmax - xmin);
        gy = (y - ymin) * (double)(ny - 1) / (ymax - ymin);

        if (blend_polygon_set_vertex(&tmp, i, gx, gy) != SUCCESS) {
            blend_polygon_free(&tmp);
            return FAIL;
        }
    }

    blend_polygon_free(dst);
    dst->n_vertices = tmp.n_vertices;
    dst->vertices = tmp.vertices;
    return SUCCESS;
}

int blend_polygon_is_xy_monotone(const polygon *poly, int *is_xy_monotone)
{
    if (is_xy_monotone == NULL) {
        return FAIL;
    }
    if (blend_polygon_validate(poly) != SUCCESS) {
        return FAIL;
    }

    *is_xy_monotone = blend_polygon_axis_is_monotone(poly, 0) && blend_polygon_axis_is_monotone(poly, 1);
    return SUCCESS;
}

int blend_polygon_is_xy_monotone_strict(const polygon *poly, int *is_xy_monotone)
{
    vertex *bottom = NULL, *left = NULL, *top = NULL, *right = NULL;
    blend_polygon_path_candidate candidates[4][2];
    double xmin, xmax, ymin, ymax;
    double tolerance;
    size_t n_input;
    size_t bottom_count = 0, left_count = 0, top_count = 0, right_count = 0;
    size_t i;
    int is_regular_xy_monotone = 0;
    int candidates_initialized = 0;
    int selected[4] = {-1, -1, -1, -1};

    if (is_xy_monotone == NULL) {
        return FAIL;
    }
    *is_xy_monotone = 0;

    if (blend_polygon_validate(poly) != SUCCESS) {
        return FAIL;
    }
    if (blend_polygon_is_xy_monotone(poly, &is_regular_xy_monotone) != SUCCESS) {
        return FAIL;
    }
    if (!is_regular_xy_monotone) {
        return SUCCESS;
    }
    if (blend_polygon_bounds(poly, &xmin, &xmax, &ymin, &ymax) != SUCCESS) {
        return FAIL;
    }

    tolerance = blend_polygon_coordinate_tolerance(poly);
    n_input = poly->n_vertices;
    if (n_input > 3 && blend_polygon_vertices_equal(poly->vertices[0], poly->vertices[n_input - 1], tolerance)) {
        n_input--;
    }

    bottom = (vertex *)calloc(n_input, sizeof(*bottom));
    left = (vertex *)calloc(n_input, sizeof(*left));
    top = (vertex *)calloc(n_input, sizeof(*top));
    right = (vertex *)calloc(n_input, sizeof(*right));
    if (bottom == NULL || left == NULL || top == NULL || right == NULL) {
        BLEND_Report(BLEND_MSG_ERROR, "polygon: could not allocate strict xy-monotone work arrays\n");
        free(bottom);
        free(left);
        free(top);
        free(right);
        return FAIL;
    }

    memset(candidates, 0, sizeof(candidates));
    for (i = 0; i < 8; i++) {
        if (blend_polygon_path_candidate_alloc(&candidates[i / 2][i % 2], n_input) != SUCCESS) {
            goto cleanup_fail;
        }
    }
    candidates_initialized = 1;

    for (i = 0; i < n_input; i++) {
        vertex point = poly->vertices[i];

        if (fabs(point.y - ymin) <= tolerance) bottom[bottom_count++] = point;
        if (fabs(point.x - xmin) <= tolerance) left[left_count++] = point;
        if (fabs(point.y - ymax) <= tolerance) top[top_count++] = point;
        if (fabs(point.x - xmax) <= tolerance) right[right_count++] = point;
    }

    if (bottom_count > 0 && left_count > 0 && top_count > 0 && right_count > 0) {
        vertex bottom_left, bottom_right, left_bottom, left_top;
        vertex top_left, top_right, right_bottom, right_top;
        size_t bottom_left_index = 0, bottom_right_index = 0;
        size_t left_bottom_index = 0, left_top_index = 0;
        size_t top_left_index = 0, top_right_index = 0;
        size_t right_bottom_index = 0, right_top_index = 0;

        qsort(bottom, bottom_count, sizeof(*bottom), blend_polygon_compare_vertex);
        qsort(left, left_count, sizeof(*left), blend_polygon_compare_y_asc_x_asc);
        qsort(top, top_count, sizeof(*top), blend_polygon_compare_vertex);
        qsort(right, right_count, sizeof(*right), blend_polygon_compare_y_asc_x_asc);

        bottom_left = bottom[0];
        bottom_right = bottom[bottom_count - 1];
        left_bottom = left[0];
        left_top = left[left_count - 1];
        top_left = top[0];
        top_right = top[top_count - 1];
        right_bottom = right[0];
        right_top = right[right_count - 1];

        if (blend_polygon_find_vertex(poly->vertices, n_input, bottom_left, tolerance,
                                      &bottom_left_index) != SUCCESS ||
            blend_polygon_find_vertex(poly->vertices, n_input, bottom_right, tolerance,
                                      &bottom_right_index) != SUCCESS ||
            blend_polygon_find_vertex(poly->vertices, n_input, left_bottom, tolerance,
                                      &left_bottom_index) != SUCCESS ||
            blend_polygon_find_vertex(poly->vertices, n_input, left_top, tolerance,
                                      &left_top_index) != SUCCESS ||
            blend_polygon_find_vertex(poly->vertices, n_input, top_left, tolerance,
                                      &top_left_index) != SUCCESS ||
            blend_polygon_find_vertex(poly->vertices, n_input, top_right, tolerance,
                                      &top_right_index) != SUCCESS ||
            blend_polygon_find_vertex(poly->vertices, n_input, right_bottom, tolerance,
                                      &right_bottom_index) != SUCCESS ||
            blend_polygon_find_vertex(poly->vertices, n_input, right_top, tolerance,
                                      &right_top_index) != SUCCESS) {
            goto cleanup;
        }

        if (blend_polygon_collect_strict_path_candidate(poly, n_input, bottom_left_index, left_bottom_index,
                                                        1, bottom_left, left_bottom, xmin, xmax, ymin, ymax,
                                                        -1, 1, tolerance, &candidates[0][0]) != SUCCESS ||
            blend_polygon_collect_strict_path_candidate(poly, n_input, bottom_left_index, left_bottom_index,
                                                        -1, bottom_left, left_bottom, xmin, xmax, ymin, ymax,
                                                        -1, 1, tolerance, &candidates[0][1]) != SUCCESS ||
            blend_polygon_collect_strict_path_candidate(poly, n_input, left_top_index, top_left_index,
                                                        1, left_top, top_left, xmin, xmax, ymin, ymax,
                                                        1, 1, tolerance, &candidates[1][0]) != SUCCESS ||
            blend_polygon_collect_strict_path_candidate(poly, n_input, left_top_index, top_left_index,
                                                        -1, left_top, top_left, xmin, xmax, ymin, ymax,
                                                        1, 1, tolerance, &candidates[1][1]) != SUCCESS ||
            blend_polygon_collect_strict_path_candidate(poly, n_input, right_top_index, top_right_index,
                                                        1, right_top, top_right, xmin, xmax, ymin, ymax,
                                                        -1, 1, tolerance, &candidates[2][0]) != SUCCESS ||
            blend_polygon_collect_strict_path_candidate(poly, n_input, right_top_index, top_right_index,
                                                        -1, right_top, top_right, xmin, xmax, ymin, ymax,
                                                        -1, 1, tolerance, &candidates[2][1]) != SUCCESS ||
            blend_polygon_collect_strict_path_candidate(poly, n_input, bottom_right_index, right_bottom_index,
                                                        1, bottom_right, right_bottom, xmin, xmax, ymin, ymax,
                                                        1, 1, tolerance, &candidates[3][0]) != SUCCESS ||
            blend_polygon_collect_strict_path_candidate(poly, n_input, bottom_right_index, right_bottom_index,
                                                        -1, bottom_right, right_bottom, xmin, xmax, ymin, ymax,
                                                        1, 1, tolerance, &candidates[3][1]) != SUCCESS) {
            goto cleanup_fail;
        }

        if (blend_polygon_choose_strict_path_candidates(poly, n_input, xmin, xmax, ymin, ymax,
                                                        tolerance, candidates, selected) == SUCCESS) {
            *is_xy_monotone = 1;
        }
    }

cleanup:
    if (candidates_initialized) {
        for (i = 0; i < 8; i++) {
            blend_polygon_path_candidate_free(&candidates[i / 2][i % 2]);
        }
    }
    free(bottom);
    free(left);
    free(top);
    free(right);
    return SUCCESS;

cleanup_fail:
    if (candidates_initialized) {
        for (i = 0; i < 8; i++) {
            blend_polygon_path_candidate_free(&candidates[i / 2][i % 2]);
        }
    }
    free(bottom);
    free(left);
    free(top);
    free(right);
    return FAIL;
}

static int blend_polygon_copy_candidate_in_source_order(const polygon *src, const polygon *candidate, polygon *dst)
{
    vertex *vertices = NULL;
    size_t count = 0;
    size_t capacity = 0;
    size_t n_src, n_candidate;
    double tolerance;
    size_t i, j;
    polygon tmp;

    if (src == NULL || candidate == NULL || dst == NULL ||
        src->vertices == NULL || candidate->vertices == NULL) {
        return FAIL;
    }

    tolerance = fmax(blend_polygon_coordinate_tolerance(src),
                     blend_polygon_coordinate_tolerance(candidate));
    n_src = src->n_vertices;
    n_candidate = candidate->n_vertices;
    if (n_src > 3 && blend_polygon_vertices_equal(src->vertices[0], src->vertices[n_src - 1], tolerance)) {
        n_src--;
    }
    if (n_candidate > 3 &&
        blend_polygon_vertices_equal(candidate->vertices[0], candidate->vertices[n_candidate - 1], tolerance)) {
        n_candidate--;
    }

    for (i = 0; i < n_src; i++) {
        for (j = 0; j < n_candidate; j++) {
            if (blend_polygon_vertices_equal(src->vertices[i], candidate->vertices[j], tolerance)) {
                if (blend_polygon_append_clean(&vertices, &count, &capacity, src->vertices[i],
                                               tolerance) != SUCCESS) {
                    free(vertices);
                    return FAIL;
                }
                break;
            }
        }
    }

    if (count < 3) {
        free(vertices);
        return FAIL;
    }

    tmp.n_vertices = count;
    tmp.vertices = vertices;
    {
        blend_verbosity verbosity = blend_get_verbosity();
        int validate_status;

        blend_set_verbosity(BLEND_MSG_QUIET);
        validate_status = blend_polygon_validate(&tmp);
        blend_set_verbosity(verbosity);
        if (validate_status != SUCCESS) {
            free(vertices);
            return FAIL;
        }
    }

    blend_polygon_free(dst);
    dst->n_vertices = tmp.n_vertices;
    dst->vertices = tmp.vertices;
    return SUCCESS;
}

int blend_polygon_xy_monotone_envelope(const polygon *src, polygon *dst)
{
    polygon hull_poly = {0};
    vertex *points = NULL;
    vertex *hull = NULL;
    double tolerance;
    size_t n_input, n_points = 0, k = 0;
    size_t i;
    int is_xy_monotone = 0;

    if (src == NULL || dst == NULL || src == dst) {
        return FAIL;
    }
    if (blend_polygon_validate(src) != SUCCESS) {
        return FAIL;
    }

    if (blend_polygon_is_xy_monotone(src, &is_xy_monotone) != SUCCESS) {
        return FAIL;
    }
    if (is_xy_monotone) {
        blend_polygon_free(dst);
        return blend_polygon_copy(src, dst);
    }

    tolerance = blend_polygon_coordinate_tolerance(src);
    n_input = src->n_vertices;
    if (n_input > 3 && blend_polygon_vertices_equal(src->vertices[0], src->vertices[n_input - 1], tolerance)) {
        n_input--;
    }

    points = (vertex *)calloc(n_input, sizeof(*points));
    hull = (vertex *)calloc(2 * n_input, sizeof(*hull));
    if (points == NULL || hull == NULL) {
        BLEND_Report(BLEND_MSG_ERROR, "polygon: could not allocate xy-monotone envelope\n");
        free(points);
        free(hull);
        return FAIL;
    }

    for (i = 0; i < n_input; i++) {
        points[i] = src->vertices[i];
    }
    qsort(points, n_input, sizeof(*points), blend_polygon_compare_vertex);

    for (i = 0; i < n_input; i++) {
        if (n_points == 0 || !blend_polygon_vertices_equal(points[i], points[n_points - 1], tolerance)) {
            points[n_points] = points[i];
            n_points++;
        }
    }

    if (n_points < 3) {
        BLEND_Report(BLEND_MSG_ERROR, "polygon: xy-monotone envelope requires at least three unique vertices\n");
        free(points);
        free(hull);
        return FAIL;
    }

    for (i = 0; i < n_points; i++) {
        while (k >= 2 && blend_polygon_orientation(hull[k - 2], hull[k - 1], points[i]) <= tolerance) {
            k--;
        }
        hull[k] = points[i];
        k++;
    }

    {
        size_t lower_size = k;
        size_t reverse_index;

        for (reverse_index = n_points - 1; reverse_index > 0; reverse_index--) {
            vertex point = points[reverse_index - 1];

            while (k > lower_size && blend_polygon_orientation(hull[k - 2], hull[k - 1], point) <= tolerance) {
                k--;
            }
            hull[k] = point;
            k++;
        }
    }

    if (k > 1) {
        k--;
    }
    if (k < 3) {
        BLEND_Report(BLEND_MSG_ERROR, "polygon: could not build a non-degenerate xy-monotone envelope\n");
        free(points);
        free(hull);
        return FAIL;
    }

    hull_poly.n_vertices = k;
    hull_poly.vertices = hull;
    if (blend_polygon_validate(&hull_poly) != SUCCESS ||
        blend_polygon_is_xy_monotone(&hull_poly, &is_xy_monotone) != SUCCESS ||
        !is_xy_monotone) {
        free(points);
        free(hull);
        return FAIL;
    }

    {
        polygon ordered = {0};

        if (blend_polygon_copy_candidate_in_source_order(src, &hull_poly, &ordered) == SUCCESS &&
            blend_polygon_is_xy_monotone(&ordered, &is_xy_monotone) == SUCCESS &&
            is_xy_monotone) {
            free(hull);
            blend_polygon_free(dst);
            dst->n_vertices = ordered.n_vertices;
            dst->vertices = ordered.vertices;
            free(points);
            return SUCCESS;
        }
        blend_polygon_free(&ordered);
    }

    blend_polygon_free(dst);
    dst->n_vertices = hull_poly.n_vertices;
    dst->vertices = hull_poly.vertices;
    free(points);
    return SUCCESS;
}

int blend_polygon_xy_monotone_envelope_strict(const polygon *src, polygon *dst)
{
    polygon envelope = {0};
    int is_xy_monotone = 0;

    if (src == NULL || dst == NULL || src == dst) {
        return FAIL;
    }
    if (blend_polygon_validate(src) != SUCCESS) {
        return FAIL;
    }

    if (blend_polygon_is_xy_monotone_strict(src, &is_xy_monotone) != SUCCESS) {
        return FAIL;
    }
    if (is_xy_monotone) {
        blend_polygon_free(dst);
        return blend_polygon_copy(src, dst);
    }

    if (blend_polygon_xy_monotone_envelope(src, &envelope) != SUCCESS) {
        return FAIL;
    }
    if (blend_polygon_is_xy_monotone_strict(&envelope, &is_xy_monotone) != SUCCESS ||
        !is_xy_monotone) {
        blend_polygon_free(&envelope);
        return FAIL;
    }

    blend_polygon_free(dst);
    dst->n_vertices = envelope.n_vertices;
    dst->vertices = envelope.vertices;
    return SUCCESS;
}

static int blend_polygon_candidate_is_xy_monotone(const polygon *candidate, int *is_xy_monotone)
{
    blend_verbosity verbosity;
    int status;

    if (candidate == NULL || is_xy_monotone == NULL) {
        return FAIL;
    }

    verbosity = blend_get_verbosity();
    blend_set_verbosity(BLEND_MSG_QUIET);
    status = blend_polygon_is_xy_monotone(candidate, is_xy_monotone);
    blend_set_verbosity(verbosity);

    return status;
}

static int blend_polygon_candidate_is_xy_monotone_strict(const polygon *candidate, int *is_xy_monotone)
{
    blend_verbosity verbosity;
    int status;

    if (candidate == NULL || is_xy_monotone == NULL) {
        return FAIL;
    }

    verbosity = blend_get_verbosity();
    blend_set_verbosity(BLEND_MSG_QUIET);
    status = blend_polygon_is_xy_monotone_strict(candidate, is_xy_monotone);
    blend_set_verbosity(verbosity);

    return status;
}

static int blend_polygon_xy_monotone_piecewise_envelope(const polygon *src, polygon *dst)
{
    polygon tmp = {0};
    vertex *bottom = NULL, *left = NULL, *top = NULL, *right = NULL;
    vertex *bottom_to_left = NULL, *left_to_top = NULL, *right_to_top = NULL, *bottom_to_right = NULL;
    vertex *output = NULL;
    vertex *chain = NULL;
    double xmin, xmax, ymin, ymax;
    double tolerance;
    size_t bottom_count = 0, left_count = 0, top_count = 0, right_count = 0;
    size_t bottom_to_left_count = 0, left_to_top_count = 0, right_to_top_count = 0, bottom_to_right_count = 0;
    size_t output_count = 0, output_capacity = 0, chain_count = 0;
    size_t i;
    int is_xy_monotone = 0;
    int status = FAIL;

    if (src == NULL || dst == NULL || src == dst) {
        return FAIL;
    }
    if (blend_polygon_validate(src) != SUCCESS) {
        return FAIL;
    }

    if (blend_polygon_is_xy_monotone(src, &is_xy_monotone) != SUCCESS) {
        return FAIL;
    }
    if (is_xy_monotone) {
        blend_polygon_free(dst);
        return blend_polygon_copy(src, dst);
    }

    if (blend_polygon_bounds(src, &xmin, &xmax, &ymin, &ymax) != SUCCESS) {
        return FAIL;
    }
    tolerance = blend_polygon_coordinate_tolerance(src);

    bottom = (vertex *)calloc(src->n_vertices, sizeof(*bottom));
    left = (vertex *)calloc(src->n_vertices, sizeof(*left));
    top = (vertex *)calloc(src->n_vertices, sizeof(*top));
    right = (vertex *)calloc(src->n_vertices, sizeof(*right));
    bottom_to_left = (vertex *)calloc(src->n_vertices, sizeof(*bottom_to_left));
    left_to_top = (vertex *)calloc(src->n_vertices, sizeof(*left_to_top));
    right_to_top = (vertex *)calloc(src->n_vertices, sizeof(*right_to_top));
    bottom_to_right = (vertex *)calloc(src->n_vertices, sizeof(*bottom_to_right));
    if (bottom == NULL || left == NULL || top == NULL || right == NULL ||
        bottom_to_left == NULL || left_to_top == NULL || right_to_top == NULL ||
        bottom_to_right == NULL) {
        BLEND_Report(BLEND_MSG_ERROR, "polygon: could not allocate piecewise envelope work arrays\n");
        goto cleanup;
    }

    for (i = 0; i < src->n_vertices; i++) {
        vertex point = src->vertices[i];

        if (fabs(point.y - ymin) <= tolerance) bottom[bottom_count++] = point;
        if (fabs(point.x - xmin) <= tolerance) left[left_count++] = point;
        if (fabs(point.y - ymax) <= tolerance) top[top_count++] = point;
        if (fabs(point.x - xmax) <= tolerance) right[right_count++] = point;
    }

    if (bottom_count == 0 || left_count == 0 || top_count == 0 || right_count == 0) {
        BLEND_Report(BLEND_MSG_DEBUG,
                     "polygon: piecewise envelope could not find all bounding-box sides\n");
        goto cleanup;
    }

    qsort(bottom, bottom_count, sizeof(*bottom), blend_polygon_compare_vertex);
    qsort(left, left_count, sizeof(*left), blend_polygon_compare_y_asc_x_asc);
    qsort(top, top_count, sizeof(*top), blend_polygon_compare_vertex);
    qsort(right, right_count, sizeof(*right), blend_polygon_compare_y_asc_x_asc);

    {
        vertex bottom_left = bottom[0];
        vertex bottom_right = bottom[bottom_count - 1];
        vertex left_bottom = left[0];
        vertex left_top = left[left_count - 1];
        vertex top_left = top[0];
        vertex top_right = top[top_count - 1];
        vertex right_bottom = right[0];
        vertex right_top = right[right_count - 1];

        for (i = 0; i < src->n_vertices; i++) {
            vertex point = src->vertices[i];

            if (blend_polygon_point_on_bbox_side(point, xmin, xmax, ymin, ymax, tolerance)) {
                continue;
            }

            if (blend_polygon_value_between_open(point.x, bottom_left.x, left_bottom.x, tolerance) &&
                blend_polygon_value_between_open(point.y, bottom_left.y, left_bottom.y, tolerance)) {
                bottom_to_left[bottom_to_left_count++] = point;
            }
            else if (blend_polygon_value_between_open(point.x, left_top.x, top_left.x, tolerance) &&
                     blend_polygon_value_between_open(point.y, left_top.y, top_left.y, tolerance)) {
                left_to_top[left_to_top_count++] = point;
            }
            else if (blend_polygon_value_between_open(point.x, right_top.x, top_right.x, tolerance) &&
                     blend_polygon_value_between_open(point.y, right_top.y, top_right.y, tolerance)) {
                right_to_top[right_to_top_count++] = point;
            }
            else if (blend_polygon_value_between_open(point.x, bottom_right.x, right_bottom.x, tolerance) &&
                     blend_polygon_value_between_open(point.y, bottom_right.y, right_bottom.y, tolerance)) {
                bottom_to_right[bottom_to_right_count++] = point;
            }
        }

        qsort(bottom_to_left, bottom_to_left_count, sizeof(*bottom_to_left), blend_polygon_compare_y_asc_x_desc);
        qsort(left_to_top, left_to_top_count, sizeof(*left_to_top), blend_polygon_compare_y_asc_x_asc);
        qsort(right_to_top, right_to_top_count, sizeof(*right_to_top), blend_polygon_compare_y_asc_x_desc);
        qsort(bottom_to_right, bottom_to_right_count, sizeof(*bottom_to_right), blend_polygon_compare_y_asc_x_asc);

        if (blend_polygon_append_chain(&output, &output_count, &output_capacity, bottom, bottom_count,
                                       tolerance) != SUCCESS) {
            goto cleanup;
        }

        if (blend_polygon_corner_chain(bottom_right, right_bottom, bottom_to_right, bottom_to_right_count,
                                       1, 1, tolerance, &chain, &chain_count) != SUCCESS ||
            blend_polygon_append_chain(&output, &output_count, &output_capacity, chain, chain_count,
                                       tolerance) != SUCCESS) {
            goto cleanup;
        }
        free(chain);
        chain = NULL;

        if (blend_polygon_append_chain(&output, &output_count, &output_capacity, right, right_count,
                                       tolerance) != SUCCESS) {
            goto cleanup;
        }

        if (blend_polygon_corner_chain(right_top, top_right, right_to_top, right_to_top_count,
                                       -1, 1, tolerance, &chain, &chain_count) != SUCCESS ||
            blend_polygon_append_chain(&output, &output_count, &output_capacity, chain, chain_count,
                                       tolerance) != SUCCESS) {
            goto cleanup;
        }
        free(chain);
        chain = NULL;

        qsort(top, top_count, sizeof(*top), blend_polygon_compare_x_desc);
        if (blend_polygon_append_chain(&output, &output_count, &output_capacity, top, top_count,
                                       tolerance) != SUCCESS) {
            goto cleanup;
        }

        if (blend_polygon_corner_chain(left_top, top_left, left_to_top, left_to_top_count,
                                       1, 1, tolerance, &chain, &chain_count) != SUCCESS) {
            goto cleanup;
        }
        for (i = chain_count; i > 0; i--) {
            if (blend_polygon_append_clean(&output, &output_count, &output_capacity, chain[i - 1],
                                           tolerance) != SUCCESS) {
                goto cleanup;
            }
        }
        free(chain);
        chain = NULL;

        qsort(left, left_count, sizeof(*left), blend_polygon_compare_y_asc_x_asc);
        for (i = left_count; i > 0; i--) {
            if (blend_polygon_append_clean(&output, &output_count, &output_capacity, left[i - 1],
                                           tolerance) != SUCCESS) {
                goto cleanup;
            }
        }

        if (blend_polygon_corner_chain(bottom_left, left_bottom, bottom_to_left, bottom_to_left_count,
                                       -1, 1, tolerance, &chain, &chain_count) != SUCCESS) {
            goto cleanup;
        }
        for (i = chain_count; i > 0; i--) {
            if (blend_polygon_append_clean(&output, &output_count, &output_capacity, chain[i - 1],
                                           tolerance) != SUCCESS) {
                goto cleanup;
            }
        }
    }

    if (output_count > 1 && blend_polygon_vertices_equal(output[0], output[output_count - 1], tolerance)) {
        output_count--;
    }

    tmp.n_vertices = output_count;
    tmp.vertices = output;
    if (blend_polygon_candidate_is_xy_monotone(&tmp, &is_xy_monotone) != SUCCESS ||
        !is_xy_monotone) {
        BLEND_Report(BLEND_MSG_DEBUG,
                     "polygon: piecewise envelope did not produce an xy-monotone polygon\n");
        tmp.vertices = NULL;
        tmp.n_vertices = 0;
        free(output);
        output = NULL;
        goto cleanup;
    }

    blend_polygon_free(dst);
    dst->n_vertices = tmp.n_vertices;
    dst->vertices = tmp.vertices;
    output = NULL;
    status = SUCCESS;

cleanup:
    free(bottom);
    free(left);
    free(top);
    free(right);
    free(bottom_to_left);
    free(left_to_top);
    free(right_to_top);
    free(bottom_to_right);
    free(chain);
    free(output);
    return status;
}

static int blend_polygon_copy_directed_open(const polygon *src, int reverse, polygon *dst)
{
    double tolerance;
    size_t n_input;
    size_t i;

    if (src == NULL || dst == NULL || src->vertices == NULL || src->n_vertices == 0) {
        return FAIL;
    }

    tolerance = blend_polygon_coordinate_tolerance(src);
    n_input = src->n_vertices;
    if (n_input > 3 && blend_polygon_vertices_equal(src->vertices[0], src->vertices[n_input - 1], tolerance)) {
        n_input--;
    }

    if (blend_polygon_alloc(dst, n_input) != SUCCESS) {
        return FAIL;
    }

    for (i = 0; i < n_input; i++) {
        dst->vertices[i] = reverse ? src->vertices[n_input - 1 - i] : src->vertices[i];
    }

    return SUCCESS;
}

static int blend_polygon_traversal_piecewise_candidate(const polygon *src, int reverse, polygon *dst)
{
    polygon directed = {0};
    polygon tmp = {0};
    vertex *bottom = NULL, *left = NULL, *top = NULL, *right = NULL;
    vertex *sector = NULL;
    vertex *output = NULL;
    vertex *chain = NULL;
    double xmin, xmax, ymin, ymax;
    double tolerance;
    size_t bottom_count = 0, left_count = 0, top_count = 0, right_count = 0;
    size_t output_count = 0, output_capacity = 0, sector_count = 0, chain_count = 0;
    size_t bottom_left_index = 0, bottom_right_index = 0, left_bottom_index = 0, left_top_index = 0;
    size_t top_left_index = 0, top_right_index = 0, right_bottom_index = 0, right_top_index = 0;
    size_t i;
    int is_xy_monotone = 0;
    int status = FAIL;

    if (blend_polygon_copy_directed_open(src, reverse, &directed) != SUCCESS) {
        return FAIL;
    }
    if (blend_polygon_bounds(&directed, &xmin, &xmax, &ymin, &ymax) != SUCCESS) {
        goto cleanup;
    }
    tolerance = blend_polygon_coordinate_tolerance(&directed);

    bottom = (vertex *)calloc(directed.n_vertices, sizeof(*bottom));
    left = (vertex *)calloc(directed.n_vertices, sizeof(*left));
    top = (vertex *)calloc(directed.n_vertices, sizeof(*top));
    right = (vertex *)calloc(directed.n_vertices, sizeof(*right));
    sector = (vertex *)calloc(directed.n_vertices, sizeof(*sector));
    if (bottom == NULL || left == NULL || top == NULL || right == NULL || sector == NULL) {
        BLEND_Report(BLEND_MSG_ERROR, "polygon: could not allocate traversal candidate work arrays\n");
        goto cleanup;
    }

    for (i = 0; i < directed.n_vertices; i++) {
        vertex point = directed.vertices[i];

        if (fabs(point.y - ymin) <= tolerance) bottom[bottom_count++] = point;
        if (fabs(point.x - xmin) <= tolerance) left[left_count++] = point;
        if (fabs(point.y - ymax) <= tolerance) top[top_count++] = point;
        if (fabs(point.x - xmax) <= tolerance) right[right_count++] = point;
    }

    if (bottom_count == 0 || left_count == 0 || top_count == 0 || right_count == 0) {
        goto cleanup;
    }

    qsort(bottom, bottom_count, sizeof(*bottom), blend_polygon_compare_vertex);
    qsort(left, left_count, sizeof(*left), blend_polygon_compare_y_asc_x_asc);
    qsort(top, top_count, sizeof(*top), blend_polygon_compare_vertex);
    qsort(right, right_count, sizeof(*right), blend_polygon_compare_y_asc_x_asc);

    if (blend_polygon_find_vertex(directed.vertices, directed.n_vertices, bottom[0], tolerance,
                                  &bottom_left_index) != SUCCESS ||
        blend_polygon_find_vertex(directed.vertices, directed.n_vertices, bottom[bottom_count - 1], tolerance,
                                  &bottom_right_index) != SUCCESS ||
        blend_polygon_find_vertex(directed.vertices, directed.n_vertices, left[0], tolerance,
                                  &left_bottom_index) != SUCCESS ||
        blend_polygon_find_vertex(directed.vertices, directed.n_vertices, left[left_count - 1], tolerance,
                                  &left_top_index) != SUCCESS ||
        blend_polygon_find_vertex(directed.vertices, directed.n_vertices, top[0], tolerance,
                                  &top_left_index) != SUCCESS ||
        blend_polygon_find_vertex(directed.vertices, directed.n_vertices, top[top_count - 1], tolerance,
                                  &top_right_index) != SUCCESS ||
        blend_polygon_find_vertex(directed.vertices, directed.n_vertices, right[0], tolerance,
                                  &right_bottom_index) != SUCCESS ||
        blend_polygon_find_vertex(directed.vertices, directed.n_vertices, right[right_count - 1], tolerance,
                                  &right_top_index) != SUCCESS) {
        goto cleanup;
    }

    if (blend_polygon_append_chain(&output, &output_count, &output_capacity, bottom, bottom_count,
                                   tolerance) != SUCCESS) {
        goto cleanup;
    }

    if (blend_polygon_collect_traversal_sector(&directed, bottom_right_index, right_bottom_index,
                                               bottom[bottom_count - 1], right[0],
                                               xmin, xmax, ymin, ymax, tolerance,
                                               sector, &sector_count) != SUCCESS ||
        blend_polygon_corner_chain_from_order(bottom[bottom_count - 1], right[0], sector, sector_count,
                                              1, 1, tolerance, &chain, &chain_count) != SUCCESS ||
        blend_polygon_append_chain(&output, &output_count, &output_capacity, chain, chain_count,
                                   tolerance) != SUCCESS) {
        goto cleanup;
    }
    free(chain);
    chain = NULL;

    if (blend_polygon_append_chain(&output, &output_count, &output_capacity, right, right_count,
                                   tolerance) != SUCCESS) {
        goto cleanup;
    }

    if (blend_polygon_collect_traversal_sector(&directed, right_top_index, top_right_index,
                                               right[right_count - 1], top[top_count - 1],
                                               xmin, xmax, ymin, ymax, tolerance,
                                               sector, &sector_count) != SUCCESS ||
        blend_polygon_corner_chain_from_order(right[right_count - 1], top[top_count - 1], sector, sector_count,
                                              -1, 1, tolerance, &chain, &chain_count) != SUCCESS ||
        blend_polygon_append_chain(&output, &output_count, &output_capacity, chain, chain_count,
                                   tolerance) != SUCCESS) {
        goto cleanup;
    }
    free(chain);
    chain = NULL;

    qsort(top, top_count, sizeof(*top), blend_polygon_compare_x_desc);
    if (blend_polygon_append_chain(&output, &output_count, &output_capacity, top, top_count,
                                   tolerance) != SUCCESS) {
        goto cleanup;
    }

    if (blend_polygon_collect_traversal_sector(&directed, top_left_index, left_top_index,
                                               top[top_count - 1], left[left_count - 1],
                                               xmin, xmax, ymin, ymax, tolerance,
                                               sector, &sector_count) != SUCCESS ||
        blend_polygon_corner_chain_from_order(top[top_count - 1], left[left_count - 1], sector, sector_count,
                                              -1, -1, tolerance, &chain, &chain_count) != SUCCESS ||
        blend_polygon_append_chain(&output, &output_count, &output_capacity, chain, chain_count,
                                   tolerance) != SUCCESS) {
        goto cleanup;
    }
    free(chain);
    chain = NULL;

    qsort(left, left_count, sizeof(*left), blend_polygon_compare_y_asc_x_asc);
    for (i = left_count; i > 0; i--) {
        if (blend_polygon_append_clean(&output, &output_count, &output_capacity, left[i - 1],
                                       tolerance) != SUCCESS) {
            goto cleanup;
        }
    }

    if (blend_polygon_collect_traversal_sector(&directed, left_bottom_index, bottom_left_index,
                                               left[0], bottom[0],
                                               xmin, xmax, ymin, ymax, tolerance,
                                               sector, &sector_count) != SUCCESS ||
        blend_polygon_corner_chain_from_order(left[0], bottom[0], sector, sector_count,
                                              1, -1, tolerance, &chain, &chain_count) != SUCCESS ||
        blend_polygon_append_chain(&output, &output_count, &output_capacity, chain, chain_count,
                                   tolerance) != SUCCESS) {
        goto cleanup;
    }

    if (output_count > 1 && blend_polygon_vertices_equal(output[0], output[output_count - 1], tolerance)) {
        output_count--;
    }

    tmp.n_vertices = output_count;
    tmp.vertices = output;
    if (blend_polygon_candidate_is_xy_monotone(&tmp, &is_xy_monotone) != SUCCESS ||
        !is_xy_monotone) {
        tmp.vertices = NULL;
        tmp.n_vertices = 0;
        goto cleanup;
    }

    blend_polygon_free(dst);
    dst->n_vertices = tmp.n_vertices;
    dst->vertices = tmp.vertices;
    output = NULL;
    status = SUCCESS;

cleanup:
    blend_polygon_free(&directed);
    free(bottom);
    free(left);
    free(top);
    free(right);
    free(sector);
    free(chain);
    free(output);
    return status;
}

static int blend_polygon_contains_point_fast(const polygon *poly, double x, double y, int *inside)
{
    int is_inside = 0;
    size_t i, j;

    if (poly == NULL || poly->vertices == NULL || poly->n_vertices < 3 ||
        inside == NULL || !isfinite(x) || !isfinite(y)) {
        return FAIL;
    }

    for (i = 0, j = poly->n_vertices - 1; i < poly->n_vertices; j = i++) {
        const vertex *vi = &poly->vertices[i];
        const vertex *vj = &poly->vertices[j];

        if (blend_polygon_point_on_segment(x, y, vj, vi)) {
            *inside = 1;
            return SUCCESS;
        }
        if (((vi->y > y) != (vj->y > y)) &&
            (x < (vj->x - vi->x) * (y - vi->y) / (vj->y - vi->y) + vi->x)) {
            is_inside = !is_inside;
        }
    }

    *inside = is_inside;
    return SUCCESS;
}

static int blend_polygon_grid_iou(const polygon *a, const polygon *b,
                                  double xmin, double xmax, double ymin, double ymax,
                                  int nx, int ny, double *iou)
{
    int ix, iy;
    size_t intersection_count = 0;
    size_t union_count = 0;

    if (iou == NULL || a == NULL || b == NULL || nx < 1 || ny < 1 ||
        !isfinite(xmin) || !isfinite(xmax) || !isfinite(ymin) || !isfinite(ymax) ||
        xmin >= xmax || ymin >= ymax) {
        return FAIL;
    }

    for (iy = 0; iy < ny; iy++) {
        double y = ny == 1 ? ymin : ymin + (double)iy * (ymax - ymin) / (double)(ny - 1);

        for (ix = 0; ix < nx; ix++) {
            double x = nx == 1 ? xmin : xmin + (double)ix * (xmax - xmin) / (double)(nx - 1);
            int inside_a = 0;
            int inside_b = 0;

            if (blend_polygon_contains_point_fast(a, x, y, &inside_a) != SUCCESS ||
                blend_polygon_contains_point_fast(b, x, y, &inside_b) != SUCCESS) {
                return FAIL;
            }

            if (inside_a || inside_b) {
                union_count++;
                if (inside_a && inside_b) {
                    intersection_count++;
                }
            }
        }
    }

    if (union_count == 0) {
        return FAIL;
    }

    *iou = (double)intersection_count / (double)union_count;
    return SUCCESS;
}

static int blend_polygon_keep_highest_grid_iou_candidate(const polygon *src, polygon *best,
                                                         const polygon *candidate,
                                                         double xmin, double xmax,
                                                         double ymin, double ymax,
                                                         int nx, int ny, double *best_iou)
{
    polygon ordered = {0};
    double iou;
    int is_xy_monotone = 0;

    if (candidate == NULL || candidate->vertices == NULL || candidate->n_vertices < 3 ||
        best == NULL || best_iou == NULL || src == NULL) {
        return FAIL;
    }

    if (blend_polygon_copy_candidate_in_source_order(src, candidate, &ordered) != SUCCESS) {
        return SUCCESS;
    }
    if (blend_polygon_candidate_is_xy_monotone(&ordered, &is_xy_monotone) != SUCCESS ||
        !is_xy_monotone) {
        blend_polygon_free(&ordered);
        return SUCCESS;
    }

    if (blend_polygon_grid_iou(src, &ordered, xmin, xmax, ymin, ymax, nx, ny, &iou) != SUCCESS) {
        blend_polygon_free(&ordered);
        return FAIL;
    }

    if (best->vertices == NULL || iou > *best_iou ||
        (fabs(iou - *best_iou) <= 1.0e-12 &&
         ordered.n_vertices > best->n_vertices)) {
        polygon copy = {0};

        if (blend_polygon_copy(&ordered, &copy) != SUCCESS) {
            blend_polygon_free(&ordered);
            return FAIL;
        }
        blend_polygon_free(best);
        best->n_vertices = copy.n_vertices;
        best->vertices = copy.vertices;
        *best_iou = iou;
    }

    blend_polygon_free(&ordered);
    return SUCCESS;
}

static int blend_polygon_keep_highest_strict_grid_iou_candidate(const polygon *src, polygon *best,
                                                                const polygon *candidate,
                                                                double xmin, double xmax,
                                                                double ymin, double ymax,
                                                                int nx, int ny, double *best_iou)
{
    polygon ordered = {0};
    double iou;
    int is_xy_monotone = 0;

    if (candidate == NULL || candidate->vertices == NULL || candidate->n_vertices < 3 ||
        best == NULL || best_iou == NULL || src == NULL) {
        return FAIL;
    }

    if (blend_polygon_copy_candidate_in_source_order(src, candidate, &ordered) != SUCCESS) {
        return SUCCESS;
    }
    if (blend_polygon_candidate_is_xy_monotone_strict(&ordered, &is_xy_monotone) != SUCCESS ||
        !is_xy_monotone) {
        blend_polygon_free(&ordered);
        return SUCCESS;
    }

    if (blend_polygon_grid_iou(src, &ordered, xmin, xmax, ymin, ymax, nx, ny, &iou) != SUCCESS) {
        blend_polygon_free(&ordered);
        return FAIL;
    }

    if (best->vertices == NULL || iou > *best_iou ||
        (fabs(iou - *best_iou) <= 1.0e-12 &&
         ordered.n_vertices > best->n_vertices)) {
        polygon copy = {0};

        if (blend_polygon_copy(&ordered, &copy) != SUCCESS) {
            blend_polygon_free(&ordered);
            return FAIL;
        }
        blend_polygon_free(best);
        best->n_vertices = copy.n_vertices;
        best->vertices = copy.vertices;
        *best_iou = iou;
    }

    blend_polygon_free(&ordered);
    return SUCCESS;
}

int blend_polygon_xy_monotone_best_piecewise_envelope(const polygon *src, polygon *dst,
                                                      double xmin, double xmax,
                                                      double ymin, double ymax,
                                                      int nx, int ny)
{
    polygon candidate = {0};
    polygon best = {0};
    double best_iou = 0.0;
    int is_xy_monotone = 0;
    int reverse;
    int candidate_count = 0;
    double start;
    size_t sample_count;

    if (src == NULL || dst == NULL || src == dst || nx < 1 || ny < 1 ||
        !isfinite(xmin) || !isfinite(xmax) || !isfinite(ymin) || !isfinite(ymax) ||
        xmin >= xmax || ymin >= ymax) {
        return FAIL;
    }
    if (blend_polygon_validate(src) != SUCCESS) {
        return FAIL;
    }

    if (blend_polygon_is_xy_monotone(src, &is_xy_monotone) != SUCCESS) {
        return FAIL;
    }
    if (is_xy_monotone) {
        blend_polygon_free(dst);
        return blend_polygon_copy(src, dst);
    }

    start = blend_elapsed_seconds();
    sample_count = (size_t)nx * (size_t)ny;
    BLEND_Report(BLEND_MSG_TIMING,
                 "polygon: evaluating best IoU xy-monotone envelope candidates on %zu grid samples (%d x %d)\n",
                 sample_count, nx, ny);

    for (reverse = 0; reverse <= 1; reverse++) {
        if (blend_polygon_traversal_piecewise_candidate(src, reverse, &candidate) == SUCCESS) {
            candidate_count++;
            if (blend_polygon_keep_highest_grid_iou_candidate(src, &best, &candidate,
                                                              xmin, xmax, ymin, ymax, nx, ny,
                                                              &best_iou) != SUCCESS) {
                blend_polygon_free(&candidate);
                blend_polygon_free(&best);
                return FAIL;
            }
            blend_polygon_free(&candidate);
        }
    }

    if (blend_polygon_xy_monotone_piecewise_envelope(src, &candidate) == SUCCESS) {
        candidate_count++;
        if (blend_polygon_keep_highest_grid_iou_candidate(src, &best, &candidate,
                                                          xmin, xmax, ymin, ymax, nx, ny,
                                                          &best_iou) != SUCCESS) {
            blend_polygon_free(&candidate);
            blend_polygon_free(&best);
            return FAIL;
        }
        blend_polygon_free(&candidate);
    }

    if (blend_polygon_xy_monotone_envelope(src, &candidate) == SUCCESS) {
        candidate_count++;
        if (blend_polygon_keep_highest_grid_iou_candidate(src, &best, &candidate,
                                                          xmin, xmax, ymin, ymax, nx, ny,
                                                          &best_iou) != SUCCESS) {
            blend_polygon_free(&candidate);
            blend_polygon_free(&best);
            return FAIL;
        }
        blend_polygon_free(&candidate);
    }

    if (best.vertices == NULL) {
        return FAIL;
    }

    blend_polygon_free(dst);
    dst->n_vertices = best.n_vertices;
    dst->vertices = best.vertices;
    BLEND_Report(BLEND_MSG_TIMING,
                 "polygon: selected best IoU xy-monotone envelope from %d candidates in %.3f s (IoU = %.6f)\n",
                 candidate_count, blend_elapsed_seconds() - start, best_iou);
    return SUCCESS;
}

int blend_polygon_xy_monotone_best_piecewise_envelope_strict(const polygon *src, polygon *dst,
                                                             double xmin, double xmax,
                                                             double ymin, double ymax,
                                                             int nx, int ny)
{
    polygon candidate = {0};
    polygon best = {0};
    double best_iou = 0.0;
    int is_xy_monotone = 0;
    int reverse;
    int candidate_count = 0;
    double start;
    size_t sample_count;

    if (src == NULL || dst == NULL || src == dst || nx < 1 || ny < 1 ||
        !isfinite(xmin) || !isfinite(xmax) || !isfinite(ymin) || !isfinite(ymax) ||
        xmin >= xmax || ymin >= ymax) {
        return FAIL;
    }
    if (blend_polygon_validate(src) != SUCCESS) {
        return FAIL;
    }

    if (blend_polygon_is_xy_monotone_strict(src, &is_xy_monotone) != SUCCESS) {
        return FAIL;
    }
    if (is_xy_monotone) {
        blend_polygon_free(dst);
        return blend_polygon_copy(src, dst);
    }

    start = blend_elapsed_seconds();
    sample_count = (size_t)nx * (size_t)ny;
    BLEND_Report(BLEND_MSG_TIMING,
                 "polygon: evaluating strict best IoU xy-monotone envelope candidates on %zu grid samples (%d x %d)\n",
                 sample_count, nx, ny);

    for (reverse = 0; reverse <= 1; reverse++) {
        if (blend_polygon_traversal_piecewise_candidate(src, reverse, &candidate) == SUCCESS) {
            candidate_count++;
            if (blend_polygon_keep_highest_strict_grid_iou_candidate(src, &best, &candidate,
                                                                     xmin, xmax, ymin, ymax, nx, ny,
                                                                     &best_iou) != SUCCESS) {
                blend_polygon_free(&candidate);
                blend_polygon_free(&best);
                return FAIL;
            }
            blend_polygon_free(&candidate);
        }
    }

    if (blend_polygon_xy_monotone_piecewise_envelope(src, &candidate) == SUCCESS) {
        candidate_count++;
        if (blend_polygon_keep_highest_strict_grid_iou_candidate(src, &best, &candidate,
                                                                 xmin, xmax, ymin, ymax, nx, ny,
                                                                 &best_iou) != SUCCESS) {
            blend_polygon_free(&candidate);
            blend_polygon_free(&best);
            return FAIL;
        }
        blend_polygon_free(&candidate);
    }

    if (blend_polygon_xy_monotone_envelope_strict(src, &candidate) == SUCCESS) {
        candidate_count++;
        if (blend_polygon_keep_highest_strict_grid_iou_candidate(src, &best, &candidate,
                                                                 xmin, xmax, ymin, ymax, nx, ny,
                                                                 &best_iou) != SUCCESS) {
            blend_polygon_free(&candidate);
            blend_polygon_free(&best);
            return FAIL;
        }
        blend_polygon_free(&candidate);
    }

    if (best.vertices == NULL) {
        return FAIL;
    }

    blend_polygon_free(dst);
    dst->n_vertices = best.n_vertices;
    dst->vertices = best.vertices;
    BLEND_Report(BLEND_MSG_TIMING,
                 "polygon: selected strict best IoU xy-monotone envelope from %d candidates in %.3f s (IoU = %.6f)\n",
                 candidate_count, blend_elapsed_seconds() - start, best_iou);
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
