/*
* Author: Rasheed Ajala
* Contact: abioyeajala@gmail.com
*/

#ifndef BLEND_BOUNDARY_H
#define BLEND_BOUNDARY_H

#include "blend_window.h"

/* Struct for the permuted vertices */
typedef struct permuted_vertex {
    double **bottom_vertices;
    int row_size_bv;
    double **bottom_to_left_vertices;
    int row_size_btlv;
    double **left_vertices;
    int row_size_lv;
    double **left_to_top_vertices;
    int row_size_lttv;
    double **top_vertices;
    int row_size_tv;
    double **right_to_top_vertices;
    int row_size_rttv;
    double **right_vertices;
    int row_size_rv;
    double **bottom_to_right_vertices;
    int row_size_btrv;
} permuted_vertex;

/* Sort vertices array and update dimensions */
int unique_vertices(double **array, int row_length);

/* Sort array in ascending order in x */
int ascend_x(double **array, int row_length);

/* Sort array in descending order in x */
int descend_x(double **array, int row_length);

/* Sort array in ascending order in y */
int ascend_y(double **array, int row_length);

/* Sort array in descending order in y*/
int descend_y(double **array, int row_length);

/* Bounding box error checks */
int bounding_boxcheck(window *data);

/* Assemble boundary vertex groups */
int assemble_bv(window *data, permuted_vertex *vertex);
int assemble_lv(window *data, permuted_vertex *vertex);
int assemble_tv(window *data, permuted_vertex *vertex);
int assemble_rv(window *data, permuted_vertex *vertex);
int assemble_btlv(window *data, permuted_vertex *vertex);
int assemble_lttv(window *data, permuted_vertex *vertex);
int assemble_rttv(window *data, permuted_vertex *vertex);
int assemble_btrv(window *data, permuted_vertex *vertex);

/* Hanging sweep */
int hanging_sweep(int j, double *vector, double **prev_array, int prev_row, double **sub_array, int sub_row,
                  double **hang_array, int hang_len, assembly_type step_type);

/* Assemble row and column stepping vectors */
int minimum_rowstep(window *data, permuted_vertex *vertex);
int maximum_rowstep(window *data, permuted_vertex *vertex);
int minimum_columnstep(window *data, permuted_vertex *vertex);
int maximum_columnstep(window *data, permuted_vertex *vertex);

/* Assemble the boundary and stepping vectors for window function */
int boundary_assembly(window *data, permuted_vertex *vertex);

/* Clear boundary classification vertices */
void blend_permuted_vertex_free(permuted_vertex *vertex);

/* Clear polygon and stepping vectors assigned to a window */
void blend_window_boundary_clear(window *data);

#endif /* BLEND_BOUNDARY_H */
