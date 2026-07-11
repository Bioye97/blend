/*
* Author: Rasheed Ajala
* Contact: abioyeajala@gmail.com
*/

#ifndef BLEND_H
#define BLEND_H

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>

#include "blend_polygon.h"

/* Constants */
#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

/* Error codes */
typedef enum blend_error_code {
    SUCCESS = 0,
    FAIL = 1,
    WFUNC_ERROR = 2
} blend_error_code;

/* Window functions */
typedef enum blend_window_function {
    WFUNC_BOXCAR = 0,
    WFUNC_COSINE,
    WFUNC_TRAPEZOID,
    WFUNC_INVALID
} blend_window_function;

/* Message verbosity */
typedef enum blend_verbosity {
    BLEND_MSG_QUIET = 0,
    BLEND_MSG_ERROR,
    BLEND_MSG_WARNING,
    BLEND_MSG_TIMING,
    BLEND_MSG_INFORMATION,
    BLEND_MSG_COMPAT,
    BLEND_MSG_DEBUG
} blend_verbosity;

/* Supported assembly type */
typedef enum { ROW_STEP = 0,
               COLUMN_STEP } assembly_type;

/* Struct for the permuted vertices */
typedef struct permuted_vertex {
    /* Bottom vertices (2D array) */
    double **bottom_vertices;
    /* Dimension of bottom vertices array */
    int row_size_bv;
    /* Bottom-to-left vertices (2D array) */
    double **bottom_to_left_vertices;
    /* Dimension of bottom-to-left vertices array */
    int row_size_btlv;
    /* Left vertices (2D array) */
    double **left_vertices;
    /* Dimension of left vertices array */
    int row_size_lv;
    /* Left-to_top vertices (2D array) */
    double **left_to_top_vertices;
    /* Dimension of left-to-top vertices array */
    int row_size_lttv;
    /* Top vertices (2D array) */
    double **top_vertices;
    /* Dimension of top vertices array */
    int row_size_tv;
    /* Right-to_top vertices (2D array) */
    double **right_to_top_vertices;
    /* Dimension of right-to-top vertices array */
    int row_size_rttv;
    /* Right vertices (2D array) */
    double **right_vertices;
    /* Dimension of right vertices array */
    int row_size_rv;
    /* Bottom-to-right vertices (2D array) */
    double **bottom_to_right_vertices;
    /* Dimension of bottom-to-right vertices array */
    int row_size_btrv;
} permuted_vertex;

/* Data structure for boundary and vector asssembly */
typedef struct window {
    /* row length of vertices array */
    int row_size;
    /* Polygon boundaries (2D array) - assigned from the model */
    double **vertices;
    /* Assembly mode - Needed for hanging sweep */
    assembly_type amode;
    /* Model dimensions */
    int nx;
    int ny;
    int nz;
    /* Row-stepping vectors (1D array): nnx2 >= nnx1 */
    int *nnx1;
    int *nnx2;
    /* Column-stepping vectors (1D array): nny2 >= nny1 */
    int *nny1;
    int *nny2;
    /* Minimum length of all row-stepping vectors, i.e., support length along x-axis: left-to-right */
    int minnx;
    /* Minimum length of all column-stepping vectors, i.e., support length along y-axis: bottom-to-top */
    int minny;
    /* Window function ratios in x: user defined */
    double ratio_x1;
    double ratio_x2;
    /* Window function ratios in y: user defined */
    double ratio_y1;
    double ratio_y2;
    /* Window function ratios in z: user defined */
    double ratio_z1;
    double ratio_z2;
    /* Window function in x: user defined */
    blend_window_function x_function;
    /* Window function in y: user defined */
    blend_window_function y_function;
    /* Window function in z: user defined */
    blend_window_function z_function;
    /* Contribution: blending weights */
    double contribution;
} window;

/* API version string */
const char *blend_api_version(void);

/* Error code description */
const char *blend_error_message(blend_error_code code);

/* Message verbosity name */
const char *blend_verbosity_name(blend_verbosity level);

/* Message verbosity lookup */
int blend_verbosity_from_name(const char *name, blend_verbosity *level);

/* Set active message verbosity */
int blend_set_verbosity(blend_verbosity level);

/* Get active message verbosity */
blend_verbosity blend_get_verbosity(void);

/* Report a message according to the active verbosity */
void BLEND_Report(blend_verbosity level, const char *format, ...);

/* Linear interpolation */
int interpolate_linear(double x0, double y0, double x1, double y1, double x, double *y);

/* Bilinear interpolation */
int interpolate_bilinear(double x0, double x1, double y0, double y1,
                         double q00, double q10, double q01, double q11,
                         double x, double y, double *value);

/* Window function name */
const char *blend_window_function_name(blend_window_function function);

/* Window function lookup */
int blend_window_function_from_name(const char *name, blend_window_function *function);

/* Print supported window function names */
void blend_print_window_function_names(FILE *fp);

/* Assign polygon support to a window */
int blend_window_set_polygon(window *data, const polygon *poly);

/* Clear polygon support assigned to a window */
void blend_window_clear_polygon(window *data);

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

/* Assemble bottom vertices */
int assemble_bv(window *data, permuted_vertex *vertex);

/* Assemble left vertices */
int assemble_lv(window *data, permuted_vertex *vertex);

/* Assemble top vertices */
int assemble_tv(window *data, permuted_vertex *vertex);

/* Assemble right vertices */
int assemble_rv(window *data, permuted_vertex *vertex);

/* Assemble bottom-to-left-vertices */
int assemble_btlv(window *data, permuted_vertex *vertex);

/* Assemble left-to-top vertices */
int assemble_lttv(window *data, permuted_vertex *vertex);

/* Assemble right-to-top vertices */
int assemble_rttv(window *data, permuted_vertex *vertex);

/* Assemble bottom-to-right vertices */
int assemble_btrv(window *data, permuted_vertex *vertex);

/* Hanging sweep */
int hanging_sweep(int j, double *vector, double **prev_array, int prev_row, double **sub_array, int sub_row,
                  double **hang_array, int hang_len, assembly_type step_type);

/* Assemble minimum row-stepping vectors: nnx1 */
int minimum_rowstep(window *data, permuted_vertex *vertex);

/* Assemble maximum row-stepping vectors: nnx2 */
int maximum_rowstep(window *data, permuted_vertex *vertex);

/* Assemble minimum column-stepping vectors: nny1 */
int minimum_columnstep(window *data, permuted_vertex *vertex);

/* Assemble maximum column-stepping vectors: nny2 */
int maximum_columnstep(window *data, permuted_vertex *vertex);

/* Assemble the boundary and stepping vectors for window function */
int boundary_assembly(window *data, permuted_vertex *vertex);

/* Clear boundary classification vertices */
void blend_permuted_vertex_free(permuted_vertex *vertex);

/* Clear polygon and stepping vectors assigned to a window */
void blend_window_boundary_clear(window *data);

/* Window/Taper function */
double window_function(int x, int x1, int x2, int n, int nmin,
                       double r1, double r2, blend_window_function func);

/* boxcar function */
double boxcar(void);

/* cosine function */
double cosine(int x, int n, int nmax, int nmin, double r1, double r2);

/* trapezoid function */
double trapezoid(int x, int n, int nmax, int nmin, double r1, double r2);

/* Embedding contribution for 1D case: Blending weight */
int embedding_contribution1d(int x, window *data);

/* Embedding contribution for 2D case: Blending weight */
int embedding_contribution2d(int x, int y, window *data);

/* Embedding contribution for 3D case: Blending weight */
int embedding_contribution3d(int x, int y, int z, window *data);

#endif /* BLEND_H */
