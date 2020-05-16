#ifndef BLEND_H
#define BLEND_H

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>
#include <math.h>

/* Constants */
#ifndef M_PI
	/* Defines pi */
	#define M_PI 3.14159265358979323846
#endif

/* Maximum window function length */
#define MAX_WINDOWFUNC_LEN 64

/* Defines a return value of success */
#define SUCCESS 0

/* Defines a return value of failure */
#define FAIL 1

/* Defines function error */
#define WFUNC_ERROR 2.0

/* Currently defined window functions */
#define WFUNC_BOXCAR "boxcar"
#define WFUNC_COSINE "cosine"
#define WFUNC_COSINE_START "scosine"
#define WFUNC_COSINE_END "ecosine"
#define WFUNC_TRAPEZOID "trapezoid"
#define WFUNC_TRAPEZOID_START "strapezoid"
#define WFUNC_TRAPEZOID_END "etrapezoid"

/* Enumerated data type */

/* Supported assembly type */
typedef enum { ROW_STEP = 0, 
	       COLUMN_STEP } assembly_type_t;

/* Data structures */

/* Struct for the permuted vertices */
typedef struct permuted_vertex_t {
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
} permuted_vertex_t;

/* Data structure for boundary and vector asssembly */
typedef struct window_t {
    /* row length of vertices array */
    int row_size;
    /* Polygon boundaries (2D array) - assigned from the model */
    double **vertices;
    /* Assembly mode - Needed for hanging sweep! */
    assembly_type_t amode;
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
    /* Minimum length of row-stepping vectors */
    int minnx;
    /* Minimum length of column-stepping vectors */
    int minny;
    /* Taper ratio: user defined */
    double ratio;
    /* Window function in x: user defined */
    char x_function[MAX_WINDOWFUNC_LEN];
    /* Window function in y: user defined */
    char y_function[MAX_WINDOWFUNC_LEN];
    /* Window function in z: user defined */
    char z_function[MAX_WINDOWFUNC_LEN];
    /* Contribution: blending weights */
    double contribution;
} window_t;

/* Function Declarations */

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
int bounding_boxcheck(window_t *data);

/* Assemble bottom vertices */
int assemble_bv(window_t *data, permuted_vertex_t *vertex);

/* Assemble left vertices */
int assemble_lv(window_t *data, permuted_vertex_t *vertex);

/* Assemble top vertices */
int assemble_tv(window_t *data, permuted_vertex_t *vertex);

/* Assemble right vertices */
int assemble_rv(window_t *data, permuted_vertex_t *vertex);

/* Assemble bottom-to-left-vertices */
int assemble_btlv(window_t *data, permuted_vertex_t *vertex);

/* Assemble left-to-top vertices */
int assemble_lttv(window_t *data, permuted_vertex_t *vertex);

/* Assemble right-to-top vertices */
int assemble_rttv(window_t *data, permuted_vertex_t *vertex);

/* Assemble bottom-to-right vertices */
int assemble_btrv(window_t *data, permuted_vertex_t *vertex);


/* Hanging sweep */
int hanging_sweep(int j, double *vector, double **prev_array, int prev_row, double **sub_array, int sub_row,
                  double **hang_array, int hang_len,  assembly_type_t step_type);

/* Assemble minimum row-stepping vectors: nnx1 */
int minimum_rowstep(window_t *data, permuted_vertex_t *vertex);

/* Assemble maximum row-stepping vectors: nnx2 */
int maximum_rowstep(window_t *data, permuted_vertex_t *vertex);

/* Assemble minimum column-stepping vectors: nny1 */
int minimum_columnstep(window_t *data, permuted_vertex_t *vertex);

/* Assemble maximum column-stepping vectors: nny2 */
int maximum_columnstep(window_t *data, permuted_vertex_t *vertex);

/* Assemble the boundary and stepping vectors for window function */
int boundary_assembly(window_t *data, permuted_vertex_t *vertex);

/* Window/Taper function */
double window_function(int x, int x1, int x2, int n, int nmin, double r, char func[MAX_WINDOWFUNC_LEN]);

/* boxcar function */
double boxcar(void);

/* cosine function */
double cosine(int x, int n, int nmax, int nmin, double r1, double r2);

/* scosine function */
double scosine(int x, int n, int nmax, int nmin, double r);

/* ecosine function */
double ecosine(int x, int n, int nmax, int nmin, double r);

/* trapezoid function */
double trapezoid(int x, int n, int nmax, int nmin, double r1, double r2);

/* strapezoid function */
double strapezoid(int x, int n, int nmax, int nmin, double r);

/* etrapezoid function */
double etrapezoid(int x, int n, int nmax, int nmin, double r);

/* Embedding contribution: Blending weight */
int embedding_contribution(int x, int y, int z, window_t *data);

#endif /* BLEND_H */
