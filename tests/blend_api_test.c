#include "blend.h"

#include <float.h>

#define ASSERT_TRUE(expr) do { \
    if (!(expr)) { \
        fprintf(stderr, "Assertion failed at %s:%d: %s\n", __FILE__, __LINE__, #expr); \
        return FAIL; \
    } \
} while (0)

#define ASSERT_EQ_INT(actual, expected) do { \
    int actual_value = (actual); \
    int expected_value = (expected); \
    if (actual_value != expected_value) { \
        fprintf(stderr, "Assertion failed at %s:%d: expected %d, got %d\n", \
                __FILE__, __LINE__, expected_value, actual_value); \
        return FAIL; \
    } \
} while (0)

#define ASSERT_NEAR(actual, expected, tolerance) do { \
    double actual_value = (actual); \
    double expected_value = (expected); \
    if (fabs(actual_value - expected_value) > (tolerance)) { \
        fprintf(stderr, "Assertion failed at %s:%d: expected %.12f, got %.12f\n", \
                __FILE__, __LINE__, expected_value, actual_value); \
        return FAIL; \
    } \
} while (0)

static double **alloc_points(int rows)
{
    int i;
    double **points = (double **)calloc((size_t)rows, sizeof(double *));

    if (points == NULL) return NULL;

    for (i = 0; i < rows; i++) {
        points[i] = (double *)calloc(2, sizeof(double));
        if (points[i] == NULL) return NULL;
    }

    return points;
}

static void set_point(double **points, int row, double x, double y)
{
    points[row][0] = x;
    points[row][1] = y;
}

static void free_points(double **points, int rows)
{
    int i;

    if (points == NULL) return;
    for (i = 0; i < rows; i++) {
        free(points[i]);
    }
    free(points);
}

static void init_rectangle(window_t *data, int nx, int ny)
{
    data->row_size = 4;
    data->nx = nx;
    data->ny = ny;
    data->nz = 5;
    data->ratio = 0.2;
    strcpy(data->x_function, WFUNC_COSINE);
    strcpy(data->y_function, WFUNC_COSINE);
    strcpy(data->z_function, WFUNC_COSINE_START);
    data->vertices = alloc_points(data->row_size);
    set_point(data->vertices, 0, 0, 0);
    set_point(data->vertices, 1, nx - 1, 0);
    set_point(data->vertices, 2, nx - 1, ny - 1);
    set_point(data->vertices, 3, 0, ny - 1);
}

static void init_hexagon(window_t *data)
{
    data->row_size = 6;
    data->nx = 10;
    data->ny = 12;
    data->nz = 5;
    data->ratio = 0.2;
    strcpy(data->x_function, WFUNC_COSINE);
    strcpy(data->y_function, WFUNC_COSINE);
    strcpy(data->z_function, WFUNC_COSINE_START);
    data->vertices = alloc_points(data->row_size);
    set_point(data->vertices, 0, 2, 0);
    set_point(data->vertices, 1, 7, 0);
    set_point(data->vertices, 2, 9, 4);
    set_point(data->vertices, 3, 9, 11);
    set_point(data->vertices, 4, 0, 11);
    set_point(data->vertices, 5, 0, 3);
}

static void init_octagon_with_quadrant_vertices(window_t *data)
{
    data->row_size = 12;
    data->nx = 10;
    data->ny = 10;
    data->nz = 5;
    data->ratio = 0.2;
    strcpy(data->x_function, WFUNC_COSINE);
    strcpy(data->y_function, WFUNC_COSINE);
    strcpy(data->z_function, WFUNC_COSINE_START);
    data->vertices = alloc_points(data->row_size);
    set_point(data->vertices, 0, 3, 0);
    set_point(data->vertices, 1, 6, 0);
    set_point(data->vertices, 2, 8, 1);
    set_point(data->vertices, 3, 9, 3);
    set_point(data->vertices, 4, 9, 6);
    set_point(data->vertices, 5, 8, 8);
    set_point(data->vertices, 6, 6, 9);
    set_point(data->vertices, 7, 3, 9);
    set_point(data->vertices, 8, 1, 8);
    set_point(data->vertices, 9, 0, 6);
    set_point(data->vertices, 10, 0, 3);
    set_point(data->vertices, 11, 1, 1);
}

static int test_sorting_and_unique_vertices(void)
{
    double **points = alloc_points(5);
    int rows;

    ASSERT_TRUE(points != NULL);
    set_point(points, 0, 3, 2);
    set_point(points, 1, 1, 4);
    set_point(points, 2, 2, 1);
    set_point(points, 3, 1, 4);
    set_point(points, 4, 0, 0);

    rows = unique_vertices(points, 5);
    ASSERT_EQ_INT(rows, 4);

    ASSERT_EQ_INT(ascend_x(points, rows), SUCCESS);
    ASSERT_NEAR(points[0][0], 0, DBL_EPSILON);
    ASSERT_NEAR(points[rows - 1][0], 3, DBL_EPSILON);

    ASSERT_EQ_INT(descend_x(points, rows), SUCCESS);
    ASSERT_NEAR(points[0][0], 3, DBL_EPSILON);
    ASSERT_NEAR(points[rows - 1][0], 0, DBL_EPSILON);

    ASSERT_EQ_INT(ascend_y(points, rows), SUCCESS);
    ASSERT_NEAR(points[0][1], 0, DBL_EPSILON);
    ASSERT_NEAR(points[rows - 1][1], 4, DBL_EPSILON);

    ASSERT_EQ_INT(descend_y(points, rows), SUCCESS);
    ASSERT_NEAR(points[0][1], 4, DBL_EPSILON);
    ASSERT_NEAR(points[rows - 1][1], 0, DBL_EPSILON);

    free_points(points, 5);
    return SUCCESS;
}

static int test_window_functions(void)
{
    char invalid[MAX_WINDOWFUNC_LEN] = "invalid";
    char boxcar_name[MAX_WINDOWFUNC_LEN] = WFUNC_BOXCAR;
    char cosine_name[MAX_WINDOWFUNC_LEN] = WFUNC_COSINE;
    char scosine_name[MAX_WINDOWFUNC_LEN] = WFUNC_COSINE_START;
    char ecosine_name[MAX_WINDOWFUNC_LEN] = WFUNC_COSINE_END;
    char trapezoid_name[MAX_WINDOWFUNC_LEN] = WFUNC_TRAPEZOID;
    char strapezoid_name[MAX_WINDOWFUNC_LEN] = WFUNC_TRAPEZOID_START;
    char etrapezoid_name[MAX_WINDOWFUNC_LEN] = WFUNC_TRAPEZOID_END;

    ASSERT_NEAR(boxcar(), 1.0, DBL_EPSILON);
    ASSERT_NEAR(cosine(1, 10, 10, 10, 0.2, 0.2), 0.146446609407, 1e-12);
    ASSERT_NEAR(cosine(5, 10, 10, 10, 0.2, 0.2), 1.0, DBL_EPSILON);
    ASSERT_NEAR(scosine(1, 10, 10, 10, 0.2), 0.146446609407, 1e-12);
    ASSERT_NEAR(scosine(5, 10, 10, 10, 0.2), 1.0, DBL_EPSILON);
    ASSERT_NEAR(ecosine(9, 10, 10, 10, 0.2), 0.853553390593, 1e-12);
    ASSERT_NEAR(ecosine(5, 10, 10, 10, 0.2), 1.0, DBL_EPSILON);
    ASSERT_NEAR(trapezoid(1, 10, 10, 10, 0.2, 0.2), 0.25, DBL_EPSILON);
    ASSERT_NEAR(trapezoid(5, 10, 10, 10, 0.2, 0.2), 1.0, DBL_EPSILON);
    ASSERT_NEAR(strapezoid(1, 10, 10, 10, 0.2), 0.25, DBL_EPSILON);
    ASSERT_NEAR(strapezoid(5, 10, 10, 10, 0.2), 1.0, DBL_EPSILON);
    ASSERT_NEAR(etrapezoid(9, 10, 10, 10, 0.2), 0.75, DBL_EPSILON);
    ASSERT_NEAR(etrapezoid(5, 10, 10, 10, 0.2), 1.0, DBL_EPSILON);

    ASSERT_NEAR(window_function(5, 1, 10, 10, 10, 0.2, boxcar_name), 1.0, DBL_EPSILON);
    ASSERT_NEAR(window_function(5, 1, 10, 10, 10, 0.2, cosine_name), 1.0, DBL_EPSILON);
    ASSERT_NEAR(window_function(1, 1, 10, 10, 10, 0.2, scosine_name), 0.146446609407, 1e-12);
    ASSERT_NEAR(window_function(9, 1, 10, 10, 10, 0.2, ecosine_name), 0.853553390593, 1e-12);
    ASSERT_NEAR(window_function(1, 1, 10, 10, 10, 0.2, trapezoid_name), 0.25, DBL_EPSILON);
    ASSERT_NEAR(window_function(1, 1, 10, 10, 10, 0.2, strapezoid_name), 0.25, DBL_EPSILON);
    ASSERT_NEAR(window_function(9, 1, 10, 10, 10, 0.2, etrapezoid_name), 0.75, DBL_EPSILON);
    ASSERT_NEAR(window_function(0, 1, 10, 10, 10, 0.2, boxcar_name), WFUNC_ERROR, DBL_EPSILON);
    ASSERT_NEAR(window_function(11, 1, 10, 10, 10, 0.2, boxcar_name), WFUNC_ERROR, DBL_EPSILON);
    ASSERT_NEAR(window_function(5, 1, 10, 10, 10, 0.5, boxcar_name), WFUNC_ERROR, DBL_EPSILON);
    ASSERT_NEAR(window_function(5, 1, 10, 10, 10, 0.2, invalid), WFUNC_ERROR, DBL_EPSILON);

    return SUCCESS;
}

static int test_boundary_box_and_assembly_helpers(void)
{
    window_t data = {0};
    permuted_vertex_t vertex = {0};

    init_rectangle(&data, 5, 4);

    ASSERT_EQ_INT(bounding_boxcheck(&data), SUCCESS);
    ASSERT_EQ_INT(assemble_bv(&data, &vertex), SUCCESS);
    ASSERT_EQ_INT(vertex.row_size_bv, 2);
    ASSERT_NEAR(vertex.bottom_vertices[0][0], 0, DBL_EPSILON);
    ASSERT_NEAR(vertex.bottom_vertices[1][0], 4, DBL_EPSILON);

    ASSERT_EQ_INT(assemble_lv(&data, &vertex), SUCCESS);
    ASSERT_EQ_INT(vertex.row_size_lv, 2);
    ASSERT_NEAR(vertex.left_vertices[0][1], 0, DBL_EPSILON);
    ASSERT_NEAR(vertex.left_vertices[1][1], 3, DBL_EPSILON);

    ASSERT_EQ_INT(assemble_tv(&data, &vertex), SUCCESS);
    ASSERT_EQ_INT(vertex.row_size_tv, 2);
    ASSERT_NEAR(vertex.top_vertices[0][0], 0, DBL_EPSILON);
    ASSERT_NEAR(vertex.top_vertices[1][0], 4, DBL_EPSILON);

    ASSERT_EQ_INT(assemble_rv(&data, &vertex), SUCCESS);
    ASSERT_EQ_INT(vertex.row_size_rv, 2);
    ASSERT_NEAR(vertex.right_vertices[0][1], 0, DBL_EPSILON);
    ASSERT_NEAR(vertex.right_vertices[1][1], 3, DBL_EPSILON);

    ASSERT_EQ_INT(assemble_btlv(&data, &vertex), SUCCESS);
    ASSERT_EQ_INT(assemble_lttv(&data, &vertex), SUCCESS);
    ASSERT_EQ_INT(assemble_rttv(&data, &vertex), SUCCESS);
    ASSERT_EQ_INT(assemble_btrv(&data, &vertex), SUCCESS);
    ASSERT_EQ_INT(vertex.row_size_btlv, 0);
    ASSERT_EQ_INT(vertex.row_size_lttv, 0);
    ASSERT_EQ_INT(vertex.row_size_rttv, 0);
    ASSERT_EQ_INT(vertex.row_size_btrv, 0);

    return SUCCESS;
}

static int test_assembly_helpers_with_quadrant_vertices(void)
{
    window_t data = {0};
    permuted_vertex_t vertex = {0};

    init_octagon_with_quadrant_vertices(&data);

    ASSERT_EQ_INT(bounding_boxcheck(&data), SUCCESS);
    ASSERT_EQ_INT(assemble_bv(&data, &vertex), SUCCESS);
    ASSERT_EQ_INT(assemble_lv(&data, &vertex), SUCCESS);
    ASSERT_EQ_INT(assemble_tv(&data, &vertex), SUCCESS);
    ASSERT_EQ_INT(assemble_rv(&data, &vertex), SUCCESS);
    ASSERT_EQ_INT(assemble_btlv(&data, &vertex), SUCCESS);
    ASSERT_EQ_INT(assemble_lttv(&data, &vertex), SUCCESS);
    ASSERT_EQ_INT(assemble_rttv(&data, &vertex), SUCCESS);
    ASSERT_EQ_INT(assemble_btrv(&data, &vertex), SUCCESS);

    ASSERT_EQ_INT(vertex.row_size_btlv, 1);
    ASSERT_NEAR(vertex.bottom_to_left_vertices[0][0], 1, DBL_EPSILON);
    ASSERT_NEAR(vertex.bottom_to_left_vertices[0][1], 1, DBL_EPSILON);

    ASSERT_EQ_INT(vertex.row_size_lttv, 1);
    ASSERT_NEAR(vertex.left_to_top_vertices[0][0], 1, DBL_EPSILON);
    ASSERT_NEAR(vertex.left_to_top_vertices[0][1], 8, DBL_EPSILON);

    ASSERT_EQ_INT(vertex.row_size_rttv, 1);
    ASSERT_NEAR(vertex.right_to_top_vertices[0][0], 8, DBL_EPSILON);
    ASSERT_NEAR(vertex.right_to_top_vertices[0][1], 8, DBL_EPSILON);

    ASSERT_EQ_INT(vertex.row_size_btrv, 1);
    ASSERT_NEAR(vertex.bottom_to_right_vertices[0][0], 8, DBL_EPSILON);
    ASSERT_NEAR(vertex.bottom_to_right_vertices[0][1], 1, DBL_EPSILON);

    return SUCCESS;
}

static int test_hanging_sweep(void)
{
    double vector[12] = {0};
    double **prev_array = alloc_points(1);
    double **sub_array = alloc_points(1);
    double **hang_array = alloc_points(2);
    int index;

    ASSERT_TRUE(prev_array != NULL && sub_array != NULL && hang_array != NULL);
    set_point(prev_array, 0, 0, 0);
    set_point(sub_array, 0, 6, 6);
    set_point(hang_array, 0, 2, 2);
    set_point(hang_array, 1, 4, 4);
    vector[0] = 0;

    index = hanging_sweep(1, vector, prev_array, 0, sub_array, 0, hang_array, 2, ROW_STEP);
    ASSERT_EQ_INT(index, 5);
    ASSERT_NEAR(vector[2], 2, DBL_EPSILON);
    ASSERT_NEAR(vector[4], 4, DBL_EPSILON);

    memset(vector, 0, sizeof(vector));
    index = hanging_sweep(1, vector, prev_array, 0, sub_array, 0, hang_array, 2, COLUMN_STEP);
    ASSERT_EQ_INT(index, 5);
    ASSERT_NEAR(vector[2], 2, DBL_EPSILON);
    ASSERT_NEAR(vector[4], 4, DBL_EPSILON);

    index = hanging_sweep(7, vector, prev_array, 0, sub_array, 0, hang_array, 0, ROW_STEP);
    ASSERT_EQ_INT(index, 7);

    return SUCCESS;
}

static int test_boundary_assembly_and_contribution(void)
{
    window_t data = {0};
    permuted_vertex_t vertex = {0};
    int i;

    init_rectangle(&data, 5, 4);
    ASSERT_EQ_INT(boundary_assembly(&data, &vertex), SUCCESS);

    ASSERT_EQ_INT(data.minnx, 5);
    ASSERT_EQ_INT(data.minny, 4);
    for (i = 0; i < data.ny; i++) {
        ASSERT_EQ_INT(data.nnx1[i], 0);
        ASSERT_EQ_INT(data.nnx2[i], 4);
    }
    for (i = 0; i < data.nx; i++) {
        ASSERT_EQ_INT(data.nny1[i], 0);
        ASSERT_EQ_INT(data.nny2[i], 3);
    }

    ASSERT_EQ_INT(embedding_contribution2d(2, 2, &data), SUCCESS);
    ASSERT_NEAR(data.contribution, 1.0, DBL_EPSILON);
    ASSERT_EQ_INT(embedding_contribution(2, 2, 3, &data), SUCCESS);
    ASSERT_NEAR(data.contribution, 1.0, DBL_EPSILON);

    ASSERT_EQ_INT(embedding_contribution2d(0, 0, &data), SUCCESS);
    ASSERT_TRUE(data.contribution > 0.0);
    ASSERT_TRUE(data.contribution < 1.0);

    return SUCCESS;
}

static int test_step_vector_functions(void)
{
    window_t data = {0};
    permuted_vertex_t vertex = {0};
    int i;

    init_hexagon(&data);
    ASSERT_EQ_INT(bounding_boxcheck(&data), SUCCESS);
    ASSERT_EQ_INT(assemble_bv(&data, &vertex), SUCCESS);
    ASSERT_EQ_INT(assemble_lv(&data, &vertex), SUCCESS);
    ASSERT_EQ_INT(assemble_tv(&data, &vertex), SUCCESS);
    ASSERT_EQ_INT(assemble_rv(&data, &vertex), SUCCESS);
    ASSERT_EQ_INT(assemble_btlv(&data, &vertex), SUCCESS);
    ASSERT_EQ_INT(assemble_lttv(&data, &vertex), SUCCESS);
    ASSERT_EQ_INT(assemble_rttv(&data, &vertex), SUCCESS);
    ASSERT_EQ_INT(assemble_btrv(&data, &vertex), SUCCESS);

    ASSERT_EQ_INT(vertex.row_size_btrv, 0);
    ASSERT_EQ_INT(vertex.row_size_btlv, 0);
    ASSERT_EQ_INT(vertex.row_size_lttv, 0);
    ASSERT_EQ_INT(vertex.row_size_rttv, 0);

    ASSERT_EQ_INT(minimum_rowstep(&data, &vertex), SUCCESS);
    ASSERT_EQ_INT(maximum_rowstep(&data, &vertex), SUCCESS);
    ASSERT_EQ_INT(minimum_columnstep(&data, &vertex), SUCCESS);
    ASSERT_EQ_INT(maximum_columnstep(&data, &vertex), SUCCESS);

    for (i = 0; i < data.ny; i++) {
        ASSERT_TRUE(data.nnx1[i] <= data.nnx2[i]);
    }
    for (i = 0; i < data.nx; i++) {
        ASSERT_TRUE(data.nny1[i] <= data.nny2[i]);
    }

    return SUCCESS;
}

int main(void)
{
    ASSERT_EQ_INT(test_sorting_and_unique_vertices(), SUCCESS);
    ASSERT_EQ_INT(test_window_functions(), SUCCESS);
    ASSERT_EQ_INT(test_boundary_box_and_assembly_helpers(), SUCCESS);
    ASSERT_EQ_INT(test_assembly_helpers_with_quadrant_vertices(), SUCCESS);
    ASSERT_EQ_INT(test_hanging_sweep(), SUCCESS);
    ASSERT_EQ_INT(test_boundary_assembly_and_contribution(), SUCCESS);
    ASSERT_EQ_INT(test_step_vector_functions(), SUCCESS);

    return SUCCESS;
}
