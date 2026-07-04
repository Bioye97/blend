#include "blend.h"

static double **alloc_vertices(int rows)
{
    int i;
    double **vertices = (double **)calloc((size_t)rows, sizeof(double *));

    if (vertices == NULL) return NULL;

    for (i = 0; i < rows; i++) {
        vertices[i] = (double *)calloc(2, sizeof(double));
        if (vertices[i] == NULL) return NULL;
    }

    return vertices;
}

static void free_vertices(double **vertices, int rows)
{
    int i;

    if (vertices == NULL) return;
    for (i = 0; i < rows; i++) {
        free(vertices[i]);
    }
    free(vertices);
}

static void set_vertex(window_t *data, int row, double x, double y)
{
    data->vertices[row][0] = x;
    data->vertices[row][1] = y;
}

int main(void)
{
    const char *filename = "ex09.txt";
    window_t data = {0};
    permuted_vertex_t polygon = {0};
    FILE *fp;
    int x, y;

    data.nx = 100;
    data.ny = 100;
    data.row_size = 9;
    data.ratio = 0.2;
    strcpy(data.x_function, WFUNC_COSINE);
    strcpy(data.y_function, WFUNC_COSINE);

    data.vertices = alloc_vertices(data.row_size);
    if (data.vertices == NULL) {
        fprintf(stderr, "Could not allocate nonagon vertices.\n");
        return FAIL;
    }

    set_vertex(&data, 0, 30, 0);
    set_vertex(&data, 1, 60, 0);
    set_vertex(&data, 2, 90, 20);
    set_vertex(&data, 3, 99, 50);
    set_vertex(&data, 4, 90, 80);
    set_vertex(&data, 5, 60, 99);
    set_vertex(&data, 6, 30, 99);
    set_vertex(&data, 7, 0, 70);
    set_vertex(&data, 8, 0, 30);

    if (boundary_assembly(&data, &polygon) != SUCCESS) {
        free_vertices(data.vertices, data.row_size);
        return FAIL;
    }

    fp = fopen(filename, "w");
    if (fp == NULL) {
        perror(filename);
        free_vertices(data.vertices, data.row_size);
        return FAIL;
    }

    for (y = 0; y < data.ny; y++) {
        for (x = 0; x < data.nx; x++) {
            if (embedding_contribution2d(x, y, &data) != SUCCESS) {
                fclose(fp);
                free_vertices(data.vertices, data.row_size);
                return FAIL;
            }
            fprintf(fp, "%d %d %.12f\n", x, y, data.contribution);
        }
        fprintf(fp, "\n");
    }

    if (fclose(fp) != 0) {
        perror(filename);
        free_vertices(data.vertices, data.row_size);
        return FAIL;
    }

    free_vertices(data.vertices, data.row_size);

    return SUCCESS;
}
