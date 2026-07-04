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

int main(void)
{
    const char *filename = "ex03.txt";
    window_t data = {0};
    permuted_vertex_t polygon = {0};
    FILE *fp;
    int x, y;

    data.nx = 100;
    data.ny = 100;
    data.row_size = 3;
    data.ratio = 0.2;
    strcpy(data.x_function, WFUNC_COSINE);
    strcpy(data.y_function, WFUNC_COSINE);

    data.vertices = alloc_vertices(data.row_size);
    if (data.vertices == NULL) {
        fprintf(stderr, "Could not allocate triangle vertices.\n");
        return FAIL;
    }

    data.vertices[0][0] = 0;
    data.vertices[0][1] = 0;
    data.vertices[1][0] = data.nx - 1;
    data.vertices[1][1] = 0;
    data.vertices[2][0] = 0;
    data.vertices[2][1] = data.ny - 1;

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
