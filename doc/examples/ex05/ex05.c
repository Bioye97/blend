#include "blend.h"

int main(void)
{
    const char *filename = "ex05.txt";
    window data = {0};
    permuted_vertex boundary = {0};
    polygon support = {0};
    FILE *fp;
    int x, y;

    data.nx = 100;
    data.ny = 100;
    data.ratio_x1 = 0.2;
    data.ratio_x2 = 0.2;
    data.ratio_y1 = 0.2;
    data.ratio_y2 = 0.2;
    data.ratio_z1 = 0.2;
    data.ratio_z2 = 0.2;
    data.x_function = WFUNC_COSINE;
    data.y_function = WFUNC_COSINE;

    if (blend_polygon_alloc(&support, 5) != SUCCESS) {
        fprintf(stderr, "Could not allocate pentagon vertices.\n");
        return FAIL;
    }

    blend_polygon_set_vertex(&support, 0, 50, 0);
    blend_polygon_set_vertex(&support, 1, 99, 45);
    blend_polygon_set_vertex(&support, 2, 75, 99);
    blend_polygon_set_vertex(&support, 3, 0, 99);
    blend_polygon_set_vertex(&support, 4, 0, 35);

    if (blend_window_set_polygon(&data, &support) != SUCCESS) {
        fprintf(stderr, "Could not assign pentagon vertices.\n");
        blend_polygon_free(&support);
        return FAIL;
    }

    if (boundary_assembly(&data, &boundary) != SUCCESS) {
        blend_window_clear_polygon(&data);
        blend_polygon_free(&support);
        return FAIL;
    }

    fp = fopen(filename, "w");
    if (fp == NULL) {
        perror(filename);
        blend_window_clear_polygon(&data);
        blend_polygon_free(&support);
        return FAIL;
    }

    for (y = 0; y < data.ny; y++) {
        for (x = 0; x < data.nx; x++) {
            if (embedding_contribution2d(x, y, &data) != SUCCESS) {
                fclose(fp);
                blend_window_clear_polygon(&data);
                blend_polygon_free(&support);
                return FAIL;
            }
            fprintf(fp, "%d %d %.12f\n", x, y, data.contribution);
        }
        fprintf(fp, "\n");
    }

    if (fclose(fp) != 0) {
        perror(filename);
        blend_window_clear_polygon(&data);
        blend_polygon_free(&support);
        return FAIL;
    }

    blend_window_clear_polygon(&data);
    blend_polygon_free(&support);

    return SUCCESS;
}
