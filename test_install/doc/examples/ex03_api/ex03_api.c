#include "blend.h"

int main(void)
{
    const char *filename = "ex03_api.txt";
    window data = {0};
    permuted_vertex boundary = {0};
    polygon support = {0};
    FILE *fp;
    int x, y, z;

    data.nx = 100;
    data.ny = 100;
    data.nz = 25;
    data.ratio_x1 = 0.2;
    data.ratio_x2 = 0.2;
    data.ratio_y1 = 0.2;
    data.ratio_y2 = 0.2;
    data.ratio_z1 = 0.2;
    data.ratio_z2 = 0.2;
    data.x_function = WFUNC_COSINE;
    data.y_function = WFUNC_COSINE;
    data.z_function = WFUNC_COSINE;

    if (blend_polygon_alloc(&support, 8) != SUCCESS) {
        fprintf(stderr, "Could not allocate isotoxal-star vertices.\n");
        return FAIL;
    }

    blend_polygon_set_vertex(&support, 0, 50, 0);
    blend_polygon_set_vertex(&support, 1, 65, 35);
    blend_polygon_set_vertex(&support, 2, 99, 50);
    blend_polygon_set_vertex(&support, 3, 65, 65);
    blend_polygon_set_vertex(&support, 4, 50, 99);
    blend_polygon_set_vertex(&support, 5, 35, 65);
    blend_polygon_set_vertex(&support, 6, 0, 50);
    blend_polygon_set_vertex(&support, 7, 35, 35);

    if (blend_window_set_polygon(&data, &support) != SUCCESS) {
        fprintf(stderr, "Could not assign isotoxal-star vertices.\n");
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

    for (z = 0; z < data.nz; z++) {
        for (y = 0; y < data.ny; y++) {
            for (x = 0; x < data.nx; x++) {
                if (embedding_contribution3d(x, y, z, &data) != SUCCESS) {
                    fclose(fp);
                    blend_window_clear_polygon(&data);
                    blend_polygon_free(&support);
                    return FAIL;
                }
                fprintf(fp, "%d %d %d %.12f\n", x, y, z, data.contribution);
            }
            fprintf(fp, "\n");
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
