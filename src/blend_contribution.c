/*
* Author: Rasheed Ajala
* Contact: abioyeajala@gmail.com
*
*/

#include "blend.h"

/* Embedding contribution 1D */
int embedding_contribution1d(int x, window *data) {

    double taper_x;

    /* The tapers are defined on [1 nx], whereas the boundaries are defined on [0, nx-1] */
    taper_x = window_function(x+1, 1, data->nx, data->nx, data->nx,
                              data->ratio_x1, data->ratio_x2, data->x_function);
    /* Error check */
    if (taper_x == WFUNC_ERROR) {
        BLEND_Report(BLEND_MSG_ERROR, "contribution: could not compute weight from x_function\n");
        return FAIL;
    }

    /* Assemble output */
    data->contribution = taper_x;

    return SUCCESS;
}

/* Embedding contribution 2D */
int embedding_contribution2d(int x, int y, window *data) {

    int nny1, nny2, nnx1, nnx2;
    double taper_x, taper_y;
    /* First define the values of these vectors */
    nny1 = data->nny1[x];
    nny2 = data->nny2[x];
    nnx1 = data->nnx1[y];
    nnx2 = data->nnx2[y];

    /* The tapers are defined on [1 nx], whereas the boundaries are defined on [0, nx-1] */
    taper_x = window_function(y+1, nny1+1, nny2+1, data->ny, data->minny,
                              data->ratio_x1, data->ratio_x2, data->x_function);
    /* Error check */
    if (taper_x == WFUNC_ERROR) {
        BLEND_Report(BLEND_MSG_ERROR, "contribution: could not compute weight from x_function\n");
        return FAIL;
    }


    /* The tapers are defined on [1 nx], whereas the boundaries are defined on [0, nx-1] */
    taper_y = window_function(x+1, nnx1+1, nnx2+1, data->nx, data->minnx,
                              data->ratio_y1, data->ratio_y2, data->y_function);
    /* Error check */
    if (taper_y == WFUNC_ERROR) {
        BLEND_Report(BLEND_MSG_ERROR, "contribution: could not compute weight from y_function\n");
        return FAIL;
    }

    /* Assemble output */
    data->contribution = taper_x * taper_y;

    return SUCCESS;
}

/* Embedding contribution 3D */
int embedding_contribution3d(int x, int y, int z, window *data) {

    int nny1, nny2, nnx1, nnx2;
    double taper_x, taper_y, taper_z;
    /* First define the values of these vectors */
    nny1 = data->nny1[x];
    nny2 = data->nny2[x];
    nnx1 = data->nnx1[y];
    nnx2 = data->nnx2[y];

     /* The tapers are defined on [1 nx], whereas the boundaries are defined on [0, nx-1] */
    taper_x = window_function(y+1, nny1+1, nny2+1, data->ny, data->minny,
                              data->ratio_x1, data->ratio_x2, data->x_function);
    /* Error check */
    if (taper_x == WFUNC_ERROR) {
        BLEND_Report(BLEND_MSG_ERROR, "contribution: could not compute weight from x_function\n");
        return FAIL;
    }


     /* The tapers are defined on [1 nx], whereas the boundaries are defined on [0, nx-1] */
    taper_y = window_function(x+1, nnx1+1, nnx2+1, data->nx, data->minnx,
                              data->ratio_y1, data->ratio_y2, data->y_function);
    /* Error check */
    if (taper_y == WFUNC_ERROR) {
        BLEND_Report(BLEND_MSG_ERROR, "contribution: could not compute weight from y_function\n");
        return FAIL;
    }

     /* The tapers are defined on [1 nx], whereas the boundaries are defined on [0, nx-1] */
    taper_z = window_function(z+1, 1, data->nz, data->nz, data->nz,
                              data->ratio_z1, data->ratio_z2, data->z_function);
    /* Error check */
    if (taper_z == WFUNC_ERROR) {
        BLEND_Report(BLEND_MSG_ERROR, "contribution: could not compute weight from z_function\n");
        return FAIL;
    }

    /* Assemble output */
    data->contribution = taper_x * taper_y * taper_z;

    return SUCCESS;
}
