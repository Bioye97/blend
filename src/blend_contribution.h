/*
* Author: Rasheed Ajala
* Contact: abioyeajala@gmail.com
*/

#ifndef BLEND_CONTRIBUTION_H
#define BLEND_CONTRIBUTION_H

#include "blend_window.h"

/* Embedding contribution for 1D case: Blending weight */
int embedding_contribution1d(int x, window *data);

/* Embedding contribution for 2D case: Blending weight */
int embedding_contribution2d(int x, int y, window *data);

/* Embedding contribution for 3D case: Blending weight */
int embedding_contribution3d(int x, int y, int z, window *data);

#endif /* BLEND_CONTRIBUTION_H */
