/*
* Author: Rasheed Ajala
* Contact: abioyeajala@gmail.com
*/

#ifndef BLEND_UTILS_H
#define BLEND_UTILS_H

/* Linear interpolation */
int interpolate_linear(double x0, double y0, double x1, double y1, double x, double *y);

/* Bilinear interpolation */
int interpolate_bilinear(double x0, double x1, double y0, double y1,
                         double q00, double q10, double q01, double q11,
                         double x, double y, double *value);

/* Trilinear interpolation */
int interpolate_trilinear(double x0, double x1, double y0, double y1, double z0, double z1,
                          double q000, double q100, double q010, double q110,
                          double q001, double q101, double q011, double q111,
                          double x, double y, double z, double *value);

#endif /* BLEND_UTILS_H */
