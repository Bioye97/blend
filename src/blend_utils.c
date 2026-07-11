/*
* Author: Rasheed Ajala
* Contact: abioyeajala@gmail.com
*
*/

#include "blend.h"

int interpolate_linear(double x0, double y0, double x1, double y1, double x, double *y)
{
    double t;

    if (y == NULL) {
        BLEND_Report(BLEND_MSG_ERROR, "utils: output pointer for linear interpolation is NULL\n");
        return FAIL;
    }

    if (!isfinite(x0) || !isfinite(y0) || !isfinite(x1) || !isfinite(y1) || !isfinite(x)) {
        BLEND_Report(BLEND_MSG_ERROR, "utils: linear interpolation requires finite inputs\n");
        return FAIL;
    }

    if (x0 == x1) {
        BLEND_Report(BLEND_MSG_ERROR, "utils: linear interpolation requires distinct x coordinates\n");
        return FAIL;
    }

    t = (x - x0) / (x1 - x0);
    *y = y0 + t * (y1 - y0);

    return SUCCESS;
}

int interpolate_bilinear(double x0, double x1, double y0, double y1,
                         double q00, double q10, double q01, double q11,
                         double x, double y, double *value)
{
    double tx, ty;

    if (value == NULL) {
        BLEND_Report(BLEND_MSG_ERROR, "utils: output pointer for bilinear interpolation is NULL\n");
        return FAIL;
    }

    if (!isfinite(x0) || !isfinite(x1) || !isfinite(y0) || !isfinite(y1) ||
        !isfinite(q00) || !isfinite(q10) || !isfinite(q01) || !isfinite(q11) ||
        !isfinite(x) || !isfinite(y)) {
        BLEND_Report(BLEND_MSG_ERROR, "utils: bilinear interpolation requires finite inputs\n");
        return FAIL;
    }

    if (x0 == x1) {
        BLEND_Report(BLEND_MSG_ERROR, "utils: bilinear interpolation requires distinct x coordinates\n");
        return FAIL;
    }
    if (y0 == y1) {
        BLEND_Report(BLEND_MSG_ERROR, "utils: bilinear interpolation requires distinct y coordinates\n");
        return FAIL;
    }

    tx = (x - x0) / (x1 - x0);
    ty = (y - y0) / (y1 - y0);

    *value = (1.0 - tx) * (1.0 - ty) * q00 +
             tx * (1.0 - ty) * q10 +
             (1.0 - tx) * ty * q01 +
             tx * ty * q11;

    return SUCCESS;
}
