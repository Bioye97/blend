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
