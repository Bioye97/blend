#include "blend_api.h"
#include "blend_boundary.h"
#include "blend_contribution.h"
#include "blend_error.h"
#include "blend_polygon.h"
#include "blend_report.h"
#include "blend_utils.h"
#include "blend_window.h"

int main(void)
{
    blend_error_code code = SUCCESS;
    blend_verbosity verbosity = BLEND_MSG_WARNING;
    blend_window_function function = WFUNC_COSINE;
    vertex point = {0.0, 0.0};
    polygon poly = {0};
    window data = {0};
    permuted_vertex boundary = {0};

    (void)code;
    (void)verbosity;
    (void)function;
    (void)point;
    (void)poly;
    (void)data;
    (void)boundary;

    return 0;
}
