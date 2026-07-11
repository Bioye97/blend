/*
* Author: Rasheed Ajala
* Contact: abioyeajala@gmail.com
*
*/

#include "blend.h"

const char *blend_error_message(blend_error_code code)
{
    switch (code) {
        case SUCCESS:
            return "success";
        case FAIL:
            return "failure";
        case WFUNC_ERROR:
            return "window function error";
        default:
            return "unknown error";
    }
}
