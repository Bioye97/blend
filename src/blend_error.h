/*
* Author: Rasheed Ajala
* Contact: abioyeajala@gmail.com
*/

#ifndef BLEND_ERROR_H
#define BLEND_ERROR_H

/* Error codes */
typedef enum blend_error_code {
    SUCCESS = 0,
    FAIL = 1,
    WFUNC_ERROR = 2
} blend_error_code;

/* Error code description */
const char *blend_error_message(blend_error_code code);

#endif /* BLEND_ERROR_H */
