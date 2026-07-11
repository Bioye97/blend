/*
* Author: Rasheed Ajala
* Contact: abioyeajala@gmail.com
*/

#ifndef BLEND_REPORT_H
#define BLEND_REPORT_H

#include "blend_error.h"

/* Message verbosity */
typedef enum blend_verbosity {
    BLEND_MSG_QUIET = 0,
    BLEND_MSG_ERROR,
    BLEND_MSG_WARNING,
    BLEND_MSG_TIMING,
    BLEND_MSG_INFORMATION,
    BLEND_MSG_COMPAT,
    BLEND_MSG_DEBUG
} blend_verbosity;

/* Message verbosity name */
const char *blend_verbosity_name(blend_verbosity level);

/* Message verbosity lookup */
int blend_verbosity_from_name(const char *name, blend_verbosity *level);

/* Set active message verbosity */
int blend_set_verbosity(blend_verbosity level);

/* Get active message verbosity */
blend_verbosity blend_get_verbosity(void);

/* Report a message according to the active verbosity */
void BLEND_Report(blend_verbosity level, const char *format, ...);

#endif /* BLEND_REPORT_H */
