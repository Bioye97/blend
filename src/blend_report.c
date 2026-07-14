/*
* Author: Rasheed Ajala
* Contact: abioyeajala@gmail.com
*
*/

#include "blend.h"

#include <stdarg.h>
#include <sys/time.h>
#include <time.h>

static int blend_verbosity_initialized = 0;
static blend_verbosity blend_current_verbosity = BLEND_MSG_WARNING;

const char *blend_verbosity_name(blend_verbosity level)
{
    switch (level) {
        case BLEND_MSG_QUIET:
            return "quiet";
        case BLEND_MSG_ERROR:
            return "error";
        case BLEND_MSG_WARNING:
            return "warning";
        case BLEND_MSG_TIMING:
            return "timing";
        case BLEND_MSG_INFORMATION:
            return "information";
        case BLEND_MSG_COMPAT:
            return "compat";
        case BLEND_MSG_DEBUG:
            return "debug";
        default:
            return "unknown";
    }
}

int blend_verbosity_from_name(const char *name, blend_verbosity *level)
{
    if (name == NULL || level == NULL) {
        return FAIL;
    }

    if (strcmp(name, "q") == 0 || strcmp(name, "quiet") == 0) {
        *level = BLEND_MSG_QUIET;
    }
    else if (strcmp(name, "e") == 0 || strcmp(name, "error") == 0) {
        *level = BLEND_MSG_ERROR;
    }
    else if (strcmp(name, "w") == 0 || strcmp(name, "warning") == 0) {
        *level = BLEND_MSG_WARNING;
    }
    else if (strcmp(name, "t") == 0 || strcmp(name, "timing") == 0) {
        *level = BLEND_MSG_TIMING;
    }
    else if (strcmp(name, "i") == 0 || strcmp(name, "information") == 0 || strcmp(name, "info") == 0) {
        *level = BLEND_MSG_INFORMATION;
    }
    else if (strcmp(name, "c") == 0 || strcmp(name, "compat") == 0 || strcmp(name, "compatibility") == 0) {
        *level = BLEND_MSG_COMPAT;
    }
    else if (strcmp(name, "d") == 0 || strcmp(name, "debug") == 0) {
        *level = BLEND_MSG_DEBUG;
    }
    else {
        return FAIL;
    }

    return SUCCESS;
}

static void blend_init_verbosity(void)
{
    const char *value;
    blend_verbosity level;

    if (blend_verbosity_initialized) {
        return;
    }

    value = getenv("BLEND_VERBOSE");
    if (value != NULL && blend_verbosity_from_name(value, &level) == SUCCESS) {
        blend_current_verbosity = level;
    }

    blend_verbosity_initialized = 1;
}

int blend_set_verbosity(blend_verbosity level)
{
    if (level < BLEND_MSG_QUIET || level > BLEND_MSG_DEBUG) {
        return FAIL;
    }

    blend_current_verbosity = level;
    blend_verbosity_initialized = 1;
    return SUCCESS;
}

blend_verbosity blend_get_verbosity(void)
{
    blend_init_verbosity();
    return blend_current_verbosity;
}

double blend_elapsed_seconds(void)
{
#if defined(CLOCK_MONOTONIC)
    struct timespec ts;

    if (clock_gettime(CLOCK_MONOTONIC, &ts) == 0) {
        return (double)ts.tv_sec + (double)ts.tv_nsec / 1.0e9;
    }
#endif

    {
        struct timeval tv;

        if (gettimeofday(&tv, NULL) == 0) {
            return (double)tv.tv_sec + (double)tv.tv_usec / 1.0e6;
        }
    }

    return (double)clock() / (double)CLOCKS_PER_SEC;
}

void BLEND_Report(blend_verbosity level, const char *format, ...)
{
    va_list args;
    blend_verbosity verbosity = blend_get_verbosity();

    if (verbosity == BLEND_MSG_QUIET || level > verbosity) {
        return;
    }

    fprintf(stderr, "blend [%s]: ", blend_verbosity_name(level));
    va_start(args, format);
    vfprintf(stderr, format, args);
    va_end(args);
}
