/*
* Author: Rasheed Ajala
* Contact: abioyeajala@gmail.com
*
*/

#include "blend_module.h"

#include <stdio.h>
#include <string.h>

static const blend_module modules[] = {
    {"docs", "Show documentation for a BLEND module", blend_docs_module},
    {"simplify", "Simplify polygon vertices", blend_simplify_module},
    {"monotone", "Check or enforce xy-monotone polygon boundaries", blend_monotone_module},
    {"window1d", "Generate 1D window weights", blend_window1d_module},
    {"window2d", "Generate 2D window weights", blend_window2d_module},
    {"window3d", "Generate 3D window weights", blend_window3d_module}
};

const blend_module *blend_module_registry(int *count)
{
    if (count != NULL) {
        *count = (int)(sizeof(modules) / sizeof(modules[0]));
    }
    return modules;
}

const blend_module *blend_find_module(const char *name)
{
    int i, count;
    const blend_module *registry = blend_module_registry(&count);

    if (name == NULL) {
        return NULL;
    }

    for (i = 0; i < count; i++) {
        if (strcmp(name, registry[i].name) == 0) {
            return &registry[i];
        }
    }

    return NULL;
}

void blend_print_module_list(void)
{
    int i, count;
    const blend_module *registry = blend_module_registry(&count);

    for (i = 0; i < count; i++) {
        printf("  %-10s %s\n", registry[i].name, registry[i].purpose);
    }
}

void blend_print_module_names(void)
{
    int i, count;
    const blend_module *registry = blend_module_registry(&count);

    for (i = 0; i < count; i++) {
        printf("%s\n", registry[i].name);
    }
}
