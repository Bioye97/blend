/*
* Author: Rasheed Ajala
* Contact: abioyeajala@gmail.com
*
*/

#include "blend_module.h"

#include <stdio.h>

int blend_docs_module(int argc, char **argv)
{
    const blend_module *module = NULL;

    if (argc > 0) {
        module = blend_find_module(argv[0]);
    }

    if (module != NULL) {
        printf("%s - %s\n", module->name, module->purpose);
        return 0;
    }

    printf("BLEND modules:\n");
    blend_print_module_list();
    return 0;
}
