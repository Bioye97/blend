/*
* Author: Rasheed Ajala
* Contact: abioyeajala@gmail.com
*
*/

#ifndef BLEND_MODULE_H
#define BLEND_MODULE_H

typedef struct blend_module {
    const char *name;
    const char *purpose;
    int (*run)(int argc, char **argv);
} blend_module;

int blend_docs_module(int argc, char **argv);
int blend_monotone_module(int argc, char **argv);
int blend_window1d_module(int argc, char **argv);
int blend_window2d_module(int argc, char **argv);
int blend_window3d_module(int argc, char **argv);

const blend_module *blend_module_registry(int *count);
const blend_module *blend_find_module(const char *name);
void blend_print_module_list(void);
void blend_print_module_names(void);

#endif /* BLEND_MODULE_H */
