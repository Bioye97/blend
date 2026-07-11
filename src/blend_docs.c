/*
* Author: Rasheed Ajala
* Contact: abioyeajala@gmail.com
*
*/

#include "blend_module.h"

#include <stdio.h>
#include <string.h>

static void docs_usage(FILE *fp)
{
    fprintf(fp, "blend docs - Show documentation for a BLEND module\n\n");
    fprintf(fp, "usage: blend docs [-Q] [-V[q|e|w|t|i|c|d]] <module-name> [<-option>]\n\n");

    fprintf(fp, "REQUIRED ARGUMENTS:\n\n");
    fprintf(fp, "  <module-name>\n");
    fprintf(fp, "     One of the BLEND modules. Run blend --show-modules to list the available\n");
    fprintf(fp, "     module names.\n\n");

    fprintf(fp, "OPTIONAL ARGUMENTS:\n\n");
    fprintf(fp, "  -Q\n");
    fprintf(fp, "     Only display the documentation target and do not open it in a viewer.\n");
    fprintf(fp, "     This is currently the default behavior for BLEND.\n\n");

    fprintf(fp, "  <-option>\n");
    fprintf(fp, "     The one-letter option of the module in question, for example -R. Option\n");
    fprintf(fp, "     specific documentation is reserved for a future documentation backend.\n\n");

    fprintf(fp, "  -V[q|e|w|t|i|c|d]\n");
    fprintf(fp, "     Change the verbosity level. Choose among q, e, w, t, i, c, and d.\n\n");

    fprintf(fp, "  -^ (or -)\n");
    fprintf(fp, "     Print short synopsis message.\n\n");

    fprintf(fp, "  -+ (or +)\n");
    fprintf(fp, "     Print longer synopsis message.\n\n");

    fprintf(fp, "  -? (or no arguments)\n");
    fprintf(fp, "     Print this usage message.\n");
}

static void docs_short_synopsis(FILE *fp)
{
    fprintf(fp, "usage: blend docs [-Q] [-V[q|e|w|t|i|c|d]] <module-name> [<-option>]\n");
}

static void docs_long_synopsis(FILE *fp)
{
    docs_short_synopsis(fp);
    fprintf(fp, "\nShow documentation for a BLEND module. Run blend --show-modules to list module names.\n");
}

int blend_docs_module(int argc, char **argv)
{
    const blend_module *module = NULL;
    const char *module_name = NULL;
    const char *option = NULL;
    int quiet = 0;
    int i;

    if (argc == 0) {
        docs_usage(stdout);
        return 0;
    }

    if (strcmp(argv[0], "-?") == 0 || strcmp(argv[0], "--help") == 0) {
        docs_usage(stdout);
        return 0;
    }
    if (strcmp(argv[0], "-^") == 0 || strcmp(argv[0], "-") == 0) {
        docs_short_synopsis(stdout);
        return 0;
    }
    if (strcmp(argv[0], "-+") == 0 || strcmp(argv[0], "+") == 0) {
        docs_long_synopsis(stdout);
        return 0;
    }

    for (i = 0; i < argc; i++) {
        if (strcmp(argv[i], "-Q") == 0) {
            quiet = 1;
        }
        else if (argv[i][0] == '-' && argv[i][1] != '\0' && module_name == NULL) {
            docs_usage(stderr);
            return 1;
        }
        else if (module_name == NULL) {
            module_name = argv[i];
        }
        else if (option == NULL) {
            option = argv[i];
        }
        else {
            docs_usage(stderr);
            return 1;
        }
    }

    if (module_name == NULL) {
        docs_usage(stdout);
        return 0;
    }

    module = blend_find_module(module_name);
    if (module == NULL) {
        fprintf(stderr, "blend docs: unknown module %s\n", module_name);
        fprintf(stderr, "Run blend --show-modules to list available modules.\n");
        return 1;
    }

    if (quiet) {
        printf("blend docs %s\n", module->name);
    }
    else {
        printf("%s - %s\n", module->name, module->purpose);
    }
    if (option != NULL) {
        printf("option: %s\n", option);
        printf("Option-specific documentation is not yet available.\n");
    }

    return 0;
}
