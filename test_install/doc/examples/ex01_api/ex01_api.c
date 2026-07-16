#include "blend.h"

int main(void)
{
    const int n = 100;
    const int nmin = n;
    const double ratio = 0.2;
    const char *filename = "ex01_api.txt";
    FILE *fp;
    int i;

    fp = fopen(filename, "w");
    if (fp == NULL) {
        perror(filename);
        return FAIL;
    }

    for (i = 1; i <= n; i++) {
        double weight = cosine(i, n, n, nmin, ratio, ratio);
        fprintf(fp, "%d %.12f\n", i, weight);
    }

    if (fclose(fp) != 0) {
        perror(filename);
        return FAIL;
    }

    return SUCCESS;
}
