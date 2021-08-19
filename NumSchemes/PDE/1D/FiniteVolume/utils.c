#include "utils.h"

void print_vec(double *vec, int n, char *vec_name)
{
    int i, j;
    int step = n / 5;
    printf("%s\n", vec_name);
    for (i = 0; i < n; i += step)
    {
        printf("%.3e ", vec[i]);
    }
    printf("\n");
}

void write_vecs(double **list_vec, int n, int m, char **vec_names, char *filename)
{
    int i, j;
    FILE *fptr;

    // use appropriate location if you are using MacOS or Linux
    printf("Print results to %s file\n", filename);
    fptr = fopen(filename, "w");

    if (fptr == NULL)
    {
        printf("Error while opening file");
        exit(1);
    }
    // Write header
    for (j = 0; j < m; j++) {
        fprintf(fptr, "%13s", vec_names[j]);
    }
    fprintf(fptr, "\n");

    // Write vector contents
    for (i = 0; i < n; i++) {
        for (j = 0; j < m; j++) {
            fprintf(fptr, "%13.4e", list_vec[j][i]);
        }
        fprintf(fptr, "\n");
    }

    fclose(fptr);
}

double *alloc_vec(int n)
{
    int i;
    double *vec;
    vec = (double *)malloc(n * sizeof(double));
    for (i = 0; i < n; i++)
    {
        vec[i] = 0.0;
    }
    return vec;
}

double *copy_vec(double *vec, int n)
{
    int i;
    double *copied_vec = alloc_vec(n);
    for (i = 0; i < n; i++)
    {
        copied_vec[i] = vec[i];
    }
    return copied_vec;
}

double *linspace(double begin, double end, int npoints)
{
    double *array = alloc_vec(npoints);
    int i;
    double h = (end - begin) / (double)(npoints - 1);
    for (i = 0; i < npoints; i++)
    {
        array[i] = begin + h * i;
    }
    return array;
}
