#include "utils.h"

void print_matrix(double **mat, int n, char *mat_name) {
    int i, j;
    int step = n / 5;
    printf("%s\n", mat_name);
    for (i = 0; i < n; i+=step) {
        for (j = 0; j < n; j+=step) {
            printf("%.3e ", mat[i][j]);
        }
        printf("\n");
    }
    printf("\n");
}

void print_vec(double *vec, int n, char *vec_name) {
    int i, j;
    int step = n / 5;
    printf("%s\n", vec_name);
    for (i = 0; i < n; i+=step) {
        printf("%.3e ", vec[i]);
    }
    printf("\n");
}

double** alloc_mat(int n) {
    int i, j;
    double **mat;
    mat = (double**) malloc(n * n * sizeof(double));
    for (i = 0; i < n; i++) mat[i] = (double*) malloc(n * sizeof(double));
    for (i = 0; i < n; i++) {
        for (j = 0; j < n; j++) {
            mat[i][j] = 0.0;
        }
    }
    return mat;
}

double** copy_mat(double **mat, int n) {
    int i, j;
    double **copied_mat = alloc_mat(n);
    for (i = 0; i < n; i++) {
        for (j = 0; j < n; j++) {
            copied_mat[i][j] = mat[i][j];
        }
    }
    return copied_mat;
}

double* alloc_vec(int n) {
    int i;
    double *vec;
    vec = (double*) malloc(n * sizeof(double));
    for (i = 0; i < n; i++) {
        vec[i] = 0.0;
    }
    return vec;
}

double* copy_vec(double *vec, int n) {
    int i;
    double *copied_vec = alloc_vec(n);
    for (i = 0; i < n; i++) {
        copied_vec[i] = vec[i];
    }
    return copied_vec;
}