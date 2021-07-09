#ifndef _UTILS_
#define _UTILS_

#include <stdio.h>
#include <stdlib.h>

// Print full matrix in console
void print_matrix(double **mat, int n, char *mat_name);

// Square matrix allocation and set all the coefficients to 
// zero
double** alloc_mat(int n);

// Square matrix copy
double** copy_mat(double **mat, int n);

// Print vector in console
void print_vec(double *vec, int n, char *vec_name);

// Vector allocation and set coefficients to zero
double* alloc_vec(int n);

// Vector copy
double* copy_vec(double *vec, int n);

#endif