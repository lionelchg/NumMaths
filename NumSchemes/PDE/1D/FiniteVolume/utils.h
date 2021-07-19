#ifndef _UTILS_
#define _UTILS_

#include <stdio.h>
#include <stdlib.h>
#include "math.h"

// Print vector in console
void print_vec(double *vec, int n, char *vec_name);

// Write list of vec to file
void write_vecs(double **list_vec, int n, int m, char **vec_names, char *filename);

// Vector allocation and set coefficients to zero
double* alloc_vec(int n);

// Vector copy
double* copy_vec(double *vec, int n);

// Linearly spaced array
double* linspace(double begin, double end, int npoints);

#endif