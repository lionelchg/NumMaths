#ifndef _LINALG_
#define _LINALG_

#include "utils.h"

// Matrix multiplication 
double** matmul(double **A, double**B, int n);

// Matrix-Vector multiplication
double* matvecmul(double **A, double *x, int n);

// LU decomposition for appropriate matrix
void LU(double **A, int n);

// Creation of Poisson1D dirichlet matrix
double** poisson_mat(int n);

// Solve upper triangular system
double *solve_upper(double **A, double *b, int n);

// Solve lower triangular system
double *solve_lower(double **A, double *b, int n);

#endif