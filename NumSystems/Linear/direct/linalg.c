#include "linalg.h"

double** matmul(double **A, double**B, int n) {
    int i, j, k;
    double **C;
    C = alloc_mat(n);
    for (i = 0; i < n; i++) {
        for (j = 0; j < n; j++) {
            for (k = 0; k < n; k++) {
                C[i][j] += A[i][k] * B[k][j];
            }
        }
    }
    return C;
}

double* matvecmul(double **A, double *x, int n) {
    int i, j, k;
    double *y = alloc_vec(n);
    for (i = 0; i < n; i++) {
        y[i] = 0.0;
        for (j = 0; j < n; j++) {
            y[i] += A[i][j] * x[j];
        }
    }
    return y;
}

double** poisson_mat(int n) {
    int i, j, k;
    double **mat;
    mat = alloc_mat(n);
    mat[0][0] = 1.0;
    mat[n - 1][n - 1] = 1.0;
    for (i = 1; i < n -1; i++) {
        mat[i][i - 1] = 1.0;
        mat[i][i] = - 2.0;
        mat[i][i + 1] = 1.0;
    }
    return mat;
}

void LU(double **A, int n) {
    int i, j, k;
    for (k = 0; k < n - 1; k++) {
        for (j = k + 1; j < n; j++) A[j][k] = A[j][k] / A[k][k];
        for (i = k + 1; i < n; i++) {
            for (j = k + 1; j < n; j++) {
                A[i][j] -= A[i][k] * A[k][j];
            }
        }
    }
}

double* solve_lower(double **A, double *b, int n) {
    double *sol = alloc_vec(n);
    double tmp_sum;
    int i, j;
    for (i = 0; i < n; i++) {
        tmp_sum = 0.0;
        for (j = 0; j < i; j++) {
            tmp_sum += A[i][j] * sol[j];
        }
        // implicit A[i][i] = 1.0
        sol[i] = (b[i] - tmp_sum);
    }
    return sol;
}

double* solve_upper(double **A, double *b, int n) {
    double *sol = alloc_vec(n);
    double tmp_sum;
    int i, j;
    for (i = n - 1; i >= 0; i--) {
        tmp_sum = 0.0;
        for (j = i + 1; j < n; j++) {
            tmp_sum += A[i][j] * sol[j];
        }
        sol[i] = (b[i] - tmp_sum) / A[i][i];
    }
    return sol;
}