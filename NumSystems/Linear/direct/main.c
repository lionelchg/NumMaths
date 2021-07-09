#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include "linalg.h"

int main (int argc, char *argv[])
{
    // local
    int i, j, k, n;
    double **A, **B, **C;
    clock_t start_t, end_t;
    double total_t;

    // Get the argument for number of nodes
    n = atoi(argv[1]);

    // Allocation
    A = alloc_mat(n);
    B = alloc_mat(n);

    // Give values to A and B
    for (i = 0; i < n; i++) {
        for (j = 0; j < n; j++) {
            A[i][j] = i + j;
            B[i][j] = i * j;
        }
    }

    // C = A x B - row access first
    start_t = clock();
    C = matmul(A, B, n);
    end_t = clock();
    total_t = (double)(end_t - start_t) / CLOCKS_PER_SEC;
    printf("CPU time for multiplication: %3e s\n", total_t);
    print_matrix(C, n, "C = A x B");

    // A = poisson mat
    A = poisson_mat(n);
    print_matrix(A, n, "A = Poisson matrix");

    // Solve the zero rhs problem with 0.0 and 1.0
    // potential at boundaries
    double *rhs = alloc_vec(n);
    double *pot = alloc_vec(n);
    double phi_l = 0.0, phi_r = 1.0;
    rhs[0] = phi_l;
    rhs[n - 1] = phi_r;

    // Transform the matrix into LU form
    LU(A, n);
    print_matrix(A, n, "LU(A)");
    double *tmp_sol;
    tmp_sol = solve_lower(A, rhs, n);
    pot = solve_upper(A, tmp_sol, n);
    print_vec(pot, n, "Solution of Poisson problem:");

    return 0;
}