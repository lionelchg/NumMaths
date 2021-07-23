#include "schemes.h"

int* scheme_id(char* scheme) {
    int* ischeme = (int *)malloc(2 * sizeof(int));
    if (scheme == "LW") {
        ischeme[0] = 0;
        ischeme[1] = 0;
    }
    else if (scheme == "FOU") {
        ischeme[0] = 1;
        ischeme[1] = 0;
    }
    return ischeme;
}

void rungekutta(double *u, int nnodes, double dx, double dt, 
                double a, int *ischeme, int nt) {
    int i, t;
    double cfl = a * dt / dx;
    double *res = alloc_vec(nnodes);
    // Temporal loop
    for (t = 0; t < nt; t++) {
        // Spatial loop to compute residual
        for (i = 0; i < nnodes; i++) {
            if (ischeme[0] == 0) {
                res[i] = cfl * (fv_expl_scheme(u, nnodes, i, cfl, ischeme[1]) 
                        - fv_expl_scheme(u, nnodes, i - 1, cfl, ischeme[1]));
            }
            else if (ischeme[0] == 1) {
                res[i] = cfl * (fv_scheme(u, nnodes, i, ischeme[1]) 
                        - fv_scheme(u, nnodes, i - 1, ischeme[1]));
            }
        }
        // Apply the residual
        for (i = 0; i < nnodes; i++) {
            u[i] -= res[i];
        }
    }
}

double fv_expl_scheme(double *u, int nnodes, int i, double sigma, int ischeme1) {
    if (ischeme1 == 0) {
        return LW_flux(u, nnodes, i, sigma);
    }
}

double LW_flux(double *u, int nnodes, int i, double sigma) {
    int iperio = periodic_index(i, nnodes);
    int i1 = periodic_index(i + 1, nnodes);
    return 0.5 * (u[iperio] + u[i1]) - 0.5 * sigma * (u[i1] - u[iperio]);
}

double fv_scheme(double *u, int nnodes, int i, int ischeme1) {
    if (ischeme1 == 0) {
        return FOU_flux(u, nnodes, i);
    }
}

double FOU_flux(double *u, int nnodes, int i) {
    int iperio = periodic_index(i, nnodes);
    return u[iperio];
}

int periodic_index(int index, int nnodes) {
    if (index < 0) {
        return index - 1 + nnodes;
    }
    else if (index > nnodes - 1) {
        return index + 1 - nnodes;
    }
    else {
        return index;
    }
}
