#include "schemes.h"

void rungekutta(double *u, int nnodes, double dx, double dt, 
                double a, int *ischeme, int nt) {
    int i;
    double cfl = a * dt / dx;
    double *res = alloc_vec(nnodes);
    for (i = 0; i < nt; i++) {
        if (ischeme[0] == 0) {
            res[i] = cfl * (fv_expl_scheme(u, nnodes, i, cfl, ischeme[1]) 
                    - fv_expl_scheme(u, nnodes, i - 1, cfl, ischeme[1]));
        }
        else if (ischeme[1] == 1) {
            res[i] = cfl * (fv_scheme(u, nnodes, i, ischeme[1]) 
                    - fv_scheme(u, nnodes, i - 1, ischeme[1]));
        }
        u[i] -= res[i];
    }
}

double fv_expl_scheme(double *u, int nnodes, int i, double sigma, int ischeme) {
    if (ischeme == 0) {
        return LW_flux(u, nnodes, i, sigma);
    }
}

double LW_flux(double *u, int nnodes, int i, double sigma) {
    int i1 = periodic_index(i + 1, nnodes);
    return 0.5 * (u[i] + u[i1]) - 0.5 * sigma * (u[i1] - u[i]);
}

double fv_scheme(double *u, int nnodes, int i, int ischeme) {
    if (ischeme == 0) {
        return FOU_flux(u, nnodes, i);
    }
}

double FOU_flux(double *u, int nnodes, int i) {
    return u[i];
}

int periodic_index(int index, int nnodes) {
    if (index < 0) {
        return (index - 1) % nnodes;
    }
    else if (index > nnodes - 1) {
        return (index + 1) % nnodes;
    }
    else {
        return index;
    }
}
