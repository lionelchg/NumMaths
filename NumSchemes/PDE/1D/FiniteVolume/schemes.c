#include "schemes.h"

int* scheme_id(char* scheme, char* limiter) {
    int* ischeme = (int *)malloc(3 * sizeof(int));

    // Scheme parsing
    if (strcmp(scheme, "LW") == 0) {
        ischeme[0] = 0;
        ischeme[1] = 0;
    }
    else if (strcmp(scheme, "WB") == 0) {
        ischeme[0] = 0;
        ischeme[1] = 1;
    }
    else if (strcmp(scheme, "FOU") == 0) {
        ischeme[0] = 1;
        ischeme[1] = 0;
    }
    else if (strcmp(scheme, "Fromm") == 0) {
        ischeme[0] = 1;
        ischeme[1] = 1;
    }
    else if (strcmp(scheme, "Quick") == 0) {
        ischeme[0] = 1;
        ischeme[1] = 2;
    }
    else if (strcmp(scheme, "Third") == 0) {
        ischeme[0] = 1;
        ischeme[1] = 3;
    }
    else if (strcmp(scheme, "SOU") == 0) {
        ischeme[0] = 1;
        ischeme[1] = 4;
    }

    // Limiter parsing
    if (strcmp(limiter, "no") == 0) {
        ischeme[2] = 0;
    }
    else if (strcmp(limiter, "vanleer") == 0) {
        ischeme[2] = 1;
    }
    else if (strcmp(limiter, "superbee") == 0) {
        ischeme[2] = 2;
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
                res[i] = cfl * (fv_expl_scheme(u, nnodes, i, cfl, ischeme) 
                        - fv_expl_scheme(u, nnodes, i - 1, cfl, ischeme));
            }
            else if (ischeme[0] == 1) {
                res[i] = cfl * (fv_scheme(u, nnodes, i, ischeme) 
                        - fv_scheme(u, nnodes, i - 1, ischeme));
            }
        }
        // Apply the residual
        for (i = 0; i < nnodes; i++) {
            u[i] -= res[i];
        }
    }
}

double fv_expl_scheme(double *u, int nnodes, int i, double sigma, int* ischeme) {
    if (ischeme[1] == 0) {
        return LW_flux(u, nnodes, i, sigma, ischeme[2]);
    }
    else if (ischeme[1] == 1) {
        return WB_flux(u, nnodes, i, sigma, ischeme[2]);
    }
}

double LW_flux(double *u, int nnodes, int i, double sigma, int ilim) {
    int iperio = periodic_index(i, nnodes);
    int i1 = periodic_index(i + 1, nnodes);
    int im1 = periodic_index(i - 1, nnodes);
    double psi = 1.0;
    if (ilim == 1) 
        psi = van_leer(grad_ratio(u, iperio, im1, i1));
    else if (ilim == 2)
        psi = superbee(grad_ratio(u, iperio, im1, i1));
    return u[iperio] + 0.5 * (1.0 - sigma) * (u[i1] - u[iperio]) * psi;
}

double WB_flux(double *u, int nnodes, int i, double sigma, int ilim) {
    int iperio = periodic_index(i, nnodes);
    int i1 = periodic_index(i + 1, nnodes);
    int im1 = periodic_index(i - 1, nnodes);
    double psi = 1.0;
    if (ilim == 1) 
        psi = van_leer(grad_ratio(u, iperio, im1, i1));
    else if (ilim == 2)
        psi = superbee(grad_ratio(u, iperio, im1, i1));
    return u[iperio] + 0.5 * (1.0 - sigma) * (u[iperio] - u[im1]) * psi;
}


double fv_scheme(double *u, int nnodes, int i, int* ischeme) {
    if (ischeme[1] == 0) {
        return FOU_flux(u, nnodes, i);
    }
    else if (ischeme[1] == 1) {
        return kappa_flux(u, nnodes, i, 0.0, ischeme[2]);
    }
    else if (ischeme[1] == 2) {
        return kappa_flux(u, nnodes, i, 0.5, ischeme[2]);
    }
    else if (ischeme[1] == 3) {
        return kappa_flux(u, nnodes, i, 1.0 / 3.0, ischeme[2]);
    }
    else if (ischeme[1] == 4) {
        return kappa_flux(u, nnodes, i, -1.0, ischeme[2]);
    }
}

double FOU_flux(double *u, int nnodes, int i) {
    int iperio = periodic_index(i, nnodes);
    return u[iperio];
}

double kappa_flux(double *u, int nnodes, int i, double kappa, int ilim) {
    int im1 = periodic_index(i - 1, nnodes);
    int iperio = periodic_index(i, nnodes);
    int i1 = periodic_index(i + 1, nnodes);
    double psi = 1.0;
    if (ilim == 1)
        psi = van_leer(grad_ratio(u, iperio, i1, im1));
    else if (ilim == 2)
        psi = superbee(grad_ratio(u, iperio, i1, im1));
    return u[iperio] + 0.25 * ((1.0 + kappa) * (u[i1] - u[iperio]) 
                    + (1.0 - kappa) * (u[iperio] - u[im1])) * psi;
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

double grad_ratio(double *u, int i, int i1, int im1) {
    double r_i;
    if (u[i1] == u[i]) {
        r_i = 0.0;
    }
    else if (u[i] == u[im1]) {
        r_i = 10.0 * (u[i1] - u[i]) / fabs(u[i1] - u[i]);
    }
    else {
        r_i = (u[i1] - u[i]) / (u[i] - u[im1]);
    }
    return r_i;
}

double van_leer(double r) {
    if (r <= 0) return 0.0;
    else return (r + fabs(r)) / (1.0 + r);
}

double superbee(double r) {
    if (r <= 0) return 0.0;
    else return fmax(fmin(2 * r, 1.0), fmin(r, 2.0));
}