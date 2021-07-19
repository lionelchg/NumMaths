#include "profiles.h"

double* gaussian(double *x, double x0, double sigma_x, int n) {
    double *result = alloc_vec(n);
    int i;
    for (i = 0; i < n; i++) {
        result[i] = exp(- pow((x[i] - x0), 2) / 2 / pow(sigma_x, 2));
    }
    return result;
}

double* step(double *x, double x0, double sigma_x, int n) {
    double *result = alloc_vec(n);
    int i;
    for (i = 0; i < n; i++) {
        if (fabs(x[i] - x0) < 0.5 * sigma_x) {
            result[i] = 1.0;
        }
        else {
            result[i] = 0.0;
        }
    }
    return result;
}


double* packet_wave(double *x, double x0, double lam, double sigma_x, int n) {
    double *result = alloc_vec(n);
    int i;
    for (i = 0; i < n; i++) {
        if (fabs(x[i] - x0) < 0.5 * sigma_x) {
            result[i] = sin(2 * M_PI / lam * (x[i] - x0));
        }
        else {
            result[i] = 0.0;
        }
    }
    return result;
}

