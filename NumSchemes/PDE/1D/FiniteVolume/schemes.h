#ifndef _SCHEMES_
#define _SCHEMES_
#include "utils.h"
// #include "defs.h"

// Retrieve the scheme identifier in the code
int* scheme_id(char* scheme);

// Apply the temporal scheme (temporal iterations)
void rungekutta(double *u, int nnodes, double dx, double dt,
                double a, int *ischeme, int nt);

// The numerical scheme (euler time integration)
double fv_expl_scheme(double *u, int nnodes, int i, double sigma, int ischeme);

// General numerical scheme
double fv_scheme(double *u, int nnodes, int i, int ischeme);

// Return the generalized index for periodic boundary conditions in 1D
int periodic_index(int index, int nnodes);

// List of schemes
double LW_flux(double *u, int nnodes, int i, double sigma);
double FOU_flux(double *u, int nnodes, int i);

#endif