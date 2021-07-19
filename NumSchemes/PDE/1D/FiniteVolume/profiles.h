#ifndef _PROFILES_
#define _PROFILES_

#include "utils.h"

double* gaussian(double *x, double x0, double sigma_x, int n);

double* step(double *x, double x0, double sigma_x, int n);

double* packet_wave(double *x, double x0, double lam, double sigma_x, int n);

#endif