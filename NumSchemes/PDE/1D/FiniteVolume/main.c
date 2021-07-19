#include "stdio.h"
#include "stdlib.h"
#include "profiles.h"
#include "sys/stat.h"

int main(int argc, char* argv[]) {
    int nschemes = 2;
    char *schemes[2] = {"LW", "SOU"};
    char *data_dir = "data";
    mkdir(data_dir, 0777);

    // Physical properties
    double xmin = -1.0, xmax = 1.0;
    double conv_speed = 1.0;

    // CFLs and mesh
    double cfls[5] = {0.1, 0.3, 0.5, 0.7, 0.9};
    int nnx = 201;

    // Running simulation for a given cfl
    double cfl = cfls[0];
    double Lx = xmax - xmin;
    int ncx = nnx - 1;
    double dx = (xmax - xmin) / (double)(ncx);
    double x0 = (xmax + xmin) / 2;
    double *x = linspace(xmin, xmax, nnx);
    double dt = dx * cfl / conv_speed;

    // Create 4 profiles
    double *u_gauss = gaussian(x, x0, 0.3, nnx);
    double *u_step = step(x, x0, 1.0, nnx);
    double *u_2pw = packet_wave(x, x0, 0.5, 1.0, nnx);
    double *u_4pw = packet_wave(x, x0, 0.25, 1.0, nnx);

    // Iterate using scheme

    // Print results
    double *results[5];
    results[0] = x;
    results[1] = u_gauss;
    results[2] = u_step;
    results[3] = u_2pw;
    results[4] = u_4pw;
    char *vec_names[5] = {"position", "Gaussian", "Step", "Sin_2", "Sin_4"};
    write_vecs(results, nnx, 5, vec_names, "data/results.dat");
    
    return 0;
}