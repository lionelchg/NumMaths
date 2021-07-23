#include "stdio.h"
#include "stdlib.h"
#include "sys/stat.h"
#include "profiles.h"
#include "schemes.h"
#include "defs.h"
#include "string.h"

int main(int argc, char** argv) {
    // All arguments of CLI are schemes
    int nschemes = argc - 1;
    char **schemes = argv + 1;
    char *data_dir = "data/";
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

    // Number of iterations
    double n_periods = 1.0;
    int nt = (int)(n_periods * Lx / conv_speed / dt);

    // Print information
    printf("cfl = %.2f - nits = %d - dx = %.2e - dt = %.2e\n", cfl, nt, dx, dt);

    // Iterate on schemes
    int index_scheme;
    for (index_scheme = 0; index_scheme < nschemes; index_scheme++) {
        // printf("Running %s schemen", schemes[index_scheme]);
        // Create 4 profiles
        double *u_gauss = gaussian(x, x0, 0.3, nnx);
        double *u_step = step(x, x0, 1.0, nnx);
        double *u_2pw = packet_wave(x, x0, 0.5, 1.0, nnx);
        double *u_4pw = packet_wave(x, x0, 0.25, 1.0, nnx);

        // Iterate using scheme
        int* ischeme = scheme_id(schemes[index_scheme]); 
        rungekutta(u_gauss, nnx, dx, dt, conv_speed, ischeme, nt);
        rungekutta(u_step, nnx, dx, dt, conv_speed, ischeme, nt);
        rungekutta(u_2pw, nnx, dx, dt, conv_speed, ischeme, nt);
        rungekutta(u_4pw, nnx, dx, dt, conv_speed, ischeme, nt);

        // Print results
        double *results[5];
        results[0] = x;
        results[1] = u_gauss;
        results[2] = u_step;
        results[3] = u_2pw;
        results[4] = u_4pw;
        char *vec_names[5] = {"position", "Gaussian", "Step", "Sin_2", "Sin_4"};
        
        char data_fn[lenstr];
        strcpy(data_fn, data_dir);
        strcat(data_fn, schemes[index_scheme]);
        strcat(data_fn, ".dat");
        write_vecs(results, nnx, 5, vec_names, data_fn);
    }
    return 0;
}