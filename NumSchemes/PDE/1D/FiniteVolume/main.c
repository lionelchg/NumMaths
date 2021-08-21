#include "stdio.h"
#include "stdlib.h"
#include "sys/stat.h"
#include "profiles.h"
#include "schemes.h"
#include "defs.h"
#include "string.h"
#include "parsing.h"

int main(int argc, char** argv) {
    // All arguments of CLI are schemes
    int index_scheme, nschemes;
    char **schemes;
    char *data_dir = "data/";
    mkdir(data_dir, 0777);

    // Parsing CLI
    struct arguments arguments;
    arguments.nschemes = 1;
    arguments.outfile = "solut.dat";
    arguments.verbose = 0;
    argp_parse(&argp, argc, argv, 0, 0, &arguments);

    // Copy to local variables
    nschemes = arguments.nschemes;
    schemes = arguments.args;
    printf("Number of schemes: %d\n", nschemes);
    for (index_scheme = 0; index_scheme < arguments.nschemes; index_scheme++) {
        printf("%s\n", schemes[index_scheme]);
    } 
    printf("Output filename: %s\n", arguments.outfile);

    // Physical properties
    double xmin = -1.0, xmax = 1.0;
    double x0 = (xmax + xmin) / 2;
    double Lx = xmax - xmin;
    double conv_speed = 1.0;

    // CFLs and mesh
    int nnx = 201;
    int ncx = nnx - 1;
    double dx = (xmax - xmin) / (double)(ncx);
    double *x = linspace(xmin, xmax, nnx);
    int icfl;
    double cfls[5] = {0.1, 0.3, 0.5, 0.7, 0.9};
    double n_periods = 1.0;

    // Looping on CFLs
    for (icfl = 0; icfl < 5; icfl++) {
        // Running simulation for a given cfl
        double cfl = cfls[0];
        double dt = dx * cfl / conv_speed;

        // Number of iterations
        int nt = (int)(n_periods * Lx / conv_speed / dt);

        // Print information
        printf("cfl = %.2f - nits = %d - dx = %.2e - dt = %.2e\n", cfl, nt, dx, dt);

        // Iterate on schemes
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
            sprintf(data_fn, "%s%s_cfl_%d.dat", data_dir, schemes[index_scheme], icfl);
            write_vecs(results, nnx, 5, vec_names, data_fn);
        }
    }
    return 0;
}