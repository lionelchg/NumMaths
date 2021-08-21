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
    char **schemes, **limiters;
    char *data_dir = "data/";
    mkdir(data_dir, 0777);

    // Parsing CLI
    struct arguments arguments;
    arguments.nschemes = 1;
    arguments.outfile = "data/solut.h5";
    arguments.verbose = 0;
    argp_parse(&argp, argc, argv, 0, 0, &arguments);

    // Copy to local variables
    nschemes = arguments.nschemes;
    schemes = arguments.args;
    limiters = arguments.args + nschemes;
    printf("Number of schemes: %d\n", nschemes);
    for (index_scheme = 0; index_scheme < arguments.nschemes; index_scheme++) {
        printf("%s - %s\n", schemes[index_scheme], limiters[index_scheme]);
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

    // Creation of HDF5 file and related variables
    hid_t file, group;
    herr_t status;
    char *grpname_base = "/cfl_";
    char grpname[lenstr], dsetname[lenstr], scheme_name[lenstr];
    file = H5Fcreate(arguments.outfile, H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);

    // Write common x vector
    sprintf(dsetname, "x", index_scheme);
    write_dset_1d(file, dsetname, nnx, x);

    // Looping on CFLs
    for (icfl = 0; icfl < 5; icfl++) {
        // Create HDF group
        sprintf(grpname, "%s%d", grpname_base, icfl);
        group = create_group(file, grpname, cfls[icfl]);

        // Running simulation for a given cfl
        double cfl = cfls[0];
        double dt = dx * cfl / conv_speed;

        // Number of iterations
        int nt = (int)(n_periods * Lx / conv_speed / dt);

        // Print information
        printf("cfl = %.2f - nits = %d - dx = %.2e - dt = %.2e\n", cfl, nt, dx, dt);

        // Iterate on schemes
        for (index_scheme = 0; index_scheme < nschemes; index_scheme++) {
            // Create 4 profiles
            double *u_gauss = gaussian(x, x0, 0.3, nnx);
            double *u_step = step(x, x0, 1.0, nnx);
            double *u_2pw = packet_wave(x, x0, 0.5, 1.0, nnx);
            double *u_4pw = packet_wave(x, x0, 0.25, 1.0, nnx);

            // Iterate using scheme
            int* ischeme = scheme_id(schemes[index_scheme], limiters[index_scheme]); 
            rungekutta(u_gauss, nnx, dx, dt, conv_speed, ischeme, nt);
            rungekutta(u_step, nnx, dx, dt, conv_speed, ischeme, nt);
            rungekutta(u_2pw, nnx, dx, dt, conv_speed, ischeme, nt);
            rungekutta(u_4pw, nnx, dx, dt, conv_speed, ischeme, nt);

            // Print results
            double *results[4];
            results[0] = u_gauss;
            results[1] = u_step;
            results[2] = u_2pw;
            results[3] = u_4pw;
            
            sprintf(dsetname, "scheme_%d", index_scheme);
            sprintf(scheme_name, "%s_%s", schemes[index_scheme], limiters[index_scheme]);
            write_dset_2d(group, dsetname, scheme_name, 4, nnx, results);
        }
        // Close group
        status = H5Gclose (group);
    }

    // Close HDF5 file
    status = H5Fclose (file);

    return 0;
}