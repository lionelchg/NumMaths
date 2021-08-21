#include "utils.h"

void print_vec(double *vec, int n, char *vec_name)
{
    int i, j;
    int step = n / 5;
    printf("%s\n", vec_name);
    for (i = 0; i < n; i += step)
    {
        printf("%.3e ", vec[i]);
    }
    printf("\n");
}

void write_vecs(double **list_vec, int n, int m, char **vec_names, char *filename)
{
    int i, j;
    FILE *fptr;

    // use appropriate location if you are using MacOS or Linux
    printf("Print results to %s file\n", filename);
    fptr = fopen(filename, "w");

    if (fptr == NULL)
    {
        printf("Error while opening file");
        exit(1);
    }
    // Write header
    for (j = 0; j < m; j++) {
        fprintf(fptr, "%13s", vec_names[j]);
    }
    fprintf(fptr, "\n");

    // Write vector contents
    for (i = 0; i < n; i++) {
        for (j = 0; j < m; j++) {
            fprintf(fptr, "%13.4e", list_vec[j][i]);
        }
        fprintf(fptr, "\n");
    }

    fclose(fptr);
}

double *alloc_vec(int n)
{
    int i;
    double *vec;
    vec = (double *)malloc(n * sizeof(double));
    for (i = 0; i < n; i++)
    {
        vec[i] = 0.0;
    }
    return vec;
}

double *copy_vec(double *vec, int n)
{
    int i;
    double *copied_vec = alloc_vec(n);
    for (i = 0; i < n; i++)
    {
        copied_vec[i] = vec[i];
    }
    return copied_vec;
}

double *linspace(double begin, double end, int npoints)
{
    double *array = alloc_vec(npoints);
    int i;
    double h = (end - begin) / (double)(npoints - 1);
    for (i = 0; i < npoints; i++)
    {
        array[i] = begin + h * i;
    }
    return array;
}

hid_t create_group(hid_t file, char *grpname, double cfl) {
    // HDF5 variables
    hid_t group, space, attr;
    herr_t status;

    // Local variables
    char *attrcfl = "cfl";


    group = H5Gcreate (file, grpname, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    
    /*
    * Create float attribute
    */
    space = H5Screate(H5S_SCALAR);
    attr = H5Acreate (group, attrcfl, H5T_IEEE_F64LE, space, H5P_DEFAULT,
                H5P_DEFAULT);
    status = H5Awrite (attr, H5T_NATIVE_DOUBLE, &cfl);

    /*
    * Close and release resources.
    */
    status = H5Aclose (attr);
    status = H5Sclose (space);

    return group;
}

void write_dset_1d(hid_t group, char *dsetname, int dim, double *wdata) {
    // HDF5 variables
    hid_t space, dset;
    herr_t status;
    hsize_t dims[1] = {dim};

    // Local variables
    char *attrscheme = "scheme";
    int i, j;
    double wbuffer[dim];

    // Init data
    for (i = 0; i < dim; i++)
        wbuffer[i] = wdata[i];

    // Create dataset
    space = H5Screate_simple(1, dims, NULL);
    dset = H5Dcreate(group, dsetname, H5T_IEEE_F64LE, space, H5P_DEFAULT,
                    H5P_DEFAULT, H5P_DEFAULT);
    status = H5Dwrite(dset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT,
                    wbuffer);

    // Close and release resources.
    status = H5Dclose (dset);
    status = H5Sclose (space);
}

void write_dset_2d(hid_t group, char *dsetname, char *scheme, int dim0, int dim1, double **wdata) {
    // HDF5 variables
    hid_t space, dset, atype, attr;
    herr_t status;
    hsize_t dims[2] = {dim0, dim1};

    // Local variables
    char *attrscheme = "scheme";
    int i, j;
    double wbuffer[dim0][dim1];

    // Init data
    for (i = 0; i < dim0; i++)
        for (j = 0; j < dim1; j++)
            wbuffer[i][j] = wdata[i][j];

    // Create dataset
    space = H5Screate_simple(2, dims, NULL);
    dset = H5Dcreate(group, dsetname, H5T_IEEE_F64LE, space, H5P_DEFAULT,
                    H5P_DEFAULT, H5P_DEFAULT);
    status = H5Dwrite(dset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT,
                    wbuffer[0]);

    // Create attribute for the name of the scheme
    space = H5Screate(H5S_SCALAR);
    atype = H5Tcopy(H5T_C_S1);
    H5Tset_size(atype, 10);
    attr = H5Acreate (dset, attrscheme, atype, space, H5P_DEFAULT,
                H5P_DEFAULT);
    status = H5Awrite (attr, atype, scheme);

    // Close and release resources.
    status = H5Aclose (attr);
    status = H5Dclose (dset);
    status = H5Sclose (space);
}
