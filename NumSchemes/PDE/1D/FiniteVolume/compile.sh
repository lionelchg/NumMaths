#!/bin/bash
HDF5_DIR=/usr/local/Cellar/hdf5/1.12.0_4
gcc-11 -I$HDF5_DIR/include -L$HDF5_DIR/lib -Wall -o main.x profiles.c schemes.c utils.c main.c libargp.a -lhdf5