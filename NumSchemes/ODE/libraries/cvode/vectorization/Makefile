# -*- mode: makefile -*-
# -----------------------------------------------------------------
# Programmer(s): Cody J. Balos @ LLNL
# -----------------------------------------------------------------
# SUNDIALS Copyright Start
# Copyright (c) 2002-2021, Lawrence Livermore National Security
# and Southern Methodist University.
# All rights reserved.
#
# See the top-level LICENSE and NOTICE files for details.
#
# SPDX-License-Identifier: BSD-3-Clause
# SUNDIALS Copyright End
# -----------------------------------------------------------------
# Makefile for CVODE Fortran 2003 interface serial examples
#
# This file is generated from a template using variables
# set at configuration time. It can be used as a template for
# other user Makefiles.
# -----------------------------------------------------------------

SHELL = sh

prefix     = /Users/bogopolsky/THESE/NUMERICS/sundials/instdir
includedir = ${prefix}/fortran
libdir     = ${prefix}/lib

F90      = mpif90
F90FLAGS = -O3 -g -march=native
F90LIBS  =  # -lm /usr/lib64/librt.so

# ------------------------------------------------------------------------------

INCLUDES  = -I${includedir}
LIBRARIES = -lsundials_fcvode_mod -lsundials_cvode ${F90LIBS}
LINKFLAGS = -Wl,-rpath,${prefix}/lib

# ------------------------------------------------------------------------------

EXAMPLES = implicit_solver
EXAMPLES_DEPENDENCIES =

OBJECTS = ${EXAMPLES:=.o}
OBJECTS_DEPENDENCIES = ${EXAMPLES_DEPENDENCIES:=.o}

# ------------------------------------------------------------------------------

.SUFFIXES : .o .f90

.f90.o :
	${F90} ${F90FLAGS} ${INCLUDES} -c $<

# ------------------------------------------------------------------------------

all: examples

examples: ${OBJECTS_DEPENDENCIES} ${OBJECTS}
	@for i in ${EXAMPLES} ; do \
	  echo "${F90} -o $${i} $${i}.o ${OBJECTS_DEPENDENCIES} ${F90FLAGS} ${INCLUDES} -L${libdir} ${LIBRARIES} ${LINKFLAGS}" ; \
	  ${F90} -o $${i} $${i}.o ${OBJECTS_DEPENDENCIES} ${F90FLAGS} ${INCLUDES} -L${libdir} ${LIBRARIES} ${LINKFLAGS} ; \
	done

clean:
	rm -f *.o *.mod
	rm -f ${OBJECTS}
	rm -f ${EXAMPLES}

# ------------------------------------------------------------------------------
