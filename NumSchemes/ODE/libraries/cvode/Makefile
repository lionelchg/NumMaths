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

prefix     = /softs/local_intel/sundials/5.7.0
includedir = ${prefix}/fortran
libdir     = ${prefix}/lib64

F90      = /softs/intel/compilers_and_libraries_2018.1.163/linux/bin/intel64/ifort
F90FLAGS = -O3
F90LIBS  =  -lm /usr/lib64/librt.so

# ------------------------------------------------------------------------------

INCLUDES  = -I${includedir}
LIBRARIES = -lsundials_fcvode_mod -lsundials_cvode ${F90LIBS}
LINKFLAGS = -Wl,-rpath,/softs/local_intel/sundials/5.7.0/lib64

INCLUDES_KLU  = ${INCLUDES} -I
LIBRARIES_KLU = ${LIBRARIES} \
                -lsundials_sunmatrixsparse -lsundials_fsunmatrixsparse_mod \
                -lsundials_sunlinsolklu -lsundials_fsunlinsolklu_mod \

LINKFLAGS_KLU = ${LINKFLAGS}:

# ------------------------------------------------------------------------------

EXAMPLES =  fcv_analytical_dns fcv_brusselator_dns_f2003
EXAMPLES_DEPENDENCIES =

OBJECTS = ${EXAMPLES:=.o}
OBJECTS_DEPENDENCIES = ${EXAMPLES_DEPENDENCIES:=.o}

# KLU Examples
EXAMPLES_KLU =
OBJECTS_KLU  = ${EXAMPLES_KLU:=.o}

# ------------------------------------------------------------------------------

.SUFFIXES : .o .f90

.f90.o :
	${F90} ${F90FLAGS} ${INCLUDES} -c $<

# ------------------------------------------------------------------------------

all: examples examples_klu

examples: ${OBJECTS_DEPENDENCIES} ${OBJECTS}
	@for i in ${EXAMPLES} ; do \
	  echo "${F90} -o $${i} $${i}.o ${OBJECTS_DEPENDENCIES} ${F90FLAGS} ${INCLUDES} -L${libdir} ${LIBRARIES} ${LINKFLAGS}" ; \
	  ${F90} -o $${i} $${i}.o ${OBJECTS_DEPENDENCIES} ${F90FLAGS} ${INCLUDES} -L${libdir} ${LIBRARIES} ${LINKFLAGS} ; \
	done

examples_klu: ${OBJECTS_KLU}
	@for i in ${EXAMPLES_KLU} ; do \
	  echo "${F90} -o $${i} $${i}.o ${F90FLAGS} ${INCLUDES_KLU} -L${libdir} ${LIBRARIES_KLU} ${LINKFLAGS_KLU}" ; \
	  ${F90} -o $${i} $${i}.o ${F90FLAGS} ${INCLUDES_KLU} -L${libdir} ${LIBRARIES_KLU} ${LINKFLAGS_KLU} ; \
	done

clean:
	rm -f *.o *.mod
	rm -f ${OBJECTS} ${OBJECTS_KLU}
	rm -f ${EXAMPLES} ${EXAMPLES_KLU}

# ------------------------------------------------------------------------------
