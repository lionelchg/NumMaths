# -*- mode: makefile -*-
# -----------------------------------------------------------------
# Programmer(s): Radu Serban @ LLNL
#                David J. Gardner @ LLNL
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
# Makefile for CVODE serial examples
#
# This file is generated from a template using variables
# set at configuration time. It can be used as a template for
# other user Makefiles.
# -----------------------------------------------------------------

SHELL = sh

prefix       = /softs/local_intel/sundials/5.7.0
exec_prefix  = /softs/local_intel/sundials/5.7.0
includedir   = /softs/local_intel/sundials/5.7.0/include
libdir       = /softs/local_intel/sundials/5.7.0/lib64

CPP      = /softs/intel/compilers_and_libraries_2018.1.163/linux/bin/intel64/icc
CPPFLAGS = -O3 -DNDEBUG
CC       = /softs/intel/compilers_and_libraries_2018.1.163/linux/bin/intel64/icc
CFLAGS   = -O3 -DNDEBUG
LDFLAGS  = 
LIBS     =  -lm /usr/lib64/librt.so

LINKFLAGS = -Wl,-rpath,/softs/local_intel/sundials/5.7.0/lib64

# -----------------------------------------------------------------------------------------

LIBRARIES_LAPACK = -lsundials_sunlinsollapackdense -lsundials_sunlinsollapackband 
LINKFLAGS_LAPACK = ${LINKFLAGS}:

INCLUDES_KLU  = 
LIBRARIES_KLU = -lsundials_sunlinsolklu 
LINKFLAGS_KLU = ${LINKFLAGS}:

INCLUDES_SLUMT  = 
LIBRARIES_SLUMT = -lsundials_sunlinsolsuperlumt   
LINKFLAGS_SLUMT = ${LINKFLAGS}:

TMP_INCS  = ${includedir} ${INCLUDES_SLUMT} ${INCLUDES_KLU}
INCLUDES  = $(addprefix -I, ${TMP_INCS})
LIBRARIES = -lsundials_cvode -lsundials_nvecserial -lsundials_nvecmanyvector ${LIBS}

# -----------------------------------------------------------------------------------------

EXAMPLES =  cvAdvDiff_bnd cvDirectDemo_ls cvDiurnal_kry_bp cvDiurnal_kry cvDisc_dns cvKrylovDemo_ls cvKrylovDemo_prec cvParticle_dns cvPendulum_dns cvRoberts_dns cvRoberts_dns_uw cvRoberts_dns_negsol cvRoberts_dns_constraints
EXAMPLES_DEPENDENCIES = 

OBJECTS = ${EXAMPLES:=.o}
OBJECTS_DEPENDENCIES = ${EXAMPLES_DEPENDENCIES:=.o}

# LAPACK Examples
EXAMPLES_BL = 
EXAMPLES_DEPENDENCIES_BL = 

OBJECTS_BL = ${EXAMPLES_BL:=.o}
OBJECTS_DEPENDENCIES_BL = ${EXAMPLES_DEPENDENCIES_BL:=.o}

# KLU Examples
EXAMPLES_KLU = 
EXAMPLES_DEPENDENCIES_KLU = 

OBJECTS_KLU = ${EXAMPLES_KLU:=.o}
OBJECTS_DEPENDENCIES_KLU = ${EXAMPLES_DEPENDENCIES_KLU:=.o}

# SuperLU_MT Examples
EXAMPLES_SLUMT = 
EXAMPLES_DEPENDENCIES_SLUMT = 

OBJECTS_SLUMT = ${EXAMPLES_SLUMT:=.o}
OBJECTS_DEPENDENCIES_SLUMT = ${EXAMPLES_DEPENDENCIES_SLUMT:=.o}

# -----------------------------------------------------------------------------------------

.SUFFIXES : .o .c

.c.o :
	${CC} ${CPPFLAGS} ${CFLAGS} ${INCLUDES} -c $<

# -----------------------------------------------------------------------------------------

all: examples examples_bl examples_klu examples_slumt

examples: ${OBJECTS}
	@for i in ${EXAMPLES} ; do \
	  echo "${CC} -o $${i} $${i}.o ${OBJECTS_DEPENDENCIES} ${CFLAGS} ${LDFLAGS} ${INCLUDES} -L${libdir} ${LIBRARIES} ${LINKFLAGS}" ; \
	  ${CC} -o $${i} $${i}.o ${OBJECTS_DEPENDENCIES} ${CFLAGS} ${LDFLAGS} ${INCLUDES} -L${libdir} ${LIBRARIES} ${LINKFLAGS} ; \
	done

examples_bl: ${OBJECTS_BL}
	@for i in ${EXAMPLES_BL} ; do \
	  echo "${CC} -o $${i} $${i}.o ${OBJECTS_DEPENDENCIES_BL} ${CFLAGS} ${LDFLAGS} ${INCLUDES} -L${libdir} ${LIBRARIES} ${LIBRARIES_LAPACK} ${LINKFLAGS_LAPACK}" ; \
	  ${CC} -o $${i} $${i}.o ${OBJECTS_DEPENDENCIES_BL} ${CFLAGS} ${LDFLAGS} ${INCLUDES} -L${libdir} ${LIBRARIES} ${LIBRARIES_LAPACK} ${LINKFLAGS_LAPACK} ; \
	done

examples_klu: ${OBJECTS_KLU}
	@for i in ${EXAMPLES_KLU} ; do \
	  echo "${CC} -o $${i} $${i}.o ${OBJECTS_DEPENDENCIES_KLU} ${CFLAGS} ${LDFLAGS} ${INCLUDES} -L${libdir} ${LIBRARIES} ${LIBRARIES_KLU} ${LINKFLAGS_KLU}" ; \
	  ${CC} -o $${i} $${i}.o ${OBJECTS_DEPENDENCIES_KLU} ${CFLAGS} ${LDFLAGS} ${INCLUDES} -L${libdir} ${LIBRARIES} ${LIBRARIES_KLU} ${LINKFLAGS_KLU} ; \
	done

examples_slumt: ${OBJECTS_SLUMT}
	@for i in ${EXAMPLES_SLUMT} ; do \
	  echo "${CC} -o $${i} $${i}.o ${OBJECTS_DEPENDENCIES_SLUMT} ${CFLAGS} ${LDFLAGS} ${INCLUDES} -L${libdir} ${LIBRARIES} ${LIBRARIES_SLUMT} ${LINKFLAGS_SLUMT}" ; \
	  ${CC} -o $${i} $${i}.o ${OBJECTS_DEPENDENCIES_SLUMT} ${CFLAGS} ${LDFLAGS} ${INCLUDES} -L${libdir} ${LIBRARIES} ${LIBRARIES_SLUMT} ${LINKFLAGS_SLUMT} ; \
	done

${OBJECTS}: ${OBJECTS_DEPENDENCIES}

${OBJECTS_BL}: ${OBJECTS_DEPENDENCIES_BL}

${OBJECTS_KLU}: ${OBJECTS_DEPENDENCIES_KLU}

${OBJECTS_SLUMT}: ${OBJECTS_DEPENDENCIES_SLUMT}

clean:
	rm -f ${OBJECTS_DEPENDENCIES} ${OBJECTS_DEPENDENCIES_BL} ${OBJECTS_DEPENDENCIES_KLU} ${OBJECTS_DEPENDENCIES_SLUMT}
	rm -f ${OBJECTS} ${OBJECTS_BL} ${OBJECTS_KLU} ${OBJECTS_SLUMT}
	rm -f ${EXAMPLES} ${EXAMPLES_BL} ${EXAMPLES_KLU} ${EXAMPLES_SLUMT}

# -----------------------------------------------------------------------------------------
