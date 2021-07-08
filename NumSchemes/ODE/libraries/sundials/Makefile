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

prefix       = /Users/cheng/opt/sundials/instdir
exec_prefix  = /Users/cheng/opt/sundials/instdir
includedir   = /Users/cheng/opt/sundials/instdir/include
libdir       = /Users/cheng/opt/sundials/instdir/lib

CPP      = /usr/local/bin/g++-10
CPPFLAGS = -O3 -DNDEBUG
CC       = /usr/local/bin/gcc-10
CFLAGS   = -O3 -DNDEBUG
LDFLAGS  = 
LIBS     =  -lm

LINKFLAGS = -Wl,-rpath,/Users/cheng/opt/sundials/instdir/lib

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

EXAMPLES =  cvRoberts_dns cvExample_dns
EXAMPLES_DEPENDENCIES = 

OBJECTS = ${EXAMPLES:=.o}
OBJECTS_DEPENDENCIES = ${EXAMPLES_DEPENDENCIES:=.o}

# -----------------------------------------------------------------------------------------

.SUFFIXES : .o .c

.c.o :
	${CC} ${CPPFLAGS} ${CFLAGS} ${INCLUDES} -c $<

# -----------------------------------------------------------------------------------------

all: examples

examples: ${OBJECTS}
	@for i in ${EXAMPLES} ; do \
	  echo "${CC} -o $${i} $${i}.o ${OBJECTS_DEPENDENCIES} ${CFLAGS} ${LDFLAGS} ${INCLUDES} -L${libdir} ${LIBRARIES} ${LINKFLAGS}" ; \
	  ${CC} -o $${i} $${i}.o ${OBJECTS_DEPENDENCIES} ${CFLAGS} ${LDFLAGS} ${INCLUDES} -L${libdir} ${LIBRARIES} ${LINKFLAGS} ; \
	done

${OBJECTS}: ${OBJECTS_DEPENDENCIES}


clean:
	rm -f ${OBJECTS_DEPENDENCIES} ${OBJECTS_DEPENDENCIES_BL} ${OBJECTS_DEPENDENCIES_KLU} ${OBJECTS_DEPENDENCIES_SLUMT}
	rm -f ${OBJECTS} ${OBJECTS_BL} ${OBJECTS_KLU} ${OBJECTS_SLUMT}
	rm -f ${EXAMPLES} ${EXAMPLES_BL} ${EXAMPLES_KLU} ${EXAMPLES_SLUMT}

# -----------------------------------------------------------------------------------------
