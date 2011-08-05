#                                               -*- Autoconf -*-
# Process this file with autoconf to produce a configure script.

AC_PREREQ([2.61])
#AC_INIT(Imogene, 1.0, [Hervé Rouault rouault@lps.ens.fr Marc Santolini marc.santolini@lps.ens.fr])
AC_INIT(Imogene, m4_esyscmd([tools/git-version-gen .tarball-version]), [rouault@lps.ens.fr marc.santolini@lps.ens.fr])
AM_INIT_AUTOMAKE()
AC_CONFIG_SRCDIR([src/imogene.cpp])
AC_CONFIG_HEADERS([config.h])
AC_CONFIG_MACRO_DIR([m4])
m4_include([m4/ax_blas.m4])


# Check for python
AM_PATH_PYTHON([2.5],,[AC_MSG_ERROR([no python version >=2.5 detected])])

# Checks for programs.
AC_PROG_CXX
AC_PROG_CC
AM_PROG_CC_C_O
AC_PROG_RANLIB

# Checks for libraries.
AC_CHECK_LIB([m],[cos],,[AC_MSG_ERROR([math functions required.])])

##AX_BLAS( [ AC_MSG_NOTICE( [ blas libs variable : x$BLAS_LIBS BLAS_LIBS ] ) ], )
AX_BLAS(LIBS="$BLAS_LIBS $LIBS $FLIBS",AC_CHECK_LIB([gslcblas],[cblas_dgemm],,AC_MSG_ERROR([blas library required.])))

#AX_BLAS([],[AC_MSG_ERROR([blas library required.])])
AC_CHECK_LIB([gsl],[gsl_blas_dgemm],,[AC_MSG_ERROR([libgsl required.])])

# Checks for header files.
AC_CHECK_HEADERS([gsl/gsl_rng.h],[],[AC_MSG_ERROR([gsl headers required.])])

# Checks for typedefs, structures, and compiler characteristics.
AC_HEADER_STDBOOL

# Checks for library functions.
AC_CHECK_FUNCS([floor pow sqrt],[],[AC_MSG_ERROR([floor,pow and sqrt functions required])])

AC_CONFIG_FILES([Makefile
                 scripts/Makefile
                 scripts/python/Makefile
                 scripts/sh/Makefile
                 data/Makefile
                 doc/Makefile
                 src/Makefile])

AC_MSG_NOTICE([ Next step instruction
               Launch the script for alignment generation...
               ])

AC_OUTPUT
