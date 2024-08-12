AC_INIT([llt],1.0.0)
AC_CONFIG_SRCDIR([src/llt.h])

LLT_MAJOR_VERSION=1
LLT_MINOR_VERSION=0
LLT_MICRO_VERSION=0
LLT_INTERFACE_AGE=0
LLT_BINARY_AGE=0
LLT_VERSION=$LLT_MAJOR_VERSION.$LLT_MINOR_VERSION.$LLT_MICRO_VERSION

AC_SUBST(LLT_MAJOR_VERSION)
AC_SUBST(LLT_MINOR_VERSION)
AC_SUBST(LLT_MICRO_VERSION)
AC_SUBST(LLT_VERSION)

# libtool versioning
LT_RELEASE=$LLT_MAJOR_VERSION.$LLT_MINOR_VERSION
LT_CURRENT=`expr $LLT_MICRO_VERSION - $LLT_INTERFACE_AGE`
LT_REVISION=$LLT_INTERFACE_AGE
LT_AGE=`expr $LLT_BINARY_AGE - $LLT_INTERFACE_AGE`
AC_SUBST(LT_RELEASE)
AC_SUBST(LT_CURRENT)
AC_SUBST(LT_REVISION)
AC_SUBST(LT_AGE)

# For automake.
VERSION=$LLT_VERSION
PACKAGE=llt

AC_SUBST(PACKAGE)
AC_SUBST(VERSION)

AM_INIT_AUTOMAKE(subdir-objects)

##$PACKAGE, $VERSION, no-define)

# Specify a configuration file
AM_CONFIG_HEADER(config.h)

AC_DEFINE_UNQUOTED(LLT_MAJOR_VERSION, $LLT_MAJOR_VERSION, [Major version])
AC_DEFINE_UNQUOTED(LLT_MINOR_VERSION, $LLT_MINOR_VERSION, [Minor version])
AC_DEFINE_UNQUOTED(LLT_MICRO_VERSION, $LLT_MICRO_VERSION, [Micro version])
AC_DEFINE_UNQUOTED(LLT_INTERFACE_AGE, $LLT_INTERFACE_AGE, [Interface age])
AC_DEFINE_UNQUOTED(LLT_BINARY_AGE, $LLT_BINARY_AGE, [Binary age])

dnl Initialize libtool
AM_PROG_LIBTOOL

dnl Initialize maintainer mode
#AM_MAINTAINER_MODE

AC_PROG_CC
AM_PROG_CC_C_O

if test x$GCC = xyes ; then
  CFLAGS="$CFLAGS -Wall -Werror-implicit-function-declaration -Wstrict-prototypes -Wmissing-prototypes -Wmissing-declarations"
fi

AC_PROG_AWK
AC_ISC_POSIX
AC_HEADER_STDC

AC_PATH_PROG(XFOIL, xfoil, no)
if test x$XFOIL = xno ; then
  echo \#################################################################
  echo \#
  echo \# XFOIL aerofoil design and analysis program not on the path
  echo \#
  echo \# It is not required for LLT, but is recommended. You can get it
  echo \# from
  echo \#
  echo \# https://web.mit.edu/drela/Public/web/xfoil/
  echo \#
  echo \# or \(more convenient for installation\)
  echo \#
  echo \# https://github.com/RobotLocomotion/xfoil
  echo \#
  echo \#################################################################
else
   AC_DEFINE(HAVE_XFOIL, 1,
   [XFOIL aerofoil design and analysis program available])
fi

AC_CHECK_PROG(BLASWRAP, blaswrap-config, yes, no)
if test x$BLASWRAP = xyes ; then
   AC_DEFINE(HAVE_BLASWRAP, 1, [BLAS wrapper header available])
   blaswrap_cflags=`blaswrap-config --cflags`
   blaswrap_libs=`blaswrap-config --libs`
else
  echo \#################################################################
  echo \#
  echo \# BLAS wrapper library not installed or not properly installed.
  echo \# You can get it from
  echo \#
  echo \# https://github.com/mjcarley/blaswrap
  echo \#
  echo \#################################################################
  exit
fi

AC_CHECK_LIB(matheval, evaluator_create, [], [
  echo \#################################################################
  echo \# 
  echo \# libmatheval not installed, symbolic differentiation of motion
  echo \# expressions will not be available. You can install it from
  echo \#
  echo \# https://www.gnu.org/software/libmatheval/
  echo \#
  echo \#################################################################
])

PKG_CHECK_MODULES([GLIB], [glib-2.0], [],
[
  echo \#################################################################
  echo \# 
  echo \# GLIB-2.0 not installed. You can install it using your 
  echo \# distribution's package manager or from:
  echo \#
  echo \# https://gitlab.gnome.org/GNOME/glib/
  echo \#
  echo \#################################################################
  exit
])

PKG_CHECK_MODULES([AMC], [amc], [],
[
  echo \#################################################################
  echo \# 
  echo \# AMC \(Affine Motion Calculation\) not installed. You can get it
  echo \# from:
  echo \#
  echo \# https://github.com/mjcarley/amc
  echo \#
  echo \#################################################################
  exit
])

AX_GCC_X86_CPU_SUPPORTS(avx, hasavx=yes) 
AX_GCC_X86_CPU_SUPPORTS(avx2, hasavx2=yes) 
AX_GCC_X86_CPU_SUPPORTS(fma, hasfma=yes) 

if test x$hasavx = xyes; then
   SIMD_FLAGS="$SIMD_FLAGS -mavx -DWBFMM_USE_AVX"
fi
if test x$hasfma = xyes; then
   SIMD_FLAGS="$SIMD_FLAGS -mfma -DWBFMM_USE_AVX"
fi

AC_CHECK_LIB(m, cos)
AC_CHECK_LIB(matheval, evaluator_evaluate)
AC_CONFIG_MACRO_DIR([m4])

AC_OPENMP
AC_PROG_F77
AC_F77_LIBRARY_LDFLAGS
AC_F77_DUMMY_MAIN
AC_F77_WRAPPERS

#DX_DOXYGEN_FEATURE(ON)
##DX_DOT_FEATURE(ON)
#DX_HTML_FEATURE(ON)
##DX_CHM_FEATURE(OFF)
##DX_CHI_FEATURE(OFF)
##DX_MAN_FEATURE(OFF)
##DX_RTF_FEATURE(OFF)
##DX_XML_FEATURE(OFF)
##DX_PDF_FEATURE(OFF)
##DX_PS_FEATURE(OFF)
#DX_INIT_DOXYGEN(llt, doc/llt.dxy)

CPPFLAGS="$glib_cflags $blaswrap_cflags"
CFLAGS="$CFLAGS $OPENMP_CFLAGS $SIMD_FLAGS $AMC_CFLAGS $GLIB_CFLAGS $blaswrap_cflags"
LIBS="$AMC_LIBS $GLIB_LIBS $blaswrap_libs $LIBS"

AC_SUBST(CPPFLAGS)
AC_SUBST(CFLAGS)
AC_SUBST(LDFLAGS)
AC_SUBST(XFOIL)

LIBS="$glib_libs $LIBS $FLIBS"

AC_CONFIG_FILES([
Makefile
llt.pc
src/Makefile
src/llt-xfoil
tools/Makefile])

AC_OUTPUT

