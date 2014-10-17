#!/bin/bash

rm -f *.x most_* F90_*

# check for FORTRAN-90 compiler
# put your favourite compiler in front if you have more than one
# WARNING: Portland Group Compiler pgf is reported to be buggy

for MOST_F90 in ifort gfortran xlf sunf90 nagfor f90 f95 pgf90 g95 af90 NO_F90
do
   `hash 2>/dev/null $MOST_F90`
   if [ $? = 0 ] ; then break ; fi
done

if [ $MOST_F90 != "NO_F90" ] ; then
   F90_PATH=`which $MOST_F90`
   echo >  most_compiler MOST_F90=$MOST_F90 " # " $F90_PATH
   echo >> most_compiler "MPIMOD=mpimod_stub"
   echo "Found FORTRAN-90 compiler     at: $F90_PATH"
   if [ $MOST_F90 = "sunf90" ] ; then
      echo >> most_compiler "MOST_F90_OPTS=-O3"
      echo > most_debug_options "MOST_F90_OPTS=-g -C -ftrap=common"
      echo > most_precision_options "MOST_PREC=-r8"
   elif [ $MOST_F90 = "ifort" ] ; then
      echo >> most_compiler "MOST_F90_OPTS=-O"
      echo > most_debug_options "MOST_F90_OPTS=-g -C -fpe0 -traceback"
      echo > most_precision_options "MOST_PREC=-r8"
   elif [ $MOST_F90 = "nagfor" ] ; then
      echo >> most_compiler "MOST_F90_OPTS=-O -kind=byte"
      echo > most_debug_options "MOST_F90_OPTS=-g -C -fpe0 -traceback -kind=byte"
   elif [ $MOST_F90 = "gfortran" ] ; then
      # get major and minor version number
      GFMAJ=`gfortran -dumpversion | head -1 | sed -e 's/.*)//' | awk 'BEGIN {FS="."}{print $1}'`
      GFMIN=`gfortran -dumpversion | head -1 | sed -e 's/.*)//' | awk 'BEGIN {FS="."}{print $2}'`
      GFVER=`expr 100 '*' $GFMAJ + $GFMIN`
      echo "gfortran version " $GFMAJ.$GFMIN
      # flags for gfortran version >= 4.5 [ -fcheck=all -finit-real=snan ]
      if [ "$GFVER" -ge "405" ] ; then
         echo >> most_compiler "MOST_F90_OPTS=-O3"
         echo > most_debug_options "MOST_F90_OPTS=-g -O0 -ffpe-trap=invalid,zero,overflow -fcheck=all -finit-real=snan"
         echo > most_precision_options "MOST_PREC=-fdefault-real-8"
      # flags for gfortran version 4.2, 4.3 and 4.4
      elif [ "$GFVER" -ge "402" ] ; then
         echo >> most_compiler "MOST_F90_OPTS=-O3"
         echo > most_debug_options "MOST_F90_OPTS=-g -ffpe-trap=invalid,zero,overflow -fbounds-check"
      # flags for gfortran version <= 4.1 [ -frecord-marker=4 ]
      else
         echo >> most_compiler "MOST_F90_OPTS=-O3 -frecord-marker=4"
         echo > most_debug_options "MOST_F90_OPTS=-g -ffpe-trap=invalid,zero,overflow -fbounds-check -frecord-marker=4"
      fi
   else
      echo >> most_compiler "MOST_F90_OPTS=-O"
      echo > most_debug_options "MOST_F90_OPTS=-g"
   fi
   # disable assembler routines in debug mode
   echo >> most_debug_options "LEGMOD=legsym"
   echo >> most_debug_options "LEGFAST="
else
   echo "****************************************************"
   echo "* Sorry - didn't find any FORTRAN-90 compiler      *"
   echo "* Please install one (e.g. gfortran or g95)        *"
   echo "* or update list of known compilers in this script *"
   echo "****************************************************"
   exit 1
fi


#check for C compiler
# put your favourite compiler in front if you have more than one

for MOST_CC in cc gxlc gcc suncc cc NO_CC
do
   `hash 2>/dev/null $MOST_CC`
   if [ $? = 0 ] ; then break ; fi
done

if [ $MOST_CC != "NO_CC" ] ; then
   CC_PATH=`which $MOST_CC`
   echo >> most_compiler MOST_CC=$MOST_CC " # " $CC_PATH
   echo >> most_compiler "MOST_CC_OPTS=-O3"
   echo "Found C compiler              at: $CC_PATH"
   if [ $MOST_CC = "suncc" ] ; then
      echo >> most_debug_options "MOST_CC_OPTS=-g -C -ftrap=common"
   else
      echo >> most_debug_options "MOST_CC_OPTS=-g"
   fi
else
   echo "****************************************************"
   echo "* Sorry - didn't find any C compiler               *"
   echo "* Please install one (e.g. gcc or g++              *"
   echo "* or update list of known compilers in this script *"
   echo "****************************************************"
   exit 2
fi

#check for MPI-FORTRAN-90 compiler
for MPI_F90 in openmpif90 mpif90 mpf90 mpxlf90 NO_F90
do
   `hash 2>/dev/null $MPI_F90`
   if [ $? = 0 ] ; then break ; fi
done

if [ $MPI_F90 != "NO_F90" ] ; then
   F90_PATH=`which $MPI_F90`
   if [ $MPI_F90 = "openmpif90" ] ; then
      echo > most_compiler_mpi "MPI_RUN=openmpiexec"
      echo >> most_compiler_mpi "MOST_F90_OPTS=-O3"
   else
      echo > most_compiler_mpi "MPI_RUN=mpiexec"
      echo >> most_compiler_mpi "MOST_F90_OPTS=-O3"
   fi
   echo >> most_compiler_mpi MOST_F90=$MPI_F90
   echo >> most_compiler_mpi "MPIMOD=mpimod"
   echo "Found MPI FORTRAN-90 compiler at: $F90_PATH"
else
   echo "********************************************"
   echo "* No Message Passing Interface (MPI) found *"
   echo "* Models run on single CPU                 *"
   echo "********************************************"
fi


#check for MPI-C compiler
if [ $MPI_F90 != "NO_F90" ] ; then
   for MPI_CC in openmpicc gxlc mpicc mpcc g++ gcc NO_CC
   do
      `hash 2>/dev/null $MPI_CC`
      if [ $? = 0 ] ; then break ; fi
   done
   if [ $MPI_CC != "NO_CC" ] ; then
      CC_PATH=`which $MPI_CC`
      echo >> most_compiler_mpi MOST_CC=$MPI_CC
      echo >> most_compiler_mpi "MOST_CC_OPTS=-O3"
      echo "Found MPI C compiler          at: $CC_PATH"
   else
      echo "********************************************"
      echo "* No Message Passing Interface (MPI) found *"
      echo "* Models run on single CPU                 *"
      echo "********************************************"
   fi
fi

#check for 86_64 operating system
MACHINE=`uname -m`

#check for GNU assembler
for MOST_AS in as NO_AS
do
   `hash 2>/dev/null $MOST_AS`
   if [ $? = 0 ] ; then break ; fi
done

if [ $MOST_AS != "NO_AS" ] ; then
   AS_PATH=`which $MOST_AS`
   echo >> most_compiler     MOST_AS=$MOST_AS
   echo >> most_compiler_mpi MOST_AS=$MOST_AS
   echo "Found GNU assembler           at: $AS_PATH"
fi

# activate fast assembler Legendre module for OS X and Linux 64 bit
if [ $MACHINE = "x86_64" -a $MOST_AS = "as" ] ; then
   # use fast vectorized assembler routines for Legendre transformation
   echo >> most_compiler_mpi LEGMOD=legini
   echo >> most_compiler_mpi LEGFAST=legfast32.o
   echo >> most_compiler     LEGMOD=legini
   echo >> most_compiler     LEGFAST=legfast32.o
   UNAMES=`uname -s`
   # generate correct symbol names for accessing gfortran module variables
   if [ $UNAMES = "Darwin" ] ; then
      echo "Fast Legendre Transformation    : ACTIVE (Darwin)"
   else
      echo "Fast Legendre Transformation    : ACTIVE (Linux)"
   fi
else
   # use FORTRAN routines for Legendre transformation
   echo >> most_compiler_mpi LEGMOD=legsym
   echo >> most_compiler_mpi LEGFAST=
   echo >> most_compiler     LEGMOD=legsym
   echo >> most_compiler     LEGFAST=
fi

#check for Xlib
  if [ -e "/usr/X11/lib64/libX11.so" ] ; then
   XLIB_PATH="/usr/X11/lib64"
elif [ -e "/usr/lib64/libX11.so" ] ; then
   XLIB_PATH="/usr/lib64"
elif [ -e "/usr/X11/lib" ] ; then
   XLIB_PATH="/usr/X11/lib"
elif [ -e "/usr/lib/X11" ] ; then
   XLIB_PATH="/usr/lib/X11"
elif [ -e "/opt/X11" ] ; then
   XLIB_PATH="/opt/X11"
fi
if [ -n ${XLIB_PATH} ] ; then
   echo "Found Xlib (X11)              at: $XLIB_PATH"
   GUILIB="-L$XLIB_PATH -lX11"
else
   echo "**********************************************"
   echo "* Didn't find Xlib (X11) at standard paths   *"
   echo "* Hopefully the compiler will find it itself *"
   echo "**********************************************"
   GUILIB=" -lX11"
fi

echo >> most_compiler     "GUILIB=$GUILIB"
echo >> most_compiler_mpi "GUILIB=$GUILIB"
echo >> most_compiler     "GUIMOD=guimod"
echo >> most_compiler_mpi "GUIMOD=guimod"
echo >> most_compiler     "PUMAX=pumax"
echo >> most_compiler_mpi "PUMAX=pumax"
echo  > makefile "most.x: most.c"
echo >> makefile "	$MOST_CC -o most.x most.c -lm $GUILIB"
echo >> makefile "clean:"
echo >> makefile "	rm -f *.o *.x F90* most_*"

#create directories

mkdir -p puma/bld
mkdir -p puma/bin
mkdir -p puma/run
mkdir -p sam/bld
mkdir -p sam/bin
mkdir -p sam/run
mkdir -p plasim/bld
mkdir -p plasim/bin
mkdir -p plasim/run
mkdir -p plasim/dat/T21_exo
mkdir -p plasim/dat/T42_exo

#check FORTRAN I/O

export MOST_F90
export MOST_CC
HOSTNAME=`hostname`
export HOSTNAME
MOSTARCH=`uname -a`
export MOSTARCH
make -e -f makecheck
./f90check.x
./cc_check.x
echo >> most_info.txt "FORTRAN Compiler: $MOST_F90"
echo >> most_info.txt "C       Compiler: $MOST_CC"
cat most_info.txt
make

echo
echo "configuration complete - run <most.x>"
echo
#
#end of configure.sh
