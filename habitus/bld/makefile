MOST_F90=gfortran
MOST_F90_OPTS=-O3 -ffpe-trap=invalid,zero,overflow -ffpe-summary=none -finit-real=snan
MPIMOD=mpimod_stub
MOST_CC=cc  #  /usr/bin/cc
MOST_CC_OPTS=-O3
MOST_AS=as
LEGMOD=legini
LEGFAST=legfast32.o
GUILIB=-L/opt/X11 -lX11
GUIMOD=guimod
PUMAX=pumax
# Makefile template for MoSt
# Compiler names and options are inserted by MoSt

%.o : %.cpp
	c++ -c $(MOST_CC_OPTS) $<

%.o : %.c
	$(MOST_CC) -c $(MOST_CC_OPTS) $<

%.o : %.f90
	$(MOST_F90) -c $(MOST_F90_OPTS) $(MOST_PREC) $<

%.o : %.s
	$(MOST_AS) -o $@ $<

# MPIMOD=mpimod     # single (mpimod_stub) or parallel (mpimod)
# GUIMOD=guimod     # with GUI or without (guimod_stub)
# LEGMOD=legsym     # Legendre module (FORTRAN or Assembler)
# PUMAX=pumax       # with GUI module or without (pumax_stub)
# GUILIB=${MOST_XLIB} -lX11
# MOST_F90=f90      # name of the FORTRAN compiler
# MOST_F90_OPTS=-O  # FORTRAN options
# MOST_PREC         # May be used to set double precision
# MOST_CC=cc        # name of the C compiler for GUI
# MOST_CC_OPTS=-O   # C options
# MOST_AS=as        # name of GNU assembler

# fftmod is faster than fft991mod, but fft991mod supports
# more resolutions, e.g. T63 and T106

# legfast is the assembler version of legsym and requires currently:
# 64-bit operating system, GNU assembler (as), AMD or Intel CPU
# and PUMA compiled for single precision (default).
# PUMA with legfast runs faster at the factor of 1.3 (T21), 1.7 (T42), 2.3 (T85)
# and factor 2.6 for all higher resolutions

ifdef LEGFAST
   ifdef MOST_PREC
      LEGFAST=legfast64.o
   endif
endif

FFTMOD=fftmod

OBJ=habitus.o

habitus.x:	$(OBJ)
	c++ -o habitus.x $(MOST_CC_OPTS) $(OBJ)

habitus.o:		habitus.cpp

clean:
	rm -f *.o *.mod
