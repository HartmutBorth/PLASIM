NGX=64

F90=gfortran
OPTS=-O2 -ffpe-trap=invalid,zero,overflow  -finit-real=snan
DEBOPTS=-g -ffpe-trap=invalid,zero,overflow -ffpe-summary=none -fcheck=all -finit-real=snan

%.o : %.c
	$(CC) -O3 -c $<

%.o : %.f90
	$(F90) -c $(OPTS) $<

cat.x:	cat.o fftfor.o fftasm.o rstmod.o simmod.o usermod.o guimod.o pumax.o
	$(F90) -o cat.x $(OPTS) cat.o fftfor.o fftasm.o rstmod.o simmod.o usermod.o guimod.o pumax.o -L/opt/X11 -lX11
cat.o:	cat.f90
fftfor.o:	fftfor.f90

fftasm.o:	fftasm.s
	as -o fftasm.o fftasm.s

rstmod.o:	rstmod.f90

simmod.o:	simmod.f90

usermod.o:	usermod.f90

guimod.o:	guimod.f90

pumax.o:	pumax.c

clean:
	rm -f cat.x *.o *.mod

cleanall:
	rm -f cat.x *.o *.mod cat_gp cat_gp.nc cat_sp cat_sp.nc cat_tseri cat_cfl cat_diag cat_rstini cat_rstfin
	make cleanrsttest

run:
	time ./cat.x $(NGX)

nc:	
	cdo -f nc copy cat_gp cat_gp.nc

view:
	ncview cat_gp.nc
all:
	make cleanall
	make
	make run
	make nc
	make view
rsttest:
	cp ../../tools/cat_rst_test.csh .       
	csh cat_rst_test.csh
	rm -f cat_rst_test.csh
cleanrsttest:
	cp ../../tools/clean_cat_rst_test.csh .       
	csh clean_cat_rst_test.csh
	rm -f clean_cat_rst_test.csh
