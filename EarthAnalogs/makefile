# This overwrites and passes on options:
# make FC=pgfortran FFLAGS=-fast LFLAGS=-fast
# make FC=ifort FFLAGS='-O2 -warn all' LFLAGS='-O2 -warn all'

FC = gfortran
FFLAGS = -O2 -Wall
LFLAGS = -O2 -Wall


insol3d_earth: filemanager.o topo3d_modules.o sunpos.o topo3d_subs.o \
	insol3d_earth.o mk_atmosphere.o topo3d_common.o flux_noatm.o \
	makefile
	$(FC) $(LFLAGS) insol3d_earth.o mk_atmosphere.o sunpos.o \
	topo3d_common.o topo3d_subs.o flux_noatm.o 

insol_flat: sunpos.o insol_flat.o mk_atmosphere.o flux_noatm.o makefile
	$(FC) $(LFLAGS) sunpos.o insol_flat.o mk_atmosphere.o flux_noatm.o


%.o:	%.for makefile
	$(FC) $(FFLAGS) -c $<

%.o:	%.f90 makefile
	$(FC) $(FFLAGS) -c $<


.PHONY: clean 

clean:
	rm -f  *.o  a.out  *.mod 



