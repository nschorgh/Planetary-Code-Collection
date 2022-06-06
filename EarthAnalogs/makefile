# This overwrites and passes on options:
# make CC=pgfortran CFLAGS=-fast LFLAGS=-fast
# make CC=ifort CFLAGS='-O2 -warn all' LFLAGS='-O2 -warn all'

CC = gfortran
CFLAGS = -O2 -Wall
LFLAGS = -O2 -Wall


insol3d_earth: filemanager.o topo3d_modules.o sunpos.o topo3d_subs.o \
	insol3d_earth.o mk_atmosphere.o topo3d_common.o flux_noatm.o \
	makefile
	$(CC) $(LFLAGS) insol3d_earth.o mk_atmosphere.o sunpos.o \
	topo3d_common.o topo3d_subs.o flux_noatm.o 

insol_flat: sunpos.o insol_flat.o mk_atmosphere.o flux_noatm.o makefile
	$(CC) $(LFLAGS) sunpos.o insol_flat.o mk_atmosphere.o flux_noatm.o


%.o:	%.for makefile
	$(CC) $(CFLAGS) -c $<

%.o:	%.f90 makefile
	$(CC) $(CFLAGS) -c $< 


.PHONY: clean 

clean:
	rm -f  *.o  a.out  *.mod 



