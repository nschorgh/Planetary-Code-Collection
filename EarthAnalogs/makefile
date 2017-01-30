CC = gfortran
CFLAGS = -O3 -Wall
LFLAGS = -O3 -Wall

#CC = ifort
#CFLAGS = -O3 -warn all  -warn notruncated_source
#LFLAGS = -O3 -warn all

SR = ../Common/

insol_earth: topos.o crater_modules.o sunpos.o insol_earth.o \
	mk_atmosphere.o crater_common.o model_subs.o flux_wshad.o makefile
	$(CC) $(LFLAGS) topos.o insol_earth.o mk_atmosphere.o \
	sunpos.o crater_common.o model_subs.o flux_wshad.o 

insol_flat: sunpos.o insol_flat.o mk_atmosphere.o flux_wshad.o makefile
	$(CC) $(LFLAGS) sunpos.o insol_flat.o mk_atmosphere.o flux_wshad.o

cratersQ_earth: topos.o crater_modules.o sunpos.o cratersQ_earth.o sunpos.o \
	model_subs.o crater_common.o mk_atmosphere.o flux_wshad.o makefile
	make -C $(SR) common
	$(CC) $(LFLAGS) cratersQ_earth.o $(SR)generalorbit.o topos.o \
	model_subs.o crater_common.o sunpos.o mk_atmosphere.o \
	$(SR)grids.o $(SR)conductionQ.o $(SR)tridag.o flux_wshad.o 

cratersQ_earth_static: topos.o crater_modules.o sunpos.o crater_common.o \
	cratersQ_earth_static.o mk_atmosphere.o model_subs.o flux_wshad.o \
	makefile
	make -C $(SR) common
	$(CC) $(LFLAGS) cratersQ_earth_static.o topos.o model_subs.o \
	crater_common.o sunpos.o mk_atmosphere.o flux_wshad.o $(SR)grids.o \
	$(SR)conductionQ.o $(SR)tridag.o 


%.o:	%.f makefile
	$(CC) $(CFLAGS) -c $<

%.o:	%.for makefile
	$(CC) $(CFLAGS) -c $<

%.o:	%.f90 makefile
	$(CC) $(CFLAGS) -c $< 


.PHONY: clean 

clean:
	rm -f  *.o  a.out  *.o  *.mod 



