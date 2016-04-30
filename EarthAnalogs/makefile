CC = gfortran
CFLAGS = -O3 -Wall
LFLAGS = -O3 -Wall

#CC = ifort
#CFLAGS = -O3 -warn all  -warn notruncated_source
#LFLAGS = -O3 -warn all

SR = ../Common/

insol_earth: topos.o crater_modules.o sunpos.o insol_earth.o \
	mk_atmosphere.o crater_common.o model_subs.o makefile
	make -C $(SR) common
	$(CC) $(LFLAGS) topos.o insol_earth.o mk_atmosphere.o \
	sunpos.o crater_common.o model_subs.o

%.o:	%.f makefile
	$(CC) $(CFLAGS) -c $<

%.o:	%.for makefile
	$(CC) $(CFLAGS) -c $<

%.o:	%.f90 makefile
	$(CC) $(CFLAGS) -c $< 


.PHONY: clean 

clean:
	rm -f  *.o  a.out  *.o  *.mod 



