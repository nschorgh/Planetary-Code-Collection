CC = gfortran
CFLAGS = -O2 -Wall
LFLAGS = -O2 -Wall 

SR = ../Common/

oscidea1: oscidea_subs.o oscidea1.o subl_h2o.o makefile
	$(CC) $(LFLAGS) oscidea1.o subl_h2o.o oscidea_subs.o

subsdiff moon_subsdiff: moon_subsdiff.o subsdiff_subs.o bet.o makefile
	make -C $(SR) common
	$(CC) $(LFLAGS) moon_subsdiff.o bet.o $(SR)ran2.o $(SR)psv.o \
	subsdiff_subs.o

equilbr moon_subsdiff_equilbr: moon_subsdiff_equilbr.o bet.o subl_h2o.o \
	makefile
	make -C $(SR) common
	$(CC) $(LFLAGS) moon_subsdiff_equilbr.o bet.o $(SR)psv.o subl_h2o.o


common-lunar: subl_h2o.o subl_species.o soilthprop_moon.o bet.o \
	subsdiff_subs.o makefile


%.o:	%.f90 makefile
	$(CC) $(CFLAGS) -c $<  


clean:
	rm -f  *.o  a.out *.mod





