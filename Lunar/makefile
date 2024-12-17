FC = gfortran
FFLAGS = -O2 -Wall
LFLAGS = -O2 -Wall 

SR = ../Common/

oscidea1: oscidea_subs.o oscidea1.o subl_h2o.o makefile
	$(FC) $(LFLAGS) oscidea1.o subl_h2o.o oscidea_subs.o

subsdiff moon_subsdiff: moon_subsdiff.o subsdiff_subs.o bet.o makefile
	make -C $(SR) common
	$(FC) $(LFLAGS) moon_subsdiff.o bet.o $(SR)ran2.o $(SR)psv.o \
	subsdiff_subs.o

equilbr moon_subsdiff_equilbr: moon_subsdiff_equilbr.o bet.o subl_h2o.o \
	makefile
	make -C $(SR) common
	$(FC) $(LFLAGS) moon_subsdiff_equilbr.o bet.o $(SR)psv.o subl_h2o.o


ensemble: ensemble.o flux_moon.o soilthprop_moon.o makefile
	make -C $(SR) common
	$(FC) $(LFLAGS) ensemble.o $(SR)conductionQ.o $(SR)tridag.o \
	$(SR)ran2.o $(SR)grids.o flux_moon.o soilthprop_moon.o


common-lunar: subl_h2o.o subl_species.o soilthprop_moon.o flux_moon.o bet.o \
	subsdiff_subs.o makefile


%.o:	%.f90 makefile
	$(FC) $(FFLAGS) -c $<


clean:
	rm -f  *.o  a.out *.mod
