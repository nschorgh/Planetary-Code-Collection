FC = gfortran
FFLAGS = -O3 -Wall  # -fopenmp
LFLAGS = -O3 -Wall  # -fopenmp

SR = ../Common/

# geogrid.f90 contains module grid used by *_exo.f90
# body.f90 contains module body used by *_exo.f90 and montecarlo.f90


moon_exo: geogrid.o body.o moon_exo.o montecarlo.o gasdev.o subl_subs.o \
	makefile

	make -C $(SR) common

	$(FC) $(LFLAGS) moon_exo.o montecarlo.o geogrid.o $(SR)flux_noatm.o \
	$(SR)conductionQ.o $(SR)tridag.o $(SR)grids.o $(SR)psv.o subl_subs.o \
	gasdev.o $(SR)ran2.o

ceres_exo: geogrid.o body.o ceres_exo.o montecarlo.o gasdev.o subl_subs.o \
	makefile

	make -C $(SR) common

	$(FC) $(LFLAGS) ceres_exo.o montecarlo.o geogrid.o $(SR)flux_noatm.o \
	$(SR)conductionQ.o $(SR)tridag.o $(SR)grids.o $(SR)psv.o subl_subs.o \
	gasdev.o $(SR)ran2.o $(SR)generalorbit.o


%.o:	%.for makefile
	$(FC) $(FFLAGS) -c $<

%.o:	%.f90 body.f90 makefile
	$(FC) $(FFLAGS) -c $<


.PHONY: clean 

clean:
	rm -f  *.o  a.out  *.mod

