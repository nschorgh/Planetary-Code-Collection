FC = gfortran
FFLAGS = -O2 -Wall
LFLAGS = -O2 -Wall

#FC = ifort
#FFLAGS = -O2 -warn all -warn notruncated_source
#LFLAGS = -O2 -warn all

#FC = pgfortran
#FFLAGS = -fast
#LFLAGS = -fast


SR =  ../Common/


asteroid_thermal: asteroid_thermal.f90 asteroid_thermal.o \
	oneasteroid_thermal1d.o insolonly.o common_subs.o makefile

	make -C $(SR) common

	$(FC) $(LFLAGS) asteroid_thermal.f90 oneasteroid_thermal1d.o \
	$(SR)flux_noatm.o $(SR)conductionQ.o $(SR)tridag.o $(SR)grids.o \
	$(SR)generalorbit.o insolonly.o common_subs.o -o $@


asteroid_fast2: fast_modules_asteroid.f90 fast_modules_asteroid.o \
	fast_subs_asteroid2.o asteroid_fast2.o impactstirring.o \
	common_subs.o makefile 

	make -C $(SR) common

	$(FC) $(LFLAGS) asteroid_fast2.f90 fast_subs_asteroid2.f90 \
	$(SR)flux_noatm.o $(SR)conductionQ.o $(SR)tridag.o $(SR)grids.o \
	$(SR)generalorbit.o $(SR)psv.o $(SR)derivs.o impactstirring.o \
	$(SR)ran2.o common_subs.o -o fast


asteroid_fast1: fast_modules_asteroid.f90 fast_modules_asteroid.o \
	asteroid_fast1.o fast_subs_asteroid1.o common_subs.o makefile

	make -C $(SR) common

	$(FC) $(LFLAGS) asteroid_fast1.o fast_subs_asteroid1.o \
	$(SR)flux_noatm.o $(SR)conductionQ.o $(SR)tridag.o $(SR)grids.o \
	$(SR)generalorbit.o $(SR)psv.o common_subs.o  -o fast


test_gardening1: test_gardening1.o impactstirring.o makefile

	make -C $(SR) common

	$(FC) $(LFLAGS) impactstirring.o test_gardening1.o $(SR)ran2.o \
	$(SR)grids.o


sphere1d sphere1d_implicit: sphere1d_implicit.o common_subs.o makefile

	make -C $(SR) common

	$(FC) $(LFLAGS) sphere1d_implicit.o $(SR)flux_noatm.o $(SR)tridag.o \
	$(SR)generalorbit.o common_subs.o $(SR)psv.o

# Modules require certain build order
#    body.mod and constants.mod in asteroid_thermal.f90
#    body.mod, constants.mod, and allinterfaces.mod in fast_modules_asteroid.f90


%.o:	%.for makefile
	$(FC) $(FFLAGS) -c $<

%.o:	%.f90 makefile
	$(FC) $(FFLAGS) -c $<


.PHONY: clean 

clean:
	rm -f  *.o a.out  *.mod asteroid_thermal fast

