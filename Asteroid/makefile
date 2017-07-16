CC = gfortran
CFLAGS = -O3 -Wall
LFLAGS = -O3 -Wall

#CC = ifort
#CFLAGS = -O3 -warn all  -warn notruncated_source
#LFLAGS = -O3 -warn all

SR = ../Common/

asteroid_thermal: asteroid_thermal.f90 asteroid_thermal.o \
	oneasteroid_thermal1d.o insolonly.o makefile 
	make -C $(SR) common
	$(CC) $(LFLAGS) asteroid_thermal.f90 oneasteroid_thermal1d.o \
	$(SR)flux_noatm.o $(SR)conductionQ.o $(SR)tridag.o \
	$(SR)grids.o $(SR)generalorbit.o insolonly.o -o $@

# Modules: body.mod miscparameters.mod allinterfaces.mod
# The following depend on modules:
#       asteroid_fast2, fast_subs_asteroid2

asteroid_fast2: fast_modules_asteroid2.f90 fast_modules_asteroid2.o \
	asteroid_fast2.o fast_subs_asteroid2.o impactstirring.o makefile 
	make -C $(SR) common
	$(CC) $(LFLAGS) asteroid_fast2.f90 fast_subs_asteroid2.f90 \
	$(SR)flux_noatm.o $(SR)conductionQ.o $(SR)tridag.o \
	$(SR)grids.o $(SR)generalorbit.o $(SR)psv.o \
	$(SR)derivs.o impactstirring.o $(SR)ran2.o -o fast

asteroid_fast1: fast_modules_asteroid1.f90 fast_modules_asteroid1.o \
	asteroid_fast1.o fast_subs_asteroid1.o makefile 
	make -C $(SR) common
	$(CC) $(LFLAGS) asteroid_fast1.o fast_subs_asteroid1.o \
	$(SR)flux_noatm.o $(SR)conductionQ.o $(SR)tridag.o \
	$(SR)grids.o $(SR)generalorbit.o $(SR)psv.o $(SR)subl_subs.o -o fast

test1: testgardening1.o impactstirring.o makefile
	make -C $(SR) common
	$(CC) $(LFLAGS) impactstirring.o testgardening1.o \
	$(SR)ran2.o $(SR)grids.o 

%.o:	%.f makefile
	$(CC) $(CFLAGS) -c $<

%.o:	%.for makefile
	$(CC) $(CFLAGS) -c $<

%.o:	%.f90 makefile
	$(CC) $(CFLAGS) -c $< 


.PHONY: clean 

clean:
	rm -f  *.o a.out  *.mod asteroid_thermal fast



