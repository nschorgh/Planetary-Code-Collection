CC = gfortran
CFLAGS = -O3 -Wall  # -fopenmp
LFLAGS = -O3 -Wall  # -fopenmp

SR = ../Common/

moon_exo: geogrid.o body.o moon_exo.o montecarlo.o makefile 
	make -C $(SR) common more
	$(CC) $(LFLAGS) moon_exo.o montecarlo.o geogrid.o \
	$(SR)flux_noatm.o $(SR)conductionQ.o $(SR)tridag.o \
	$(SR)grids.o $(SR)psv.o $(SR)subl_subs.o \
	$(SR)gasdev.o $(SR)ran2.o

ceres_exo: geogrid.o body.o ceres_exo.o montecarlo.o makefile 
	make -C $(SR) common more
	$(CC) $(LFLAGS) ceres_exo.o montecarlo.o geogrid.o \
	$(SR)flux_noatm.o $(SR)conductionQ.o $(SR)tridag.o \
	$(SR)grids.o $(SR)psv.o $(SR)subl_subs.o \
	$(SR)gasdev.o $(SR)ran2.o $(SR)generalorbit.o

# module in geogrid.f90 

%.o:	%.f makefile
	$(CC) $(CFLAGS) -c $<

%.o:	%.for makefile
	$(CC) $(CFLAGS) -c $<

%.o:	%.f90 body.f90 makefile
	$(CC) $(CFLAGS) -c $< 


.PHONY: clean 

clean:
	rm -f  *.o  a.out  *.mod



