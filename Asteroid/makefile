CC = gfortran
CFLAGS = -O3 -Wall
LFLAGS = -O3 -Wall

#CC = ifort
#CFLAGS = -O3 -warn all  -warn notruncated_source
#LFLAGS = -O3 -warn all

SR = ~/Thermal/ThermalLibF/Common/

asteroid_thermal: asteroid_thermal.f90 asteroid_thermal.o \
	oneasteroid_thermal1d.o insolonly.o makefile 
	make -C $(SR) common
	$(CC) $(LFLAGS) asteroid_thermal.f90 oneasteroid_thermal1d.o \
	$(SR)flux_noatm.o $(SR)conductionQ.o $(SR)tridag.o \
	$(SR)grids.o $(SR)generalorbit.o insolonly.o -o $@

%.o:	%.f makefile
	$(CC) $(CFLAGS) -c $<

%.o:	%.for makefile
	$(CC) $(CFLAGS) -c $<

%.o:	%.f90 makefile
	$(CC) $(CFLAGS) -c $< 


.PHONY: clean 

clean:
	rm -f  *.o  a.out  *.o90  *.mod  asteroid_thermal



