# This overwrites and passes on options:
# make FC=pgfortran FFLAGS=-fast LFLAGS=-fast
# make FC=ifort FFLAGS='-O3 -warn all' LFLAGS='-O3 -warn all'

FC = gfortran
FFLAGS = -O3 -Wall
LFLAGS = -O3 -Wall


SR = ../Common/



common-mars: marsorbit.o flux_mars.o


%.o:	%.f90 makefile
	$(FC) $(FFLAGS) -c $<


.PHONY: clean 

clean:
	rm -f  *.o 

