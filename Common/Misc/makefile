FC = gfortran
FFLAGS = -O2 -Wall
LFLAGS = -O2 -Wall

#FC = ifort
#FFLAGS = -O2 #-warn all -warn notruncated_source -nogen-interface
#LFLAGS = -O2 -warn all

#FC = pgfortran
#FFLAGS = -fast
#LFLAGS = -fast


all: conductionQ+iter.o conductionQ+smooth.o conductionQ+volt.o makefile


%.o:	%.f90 makefile
	$(FC) $(FFLAGS) -c $<



.PHONY: clean

clean:
	rm -f  *.o  *.mod 



