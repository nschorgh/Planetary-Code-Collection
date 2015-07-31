CC = gfortran
CFLAGS = -O3 -Wall
LFLAGS = -O3 -Wall

#CC = ifort
#C90 = ifort
#CFLAGS = -O3 -warn all  -warn notruncated_source
#C90FLAGS = -O3 -warn all
#LFLAGS = -O3 -warn all


common: generalorbit.o flux_noatm.o conductionQ.o conductionT.o \
	tridag.o grids.o psv.o derivs.o psvco2.o subl_subs.o \
	vapordiffusioni.o adsorption.o julday.o gasdev.o ran2.o \
	hpsort.o makefile

julday.o: julday.for makefile   # NR routine without pause statement
	$(CC) $(CFLAGS) -c julday.for

%.o:	%.f makefile
	$(CC) $(CFLAGS) -c $<

%.o:	%.for makefile
	$(CC) $(CFLAGS) -c $<

%.o:	%.f90 makefile
	$(CC) $(CFLAGS) -c $< 



.PHONY: clean 

clean:
	rm -f  *.o 



