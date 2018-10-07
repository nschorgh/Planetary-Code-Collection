CC = gfortran
CFLAGS = -O2 -Wall -Wno-integer-division
LFLAGS = -O2 -Wall

#CC = ifort
#CFLAGS = -O2 #-warn all -warn notruncated_source -nogen-interface
#LFLAGS = -O2 -warn all

#CC = pgfortran
#CFLAGS = -fast
#LFLAGS = -fast


common: generalorbit.o flux_noatm.o conductionQ.o conductionT.o \
	tridag.o grids.o psv.o derivs.o psvco2.o julday.o \
	gasdev.o ran2.o makefile

more:	subl_subs.o psv3.o conductionT2.o conductionQ2.o

%.o:	%.f makefile
	$(CC) $(CFLAGS) -c $<

%.o:	%.for makefile
	$(CC) $(CFLAGS) -c $<

%.o:	%.f90 makefile
	$(CC) $(CFLAGS) -c $< 



.PHONY: clean more

clean:
	rm -f  *.o  *.mod 



