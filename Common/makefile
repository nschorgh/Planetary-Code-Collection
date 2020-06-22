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
	tridag.o grids.o derivs.o psv.o psv2.o psv3.o psvco2.o \
	julday.o ran2.o subl_h2o.o makefile 

%.o:	%.f makefile
	$(CC) $(CFLAGS) -c $<

%.o:	%.for makefile
	$(CC) $(CFLAGS) -c $<

%.o:	%.f90 makefile
	$(CC) $(CFLAGS) -c $< 



.PHONY: clean

clean:
	rm -f  *.o  *.mod 



