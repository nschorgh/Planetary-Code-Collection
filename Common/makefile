CC = gfortran
CFLAGS = -O3 -Wall
LFLAGS = -O3 -Wall

#CC = ifort
#CFLAGS = -O3 -warn all  -warn notruncated_source
#LFLAGS = -O3 -warn all

#CC = pgfortran
#CFLAGS = -fast
#LFLAGS = -fast


common: generalorbit.o flux_noatm.o conductionQ.o conductionT.o \
	tridag.o grids.o psv3.o derivs.o psvco2.o subl_subs.o \
	julday.o gasdev.o ran2.o makefile

%.o:	%.f makefile
	$(CC) $(CFLAGS) -c $<

%.o:	%.for makefile
	$(CC) $(CFLAGS) -c $<

%.o:	%.f90 makefile
	$(CC) $(CFLAGS) -c $< 



.PHONY: clean 

clean:
	rm -f  *.o 



