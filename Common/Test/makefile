CC = gfortran
CFLAGS = -O2 -Wall
LFLAGS = -O2 -Wall

#CC = ifort
#CFLAGS = -O2 -warn all  -warn notruncated_source  -nogen-interface
#LFLAGS = -O2 -warn all

SR = ../



testcrankQ: testcrankQ.o makefile
	make CC='$(CC)' CFLAGS='$(CFLAGS)' -C $(SR) common
	$(CC) -O1 testcrankQ.o $(SR)flux_noatm.o $(SR)conductionQ.o \
	$(SR)tridag.o $(SR)grids.o

testcrankT: testcrankT.o makefile
	make CC='$(CC)' CFLAGS='$(CFLAGS)' -C $(SR) common
	$(CC) -O1 testcrankT.o $(SR)conductionT.o $(SR)tridag.o $(SR)grids.o



%.o:	%.f makefile
	$(CC) $(CFLAGS) -c $<

%.o:	%.f90 makefile
	$(CC) $(CFLAGS) -c $<  


.PHONY: clean cleanrun

clean:
	rm -f  *.o  a.out  *.mod
	rm -f Tprofile z Tsurface



