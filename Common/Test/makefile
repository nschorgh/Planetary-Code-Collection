CC = gfortran
CFLAGS = -O2 -Wall
LFLAGS = -O2 -Wall

#CC = ifort
#CFLAGS = -O2 -warn all  -warn notruncated_source  -nogen-interface
#LFLAGS = -O2 -warn all

SR = ../



testcrankQ: testcrankQ.o makefile
	make CC='$(CC)' CFLAGS='$(CFLAGS)' -C $(SR) common
	$(CC) $(LFLAGS) testcrankQ.o $(SR)flux_noatm.o $(SR)conductionQ.o \
	$(SR)tridag.o $(SR)grids.o

testcrankT: testcrankT.o makefile
	make CC='$(CC)' CFLAGS='$(CFLAGS)' -C $(SR) common
	$(CC) $(LFLAGS) testcrankT.o $(SR)conductionT.o $(SR)tridag.o \
	$(SR)grids.o

testcrankQ_asymp asymp: testcrankQ_asymp.o makefile
# 	use pure version of cranknQ
	make CC='$(CC)' CFLAGS='$(CFLAGS)' -C $(SR) common
	$(CC) $(CFLAGS) -c ../Misc/conductionQ+smooth.f90
	$(CC) $(LFLAGS) testcrankQ_asymp.o conductionQ+smooth.o $(SR)tridag.o \
	$(SR)grids.o

testcrankQ_conv conv: testcrankQ_conv.o makefile
	make -C $(SR) common
	$(CC) -O1 testcrankQ_conv.o $(SR)flux_noatm.o $(SR)conductionQ.o \
	$(SR)tridag.o $(SR)grids.o


%.o:	%.f90 makefile
	$(CC) $(CFLAGS) -c $<  


.PHONY: clean cleanrun

clean:
	rm -f  *.o  a.out  *.mod
	rm -f Tprofile z Tsurface



