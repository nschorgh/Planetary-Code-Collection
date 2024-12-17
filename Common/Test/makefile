FC = gfortran
FFLAGS = -O2 -Wall
LFLAGS = -O2 -Wall

#FC = ifort
#FFLAGS = -O2 -warn all  -warn notruncated_source  -nogen-interface
#LFLAGS = -O2 -warn all

SR = ../



testcrankQ: testcrankQ.o makefile
	make FC='$(FC)' FFLAGS='$(FFLAGS)' -C $(SR) common
	$(FC) $(LFLAGS) testcrankQ.o $(SR)flux_noatm.o $(SR)conductionQ.o \
	$(SR)tridag.o $(SR)grids.o

testcrankT: testcrankT.o makefile
	make FC='$(FC)' FFLAGS='$(FFLAGS)' -C $(SR) common
	$(FC) $(LFLAGS) testcrankT.o $(SR)conductionT.o $(SR)tridag.o \
	$(SR)grids.o

testcrankQ_asymp asymp: testcrankQ_asymp.o makefile
	make FC='$(FC)' FFLAGS='$(FFLAGS)' -C $(SR) common
# 	use legacy version of cranknQ
#	which happens to be in Misc/conductionQ+smooth.f90
	$(FC) $(FFLAGS) -c ../Misc/conductionQ+smooth.f90
	$(FC) $(LFLAGS) testcrankQ_asymp.o conductionQ+smooth.o $(SR)tridag.o \
	$(SR)grids.o

testcrankQ_conv conv: testcrankQ_conv.o makefile
	make -C $(SR) common
# 	use legacy version of cranknQ
	$(FC) $(FFLAGS) -c ../Misc/conductionQ+smooth.f90
	$(FC) -O1 testcrankQ_conv.o $(SR)flux_noatm.o conductionQ+smooth.o \
	$(SR)tridag.o $(SR)grids.o


%.o:	%.f90 makefile
	$(FC) $(FFLAGS) -c $<


.PHONY: clean cleanrun

clean:
	rm -f  *.o  a.out  *.mod
	rm -f Tprofile z Tsurface



