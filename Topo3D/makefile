CC = gfortran
CFLAGS = -O3 -Wall
LFLAGS = -O3 -Wall

#CC = ifort
#CFLAGS = -O3 -warn all  -warn notruncated_source
#LFLAGS = -O3 -warn all

SR = ../Common/

shadows: topos.o crater_modules.o shadows.o shadow_subs.o  \
	crater_common.o hpsort.o makefile
	make -C $(SR) common
	$(CC) $(LFLAGS) topos.o shadows.o shadow_subs.o hpsort.o \
	crater_common.o -o $@

fieldofviews: topos.o crater_modules.o fieldofviews.o crater_common.o \
	shadow_subs.o fieldofview_subs.o hpsort.o makefile
	make -C $(SR) common
	$(CC) $(LFLAGS) fieldofviews.o topos.o crater_common.o \
	hpsort.o fieldofview_subs.o shadow_subs.o -o $@

cratersQ_snapshot: topos.o crater_modules.o cratersQ_snapshot.o \
	model_subs.o crater_common.o makefile
	make -C $(SR) common
	$(CC) $(LFLAGS) cratersQ_snapshot.o topos.o model_subs.o \
	crater_common.o $(SR)grids.o $(SR)conductionQ.o $(SR)tridag.o -o $@

cratersQ_moon: topos.o crater_modules.o cratersQ_moon.o model_subs.o \
	crater_common.o makefile
	make -C $(SR) common
	$(CC) $(LFLAGS) cratersQ_moon.o crater_common.o $(SR)generalorbit.o \
	topos.o model_subs.o $(SR)grids.o $(SR)conductionQ.o \
	$(SR)tridag.o -o $@

%.o:	%.f makefile
	$(CC) $(CFLAGS) -c $<

%.o:	%.for makefile
	$(CC) $(CFLAGS) -c $<

%.o:	%.f90 makefile
	$(CC) $(CFLAGS) -c $< 


.PHONY: clean 

clean:
	rm -f  *.o  a.out  *.mod \
	shadows fieldofviews cratersQ_snapshot cratersQ_moon


