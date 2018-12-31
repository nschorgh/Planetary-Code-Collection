CC = gfortran
CFLAGS = -O3 -Wall -Wno-integer-division
LFLAGS = -O3 -Wall

#CC = ifort
#CFLAGS = -O3 -warn all  -warn notruncated_source
#LFLAGS = -O3 -warn all

SR = ../Common/
SRM = ../Mars/

shadows: topos.o crater_modules.o multigrid.o shadows.o shadow_subs.o \
	crater_common.o hpsort.o multigrid.o makefile
	make -C $(SR) common
	$(CC) $(LFLAGS) topos.o shadows.o shadow_subs.o hpsort.o \
	crater_common.o multigrid.o -o $@

fieldofviews: topos.o crater_modules.o fieldofviews.o crater_common.o \
	shadow_subs.o fieldofview_subs.o hpsort.o makefile
	make -C $(SR) common
	$(CC) $(LFLAGS) fieldofviews.o topos.o crater_common.o hpsort.o \
	fieldofview_subs.o shadow_subs.o -o $@

cratersQ_snapshot: topos.o crater_modules.o model_subs.o crater_common.o \
	cratersQ_snapshot.o makefile
	$(CC) $(LFLAGS) cratersQ_snapshot.o topos.o model_subs.o \
	crater_common.o $(SR)flux_noatm.o -o $@

cratersQ_moon: topos.o crater_modules.o model_subs.o cratersQ_moon.o \
	crater_common.o makefile
	make -C $(SR) common
	$(CC) $(LFLAGS) cratersQ_moon.o crater_common.o $(SR)generalorbit.o \
	topos.o model_subs.o $(SR)grids.o $(SR)conductionQ.o \
	$(SR)flux_noatm.o $(SR)tridag.o -o $@

cratersQ_mars: topos.o crater_modules.o model_subs.o conductionT2.o \
	conductionQ2.o cratersQ_mars.o crater_common.o flux_mars.o \
	makefile
	make -C $(SR) common more
	make -C $(SRM) common-mars
	$(CC) $(LFLAGS) cratersQ_mars.o crater_common.o topos.o \
	model_subs.o flux_mars.o $(SR)grids.o $(SR)conductionQ2.o \
	$(SR)tridag.o $(SR)conductionT2.o $(SR)julday.o \
	$(SRM)marsorbit.o $(SR)psv.o  -o $@

cratersQ_mars_shadowsonly: topos.o crater_modules.o model_subs.o \
	crater_common.o cratersQ_mars_shadowsonly.o makefile
	make -C $(SR) common
	make -C $(SRM) common-mars
	$(CC) $(LFLAGS) cratersQ_mars_shadowsonly.o crater_common.o \
	topos.o model_subs.o $(SR)flux_noatm.o $(SR)julday.o \
	$(SRM)marsorbit.o

test: topos.o crater_modules.o test_shadows1pt.o shadow_subs.o
	crater_common.o topos.o hpsort.o multigrid.o makefile
	make -C $(SR) common
	$(CC) $(LFLAGS) topos.o test_shadows1pt.o shadow_subs.o \
	hpsort.o crater_common.o multigrid.o 

fieldproperties: topos.o crater_modules.o crater_common.o model_subs.o \
	fieldproperties.o makefile
	$(CC) $(LFLAGS) fieldproperties.o topos.o crater_common.o \
	model_subs.o 

%.o:	%.f makefile
	$(CC) $(CFLAGS) -c $<

%.o:	%.for makefile
	$(CC) $(CFLAGS) -c $<

%.o:	%.f90 makefile
	$(CC) $(CFLAGS) -c $< 


.PHONY: clean 

clean:
	rm -f  *.o  a.out  *.mod \
	shadows fieldofviews cratersQ_snapshot cratersQ_moon cratersQ_mars



