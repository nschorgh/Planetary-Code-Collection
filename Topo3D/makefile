CC = gfortran
CFLAGS = -O2 -Wall -Wno-integer-division
LFLAGS = -O2 -Wall

#CC = ifort
#CFLAGS = -O2 -warn all  -warn notruncated_source
#LFLAGS = -O2 -warn all

SR = ../Common/
SRM = ../Mars/
# topos.f90 has to be compiled first for module filemanager
# topo3d_modules.f90 has to be compiled second for module allinterfaces
# model_subs.f90 contains module newhorizons and should be compiled third
# conductionT2/Q2 contain modules conductionT/Q (therefore the symlink)

shadows: topos.o topo3d_modules.o multigrid.o shadows.o shadow_subs.o \
	topo3d_common.o hpsort.o multigrid.o makefile
	make -C $(SR) common
	$(CC) $(LFLAGS) shadows.o shadow_subs.o hpsort.o \
	topo3d_common.o multigrid.o -o $@

fieldofviews: topos.o topo3d_modules.o fieldofviews.o topo3d_common.o \
	shadow_subs.o fieldofview_subs.o hpsort.o makefile
	make -C $(SR) common
	$(CC) $(LFLAGS) fieldofviews.o topo3d_common.o hpsort.o \
	fieldofview_subs.o shadow_subs.o -o $@

cratersQ_snapshot: topos.o topo3d_modules.o model_subs.o topo3d_common.o \
	cratersQ_snapshot.o makefile
	$(CC) $(LFLAGS) cratersQ_snapshot.o model_subs.o \
	topo3d_common.o $(SR)flux_noatm.o -o $@

cratersQ_moon: topos.o topo3d_modules.o model_subs.o cratersQ_moon.o \
	topo3d_common.o makefile
	make -C $(SR) common
	$(CC) $(LFLAGS) cratersQ_moon.o topo3d_common.o $(SR)generalorbit.o \
	model_subs.o $(SR)grids.o $(SR)conductionQ.o $(SR)flux_noatm.o \
	$(SR)tridag.o -o $@

cratersQ_mars: topos.o topo3d_modules.o model_subs.o conductionT2.o \
	conductionQ2.o cratersQ_mars.o topo3d_common.o flux_mars.o \
	makefile
	make -C $(SR) common
	make -C $(SRM) common-mars
	$(CC) $(LFLAGS) cratersQ_mars.o topo3d_common.o model_subs.o \
	flux_mars.o $(SR)grids.o conductionQ2.o $(SR)tridag.o \
	conductionT2.o $(SR)julday.o $(SRM)marsorbit.o $(SR)psv.o -o $@

cratersQ_mars_shadowsonly: topos.o topo3d_modules.o model_subs.o \
	topo3d_common.o cratersQ_mars_shadowsonly.o makefile
	make -C $(SR) common
	make -C $(SRM) common-mars
	$(CC) $(LFLAGS) cratersQ_mars_shadowsonly.o topo3d_common.o \
	model_subs.o $(SR)flux_noatm.o $(SR)julday.o $(SRM)marsorbit.o

test: topos.o topo3d_modules.o test_shadows1pt.o shadow_subs.o \
	topo3d_common.o hpsort.o multigrid.o makefile
	make -C $(SR) common
	$(CC) $(LFLAGS) test_shadows1pt.o shadow_subs.o hpsort.o \
	topo3d_common.o multigrid.o 

fieldproperties: topos.o topo3d_modules.o topo3d_common.o model_subs.o \
	fieldproperties.o makefile
	$(CC) $(LFLAGS) fieldproperties.o topo3d_common.o model_subs.o 

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



