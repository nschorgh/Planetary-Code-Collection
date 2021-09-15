CC = gfortran
CFLAGS = -O2 -Wall -Wno-integer-division
LFLAGS = -O2 -Wall

#CC = ifort
#CFLAGS = -O2 -warn all  -warn notruncated_source
#LFLAGS = -O2 -warn all

SR = ../Common/
SRM = ../Mars/

OBJFLS = filemanager.o topo3d_modules.o


shadows: $(OBJFLS) shadow_subs.o multigrid.o topo3d_common.o shadows.o \
	topo3d_geometry.o hpsort.o makefile

	$(CC) $(LFLAGS) shadows.o shadow_subs.o hpsort.o topo3d_common.o \
	multigrid.o topo3d_geometry.o -o $@

fieldofviews: $(OBJFLS) fieldofview_subs.o topo3d_common.o topo3d_geometry.o \
	fieldofviews.o hpsort.o makefile

	$(CC) $(LFLAGS) fieldofviews.o topo3d_common.o hpsort.o \
	fieldofview_subs.o topo3d_geometry.o -o $@

cratersQ_equilbr equilbr: $(OBJFLS) topo3d_subs.o topo3d_common.o \
	cratersQ_equilbr.o makefile
	make -C $(SR) common

	$(CC) $(LFLAGS) cratersQ_equilbr.o topo3d_subs.o topo3d_common.o \
	$(SR)flux_noatm.o 

cratersQ_moon: $(OBJFLS) topo3d_subs.o cratersQ_moon.o topo3d_common.o \
	makefile
	make -C $(SR) common
	$(CC) $(LFLAGS) cratersQ_moon.o topo3d_common.o $(SR)generalorbit.o \
	topo3d_subs.o $(SR)grids.o $(SR)conductionQ.o $(SR)flux_noatm.o \
	$(SR)tridag.o -o $@

insol3d_mars: $(OBJFLS) topo3d_subs.o topo3d_common.o insol3d_mars.o makefile

	make -C $(SR) common
	make -C $(SRM) common-mars

	$(CC) $(LFLAGS) insol3d_mars.o topo3d_common.o topo3d_subs.o \
	$(SR)flux_noatm.o $(SR)julday.o $(SRM)marsorbit.o

cratersQ_mars: $(OBJFLS) topo3d_subs.o conductionT2.o conductionQ2.o \
	cratersQ_mars.o topo3d_common.o topo3d_subs_mars.o makefile

	make -C $(SR) common
	make -C $(SRM) common-mars

	$(CC) $(LFLAGS) cratersQ_mars.o topo3d_common.o topo3d_subs.o \
	$(SRM)flux_mars.o $(SR)grids.o conductionQ2.o $(SR)tridag.o \
	conductionT2.o $(SR)julday.o $(SRM)marsorbit.o $(SR)psv.o \
	topo3d_subs_mars.o -o $@

cratersQ_mars_parallel: $(OBJFLS) topo3d_subs.o conductionT2.o conductionQ2.o \
	cratersQ_mars_parallel.o topo3d_common.o topo3d_subs_mars.o makefile

	make -C $(SR) common
	make -C $(SRM) common-mars

	$(CC) $(LFLAGS) cratersQ_mars_parallel.o topo3d_common.o topo3d_subs.o \
	$(SRM)flux_mars.o $(SR)grids.o conductionQ2.o $(SR)tridag.o \
	conductionT2.o $(SR)julday.o $(SRM)marsorbit.o $(SR)psv.o \
	topo3d_subs_mars.o 

cratersQ_mars_full: $(OBJFLS) topo3d_subs.o conductionT2.o conductionQ2.o \
	cratersQ_mars_full.o topo3d_common.o topo3d_subs_mars.o makefile

	make -C $(SR) common
	make -C $(SRM) common-mars

	$(CC) $(LFLAGS) cratersQ_mars_full.o topo3d_common.o topo3d_subs.o \
	$(SRM)flux_mars.o $(SR)grids.o conductionQ2.o $(SR)tridag.o \
	conductionT2.o $(SR)julday.o $(SRM)marsorbit.o $(SR)psv.o \
	topo3d_subs_mars.o 

fieldproperties: $(OBJFLS) topo3d_common.o topo3d_subs.o fieldproperties.o \
	makefile
	$(CC) $(LFLAGS) fieldproperties.o topo3d_common.o topo3d_subs.o 

test_shadows1pt test: $(OBJFLS) shadow_subs.o multigrid.o fieldofview_subs.o \
	topo3d_common.o topo3d_geometry.o hpsort.o test_shadows1pt.o makefile

	$(CC) $(LFLAGS) test_shadows1pt.o shadow_subs.o topo3d_common.o \
	multigrid.o topo3d_geometry.o hpsort.o fieldofview_subs.o

test_fov1pt test_fov: $(OBJFLS) fieldofview_subs.o test_fov1pt.o \
	topo3d_common.o topo3d_geometry.o hpsort.o makefile
	$(CC) $(LFLAGS) test_fov1pt.o topo3d_common.o topo3d_geometry.o \
	fieldofview_subs.o hpsort.o 

makegaussian: makegaussian.o ran1.o rlft3.o fourn.o nrutil.o makefile
	gcc -O2 makegaussian.o ran1.o rlft3.o fourn.o nrutil.o -lm	

# filemanager.f90 has to be compiled first for module filemanager
# topo3d_modules.f90 has to be compiled second for module allinterfaces
# topo3d_subs.f90 contains module newhorizons and should be compiled third
# conductionT2/Q2 contain modules conductionT/Q (therefore the symlink)
# multigrid.f90 has module multigrid needed by shadows.f90 and test_shadows1pt.f90
# shadow_subs.f90 has module azRays needed by several programs
# fieldofview_subs.f90 has module findvisibletopo needed by fieldofviews.f90 and test_shadows1pt.f90

%.o:	%.for makefile
	$(CC) $(CFLAGS) -c $<

%.o:	%.f90 filemanager.f90 makefile
	$(CC) $(CFLAGS) -c $< 



.PHONY: clean 

clean:
	rm -f  *.o  a.out  *.mod \
	shadows fieldofviews cratersQ_equilbr cratersQ_moon cratersQ_mars \
	horizons.dat viewfactors.dat 

