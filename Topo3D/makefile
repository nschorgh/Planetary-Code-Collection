CC = gfortran
CFLAGS = -O2 -Wall -Wno-integer-division
LFLAGS = -O2 -Wall

#CC = ifort
#CFLAGS = -O2 -warn all  -warn notruncated_source
#LFLAGS = -O2 -warn all

SR = ../Common/
SRM = ../Mars/
# filemanager.f90 has to be compiled first for module filemanager
# topo3d_modules.f90 has to be compiled second for module allinterfaces
# topo3d_subs.f90 contains module newhorizons and should be compiled third
# conductionT2/Q2 contain modules conductionT/Q (therefore the symlink)

shadows: filemanager.o topo3d_modules.o multigrid.o shadows.o shadow_subs.o \
	topo3d_common.o hpsort.o multigrid.o makefile
	$(CC) $(LFLAGS) shadows.o shadow_subs.o hpsort.o \
	topo3d_common.o multigrid.o -o $@

fieldofviews: filemanager.o topo3d_modules.o fieldofviews.o topo3d_common.o \
	shadow_subs.o fieldofview_subs.o hpsort.o makefile
	$(CC) $(LFLAGS) fieldofviews.o topo3d_common.o hpsort.o \
	fieldofview_subs.o shadow_subs.o -o $@

cratersQ_equilbr: filemanager.o topo3d_modules.o topo3d_subs.o \
	topo3d_common.o cratersQ_equilbr.o makefile
	make -C $(SR) common
	$(CC) $(LFLAGS) cratersQ_equilbr.o topo3d_subs.o topo3d_common.o \
	$(SR)flux_noatm.o -o $@

cratersQ_moon: filemanager.o topo3d_modules.o topo3d_subs.o cratersQ_moon.o \
	topo3d_common.o makefile
	make -C $(SR) common
	$(CC) $(LFLAGS) cratersQ_moon.o topo3d_common.o $(SR)generalorbit.o \
	topo3d_subs.o $(SR)grids.o $(SR)conductionQ.o $(SR)flux_noatm.o \
	$(SR)tridag.o -o $@

cratersQ_mars: filemanager.o topo3d_modules.o topo3d_subs.o conductionT2.o \
	conductionQ2.o cratersQ_mars.o topo3d_common.o makefile
	make -C $(SR) common
	make -C $(SRM) common-mars
	$(CC) $(LFLAGS) cratersQ_mars.o topo3d_common.o topo3d_subs.o \
	$(SRM)flux_mars.o $(SR)grids.o conductionQ2.o $(SR)tridag.o \
	conductionT2.o $(SR)julday.o $(SRM)marsorbit.o $(SR)psv.o -o $@

insol_mars: filemanager.o topo3d_modules.o topo3d_subs.o \
	topo3d_common.o insol_mars.o makefile
	make -C $(SR) common
	make -C $(SRM) common-mars
	$(CC) $(LFLAGS) insol_mars.o topo3d_common.o topo3d_subs.o \
	$(SR)flux_noatm.o $(SR)julday.o $(SRM)marsorbit.o

fieldproperties: filemanager.o topo3d_modules.o topo3d_common.o \
	topo3d_subs.o fieldproperties.o makefile
	$(CC) $(LFLAGS) fieldproperties.o topo3d_common.o topo3d_subs.o 

%.o:	%.f makefile
	$(CC) $(CFLAGS) -c $<

%.o:	%.for makefile
	$(CC) $(CFLAGS) -c $<

%.o:	%.f90 makefile
	$(CC) $(CFLAGS) -c $< 


.PHONY: clean 

clean:
	rm -f  *.o  a.out  *.mod \
	shadows fieldofviews cratersQ_equilbr cratersQ_moon cratersQ_mars



