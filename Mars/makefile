# This overwrites and passes on options:
# make CC=pgfortran CFLAGS=-fast LFLAGS=-fast
# make CC=ifort CFLAGS='-O2 -warn all' LFLAGS='-O2 -warn all'

CC = gfortran
CFLAGS = -O2 -Wall
LFLAGS = -O2 -Wall

#CC = ifort
#CFLAGS = -O2 #-warn all -warn notruncated_source -nogen-interface
#LFLAGS = -O2 -warn all

#CC = pgfortran
#CFLAGS = -fast
#LFLAGS = -fast


SR = ../Common/

exper_thermal: exper_thermal1d.o vapordiffusioni.o soilthprop.o makefile
	make -C $(SR) common
	$(CC) $(LFLAGS) exper_thermal1d.o $(SR)conductionT.o \
	$(SR)tridag.o vapordiffusioni.o $(SR)psv.o $(SR)grids.o \
	soilthprop.o

mars_thermal: mars_thermal1d.o marsorbit.o flux.o soilthprop.o makefile
	make -C $(SR) common
	$(CC) $(LFLAGS) mars_thermal1d.o $(SR)julday.o marsorbit.o \
	flux.o $(SR)conductionQ.o $(SR)tridag.o $(SR)conductionT.o \
	$(SR)grids.o soilthprop.o $(SR)psv.o 

mars_mapt2: mars_mapt2.o jsub.o marsorbit.o flux.o soilthprop.o makefile
	make -C $(SR) common
	$(CC) $(LFLAGS) mars_mapt2.o jsub.o $(SR)julday.o marsorbit.o \
	flux.o $(SR)conductionQ.o $(SR)tridag.o $(SR)conductionT.o \
	$(SR)grids.o soilthprop.o $(SR)psv.o $(SR)psvco2.o

mars_mapi: mars_mapi.o marsorbit.o jsub.o flux.o soilthprop.o makefile
	make -C $(SR) common
	$(CC) $(LFLAGS) mars_mapi.o $(SR)julday.o marsorbit.o \
	flux.o $(SR)conductionQ.o $(SR)tridag.o $(SR)conductionT.o \
	$(SR)grids.o soilthprop.o $(SR)psv.o jsub.o $(SR)psvco2.o 

mars_mapiqh2v: mars_mapiqh2v.o jsubv.o marsorbit.o flux.o soilthprop.o \
	makefile
	make -C $(SR) common
	$(CC) $(LFLAGS) mars_mapiqh2v.o jsubv.o $(SR)julday.o \
	marsorbit.o flux.o $(SR)conductionQ.o $(SR)tridag.o \
	$(SR)conductionT.o $(SR)grids.o $(SR)psv.o soilthprop.o

stabgrow_fast: fast_modules.o fast_subs_exper.o stabgrow_fast.o \
	fast_subs_univ.o soilthprop.o makefile        
	make -C $(SR) common
	$(CC) $(LFLAGS) stabgrow_fast.o fast_subs_univ.o \
	fast_subs_exper.o $(SR)tridag.o $(SR)conductionT.o \
	$(SR)grids.o soilthprop.o $(SR)psv.o $(SR)derivs.o -o $@

exper_fast: fast_modules.o exper_fast.o fast_subs_univ.o \
	fast_subs_exper.o soilthprop.o makefile
	make -C $(SR) common
	$(CC) $(LFLAGS) exper_fast.o $(SR)tridag.o $(SR)conductionT.o \
	$(SR)grids.o soilthprop.o $(SR)psv.o $(SR)derivs.o \
	fast_subs_univ.o fast_subs_exper.o -o $@

mars_fast: fast_modules.o mars_fast.o fast_subs_univ.o \
	fast_subs_mars.o flux.o soilthprop.o makefile
	make -C $(SR) common
	$(CC) $(LFLAGS) mars_fast.o $(SR)tridag.o $(SR)conductionT.o \
	$(SR)grids.o soilthprop.o $(SR)psv.o $(SR)derivs.o \
	fast_subs_univ.o flux.o $(SR)conductionQ.o $(SR)psvco2.o \
	$(SR)generalorbit.o fast_subs_mars.o -o $@

insol_driver: insol_driver.o $(SR)generalorbit.o $(SR)flux_noatm.o \
	makefile
	make -C $(SR) common
	$(CC) $(LFLAGS) insol_driver.o $(SR)generalorbit.o \
	$(SR)flux_noatm.o


common-mars: marsorbit.o vapordiffusioni.o flux.o flux_mars.o \
	soilthprop.o 

%.o:	%.f makefile
	$(CC) $(CFLAGS) -c $<

%.o:	%.for makefile
	$(CC) $(CFLAGS) -c $<

%.o:	%.f90 makefile
	$(CC) $(CFLAGS) -c $< 


.PHONY: clean 

clean:
	rm -f  *.o  a.out  *.o90  *.mod \
	exper_fast mars_fast stabgrow_fast



