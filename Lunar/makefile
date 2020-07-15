CC = gfortran
CFLAGS = -O2 -Wall
LFLAGS = -O2 -Wall 

SR = ../Common/

oscidea1: oscidea_subs.o oscidea1.o subl_h2o.o makefile 
	$(CC) $(LFLAGS) oscidea1.o subl_h2o.o oscidea_subs.o


common-lunar: subl_h2o.o makefile


%.o:	%.f90 makefile
	$(CC) $(CFLAGS) -c $<  


clean:
	rm -f  *.o  a.out *.mod





