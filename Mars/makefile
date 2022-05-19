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



common-mars: marsorbit.o flux_mars.o


%.o:	%.f makefile
	$(CC) $(CFLAGS) -c $<

%.o:	%.f90 makefile
	$(CC) $(CFLAGS) -c $< 


.PHONY: clean 

clean:
	rm -f  *.o 

