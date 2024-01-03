CC = gfortran
CFLAGS = -O2 -Wall
LFLAGS = -O2 -Wall

#CC = ifort
#CFLAGS = -O2 #-warn all -warn notruncated_source -nogen-interface
#LFLAGS = -O2 -warn all

#CC = pgfortran
#CFLAGS = -fast
#LFLAGS = -fast


all: conductionQ+iter.o conductionQ+smooth.o conductionQ+volt.o makefile


%.o:	%.f90 makefile
	$(CC) $(CFLAGS) -c $< 



.PHONY: clean

clean:
	rm -f  *.o  *.mod 



