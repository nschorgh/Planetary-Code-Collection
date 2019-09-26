CC = gcc
CFLAGS = -O2 -Wall
LFLAGS = -O2 -Wall 



testcrankQ: testcrankQ.o tridag.o conductionQ.o flux_noatm.o \
	grids.o makefile
	$(CC) -O2 testcrankQ.o conductionQ.o tridag.o flux_noatm.o \
	grids.o -lm



%.o:	%.c makefile
	$(CC) $(CFLAGS) -c $<


.PHONY: clean cleanrun 

clean:
	rm -f  *.o  a.out 

cleanrun:
	rm -f Tprofile z Tsurface


