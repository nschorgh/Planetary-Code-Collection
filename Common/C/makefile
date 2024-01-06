CC = gcc
CFLAGS = -O2 -Wall
LFLAGS = -O2 -Wall 



testcrankT: testcrankT.o tridag.o conductionT.o grids.o makefile
	$(CC) -O2 testcrankT.o conductionT.o tridag.o grids.o -lm

testcrankQ: testcrankQ.o tridag.o conductionQ.o flux_noatm.o grids.o \
	makefile
	$(CC) -O2 testcrankQ.o conductionQ.o tridag.o flux_noatm.o \
	grids.o -lm



%.o:	%.c makefile
	$(CC) $(CFLAGS) -c $<


.PHONY: clean 

clean:
	rm -f  *.o  a.out 
	rm -f Tprofile z Tsurface


