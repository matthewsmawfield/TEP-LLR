#  makefile for Bulla
# 01.04 RLR
# 
CFLAGS= -g -I../include
FFLAGS= -g -I../include

CC= gcc
#=====================================================================
BASRC= BullA.c
BAOBJ= BullA.o

bulla:$(BAOBJ)
	$(CC) $(FFLAGS) -o BullA $(BAOBJ) -lm
#=====================================================================
clean:
	rm -f core *.o 
#=====================================================================
all:
	make bulla

install:
	mv -f BullA ../../bin
