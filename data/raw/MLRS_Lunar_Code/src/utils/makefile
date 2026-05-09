INCLUDES = -I../include -I.
 
# This variable lists the required libraries.
CFLAGS = -g -I/usr/include -I. $(INCLUDES) 
CC= gcc

.o:	.c

#---------------------------------------------------------------------------

FRSSRC= frd_strip.c ../utilib/read_crd.c ../utilib/read_crd_MLRS.c

FRSOBJ= frd_strip.o ../utilib/read_crd.o ../utilib/read_crd_MLRS.o

frs:	 $(FRSOBJ)
	$(CC) $(CFLAGS) $(FRSOBJ) ../utilib/write_crd.o \
		-lm -o frd_strip
#---------------------------------------------------------------------------

clean:	
	rm -f *.o frd_strip
all:	
	make frs
install:
	mv frd_strip ../../bin
