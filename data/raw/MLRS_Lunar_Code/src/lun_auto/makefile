#  makefile for lunar auto-segment and data quality programs
#  RLR
CFLAGS= -g -Wall -I../include
FFLAGS= -I2 -g -I../include

CC= cc
F77= g77

#=====================================================================
LACSRC= lun_auto_crd.c ../utilib/getfield.c \
	../utilib/read_crd.c ../utilib/write_crd.c \
	../utilib/read_crd_MLRS.c ../utilib/write_crd_MLRS.c 
LACOBJ= lun_auto_crd.o ../utilib/getfield.o \
	../utilib/read_crd.o ../utilib/write_crd.o \
	../utilib/read_crd_MLRS.o ../utilib/write_crd_MLRS.o 

lac:	$(LACOBJ)
	$(CC) $(CFLAGS) -o lun_auto_crd $(LACOBJ) -lm

#=====================================================================
LSCSRC= lun_sdqi_crd.c ../utilib/getfield.c \
	../utilib/read_crd.c ../utilib/write_crd.c
LSCOBJ= lun_sdqi_crd.o ../utilib/getfield.o \
	../utilib/read_crd.o ../utilib/write_crd.o

lsc:	$(LSCOBJ)
	$(CC) $(CFLAGS) -o lun_sdqi_crd $(LSCOBJ) -lm

#=====================================================================
all:
	make lac
	make lsc

clean:
	rm -f lun_auto_crd lun_sdqi_crd *.o

install:
	mv -f lun_auto_crd ../../bin
	mv -f lun_sdqi_crd ../../bin
