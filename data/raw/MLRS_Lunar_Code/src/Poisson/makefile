#  makefile for Poisson
# 2/90 RLR

F77= f77
CC= cc
FFLAGS= -g -I../include
CFLAGS= -I../include


# CRD format Poisson (all targets):
PSNCSRC= Poisson_crd.f lsq.f ../utilib/read_crdf.f ../utilib/read_crd_MLRSf.f \
	../utilib/write_crdf.f ../utilib/write_crd_MLRSf.f
PSNCOBJ= Poisson_crd.o lsq.o ../utilib/read_crdf.o ../utilib/read_crd_MLRSf.o \
	../utilib/write_crdf.o ../utilib/write_crd_MLRSf.o

Poisson_crd:	$(PSNCOBJ)
	$(F77) $(FFLAGS) -o Poisson_crd $(PSNCOBJ) ../utilib/cexit.o

clean:
	rm -f Poisson_crd *.o

all:
	make Poisson_crd


install:
	mv -f Poisson_crd ../../bin
