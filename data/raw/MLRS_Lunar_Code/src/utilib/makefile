# Makefile for utilib
# History:
# 	11/07/11 - remove get_target. Replaced everywhere by target_tools. rlr.
#
FFLAGS = -c -g -I../include -I/usr/include -I.
F= f77
C= gcc

.o:	.f

ALLFSRC= read_crdf.f read_crd_MLRSf.f write_crdf.f write_crd_MLRSf.f \
	grtodoyf.f calc_pmmf.f

ALLFOBJ= read_crdf.f read_crd_MLRSf.f write_crdf.f write_crd_MLRSf.f \
	grtodoyf.o calc_pmmf.o

clean:
	rm -f *.o core

all:	 $(ALLFOBJ)
	$(F) $(FFLAGS) $(ALLFOBJ)
	$(C) $(FFLAGS) getfield.c cexit.c

install:
