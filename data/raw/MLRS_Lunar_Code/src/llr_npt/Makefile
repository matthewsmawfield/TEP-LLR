#  makefile for llr_npt
# 2/97 RR

F77= g77
#F77= f95
CC= gcc
FFLAGS= -g
CFLAGS= -g

# llr normalpoint:
NPTSRC= np12.f eulcomnp.f jeulpkg.f  \
                            error.f     \
                            ceror.c     \
                            grtjd.f     \
                            jdtgr01.f   \
                            kustner01.f \
                            mjreadnp.f \
                            jjreadnp.f \
                            misclibnp.f \
                            sbetidnp.f  \
   datio01.f lsq.f lunio02.f maknorm10.f mlrsgc01.f pmodrep02.f \
   npaug10.f nrept11.f parin11.f putnp13.f rdlun12.f refrmm.f  \
   saltne.f sublasr03.f wrlun12.f sncalc.f opio.f
#
NPTOBJ= np12.o eulcomnp.o jeulpkg.o  \
                            error.o     \
                            ceror.o     \
                            grtjd.o     \
                            jdtgr01.o   \
                            kustner01.o \
                            mjreadnp.o \
                            jjreadnp.o \
                            misclibnp.o \
                            sbetidnp.o  \
   datio01.o lsq.o lunio02.o maknorm10.o mlrsgc01.o pmodrep02.o \
   npaug10.o nrept11.o parin11.o putnp13.o rdlun12.o refrmm.o  \
   saltne.o sublasr03.o wrlun12.o sncalc.o opio.o
#
npt:	$(NPTOBJ)
#	g77 $(FFLAGS) -o llr_npt $(NPTOBJ)
#	f95 $(FFLAGS) $(NPTOBJ) -o llr_npt95 
	$(F77) $(FFLAGS) -o llr_npt $(NPTOBJ)


cosp:	cosp2fgnp.c
	gcc -o cosp2fgnp cosp2fgnp.c ../utilib/getfield.o \
	../utilib/grtojd.o ../utilib/grtodoy.o -lm

clean:
	rm -f *.o llr_npt cosp2fgnp
all:
	make npt
	make cosp

install:
	mv -f llr_npt /prod/bin
	mv -f cosp2fgnp /prod/bin
