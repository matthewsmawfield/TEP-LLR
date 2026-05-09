#
#    Makefile (jpl) for euler (genpred) program.
#    ----------------------------------------------.
#
F77= f77
FFLAGS= -g -I../include

NPTDIR = ../llr_npt

GPSRC= genpred.f dgetprdtim.f grtdoy.f pparin11.f icgetime.c \
eulcomnp.f djeulpkg.f error.f ceror.c dsublasr.f grtjd.f jdtgr01.f kustner01.f \
mjreadnp.f jjreadnp.f testeph.f misclibnp.f sbetidnp.f datio01.f lunio02.f pmodrep02.f \
usnortns.f opio.f
#
GPOBJ= genpred.o dgetprdtim.o grtdoy.o pparin11.o icgetime.o \
eulcomnp.o djeulpkg.o error.o ceror.o dsublasr.o grtjd.o jdtgr01.o kustner01.o \
mjreadnp.o jjreadnp.o testeph.o misclibnp.o sbetidnp.o datio01.o lunio02.o pmodrep02.o \
usnortns.o opio.o
#-----------------------------
genpred:    $(GPOBJ)
	$(F77) $(FFLAGS) -o genpred $(GPOBJ)

clean:
	rm -f *.o genpred core

all:
	make genpred

install:
	mv -f genpred ../../bin
