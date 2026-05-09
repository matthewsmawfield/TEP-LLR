	subroutine jreade(rjed,tsec,ierr)
	implicit double precision(a-h,o-z)
      INCLUDE 'LICENSE-BSD3.inc'

c  This routine is the link between the normal point/prediction programs
c  and the JPL supplied routines PLEPH etc that read the JPL ephemeris.
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C 08 Mar 91, 1720: add parameter to error for calling HP-UX exit.
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

	character*(*) sccsid
	parameter (sccsid = "@(#)jjreadnp.f	1.4\t04/29/91")

	include 'CETB1.INC'

	include 'CETB2.INC'

	include 'CETB4.INC'

	include 'CNTINF.INC'

	INCLUDE 'DATA2.INC'

	INCLUDE 'DEVCOM.INC'

	INCLUDE 'LIB.INC'

	dimension rrd(6)

c  Get position and velocity vectors.
	itarg=3
CC	call pleph(rjed,tsec,itarg,icent,rrd,*9998)
	call pleph(rjed+(tsec/86400),itarg,icent,rrd)
	do 2 i2=1,6
2	tabout(i2,3)=rrd(i2)
c . . . . . . . .
	itarg=10
CC	call pleph(rjed,tsec,itarg,icent,rrd,*9998)
	call pleph(rjed+(tsec/86400),itarg,icent,rrd)
	do 4 i4=1,6
4	tabout(i4,11)=rrd(i4)

c  Get libration.
	itarg=15
CC	call pleph(rjed,tsec,itarg,icent,rrd,*9998)
	call pleph(rjed+(tsec/86400),itarg,icent,rrd)
	do 1 i1=1,3
1	libr(i1)=rrd(i1)

c  Get nutation.
	itarg=14
CC	call pleph(rjed,tsec,itarg,icent,rrd,*9998)
	call pleph(rjed+(tsec/86400),itarg,icent,rrd)
	do 3 i3=1,4
3	nut(i3)=rrd(i3)

	goto 10000
9998	continue
	errmsg=' '
	write(errmsg,'(f10.1,3x,f16.9)')rjed,tsec
	call error(2,errmsg,0,0)
10000	continue
c
c	   the
	   end
