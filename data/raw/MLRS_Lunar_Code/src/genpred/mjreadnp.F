C--------------------------------------------------------------------
C
C  READE for MIT PEP ephemeris
C
C  history:
C	10/17/94 - restrict limits on input date too early or late.
C		   The previous 20 day limit was too restrictive and 
C		   unneccessary. rlr
C--------------------------------------------------------------------
	subroutine reade(jed,tsec,ierr)
	implicit double precision (a-h,o-z)

	character*(*) sccsid
	parameter (sccsid = "%W%\t%G%")

	include 'CNTINF.INC'

	INCLUDE 'DATA2.INC'

	INCLUDE 'DEVCOM.INC'

	INCLUDE 'LIB.INC'

	common /rdecmn/ ojd,bjd(3)
	save /rdecmn/

	common /mitcmn/ emptr,munptr,nutptr,libptr,jd1,jd2,meqinc,
     1		head(64)
      double precision jd1,jd2,meqinc,mithld(20)
      integer emptr,munptr,nutptr,libptr,head
      equivalence (emptr,mithld(1))
	save /mitcmn/
 
      double precision hld(6,10),hldn(2,10),hldl(3,10)
      double precision jed,tsec

#ifdef MIT
C  .5 DUE TO DEFINITION OF TSEC & JED. MIT USED FRACTION UTDAY+ JULIAN DAY 
 10   tjd = jed + tsec/86400.d0 + .5d0
CC	write(*,*) 'tjd=',tjd, "jd1=",jd1, 'jd2=',jd2

CC	if (tjd.ge.jd1+20.d0) go to 20
	if (tjd.ge.jd1+2.d0) go to 20
          ierr = 1
          goto 6790
 
20	continue
CC	if (tjd.le.jd2-20.d0) go to 30
	if (tjd.le.jd2-2.d0) go to 30
          ierr = 2
          goto 6790
 
30	if (abs(tjd-ojd).le.1.d0) go to 100
          ojd = tjd 
          dt = tjd - jd1
 
C  EARTH/MOON BARYCENTER
      call readd(2,emptr+(int(dt/4)-4)*6,hld,60,8,jerr)
      do 40 i=1,10
      do 40 j=1,3
          pos(i,j,1) = hld(j,i)
          vel(i,j,1) = hld(j+3,i)
 40   continue
      bjd(1) = jd1 + (int(dt/4)-4)*4
 
C  MOON
      call readd(2,munptr+(int(dt*2)-4)*6,hld,60,8,jerr)
      do 50 i=1,10
      do 50 j=1,3
          pos(i,j,2) = hld(j,i)
          vel(i,j,2) = hld(j+3,i)
 50   continue
      bjd(2) = jd1 + (int(dt*2)-4)/2.d0 
 
C  NUTATIONS
      call readd(2,nutptr+(int(dt*2)-4)*2,hldn,20,8,jerr)
      do 60 i=1,10
      do 60 j=1,2
 60   xnut(i,j) = hldn(j,i)
 
C  LIBRATIONS
      call readd(2,libptr+(int(dt*2)-4)*3,hldl,30,8,jerr)
      do 70 i=1,10
      do 70 j=1,3
 70   xlib(i,j) = hldl(j,i)
      bjd(3) = bjd(2)
 
C  NOW SOLVE
100	continue
	call mitsol(jed,tsec,bjd,ierr)
6790	continue
#else
      WRITE (*,*) "MIT Ephemeris is no longer supported"
      STOP
#endif
C
c	   the
      end

      SUBROUTINE mitopen

      WRITE (*,*) "MIT Ephemeris is no longer supported"
      STOP
      END
