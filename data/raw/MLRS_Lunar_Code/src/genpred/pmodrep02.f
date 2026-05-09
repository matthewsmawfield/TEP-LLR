	SUBROUTINE PMODREP
	IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      INCLUDE 'LICENSE-BSD3.inc'
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C  Wed Oct 19, 1988 15:09:48
C     1) REPLACE "COMPILER DOUBLE PRECISION" WITH
C        "IMPLICIT DOUBLE PRECISION (A-H,O-Z)".
C     2) MAKE SURE STATEMENTS DON'T EXTEND BEYOND COLUMN 72.
C     3) CHANGE OPTIONAL COMPILE STATEMENTS TO COMMENTS.
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C  Fri Oct 28, 1988 10:21:15
C     1) CHANGE HOLLERITH TO STRING; COMMENT OUT OVERLAY STATEMENT.
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C  Mon Oct 31, 1988 15:03:54
C     1) CHANGE TYPE REAL TO TYPE DOUBLE PRECISION.
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
c  3-6 Nov 89  -  tlr
c     Convert old modrep02 on HP9000 to pmodrep.  Change concerns
c     calls to jdtgr and its parameters in common /mdrep/.
c     Actually, add "pmodrep" code from euler.f to the existing
c     modrep02.f in np.dr and use variable modeop to decide whether
c     to execute prediction or normal point sections of code.
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c 30 May 90, 0939 - Fix calls to "jdtgr" to avoid "14:44:60.0" problem.
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c 19 Jul 90, 1122 - use i2.2 editing for time and date in format 1000.
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c 24 Aug 90, 1204: add ephemeris source and name to output.
c 24 Sep 90, 0928: fix output of ephemeris source (ephs)
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c 07 Feb 91, 1034: add eosrc to /mdrep/ so it can be printed.
c		   (by both prediction and normalpoint parts.)
c 30 Jul 91,     : Remove eosrc from npt output, as this is
c	only available on a data-record-basis. rlr
c 03 Nov 91,	 : Print sigma multiplier for npt fit. rlr.
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

C  PRINT MODEL PARAMETERS AND REFLECTOR COORDINATES
C
C  REVISIONS:
C	07-20-87 - HANDLE MULTIPLE REFLECTOR SETS & UT1TIDE FLAG. V01. RLR.
C
C$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$ RLR 10/85

	character*(*) sccsid
	parameter (sccsid = "@(#)pmodrep02.f	1.8\t01/06/92")
c#################################################################

	INCLUDE 'CNTINF.INC'

	INCLUDE 'COMARG.INC'

	INCLUDE 'DEFCMN.INC'

	INCLUDE 'DEVCOM.INC'

	INCLUDE 'MODEOP.INC'

	INCLUDE 'MODINF.INC'

	INCLUDE 'OBSINF.INC'

	INCLUDE 'RFLINF.INC'

	INCLUDE 'USNOC.INC'

c-----------------------------------------------------------------

	character*3 ephs
	integer d1(5),d2(5)
	double precision pcs(14)
	equivalence (pcs(1),c1)

c these put in common by euler (prd)
	DOUBLE PRECISION JDINT1,JDF1,JDINT2,JDF2
	integer rflset(5)
	character*1 eosrc
	common /mdrep/ rflset,jdint1,jdf1,jdint2,jdf2, eosrc

c#################################################################

	if(modeop .EQ. -1) then
c	this is prediction part of pmodrep
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c:print greeting & list of control & model parameters
	call jdtgr(jdint1,jdf1,d1(3),d1(1),d1(2),d1(4),d1(5),secs1)
	  if(secs1.GE.59.95d0) call jdtgr(jdint1,jdf1+1.d-9,
     1	  		d1(3),d1(1),d1(2),d1(4),d1(5),secs1)
	call jdtgr(jdint2,jdf2,d2(3),d2(1),d2(2),d2(4),d2(5),secs2)
	  if(secs2.GE.59.95d0) call jdtgr(jdint2,jdf2+1.d-9,
     1			d2(3),d2(1),d2(2),d2(4),d2(5),secs2)
        d1(3)= d1(3)+ 1900
        d2(3)= d2(3)+ 1900
        if (ioutdev .gt. 0)
     1	write(ioutdev,1000) rjd1,d1,secs1,rjd2,d2,secs2,stepjd,minelv
	write(eulog,1000) rjd1,d1,secs1,rjd2,d2,secs2,stepjd,minelv
1000	format(//
     1  1x,'INITIAL TIME: ',f13.4, 3x,2(i2.2,'/'),i4.4,2x,
     2  2(i2.2,':'),f4.1/
     3  1x,'FINAL TIME: ',  f15.4, 3x,2(i2.2,'/'),i4.4,2x,
     4  2(i2.2,':'),f4.1/
     5  1x,'EPHEMERIS POINT SEPARATION: ',f7.0/
     6  1x,'MINIMUM ELEVATION: ',f5.1/)
c
	if(ieph .EQ. 1)then
	  ephs = 'JPL'
	else
	  ephs = 'MIT'
	endif
	if (ioutdev .gt. 0) write(ioutdev,5011)ephs,fileph
	write(eulog  ,5011)ephs,fileph
5011	format(1x,a3,' ephemeris;  file name is:',/,1x,a)
c
	write(eulog,5010) twoway,polmotn,ut1tid,tolerr,itermax,
     1	vlight,aukm,emrati,k2lovi,
     2	pcs(13),pcs(14),pmtitl,(pcs(k),k=1,12),eosrc
5010	format(1x,'TWOWAY: ',l1,5x,'POLAR MOTION: ',l1,
     1	5x,'UT1 TIDE: ',L1,/
     1	1x,'TOLERANCE: ',d10.3,5x,'MAX NUMBER OF ITERATIONS: ',i4/
     2	1x,'SPEED OF LIGHT: ',d16.9,5x,'km/AU: ',d21.14/
     3	1x,'E/M: ',f11.8,5x,'K2 LOVE: ',f4.2
     4	/1x,'POLAR MOTION: K1 THROUGH K12 FOR MJD(UT1-UTC) ',f6.0,
     5	' AND MJD(A/B) ',f6.0/ 2x,'FROM ',a30/
     6	2x,'X: ',5(1x,f8.5)/ 2x,'Y: ',5(1x,f8.5)/ 2x,'UT1-UTC: ',
     7	 2(1x,f8.5),20x,'EOSRC="',a1,'"')

	DO 30 II=1,5
	  IF (.NOT.RFLPIK(II)) GO TO 30
	  WRITE(EULOG,5020) RFLSET(II),SERAD(II),SELON(II),SELAT(II)
5020      FORMAT(1X,'REFLECTOR ',I1,' COORDINATES:'/
     X	    3X, 'RADIUS:    ',F17.9/
     X	    3X, 'LONGITUDE: ',F12.8/
     X	    3X, 'LATITUDE:  ',F12.8)
30	CONTINUE
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
C:PRINT REFLECTOR & OBSERVATORY COORDINATES
      WRITE(EULOG,5030) TRAD,TLONG,TLAT,TGDLAT,RRAD,RLONG,RLAT,RGDLAT
 5030 FORMAT(
     X  1X,'TRANSMITTER COORDINATES:'/3X,'RADIUS  ',F15.9/3X,
     X      'LONGITUDE ',F12.8/
     X      3X,'LATITUDE  ',F12.8/3X,'GEODETIC LATITUDE ',F12.8/
     X      1X,'RECEIVER COORDINATES:'/3X,'RADIUS  ',F15.9/3X,
     X      'LONGITUDE ',F12.8/
     X      3X,'LATITUDE  ',F12.8/3X,'GEODETIC LATITUDE ',F12.8)
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

	else
c	this is normal point/recalc part of pmodrep
C-----------------------------------------------------------------
 1030   FORMAT(/1X,A18/
     X   1X,'TRANSMITTER COORDINATES:'
     X   /3X,'RADIUS  ',F15.9/3X,'LONGITUDE ',F13.9/
     X	 3X,'LATITUDE  ',F13.9/3X,'GEODETIC LATITUDE ',F13.9/
     X	 1X,'RECEIVER COORDINATES:'
     X   /3X,'RADIUS  ',F15.9/3X,'LONGITUDE ',F13.9/
     X	 3X,'LATITUDE  ',F13.9/3X,'GEODETIC LATITUDE ',F13.9)

C-----------------------------------------------------------------
      WRITE(NPLOGFIL,1010) CNTNAM,
     X  NAMMOD,TWOWAY,POLMOTN,UT1TID,TOLERR,ITERMAX,
     X	VLIGHT,AUKM,EMRATI,K2LOVI,
     X	CS(13),CS(14),PMTITL,(CS(K),K=1,12), sigmult
C     X	VLIGHT,AUKM,EMRATI,K2LOVI,BETAI,GAMMAI,CSNSS,EPOCH,CUTI,
 1010 FORMAT(1X,'CONTROL: ',A18/
     1  1X,'MODEL: ',A18,/1X,'TWOWAY: ',L1,5X,'POLAR MOTION: ',L1,
     2	5X,'UT1 TIDE: ',L1,/
     3	1X,'TOLERANCE: ',E10.3,5X,'MAX NUMBER OF ITERATIONS: ',I4/
     4	1X,'SPEED OF LIGHT: ',E16.9,5X,'AU/KM: ',E21.14/
     5	1X,'E/M: ',F11.8,5X,'K2 LOVE: ',F4.2,
     6	/1X,'POLAR MOTION: K1 THROUGH K12 FOR MJD(UT1-UTC) ',F6.0,
     7  ' AND MJD(A/B) ',F6.0/ 2X,'FROM ',A30/
     8	2X,'X: ',5(1X,F8.5)/ 2X,'Y: ',5(1X,F8.5)/ 2X,'UT1-UTC: ',
     9		2(1X,F8.5)/ 1X, 'N*SIGMA: ',f4.1)
C    X  5X,'BETA: ',E12.5,5X,'GAMMA: ',E12.5/
C    X	/1X,'COSINES OF LUNAR HARMONICS:'/3(3(3X,E11.4)/)
C    X	1X,'SINES OF LUNAR HARMONICS:'/3(3(3X,E11.4)/)
C    X	/1X,'UT VARIATIONS:'/1X,'EPOCH: ',F9.1/
C    X	1X,'LINEAR: ',E10.3/
C    X	1X,'SIN LONG PERIOD: ',E10.3,5X,'COS LONG PERIOD: ',E10.3/
C    X	7X,'SIN NODAL: ',E10.3,11X,'COS NODAL: ',E10.3/

C  COORDINATE SET TITLE:
      WRITE(NPLOGFIL,1015) RFSET
 1015 FORMAT(/1X,'REFLECTOR SET: ',A18)

	DO 3 II=1,5
	 NREFL= II- 1
	 IF (RFLPIK(II)) THEN
	  WRITE(NPLOGFIL,1020) NREFL,SERAD(II),SELON(II),SELAT(II)
1020	  FORMAT(1X,'REFLECTOR ',I1,' COORDINATES:'/
     X    3X, 'RADIUS:    ',F17.9/
     X    3X, 'LONGITUDE: ',F13.9/
     X    3X, 'LATITUDE:  ',F13.9)
	 ENDIF
3	CONTINUE

C  RECORD WHICH COORDINATE SETS ARE AVAILABLE
      DO 50 NOBS=1,10
C . . . . . . . . . . . . . . . . . . . . . . . . . . . BLOCK IF
        IF (OBSPIK(NOBS)) THEN
        CALL DATIO(1,1,NOBS)
C:PRINT REFLECTOR & OBSERVATORY COORDINATES
        WRITE(NPLOGFIL,1030) OBSITE,TRAD,TLONG,TLAT,TGDLAT,
     1		RRAD,RLONG,RLAT,RGDLAT
	ENDIF
C . . . . . . . . . . . . . . . . . . . . . . . . . . .END BLOCK IF
   50 CONTINUE

	endif

C
C	   THE
	   END

