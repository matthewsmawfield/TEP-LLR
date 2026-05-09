	program np12
	implicit double precision (a-h,o-z)
      INCLUDE 'LICENSE-BSD3.inc'
C *****************************************************************************
c Copyright (c) 2017, The University of Texas at Austin
c All rights reserved.
c
c Redistribution and use in source and binary forms, with or without
c modification, are permitted provided that the following conditions are met:
c
c 1. Redistributions of source code must retain the above copyright notice,
c this list of conditions and the following disclaimer.
c
c 2. Redistributions in binary form must reproduce the above copyright notice,
c this list of conditions and the following disclaimer in the documentation
c and/or other materials provided with the distribution.
c
c 3. Neither the name of the copyright holder nor the names of its contributors
c may be used to endorse or promote products derived from this software without
c specific prior written permission.
c
c THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
c AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
c IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
c ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE
c LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
c CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
c SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
c INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
c CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
c ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF
c THE POSSIBILITY OF SUCH DAMAGE.
c
C *****************************************************************************
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C  Mon Aug 29, 1988 09:50:47
C     1) REPLACE "COMPILER DOUBLE PRECISION" WITH
C        "IMPLICIT DOUBLE PRECISION (A-H,O-Z)".
C     2) MAKE SURE STATEMENTS DON'T EXTEND BEYOND COLUMN 72.
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C  Wed Oct 19, 1988 14:28:17
C     1) REMOVE IN-LINE COMMENTS.
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C  Mon Oct 31, 1988 15:29:19
C     1) CHANGE TYPE REAL TO TYPE DOUBLE PRECISION.
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C  Thu Nov 3, 1988 15:29:52
C     1) COMMENT OUT CALLS TO "SVGET", "NAMCHG", AND "UNIST".
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C  Wed Mar 29, 1989 08:09:44
C     1) KLUDGE D.P. RELATIONAL EXPRESSIONS.
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C  Wed May 3, 1989 10:14:09
C     1) ADD FAKE AZOFF/ALTOFF UNTIL IT CAN BE READ INTO
C        /SITE/ BLOCK
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C  Oct 89 - tlr
C     1) Correct handling of clock drift and of geometric
C        correction @ Mt. Fowlkes site.
c     2) fix "14:29:60.0" problem
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C  Nov 89 - tlr
C     1) Must multiply by cdr to convert azoff to radians.
C     2) Add "oobscode=obsy".
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
c	Dec 89  -  tlr
c  add command line input (normal point or recalc, file names, etc.)
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c  Tue Apr 17 09:43:45 CDT 1990 - tlr
c    Use proper geometric delay (geodt) for 71111 and 71112.
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c 05 Sep 90, 1011: remove "/us1" from path names
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c 14 Sep 90, 1529: call pmodrep instead of modrep.
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c 22 Feb 91, 1607: put "continue" before format for xdb.
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c 28 Feb 91, 1321 - use new pathnames: /data/lib, /data/llr/pred.
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c 08 Mar 91, 1738: terminate program w/call to exit through ceror.
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c 23 Apr 91, 1543 through  <:r! vidate>
c	Make output to screen like old NOVA programs.  Work on use of
c	pathnames.  Allow up to 500 normalpoints and take error exit
c	if more.  Around statement label 600, check for no data
c	records read.
c 31 Jul 91,     : correct output by writing long output to listing
c	at all times and limiting NOVA screen output to debugging.
c	rlr.
c 05 Nov 91,     : When there are no points in a npt bin, simply go on;
c	do not write and error and stop. rlr.
c 06 Nov 91,     : Allow fixed length npt bins with starting time based
c	on time of day. rlr
c 26 Nov 91,	 : Compile signal & noise stats and call s:n calculation 
c	routine. rlr
c 06 Oct 95,     : Get detector from run header and place into mini npt. rlr.
c 14 May 99,     : Add support for cstg (ilrs) normalpoints. rlr.
c 26 Oct 99,     : Y2K correction to output. rlr.
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
C:%TITLE NP -- PRODUCE LUNAR NORMALPOINTS
C:%NOSOURCE
C:%G LUNAR RANGE AND POSITION PREDICTION PROGRAM
C:%T PURPOSE
C:THIS PROGRAM IS DESIGNED TO PRODUCE LUNAR NORMALPOINTS
C:UNDER A GENERALIZED SET OF CONDITIONS.
C:FILES DETAILING MODEL PARAMETERS, REFLECTOR, TRANSMITTER, & RECEIVER
C:COORDINATES, AND CONTROL PARAMETERS ARE PRODUCED BY PROGRAM 'PREDCNTRL'.
C
C  REVISIONS:
C       10/25/84 - ADD ITERDIR TO ALLOW PREDICTIONS FORWARD OR BACKWARD
C                  IN TIME.
C                - PASS RDLUN VLIGHT FOR WAVELENGTH CALCS.
C       06/18/85 - Expand file names to 14 words; print to o/p;
C                  print o-c to crt
C       07/12/85 - Check for xcvr .ne. xmtr. rlr
C       10/29/85 - Improve modularity. rlr
C       10/31/85 - Convert to do normalpointing. rlr
C       06/09/86 - fix round-off problem. v01. rlr.
C       12/09/86 - LUNAR '86 FORMAT. V10. RLR.
C	07/21/87 - CHANGES TO ALLOW MULTIPLE REFL COORD SETS (/RFLINFO/ &
C		   /CNTINFO/), TOGGLE UT1 TIDE (/MODINFO/) AND HANDLE
C		   POINT-ANGLE-DEPENDENT MLRS TERM BEFORE & AFTER MOVE
C		   (71111. & 71112.). V11. RLR.
C	09/01/87 - GET AND PRINT SOURCE OF EARTH ORIENTATION PARAMETERS
C		   FROM SUBLASR. V12. RLR.
C
C$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$ RA & RLR

	character sccsid*(*)
	parameter (sccsid = "@(#)np12.f	1.15\t06/29/93")

C****************  INCLUDE STATEMENTS FOLLOW  ****************

	include 'ATMAUX.INC'

	INCLUDE 'AZLTRN.INC'

	INCLUDE 'BLD.INC'

	INCLUDE 'BOUNCE.INC'

	INCLUDE 'CETB1.INC'

	INCLUDE 'CETB2.INC'

	INCLUDE 'CETB4.INC'

	INCLUDE 'CNTINF.INC'

	include 'COMARG.INC'

	INCLUDE 'DATA2.INC'

	INCLUDE 'DDMU1.INC'

	INCLUDE 'DEFCMN.INC'

	INCLUDE 'DEVCOM.INC'

	INCLUDE 'DIRCTN.INC'

	INCLUDE 'EPHREC.INC'

	INCLUDE 'ENVPAR.INC'

	INCLUDE 'FUNBLK.INC'

	INCLUDE 'INTRNG.INC'

	INCLUDE 'LABELS.INC'

	INCLUDE 'LIB.INC'

	INCLUDE 'LIOCMN.INC'

	INCLUDE 'MODEOP.INC'

	INCLUDE 'MODINF.INC'

	INCLUDE 'NOISE.INC'

	INCLUDE 'NPHLD.INC'

	INCLUDE 'OBSINF.INC'

	INCLUDE 'PASSIT.INC'

	INCLUDE 'PDCMN.INC'

	INCLUDE 'PNT.INC'

	INCLUDE 'PRMAT.INC'

	INCLUDE 'PVHLD.INC'

	INCLUDE 'RFLINF.INC'

	INCLUDE 'SITE.INC'

	INCLUDE 'SRFINF.INC'

	INCLUDE 'U1UC.INC'

	INCLUDE 'USNOC.INC'

      character*16 srf
      common/labelx/srf(96),ifl(12)

C****************  END OF INCLUDE STATEMENTS  ****************

CCCCCC  COMMON STATEMENTS FROM FUNARG CCCCCCCCCCCCCCCCCCCCC
      common /falook/ fa(26)
CCCCCC  END OF FUNARG CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

cxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
	integer rename
cxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
      double precision lhatrn(6),lharcv(6),azalrcv(6)
      double precision obcdsv(10)
      double precision mlrsgc

      integer d1(6),onrefl
      character*255 filnami,filtmp, pathop
c "renamfil" is for determining if output file should replace lunar i/o file.
	character*1 eosrc, renamfil, detector
	character*9 vers

      common/twopi/ twopi,cdr,casr,ctsr

C...FOR KUSTNER TO SAVE INFO
	common/kuscmn/ frstim,tt(4),jd1,jd2,jdstep
	logical frstim
	double precision jd1,jd2,jdstep

	logical anydeb, dontread

	continue
5001	format(a)

CC	pathfil='/data/lib/'
CC	patheph='/data/llr/pred/'
CC	pathlog='/data/llr/pred/'
	pathfil='data/lib/'
	patheph='data/pred/'
	pathlog='data/pred/'
	lenfil=itrimlen(pathfil)
	leneph=itrimlen(patheph)
	lenlog=itrimlen(pathlog)
	dontread= .false.
c-----------------------------------------------------------------
	call getargs
	if(nargs .EQ. 0) then
c	  interactive input.
	  write(ioutdev,'('' Enter log file name:'')')
	  read(indev,5001)clfillog
	else
c	  command-line input: cheap way of preventing output to screen.
	  ioutdev=-1
CC	  ioutdev=24
CC	  open(24,status='scratch') This left files all over...
CC	  open(24,"/tmp/lun.scratch",status='unknown')
	endif
c-----------------------------------------------------------------
c  Data statement in eulcomnp can't properly handle implied do loop
c  with character substring.  Therefore we initialize "ip" here.
	do 161 i161=1,130
161	ip(i161:i161)=' '

CCCCCCCCCCCCCCCCCCCCCCCCCCC---------------------------
	DO 1 I1=1,20
1	DEB(I1)=.FALSE.
CCCCCCCCCCCCCCCCCCCCCCCCCCC---------------------------
c  Must multiply by cdr to convert to radians.
	azoff  = -89.d0*cdr
	altoff =   0.d0*cdr
CCCCCCCCCCCCCCCCCCCCCCCCCCC---------------------------

c  Open log file clfillog.  This file may be printed after
c  the user exits the normal point program.  We may want to give the
c  user the option of changing the output device to 6: line printer.
c  (P.S. the default value of "NPLOGFIL" is 10.)
	lenclog = itrimlen(clfillog)
	if(clfillog(1:1) .EQ. '/') then
	  pathop = clfillog(1:lenclog)
	else
	  pathop = patheph(1:leneph)//clfillog(1:lenclog)
	endif
	open(nplogfil,file=pathop,access='append')
c   position listing file for appending.
ccc use "access='append'" above instead of "gotoeof"
ccc	call gotoeof(nplogfil,nrec)
C  summary file is taken from command line ONLY.
	if (sflg) then
	    open(npsumfil,file=clfilsum,access='append')
	    write(npsumfil,1010) 
 1010	    format('Normalpoint(s):')
	endif

c:%s normalpoint driver
C	vers='30 Jul 91'
C	vers='05 Nov 91'
	vers='14 May 99'

c:read control, model, & reflector files
c slr sitefile for system configuration flag and site ID
c      call sitefilein(isch,isite)

c lunar files for everything else.
      call parin(vlite,obcdsv,vers,modeop,filnami,filtmp)

c  change the definition of dunit and open debug file.
CCCCCCCCCCCCCCCCCCCCCCCCCCC---------------------------
	anydeb = .FALSE.
	DO 2 I2=1,20
2	if(DEB(I2)) anydeb = .TRUE.
	dunit = 33
	if(anydeb) open(dunit,file='npt-debug')
c:initialize ephemeris read package
      call ephinit(aukm,aukm)

c:display model parameters
      call pmodrep
      norms= 0

c:write listing header
95	continue
	norms= norms+ 1
	if (ioutdev .gt. 0) write(ioutdev,1033) norms
	write(nplogfil,1033) norms
 1033 format(//1x,'NORMALPOINT GROUP #',i3/
     1          1X,'RECD',8X,'JD',14X,'GREG DATE',7X,'RF OBSY',4X,
     2          'OBS RANGE',4X,'REFR',6X,'O-C',5X,'X',7X,'Y',
     3          7X,'DUT',2X,'EOS',
     4          2X,'LHA',6X,'DEC',5X,'ALT')

      nrcd= 0
      oobscode= 0.d0
      nsignal= 0
      nnoise= 0

c=====================================================================

c:observation loop
c:read next observation (unless there is already one in ip buffer)
100	continue
	if (dontread) go to 102
	call rdlun(rjd,fjd,nrefl,obsy,vlight,isignal,detector,isch,isite)

c  Work to do here if we are creating normalpoints by even time 
c  intervals starting at 0 hr utc
 102  if (tflg) then
	  dontread= .false.
	  if (nrcd .eq. 0 .and. ip(1:1) .eq. '3') then	! 1st detail recd of grp
c	  start at 1 hours UTC
	      tfs= dmod(FJD+ 0.5d0,1.d0)* 86400.d0	! data time sec utc
	      tnp= (idint(tfs)/idint(nptlen)+1)*nptlen	! time of next npt break
	  elseif (ip(1:1) .eq. '3') then			! detail record
	      ltfs= tfs					! save old time
	      tfs= dmod(FJD+ 0.5d0,1.d0)* 86400.d0	! data time in sec utc
	      if (tfs .lt. ltfs) tfs= tfs+ 86400.d0	! assume sequential data
	      if (tfs .ge. tnp) then			! if data after next brk
		  tnp= tnp+ nptlen			! update to next break
		  if (nrcd .gt. 0) then			! if we have any data
		      dontread= .true.			! we just read a detail
		      go to 600				! normalpoint it
		  endif
	      endif
          elseif (ip(1:1) .eq. '1') then
	      go to 100					! read the next recd
          elseif (nrefl .lt. 0) then
	      go to 600					! all done...
	  endif
c  did we read a header set? if it's not the first one, form a norm pt
      elseif (ip(1:1) .eq. '1') then
          if (ntype .ne. 0) ntype= 1
	  if (nrcd .gt. 0) then
	      go to 600					! normalpoint it
	  else 
	      go to 100					! read the next recd
	  endif
c  are we through? if so, form a normalpoint
      elseif (nrefl .lt. 0) then
	  go to 600
      endif

c  copy data for S:N calculations
      if (isignal .eq. 1 .and. modeop .eq. 0) then
          nsignal= nsignal+ 1
          stime(nsignal)= rjd+ fjd
          somc(nsignal)= oresid
      elseif (isignal .eq. 0 .and. modeop .eq. 0) then
          nnoise= nnoise+ 1
          ntime(nnoise)= rjd+ fjd
          nomc(nnoise)= oresid
          go to 100
      endif

c  did we read a header set? if it's not the first one, form a norm pt
c      if (ip(1:1) .NE. '1') go to 104
c      if (nrcd.GT.0) go to 600
c      go to 100

c  are we through? if so, form a normalpoint
c 104  if (nrefl.LT.0) go to 600
      nrcd= nrcd+ 1

c:load observatory coordinates
c  check for different observatory and xcvr .NE. xmtr
      if (obsy.EQ.oobscode) go to 110
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
c  Set "old observatory code" to obsy so we won't keep reading
c  the same observatory file with each observation loop.
	oobscode = obsy
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      do 105 nobs=1,10
        if (obsy.NE.obcdsv(nobs)) go to 105
        call datio(1,1,nobs)
        irecvr= 1
        if (trad.NE.rrad .OR. tlat.NE.rlat .OR. tlong.NE.rlong)
     1							irecvr= 2
c...set up coords for refraction
        gdlong= tlong*cdr
        gdlat= tgdlat*cdr
        glat= tlat*cdr
        hgt= 0.d0
        go to 110
  105 continue

      WRITE(ierrout,1034) OBSY
 1034 FORMAT(1X,'OBSERVATORY NOT RECOGNIZED OR NOT SELECTED: ',F6.0)
      GO TO 100

C  THESE PARAMETERS ARE MODIFIED BY POLAR MOTION EACH PASS
110   TRANSM(1) = TRAD
      TRANSM(2) = TLONG
      TRANSM(3) = TLAT
      RECEIV(1) = RRAD
      RECEIV(2) = RLONG
      RECEIV(3) = RLAT

C:LOAD TARGET COORDINATES
      NRFP1 = NREFL+1
C				ACTS AS FILTER
      IF (RFLPIK(NRFP1).EQV. .FALSE.) GO TO 450
      REFLCT(1) = SERAD(NRFP1)
      REFLCT(2) = SELON(NRFP1) 
      REFLCT(3) = SELAT(NRFP1)

C:SET UP TIMES FOR PREDICTION
      TJDINT= RJD
      TJDUTC= RJD+ FJD
      SUTCT= FJD*86400.D0
C...ITERATE BKWRDS FROM RETURN TIME
      IF (ITERDIR.LT.0) SUTCT= SUTCT+ ORANGE
      IF (NRCD.EQ.1) JDI= TJDINT

C:WRITE TTY RECORD INFO
	CALL JDTGR(RJD,FJD,D1(3),D1(1),D1(2),D1(4),D1(5),SECOND)
	if(second .GE. 59.95)then
c	  this is to avoid times such as 14:29:60.0
	  CALL JDTGR(RJD,FJD+1.d-6,D1(3),D1(1),D1(2),D1(4),D1(5),SECOND)
	endif
      if (ioutdev .gt. 0)
     1          WRITE(IOUTDEV,1035) NRCD,TJDUTC,(D1(I),I=1,5),SECOND
 1035 FORMAT(/1X,'RECORD:',I4,
     1          2X,'TIME:',F15.6, 3X,2(I2,'/'),I2,2X,2(I2,':'),F5.2)

C:DO RANGE CALCULATION
      CALL SUBLASR(TJDINT,SUTCT,VLITE,UT1C,XPOLE,YPOLE,EOSRC,TTLTIM,
     1		LHATRN,AZALTRN,LHARCV,AZALRCV)

C:CORRECT FOR REFRACTION
      CALL REFRMM(AZALTRN(3)*CDR,TTLTIM,DALT,DTIM)
C		CONVERT FROM METERS TO NSEC
	DTIM= 2.D6*DTIM/VLIGHT
      TTLTIM= TTLTIM+ DTIM*1.D-9
      AZALTRN(3)= AZALTRN(3)+ DALT/CDR

C:CORRECT FOR 2.7 M FOLDED LIGHT PATH
	IF (ABS(OBSY-71110.D0).LT.0.1D0) TTLTIM=TTLTIM+ 121.8D0*1.D-9

C:CALCULATE FINAL O-C RESIDUAL
      OMC= (ORANGE-TTLTIM)*1.D9
      ONREFL= NREFL

C:WRITE DETAIL RECORDS BACK TO FILE
      CALL WRLUN(OMC,AZALTRN(2),AZALTRN(3))

C:SAVE CALCS AWAY FOR NORMALPOINTING
      IF (MODEOP.EQ.0)
     1	CALL NPAUG(JDI,(TJDINT-JDI)*86400.D0+SUTCT,TTLTIM,OMC)

C:PRINT REFLECTOR NUMBER, AZ, ALT, & RANGE ONTO TTY & LPT
CXXXXZZZZZZZZZZZZZZZ   REMOVE THE 2 CONTINUATION LINES BELOW.  XZXZZZZ
      if (anydeb .and. ioutdev .gt. 0) 
     1     WRITE(IOUTDEV,1050) NREFL,AZALTRN(2),AZALTRN(3),TTLTIM,OMC
CCC     1		,LHATRN(2),AZALTRN(3)
      if (anydeb) WRITE(NPLOGFIL,1050) NREFL,AZALTRN(2),AZALTRN(3),TTLTIM,OMC
CCC     1		,LHATRN(2),AZALTRN(3)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
 1050 FORMAT(5X,I1,2F10.4,F15.10,F9.2)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C 1050 FORMAT(5X,I1,2F10.4,F15.10,D13.6,2(1X,2D11.4))
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C     Y2K correction
      D1(3)= D1(3)+ 1900
      if (ioutdev .gt. 0) WRITE(IOUTDEV,1055) NRCD,TJDUTC,
     1          (D1(I),I=1,5),SECOND,NREFL,OBSY,
     2          ORANGE,DTIM,OMC,
     3          XPOLE,YPOLE,UT1C,EOSRC,
     4          LHATRN(2),LHATRN(3),AZALTRN(3)
      WRITE(NPLOGFIL,1055) NRCD,TJDUTC,(D1(I),I=1,5),SECOND,
     1          NREFL,OBSY,
     2          ORANGE,DTIM,OMC,
     3          XPOLE,YPOLE,UT1C,EOSRC,
     4          LHATRN(2),LHATRN(3),AZALTRN(3)
 1055 FORMAT(I4,1X,F18.10,1X,2(I2,'/'),I4,1X,2(I2,':'),F6.2,
     1          1X,I1,1X,F6.0,
     2          1X,F12.10,1X,F7.3,F9.2,
     3          2(1X,F7.4),1X,F8.5,1X,A1,
     4          1X,F7.4,1X,F8.4,1X,F7.4)
 
  450 GO TO 100

C========================= END OF OBS LOOP ===========================

C:NORMALPOINT PROCESSING
C...IF THIS IS RECALC-ONLY MODE, DON'T BOTHER
600	continue
C	if(nrcd .EQ. 0) then
C	  errmsg = ' no data records read.'
C	  call error(2,errmsg,0,0)
C	endif
	if(nrcd .EQ. 0) go to 699
	IF (MODEOP.NE.0) GO TO 699

	call sncalc

        CALL MAKNORM(TJDINT,SUTCT,OMC,A,B,C,SD,SE,SPAN)

c Get the basic geometric delay w/o point-angle terms from last record.
	if(abs(obsy-71111.d0) .LT. 0.1d0) then
	  geodt=geodel-
     1		mlrsgc((azaltrn(2)-24.5079)*cdr,azaltrn(3)*cdr)
	elseif(abs(obsy-71112.d0) .LT. 0.1d0) then
CCC          if (tjdint.lt.2449139.5) then
          if (tjdint.lt.2449078.5) then
            geodt=geodel-
     1          mlrsgc(azaltrn(2)*cdr+azoff,azaltrn(3)*cdr+altoff)
          else
            geodt=geodel-
     1          mlrsgct(azaltrn(2)*cdr+azoff,azaltrn(3)*cdr+altoff)
CC	WRITE(6,*) "last recd geodel, geodt",geodel,geodt
          endif
	endif

C:SET UP STATION (& REFL) COORDINATES
      TRANSM(1)= TRAD
      TRANSM(2)= TLONG
      TRANSM(3)= TLAT
      RECEIV(1)= RRAD
      RECEIV(2)= RLONG
      RECEIV(3)= RLAT

C:CALCULATE RANGE AND POINT ANGLES AT NORMALPOINT TIME
      CALL SUBLASR(TJDINT,SUTCT,VLITE,UT1C,XPOLE,YPOLE,EOSRC,
     1		TTLTIM,LHATRN,AZALTRN,LHARCV,AZALRCV)

C:CALC REFRACTION
      CALL REFRMM(AZALTRN(3)*CDR,TTLTIM,DALT,DTIM)
      DTIM= 2.D6*DTIM/VLIGHT
      TTLTIM= TTLTIM+ DTIM*1.D-9
      AZALTRN(3)= AZALTRN(3)+ DALT/CDR
      RANG= TTLTIM+ OMC*1.D-9

C:CALCULATE MOUNT GEOM OFFSET
	if(abs(obsy-71111.d0) .LT. 0.1d0) then
	  geodel=geodt+
     1		mlrsgc((AZALTRN(2)-24.5079D0)*CDR,AZALTRN(3)*CDR)
	elseif(abs(obsy-71112.d0) .LT. 0.1d0) then
CC          if (tjdint.lt.2449139.5) then
          if (tjdint.lt.2449078.5) then
            geodel= geodt+
     1          mlrsgc(AZALTRN(2)*CDR+AZOFF,AZALTRN(3)*CDR+ALTOFF)
          else
            geodel= geodt+
     1          mlrsgct(AZALTRN(2)*CDR+AZOFF,AZALTRN(3)*CDR+ALTOFF)
          endif
CC	WRITE(6,*) "npt recd geodel, geodt",geodel,geodt,AZALTRN(2),
CC     x	AZOFF,AZALTRN(3),ALTOFF
	endif

C:FORMAT AND WRITE THE STUFF
      CRANGE= RANG*(1.D0+FREQ/8.64D11) + (ELDEL+GEODEL)/1.D9
      CALL PUTNP(TJDINT,SUTCT,BH,NSHOTS,ONREFL,OBSY,
     1		RANG,CRANGE,OMC,ELDEL,GEODEL,TOFF,UNCERT,SBYS,SD,RSD,
     2		SPAN,SDQI,AZALTRN(2),AZALTRN(3),detector,isch,isite,ntype)
C Ntype '1' send out header. Only need one per run.
      if (ntype .eq. 1) ntype= 2

C:REPORT THE RESULTS
      CALL NREPT(TJDINT,SUTCT,RANG,OMC,A,B,C,SD,SE,UNCERT,SPAN,LHATRN)
      NORMS= NORMS+ 1

C...ARE WE THROUGH (LAST DATA SET -> NREFL<0)?
 699  IF (NREFL.GE.0) GO TO 95

C======================= END OF NORM POINT CREATION =======================

c-------------------------------------------------------------------
c   write the cstg normalpoints to file.
      if (ntype .eq. 2) 
     1  CALL PUTNP(TJDINT,SUTCT,BH,NSHOTS,ONREFL,OBSY,
     2		RANG,CRANGE,OMC,ELDEL,GEODEL,TOFF,UNCERT,SBYS,SD,RSD,
     3		SPAN,SDQI,AZALTRN(2),AZALTRN(3),detector,isch,isite,3)

			if(nargs .EQ. 0)then
c   (In case of interactive input, take care of I/O files.)
C:RENAME FILES TO GIVE OUTPUT FILE THE INPUT FILENAME
	renamfil = 'n'
5002	format(/,' ***********************************',/,
     1 ' Do you want to replace "',a,'" with the output file?')
	lenami=itrimlen(filnami)
	if(filnami(1:1) .EQ. '/') then
	  pathop = filnami(1:lenami)
	else
	  pathop = patheph(1:leneph)//filnami(1:lenami)
	endif
	if (ioutdev .gt. 0) write(ioutdev,5002)pathop
	read(indev,'(A1)') renamfil
c	.............................................
	if(renamfil.EQ.'y' .OR. renamfil.EQ.'Y') then
	  jgx = nullterm(filnami,255)
	  kgx = nullterm(filtmp,255)
	  irx = rename(filtmp,patheph(1:leneph)//filnami)
c         -------------------------------------------
	  if (irx .LT. 0) then
	     call perror(' np12 - call rename')
	  else
	     if (ioutdev .gt. 0) write(ioutdev,5003)
5003	     format(1x,'  lunar I/O file replaced.')
	  endif
c         -------------------------------------------
	else
	  if (ioutdev .gt. 0) write(ioutdev,5004)
5004	  format(1x,'   output is in file "npout";',
     1                ' lunar I/O file unchanged.')
	endif
c	.............................................
			endif
c-------------------------------------------------------------------
C:CLOSE FILES
	CLOSE(21,IOSTAT = IERR)
	CLOSE(22,IOSTAT = IERR)
	CLOSE(23,IOSTAT = IERR)
	CLOSE(24,IOSTAT = IERR)
c  ceror calls HP-UX "exit"
	call ceror(0)
C
C	   THE
	   END

