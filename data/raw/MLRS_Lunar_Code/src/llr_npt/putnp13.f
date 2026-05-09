	SUBROUTINE PUTNP(TJDINT,SUTCT,BH,NSHOTS,NREFL,
     1		OBSY,RANG,CRANGE,OMC,ELDEL,GEODEL,TOFF,UNCERT,
     2          calrms,sigrms,passrms,
     2		SPAN,SDQI,AZ,ALT,detector,isch,isite,ntype)
	IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      INCLUDE 'LICENSE-BSD3.inc'
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C  Wed Oct 19, 1988 15:27:09
C     1) REPLACE "COMPILER DOUBLE PRECISION" WITH
C        "IMPLICIT DOUBLE PRECISION (A-H,O-Z)".
C     2) MAKE SURE STATEMENTS DON'T EXTEND BEYOND COLUMN 72.
C     3) REMOVE IN-LINE COMMENTS.
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C  Thu Nov 3, 1988 15:28:51
C     1) COMMENT OUT CALLS TO "FLA1".
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C  Fri Nov 4, 1988 13:34:48
C     1) CHANGE FORMAL PARAMETER "RANGE" TO "RANG" BECAUSE OF
C        CONFLICT WITH /NPHLD/.
C     2) CHANGE /NPHLD/ TIME TO TIM BCS CONFLICT WITH /ATMAUX/
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C  Mon Jan 30, 1989 10:05:54
C     1) CHANGE BHSKEL, BH TO CHAR*130
C     2) REPLACE "FLA1" WITH INTERNAL READS.
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C  Mon May 1, 1989 11:18:07
C     1) TRY ROUNDING INSTEAD OF TRUNCATING WHEN PUTTING
C        THINGS IN NPSKEL ETC.
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C      Dec,   1989
C     1) Put signal to noise ratio in output normalpoint.  rlr.
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C  06 Oct 1995 : Write detector type to mini-normalpoint. rlr.
C  19 Mar 1999 : write cstg normalpoints. rlr.
C  26 Oct 1999 : Y2K fix. rlr.

C  ROUTINE TO SET UP AND WRITE LUNAR NORMALPOINTS IN THE ARCHIVAL AND
C  MINI FORMATS.
C
C  REVISIONS:
C       12-09-86 - IMPLEMENT LUNAR '86 FORMAT. V10. RLR.
C	02-09-87 - FIX NEGATIVE AZ HANDLING. V11. RLR.
C	05-05-87 - LUNAR FORMAT '86 MOD 1. V12. RLR.
C	08-22-88 - Add np time span to mini-np. v13. rlr.
C
C$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$ RLR 11/85

	character*(*) sccsid
	parameter (sccsid = "%W%\t%G%")

	INCLUDE 'NOISE.INC'

	INCLUDE 'NPHLD.INC'

	INCLUDE 'DEVCOM.INC'

CCC	INCLUDE 'PASSIT.INC'

        CHARACTER SDQI*1, SCHAR*20, BHSKEL*130, BH*130, detector*1
        CHARACTER seldel*1
        INTEGER YEAR,MON,DAY,HOUR,MIN

C  max=1500. See below and NPHLD.
        character cstgh*55, cstgn(1500)*54
        integer nnpts
        double precision calsumsq,sigsumsq
        save cstgh, cstgn, nnpts, calsumsq, sigsumsq

	INCLUDE 'ENVPAR.INC'

C  HANDLE NEGATIVE AZ
	XAZ= AZ
	IF (XAZ.LT.0.D0) XAZ= 360.D0+ XAZ

C...ARCHIVAL NORMALPOINT ( & BURST HEADER)
C----------------------------------------------------
        NPSKEL(1:1)= '4'
        CALL JDTGR(TJDINT,(SUTCT-TOFF)/86400.D0,
     1			YEAR,MON,DAY,HOUR,MIN,SEC)
	IYEAR = YEAR+1900
	ISEC = SEC*1.D7 + 0.5D0
C  PUT INTEGER VALUES STRAIGHT INTO NPSKEL. (MUST SCALE
C  SECONDS BY 1.D-7)
	WRITE(NPSKEL(3:23),5701)IYEAR,MON,DAY,HOUR,MIN, ISEC
5701	FORMAT(I4,4I2,I9)
C----------------------------------------------------
C  WRITE CRANGE TO "SCHAR" AND THEN MOVE ALL BUT RADIX TO NPSKEL.
C  (MUST READ W/F14.0 AND SCALE BY 1.D-13)
	DRANGE = CRANGE + 0.5D-13
	WRITE(SCHAR,5702)DRANGE
5702	FORMAT(F15.13)
	NPSKEL(24:24) = SCHAR(1:1)
	NPSKEL(25:37) = SCHAR(3:15)
C----------------------------------------------------
C  WRITE SCALED INTEGERS OF UNCERT, ELDEL, GEODEL.
	IELDEL  = ELDEL* 1.D4 + 0.5D0
	IGEODEL = GEODEL*1.D4 + 0.5D0
	IUNCERT = UNCERT*1.D4 + 0.5D0
	WRITE(NPSKEL(39:61),5703)IELDEL,IGEODEL,IUNCERT
5703	FORMAT(I9,I8,I6)
C----------------------------------------------------
C  WRITE UNMODIFIED INTEGER "NPTS" DIRECTLY TO NPSKEL
	WRITE(NPSKEL(75:77),5704) NPTS
5704	FORMAT(I3)
C----------------------------------------------------
      	write(npskel(78:80),5704) IDNINT(snratio*10)
C----------------------------------------------------
C  WRITE INTEGER VERSION OF "SPAN" (UNSCALED) TO NPSKEL.
	ISPAN = SPAN + 0.5D0
	WRITE(NPSKEL(81:84),5705) ISPAN
5705	FORMAT(I4)
C----------------------------------------------------

C		QUALITY INDICATOR= DATA (FOR NOW)
        NPSKEL(85:85) = '1'

C----------------------------------------------------
C  WRITE "XAZ" AND "ALT" TO TEMP ("SCHAR") AND MOVE ALL BUT RADIX
C  POINT TO NPSKEL.  READ W/F10.0 AND SCALE BY 1.D-7.
	WRITE(SCHAR,5706) XAZ
5706	FORMAT(F11.7)
	NPSKEL(91:93)  = SCHAR(1:3)
	NPSKEL(94:100) = SCHAR(5:11)
C
	WRITE(SCHAR,5706) ALT
	NPSKEL(101:103)  = SCHAR(1:3)
	NPSKEL(104:110) = SCHAR(5:11)
C----------------------------------------------------
C  WRITE SCALED INTEGER VERSION OF "OMC" TO NPSKEL.
	IOMC = OMC*1000.D0 + 0.5D0
	WRITE(NPSKEL(111:119),5707) IOMC
5707	FORMAT(I9)
C----------------------------------------------------

C...Normalpoint burst header
	BHSKEL = BH
	BHSKEL(3:14) = NPSKEL(3:14)
C   Shot-by-shot resolution is the same as cal rms in slr parlance.
	ISBYS = calrms*1.D4 + 0.5D0
	WRITE(BHSKEL(23:28),5781) ISBYS
	WRITE(BHSKEL(42:46),5782) NSHOTS
5781	FORMAT(I6)
5782	FORMAT(I5)
	BHSKEL(68:68)= 'N'
	DO 11 I=82,89
 11	BHSKEL(I:I)= ' '
	DO 12 I=91,130
 12	BHSKEL(I:I)= ' '

	WRITE(22,1010) BHSKEL
        WRITE(22,1010) NPSKEL
 1010   FORMAT(A130)

        if (ntype .eq. 0) then
C..NOW CREATE MINI NORMAL POINT
C  NOTE THAT RANGE IS CORRECTED FOR ALL EFECTS EXCEPT REFRACTION
C  AND TIME IS CORRECTED FOR THE CLOCK OFFSET, GIVING US UTC
        NPSKEL(1:1)= '5'
        CALL JDTGR(TJDINT,SUTCT/86400.D0,YEAR,MON,DAY,HOUR,MIN,SEC)
C----------------------------------------------------
	IYEAR = YEAR+1900
	ISEC = SEC*1.D7 + 0.5D0
C  PUT INTEGER VALUES STRAIGHT INTO NPSKEL. (MUST SCALE
C  SECONDS BY 1.D-7)
	WRITE(NPSKEL(3:23),5701)IYEAR,MON,DAY,HOUR,MIN,ISEC
C----------------------------------------------------
C  WRITE "RANG" TO "SCHAR" AND THEN MOVE ALL BUT RADIX TO NPSKEL.
C  (MUST READ W/F14.0 AND MULTIPLY BY 1.D-13)
	DRANG = RANG + 0.5D-13
	WRITE(SCHAR,5702)DRANG
	NPSKEL(24:24) = SCHAR(1:1)
	NPSKEL(25:37) = SCHAR(3:15)
C----------------------------------------------------
C  LOAD 4 VARIABLES INTO "NPSKEL"
	IOBSY = OBSY + 0.5D0
	IUNCERT = UNCERT*1.D4 + 0.5D0
	WRITE(NPSKEL(38:52),5783) NREFL,IOBSY,NPTS,IUNCERT
5783	FORMAT(I1,I5,I3,I6)
C----------------------------------------------------
	NPSKEL(53:55) = NPSKEL(78:80)
C----------------------------------------------------
        NPSKEL(56:56) = SDQI
C----------------------------------------------------
	IPRESS = PRESS*100.D0 + 0.5D0
CCC  TRY ROUNDING THE TEMPERATURE IN STEAD OF TRUNCATING IT.
	ICELSIUS = ((TEMP-273.15D0)*10.D0) + 0.5D0
	IHUMID = HUMID + 0.5D0
	IMICRN = WMICRN*1.D4 + 0.5D0
	WRITE(NPSKEL(57:73),5784) IPRESS,ICELSIUS,IHUMID,IMICRN
5784	FORMAT(I6,I4,I2,I5)
C----------------------------------------------------

C		DATA SET VERSION
        NPSKEL(74:74)= NPSKEL(90:90)

C		NP time span (sec)
C----------------------------------------------------
	ISPAN = SPAN + 0.5D0
	WRITE(NPSKEL(75:78),5785) ISPAN
5785	FORMAT(I4)

C		Detector type
C----------------------------------------------------
	NPSKEL(79:79)= detector

C----------------------------------------------------
	NPSKEL(80:80)= ' '

C----------------------------------------------------
C  WRITE "XAZ" AND "ALT" TO TEMP ("SCHAR") AND MOVE ALL BUT RADIX
C  POINT TO NPSKEL.  READ W/F10.0 AND SCALE BY 1.D-7.
	WRITE(SCHAR,5706) XAZ
	NPSKEL(81:83)  = SCHAR(1:3)
	NPSKEL(84:90) = SCHAR(5:11)
C
	WRITE(SCHAR,5706) ALT
	NPSKEL(91:93)  = SCHAR(1:3)
	NPSKEL(94:100) = SCHAR(5:11)
C----------------------------------------------------
C  WRITE SCALED INTEGER VERSION OF "OMC" TO NPSKEL.
	IOMC = OMC*1000.D0 + 0.5D0
	WRITE(NPSKEL(101:109),5707) IOMC
C----------------------------------------------------

        WRITE(23,1011) NPSKEL(1:109)
1011	FORMAT(A109)
        endif

C CSTG header - once per run
        if (ntype .eq. 1) then
          icalshft= 0
          CALL JDTGR(TJDINT,SUTCT/86400.D0,YEAR,MON,DAY,HOUR,MIN,SEC)
CC Y2k Fix          CALL GRTJD(mod(YEAR,100),1,0,0,0,0.0,RJD,RJDF)
          CALL GRTJD(YEAR,1,0,0,0,0.0,RJD,RJDF)
          idoy= int(tjdint+ sutct/86400.d0- rjd- rjdf)
C FIXME check occupancy #
C          if (idint(obsy+0.05) .eq. 71112) isite= 70802419
          seldel= ' '
          if (ieldel < 0) seldel= '-'
C         Assume we  are using the Varian unless told otherwise
          idet= 0
          if (detector .eq. ' ') idet= 1
          if (detector .eq. 'a' .or. detector .eq. 'A') idet= 1
          if (detector .eq. 'b' .or. detector .eq. 'B') idet= 2
          if (detector .eq. 'c' .or. detector .eq. 'C') idet= 3
          if (detector .eq. 'd' .or. detector .eq. 'D') idet= 4
          if (detector .eq. 'e' .or. detector .eq. 'E') idet= 5
          if (detector .eq. 'f' .or. detector .eq. 'F') idet= 6
          if (detector .eq. 'g' .or. detector .eq. 'G') idet= 7
          if (detector .eq. 'h' .or. detector .eq. 'H') idet= 8

C         Subjective data quality indicator: map a-e to 1-5. Blank maps to 0.
          itmp= ICHAR(sdqi)
C         ...lower case
          isdqi= itmp- 96
C         ...upper case
          if (itmp .LT. 97) isdqi= itmp- 64
          if (sdqi .EQ. ' ') isdqi = 0

C         note that the detector code suffices as the sys config indicator.
          write(cstgh(1:53),1020) nrefl+100, mod(year,100), idoy,isite,
     1      int (wmicrn*1.d4),seldel,iabs(ieldel/10),icalshft,ISBYS/10,
     1      isch,idet, idint(passrms*1.d3+ 0.05d0),isdqi
 1020     format(i7,i2,i3,i8,i4,a1,i7,i6,i4,'246',i1,i1,i4,i1)
          nnpts= 0
          calsumsq= 0.0
          sigsumsq= 0.0
        endif

C CSTG normalpoint data record.
        if (ntype .eq. 1 .or. ntype .eq. 2) then
          nnpts= nnpts+ 1
          calsumsq= calsumsq+ calrms*calrms
          sigsumsq= sigsumsq+ sigrms*sigrms
          CALL JDTGR(TJDINT,SUTCT/86400.D0,YEAR,MON,DAY,HOUR,MIN,SEC)
          isec1= hour*3600+ min*60+ int(sec)
          isec2= (sec- int(sec))*1.e7
          irang1= (rang-2.0d0)*1.0d7
          irang2= ((rang-2.0d0)*1.0d7- irang1)*1.d5
          ipress= press*10.d0+ 0.05d0
          ihumid= humid+0.05d0
          ispan= span/300.d0+ 0.5001
          if (ispan .gt. 9) ispan= 9
C          read (npskel(78:80),1024) isnr
C 1024     format(i3)
          isnr= snratio*10
          if (isnr .gt. 99) isnr= 99
          if (isnr .lt. 0) isnr= 0
          write (cstgn(nnpts),1025) isec1,isec2,irang1,irang2,
     1      idint(sigrms*1.d3+0.5),ipress,
     1      idint(temp*10.d0+0.5D0),ihumid,npts,ispan,isnr
 1025     format (i5,i7,i7,i5,i7,i5,i4,i3,i4,'02',i1,i2)
          call checksum(cstgn(nnpts))
          DO i = 1,54
            IF (cstgn(nnpts)(i:i) .EQ. ' ') cstgn(nnpts)(i:i) = '0'
          ENDDO
          if(nnpts .GT. max) then
            errmsg=' putnp: too many normalpoints in this run.  '//
     1        'Number created and max are:'
            call error(2,errmsg,nnpts,max)
          endif

        endif

C Time to write out CSTG normalpoint with RSS of range residuals
        if (ntype .eq. 3) then
          passrms= dsqrt (calsumsq/nnpts)
          write (cstgh(39:42),1030) idint(passrms*1.e3+ 0.05)
          passrms= dsqrt (sigsumsq/nnpts)
          write (cstgh(48:51),1030) idint(passrms*1.e3+ 0.05)
 1030     format(i4)
          call checksum(cstgh)
          call checksum(cstgh)
          DO i = 1,54
            IF (cstgh(i:i) .EQ. ' ') cstgh(i:i) = '0'
          ENDDO
          cstgh(55:55)='1'
          WRITE(23,1031) cstgh
1031	  FORMAT("99999"/, A55)
          do i=1,nnpts
            WRITE(23,1033) cstgn(i)
1033	    FORMAT(A54)
          enddo
        endif
C
C	   THE
	   END

        SUBROUTINE CHECKSUM(LINE)
**********************************************************************
*
*PURPOSE
*       DETERMINE THE CHECKSUM FOR ONE NORMAL POINT DATA RECORD
*       AND PLACE THAT SUM IN THE RECORD
*
*INPUT AND OUTPUT VARIABLE
*       LINE CHARACTER*54
*
*       BY BRION CONKLIN BENDIX FIELD ENGINEERING DSG/OAS 4/90
*
**********************************************************************

        CHARACTER LINE*54

        INTEGER*4 I,SUM

        SUM = 0
        DO I = 1,52
                IF(LINE(I:I) .EQ. '1') SUM = SUM + 1
                IF(LINE(I:I) .EQ. '2') SUM = SUM + 2
                IF(LINE(I:I) .EQ. '3') SUM = SUM + 3
                IF(LINE(I:I) .EQ. '4') SUM = SUM + 4
                IF(LINE(I:I) .EQ. '5') SUM = SUM + 5
                IF(LINE(I:I) .EQ. '6') SUM = SUM + 6
                IF(LINE(I:I) .EQ. '7') SUM = SUM + 7
                IF(LINE(I:I) .EQ. '8') SUM = SUM + 8
                IF(LINE(I:I) .EQ. '9') SUM = SUM + 9
        ENDDO
C       SUM = JMOD(SUM,100)
        SUM = MOD(SUM,100)
CC      ENCODE(2,100,LINE(53:54)) SUM
        write(LINE(53:54),100) SUM
 100    FORMAT(I2.2)
        RETURN
        END

********************************************************************

