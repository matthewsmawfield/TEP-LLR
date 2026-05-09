	SUBROUTINE RDLUN(RJD,FJD,NREFL,OBSY,isignal,np_break)
	IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      INCLUDE 'LICENSE-BSD3.inc'
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C  Wed Oct 19, 1988 14:18:00
C     1) REPLACE "COMPILER DOUBLE PRECISION" WITH
C        "IMPLICIT DOUBLE PRECISION (A-H,O-Z)".
C     2) MAKE SURE STATEMENTS DON'T EXTEND BEYOND COLUMN 72.
C     3) REMOVE IN-LINE COMMENTS.
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C  Fri Oct 28, 1988 09:29:15
C     1) MOVE DATA STATEMENT TO "BLKNP.FOR", AND FIX "STOP [S]".
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C  Wed Feb 1, 1989 08:25:49
C     1) REPLACE CALLS TO SCANF AND SCANB WITH INTERNAL READS.
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C 08 Mar 91, 1720: add parameter to error for calling HP-UX exit.
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C 06 Oct 95      : Copy the detector type to an output variable. rlr.
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C  RDLUN IS RESPONSIBLE FOR READING AND DECODING INFORMATION FROM INPUT 
C  LUNAR DATA RECORDS AND MAKING IT AVAILABLE TO OTHER ROUTINES.
C
C  REVISIONS:
C       10/25/84 - DECODE LASER FREQUENCY FROM Z CARD AND PASS IT TO 
C                  REFRACTION ROUTINE AS WAVELENGTH IN MICRONS.
C       02/26/85 - CONVERT FROM WRITING IP BACK TO DISK & REREADING TO
C                  CALLING SCANF.
C       10/29/85 - SET UP FOR USE IN NORMAL POINTING PGM
C       12/03/86 - LUNAR '86 FORMAT. V10. RLR.
C	05/01/87 - Lunar '86 format, mod 1 (mini-normalpt). v11. rlr.
C	06/08/88 - CORRECT HANDLING OF FREQ OFFSET! AND USE SCANB FOR
C		   MOST FIELDS. V12. RLR.
C	11/26/91 - Pick up all detail times and o-c's for S:N calcs. rlr
C	12/03/91 - Allow non-signal records to be passed back for 
C		   S:N calculations. rlr.
C	10/26/99 - Y2K fixes. rlr.
C	07/26/07 - Total rewrite for CRD format. rlr.
C$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$ rlr

	character*(*) sccsid
	parameter (sccsid = "%W%\t%G%")

	INCLUDE 'MODEOP.INC'

	INCLUDE 'PASSIT.INC'

	INCLUDE 'DEFCMN.INC'

	INCLUDE 'DEVCOM.INC'

	INCLUDE 'ENVPAR.INC'

	INCLUDE 'AZLTRN.INC'

	INCLUDE 'NOISE.INC'

	INCLUDE 'MET.INC'

	INCLUDE '../include/crd.inc'
	INCLUDE '../include/crd_MLRS.inc'

        INTEGER YEAR,MON,DAY,HOUR,MIN, doy
	integer itrimlen, lengtht
	integer prt_year, reshour, resmin, iressec, jressec
	double precision fr_day, ressec
	character line*512

	character detector*1
C       np_break*1 -> np_break 010408 for f95
	logical np_break;

	np_break= .false.
	isignal= 0
 10     READ(21,1000,END=99) LINE
 1000   FORMAT(A512)

C       Write the normalpoint headers
        IF (modeop .eq. 0 .and.
     &     (LINE(1:2) .EQ. 'h1' .OR.
     &      LINE(1:2) .EQ. 'h2' .OR.
     &      LINE(1:2) .EQ. 'h3' .OR.
     &      LINE(1:1) .EQ. 'c' .OR.
CC     &      LINE(1:2) .EQ. '12' .OR.
CC     &      LINE(1:2) .EQ. '20' .OR.
     &      LINE(1:2) .EQ. '21' .OR.
     &      LINE(1:2) .EQ. '40' .OR.
CC     &      LINE(1:2) .EQ. '50' .OR.
     &      LINE(1:2) .EQ. '60' .OR.
     &      LINE(1:2) .EQ. '00')) then
          lengtht=itrimlen(line)
          write(23,'(a)') line(1:lengtht)
        ENDIF

C       Write the fullrate records
        IF (LINE(1:2) .EQ. 'h1' .OR.
     &      LINE(1:2) .EQ. 'h2' .OR.
     &      LINE(1:2) .EQ. 'h3' .OR.
     &      LINE(1:2) .EQ. 'h4' .OR.
     &      LINE(1:1) .EQ. 'c' .OR.
     &      LINE(1:2) .EQ. '10' .OR.
     &      LINE(1:2) .EQ. '12' .OR.
     &      LINE(1:2) .EQ. '20' .OR.
     &      LINE(1:2) .EQ. '21' .OR.
     &      LINE(1:2) .EQ. '30' .OR.
     &      LINE(1:2) .EQ. '40' .OR.
     &      LINE(1:2) .EQ. '50' .OR.
     &      LINE(1:2) .EQ. '60' .OR.
     &      LINE(1:2) .EQ. '00' .OR.
     &      LINE(1:2) .EQ. '91' .OR.
     &      LINE(1:2) .EQ. '92' .OR.
     &      LINE(1:2) .EQ. '94') then
          lengtht=itrimlen(line)
          write(22,'(a)') line(1:lengtht)
        ENDIF

        IF (LINE(1:2) .EQ. 'h2') then
	  CALL read_h2 (LINE)
	  if (cdp_pad_id .eq. 7080) obsy= 71112
	  if (cdp_pad_id .eq. 7086) obsy= 71111
	  if (cdp_pad_id .eq. 0) obsy= 71110	! Replace '0' w/ pad id...

	ELSE IF (LINE(1:2) .EQ. 'h3') then
	  CALL read_h3 (LINE)
	  nrefl= sic- 100
	  wmicron= xmit_wavelength/1.e3
CC	  sdqi= data_quality_ind

	ELSE IF (LINE(1:2) .EQ. 'h4') then
	  CALL read_h4 (LINE)
C	  nshots= ?
C       Specify 'normalpoint'
	  if (modeop .eq. 0) then
            data_type= 1
            call write_h4 (line)
            lengtht= itrimlen(line)
            write(23,'(a)') line(1:lengtht)
          endif

	ELSE IF (LINE(1:2) .EQ. 'c0') then
	  CALL read_c0 (LINE)

	ELSE IF (LINE(1:2) .EQ. '10') then
	  CALL read_10 (LINE)
	  CALL grtodoy (start_year-1900, start_mon, start_day, doy)
          CALL doytomjd (start_year- 1900, DOY, rjd)
CC	write(*,*) "jd",start_year, start_mon, start_day, doy, rjd
	  rjd= rjd+ 2400000.d0
	  fjd= d10_sec_of_day/86400.d0+ 0.5d0
	  if (fjd .gt. 1.d0) then
	    rjd= rjd+ 1.d0
	    fjd= fjd- 1.d0
	  endif
	  if (filter_flag .eq. 2) isignal= 1
	  else	isignal= 0
	  ORANGE= d10_time_of_flight
        RORANGE = 0
        ELDEL  = 0
        GEODEL = 0
CC?        UNCERT = 0
CC	write(*,*)"n10:",rjd,fjd,orange,oresid

	  IF (DEB(20)) THEN
	    WRITE(DUNIT,1010) RORANGE,ELDEL,GEODEL,UNCERT,
     1		TOFF,FREQ,PRESS,ORESID,ORANGE,isignal
 1010	    FORMAT(1X,'RDLUN: ',F15.4,2F10.4,F8.4,
     1          F10.6,F8.0,F8.2,F8.4,F15.13,i2)
	  ENDIF

	ELSE IF (LINE(1:2) .EQ. '20') then
	  CALL read_20 (LINE)
CC	  SBYS = SBYS*1.D-4
	  humid= humidity
	  TEMP = temperature
	  TOFF = 0.d0
	  PRESS= pressure
C Record met, unless already too many. Need to keep #500 unused.
          if (n_met .lt. 499) then
	    n_met= n_met+ 1
            met_sod(n_met)= d20_sec_of_day
            met_h(n_met)= humidity
            met_t(n_met)= temperature
            met_p(n_met)= pressure
            met_o(n_met)= value_origin
          endif
CC	  FKSX  = FKSX *1.D-4
CC	  FKSX2 = FKSX2*1.D-2
CC	  RKSX  = RKSX *1.D-4
CC	  RKSX2 = RKSX2*1.D-2

	ELSE IF (LINE(1:2) .EQ. '40') then
	  CALL read_40 (LINE)

C  Read but don't write this record. New one written in putnp.
	ELSE IF (LINE(1:2) .EQ. '50') then
	  CALL read_50 (LINE)

	ELSE IF (LINE(1:2) .EQ. '93') then
	  CALL read_93 (LINE)
C  use the post fit (actually recalc) residual unless it isn't there
	  ORESID= range_OmC_post
	  if (dabs(ORESID- (-1.d0)) .lt. 1.d-6) ORESID= range_OmC

C  Return with info for npt or just write a bad recd
	  if (ORESID < 1.d6) then
C  Residual file record
            prt_year= mod(start_year,100)
            fr_day= d93_sec_of_day/86400.d0
            reshour= INT(fr_day * 24.d0)
            resmin=  INT((fr_day - reshour/24.d0) * 1440.d0)
            ressec=  (fr_day- reshour/24.d0- resmin/1440.d0)
     .                  * 86400.d0
            iressec= idint(ressec)
            jressec= idnint((ressec-iressec)*1.d7)

            write(24,4050) prt_year,start_mon,start_day,
     .                reshour,resmin,iressec,jressec,
     .                isignal, int(ORESID*1.d3)	! to psec
 4050       format('p',6i2,i7,1x,i1,1x,i9)
	    return
	  endif

C	  These unmatched '93' recds get written here
          lengtht= itrimlen(line)
          write(22,'(a)') line(1:lengtht)

	ELSE IF (LINE(1:2) .EQ. '94') then
	  CALL read_94 (LINE)
	  np_break= .true.
          GO TO 100
	ENDIF
	GO TO 10

C  END OF FILE
 99     NREFL= -1
        np_break= .true.
 100    RETURN
C
989	CONTINUE
	CALL ERROR(1,errmsg,IERR,IERR)
C
C	THE
	END

