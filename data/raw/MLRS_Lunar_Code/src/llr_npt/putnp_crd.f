	SUBROUTINE PUTNP(TJDINT,SUTCT,NSHOTS,NREFL,
     1		OBSY,RANG,CRANGE,OMC,ELDEL,GEODEL,TOFF,UNCERT,
     2          calrms,sigrms,passrms,sigskew,sigkurtosis,sigPmM,
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
C	??-??-?? - Add cstg normalpoint format. rlr.
c	07-27-07 - Convert to CRD format. rlr.
C
C$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$ RLR 11/85

	character*(*) sccsid
	parameter (sccsid = "%W%\t%G%")

	INCLUDE 'NOISE.INC'

	INCLUDE 'NPHLD.INC'

	INCLUDE 'DEVCOM.INC'

	INCLUDE 'MET.INC'

	INCLUDE '../include/crd.inc'
	INCLUDE '../include/crd_MLRS.inc'

        CHARACTER line*512
	CHARACTER sdqi*1

	integer itrimlen, lengtht
        integer nnpts
	integer prt_year, reshour, resmin, iressec, jressec
	double precision ressec, fr_day
        double precision sigsumsq
	double precision sigskew, sigkurtosis, sigPmM
        common /putnpstuff/ nnpts, sigsumsq
        data nnpts /0/, sigsumsq /0.d0/
     
C  Write a normalpoint and residual record
        if (ntype .eq.1 .or. ntype .eq. 2) then
	  d11_sec_of_day= sutct- 43200.d0
	  if (d11_sec_of_day .lt. 0) 
     &		d11_sec_of_day= d11_sec_of_day+ 86400.d0
	  d11_time_of_flight= RANG
	  d11_sysconfig_id= d10_sysconfig_id
	  d11_epoch_event= d10_epoch_event
          d11_detector_channel= 0
	  np_window_length= span
	  num_ranges= npts
	  bin_rms= sigrms*1.e3
	  bin_skew= sigskew
	  bin_kurtosis= sigkurtosis
	  bin_PmM= sigPmM
	  if (DABS(bin_PmM- (-1.d0)) .gt. 1.d-6) bin_PmM= sigPmM*1.d3
CC	write(*,*) "sigPmM= ",sigPmM," bin_PmM= ",bin_PmM
C	  return_rate= npts/(span*10) ! parameterize '10'
	  return_rate= snratio
C Get the appropriate met record (if there wasa  change)
CC	write(*,*) "putnp: ",n_met_next, n_met
          do i=n_met_next,n_met
            if (met_sod(i) .lt. d11_sec_of_day .and. 
     &          (met_sod(i+1) .gt. d11_sec_of_day .or. i+1 .gt. n_met))
     &      then
              d20_sec_of_day= met_sod(i)
              pressure= met_p(i)
              temperature= met_t(i)
              humidity= met_h(i)
              value_source= met_o(i)
              call write_20 (line)
              lengtht= itrimlen(line)
              write(23,'(a)') line(1:lengtht)
              n_met_next= i+1
CC	write(*,*) "putnp>: ",n_met_next, n_met
            endif
          enddo
          call write_11 (line)
          lengtht= itrimlen(line)
          write(23,'(a)') line(1:lengtht)

          prt_year= mod(start_year,100)
          fr_day= d11_sec_of_day/86400.d0
          reshour= INT(fr_day * 24.d0)
          resmin=  INT((fr_day - reshour/24.d0) * 1440.d0)
          ressec=  (fr_day- reshour/24.d0- resmin/1440.d0)
     .                  * 86400.d0
          iressec= idint(ressec)
          jressec= idnint((ressec-iressec)*1.d7)
          write(24,4050) prt_year,start_mon,start_day,
     .                reshour,resmin,iressec,jressec,
     .                INT(OMC*1.d3)       ! to psec
 4050     format('n',6i2,i7,1x,'1',1x,i9)
C  93 for npt?

          nnpts= nnpts+ 1
          sigsumsq= sigsumsq+ sigrms*sigrms
	endif

C  Wrap it up
        if (ntype .eq. 3) then
C	  d50_sysconfig_id= d10_sysconfig_id
          sess_rms= 0.d0
          if (nnpts .ge. 1)
     &      sess_rms= dsqrt (sigsumsq/nnpts)*1.e3
          sess_skew= -1
          sess_kurtosis= -1
          sess_PmM= -1
C	  sdqi is read from input fr file
          call write_50 (line)
          lengtht= itrimlen(line)
          write(23,'(a)') line(1:lengtht)
          call write_h8 (line)
          lengtht= itrimlen(line)
          write(23,'(a)') line(1:lengtht)
          call write_h9 (line)
          lengtht= itrimlen(line)
          write(23,'(a)') line(1:lengtht)
        endif
C
C	   THE
	   END

