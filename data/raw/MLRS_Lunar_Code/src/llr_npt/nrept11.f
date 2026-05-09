	SUBROUTINE NREPT(TJDINT,SUTCT,NPRANGE,OMC,A,B,C,SD,SE,UNCERT,
     X			SPAN,LHATRN)
	IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      INCLUDE 'LICENSE-BSD3.inc'
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C  Wed Oct 19, 1988 15:34:34
C     1) REPLACE "COMPILER DOUBLE PRECISION" WITH
C        "IMPLICIT DOUBLE PRECISION (A-H,O-Z)".
C     2) REMOVE IN-LINE COMMENTS.
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C  Mon Oct 31, 1988 13:18:09
C     1) REMOVE "Z" FORMAT SPECIFIER.
C     2) REMOVE "DINT" REFERENCES.
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C  Fri Nov 4, 1988 13:47:28
C     IN /NPHLD/, "TIME" <- "TIM" BCS CONFLICT WITH /ATMAUX/
C      Jul 30, 1991
C     Add Summary file handling. rlr.
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C  WRITE A REPORT ON THE NORMALPOINT FORMATION TO THE LISTING FILE.
C
C	12/18/86 - INCREASE DIGITS OF JD PRINTED W/O EXCEEDING
C		   ACCURACY OF THE SYSTEM (14-16 DIGITS). V10. RLR.
C	11/25/87 - FIX A COUPLE OF MINOR ERRORS IN HANDLING REJECTED
C		   DATA. V11. RLR.
C	10/13/97 - Printout starting "Time:" had second printed as
C		   wrong type. rlr.
C       03/02/98 - Fix handling of jd in RECD printout. lower precision output
C                  had been wrong because T was not included in JDLO. rlr.
C       10/26/99 - Y2K fix to output. rlr.
C
C$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$ RLR

	character*(*) sccsid
	parameter (sccsid = "@(#)nrept11.f	1.3\t07/31/91")

	INCLUDE 'DEFCMN.INC'

	INCLUDE 'DEVCOM.INC'

	INCLUDE 'NPHLD.INC'

	INCLUDE 'MODINF.INC'

	INCLUDE 'COMARG.INC'

	DOUBLE PRECISION NPRANGE,sec,SPAN,LHATRN(3)
	integer year,mon,day,hour,min

	if (ioutdev .gt. 0) WRITE(IOUTDEV,990)
	WRITE(NPLOGFIL,990)
 990	FORMAT(/1X,'RECD',11X,'JD',15X,'COR RANGE',8X,'PREFIT O-C',
     X		5X,'POSTFIT O-C')
	NTOT= NPTS+ NREJ
	DO 100 I=1,NTOT
	  T= ABS(TIM(I))
	  RESIN= A+ B*T+ C*T*T
	  DRES= RESID(I)- RESIN
CCCCC	  JDLO= DINT(JDI/100.D0)
	  JDLO= ((JDI+ T/86400.D0)/100.D0)
	  RJD= (JDI-JDLO*100.D0)+ T/86400.D0
          if (ioutdev .gt. 0)
     x	  WRITE(IOUTDEV,1000) I,JDLO,RJD,RANGE(I),RESID(I),DRES
	  WRITE(NPLOGFIL,1000) I,JDLO,RJD,RANGE(I),RESID(I),DRES
 1000	  FORMAT(1X,I3,5X,I5,F11.8,5X,F13.10,5X,F10.2,5X,F10.2)

	  IF (TIM(I) .LT. 0.D0) THEN
	    if (ioutdev .gt. 0) WRITE(IOUTDEV,1002)
	    WRITE(NPLOGFIL,1002)
	  ENDIF
 1002	  FORMAT(2X,'REJECTED')
	  if (ioutdev .gt. 0) WRITE(IOUTDEV,1004)
	  WRITE(NPLOGFIL,1004)
 1004	  FORMAT(1X)

 100	CONTINUE

C	calc jdlo in case np spans day
CCCCC	JDLO= DINT(TJDINT/100.)
	JDLO= (TJDINT/100.D0)
	call jdtgr(tjdint,SUTCT/86400.d0,year,mon,day,hour,min,sec)
        year= year+ 1900

	TIMEN= (TJDINT-JDLO*100.D0)+ SUTCT/86400.D0
        if (ioutdev .gt. 0)
     x	WRITE(IOUTDEV,1010) JDLO,TIMEN,year,mon,day,
     X		hour,min,sec,NPRANGE,OMC,UNCERT,SD
	WRITE(NPLOGFIL,1010) JDLO,TIMEN,year,mon,day,
     X		hour,min,sec,NPRANGE,OMC,UNCERT,SD
 1010	FORMAT( /10X,'Julian date:',10X,I5,F15.12, '(',I4,2('/',i2),
     X		     1x,2(i2,':'),f10.7,')'
     X		/10X,'Norm time delay (sec):',F16.13
     X		/10X,'Time delay res (nsec):',F11.4
     X		/10X,'Obs. uncertainty (nsec):',F9.4
     X		/10X,'Obs. consistency (nsec):',F9.4)

C	to_cm= 2.d9/vlight
	if (sflg) then
	    write(npsumfil,1020) year,mon,day,hour,min,idint(sec),
     X	        omc,sd,uncert,ntot,nrej,span/60.d0,lhatrn(2),lhatrn(3)
 1020	    format('Time: ',i4,2('/',i2),2x,2(i2,':'),i2, 
     X	    	' O-C(nsec):',f11.4,' RMS:',f9.4,' Uncert:',f9.4/
     X		'Fitted Obs: ',i3, ' Rejects: ',i3, ' Span (min): ',f5.2,
     X		' HA: ',f9.4,' Dec: ',f9.4)
	endif

C	set up for next time!
	NPTS= 0
C
C	   THE
	   END
