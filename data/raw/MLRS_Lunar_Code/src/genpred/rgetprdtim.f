	subroutine getprdtim(pjdint,pjdf,vlite,ut1c,xpole,ypole,
     1		eosrc,ttltim,lhatrn,azaltrn,lharcv,azalrcv,moreprd,
     2							done,ndays)
	IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      INCLUDE 'LICENSE-BSD3.inc'
C  This routine is called by euler to determine initial time of
c  prediction run.  "Getprdtim" is only called when automatic
c  predictions are being done.  The input values of pjdint,pjdf are
c  the current time (that program started).  The output values of
c  pjdint,pjdf are the initial time of the prediction.  The next to
c  the last parameter, moreprd, indicates whether another prediction
c  is to be done after this one.  The last parameter, "done", is true
c  if there isn't even one more prediction left to do.  The other
c  parameters are for sublasr.			    6 July 90 - tlr
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
c 23 Jul 90, 1040: add call to calsubl; catch up on today, next and
c	standard (day after tomorrow) predictions.
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C 08 Mar 91, 1109: add parameter to error for calling HP-UX exit.
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C 25 Mar 91, 1638: fix inquire to use pathnames.
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C 31 May 91 - 13 Jun 91: Add "ndays" to argument list.  If ndays is
C   not zero, C   then ndays is added to the current date, and
C   predictions are done up to the new date.
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
	character*(*) sccsid
	parameter (sccsid = "@(#)getprdtim.f	1.6\t06/17/91")
	DOUBLE PRECISION LHATRN(6),AZALTRN(6),LHARCV(6),AZALRCV(6)
	CHARACTER*1 EOSRC, tmpnam*20
	logical moreprd, done, exsts, loopup
	common /twopi/twopi,cdr,casr,ctsr

C........... INCLUDE FILES .....................................

	INCLUDE 'CNTINF.INC'

	INCLUDE 'DEFCMN.INC'

	INCLUDE 'DEVCOM.INC'
 
	INCLUDE 'INTRNG.INC'

	INCLUDE 'OBSINF.INC'

	INCLUDE 'PDCMN.INC'

	INCLUDE 'PRDARGCOM.INC'

C...............................................................
c  Set limit on number of predictions to do or look for (set limit
c  to ten thousand if we are doing a certain number of days - ndays),
c  and say initially we are not done.  In case we're predicting for a
c  specified number of days, find the end time (vjd) for the automatic
c  prediction run.  Note that getprdtim expects pjdint,pjdf to be the
c  program execution time every time this routine is called.
	limnum = 3
	if(ndays .GT. 0) limnum = 10000
	done = .FALSE.
	vjd = pjdint+ndays + pjdf

c initially set the working time to the current time and
c see where we are.
	wjdint = pjdint
	wjdf = pjdf
	moreprd = .FALSE.
	call calsubl(wjdint,wjdf,wsutc,vlite,ut1c,xpole,ypole,eosrc,
     1		     ttltim,lhatrn,azaltrn,lharcv,azalrcv,
     2		     t2rise,riset,t2set,sett)
	if(azaltrn(3) .GT. minelv) then
	  nrx=1
	else
	  nrx=2
	endif
c	go to moonrise (forward or back)
	wsutc=riset
	call addtim(wjdint,wsutc,0.d0)

c we will come back to 21 if we haven't found an unpredicted moonrise.
21	continue
c,,,,,,,,,,,,, force time to integer multiple of stepjd ,,,,,,,,,,,,,,
	rstmp=wsutc/stepjd
	qstmp=dint(rstmp)
	if(dabs(qstmp-rstmp) .GE. 0.01d0)qstmp=qstmp+1.d0
	wsutc=dint(qstmp*stepjd+0.5d0)
c'''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''
11	wjdf=wsutc/86400.d0
c	check if moon has risen yet
	call calsubl(wjdint,wjdf,wsutc,vlite,ut1c,xpole,ypole,eosrc,
     1		     ttltim,lhatrn,azaltrn,lharcv,azalrcv,
     2		     t2rise,riset,t2set,sett)
	if(azaltrn(3) .GT. minelv) goto 12
c	moon not up yet; take another step.
	wsutc=wsutc+stepjd
	goto 11

c moon is up now; see if this prediction has been done.
12	call addtim(wjdint,wsutc,0.d0)
	wjdf=wsutc/86400.d0
	call jdtgr(wjdint,wjdf,iwyr,iwmon,iwda,iwhr,iwmin,wsec)
	if(wsec .GE. 59.95d0) call jdtgr(wjdint,wjdf+1.d-9,
     1		iwyr,iwmon,iwda,iwhr,iwmin,wsec)
	call grtdoy(iwyr,iwmon,iwda,iwdoy)
	tmpnam(1:14)='LUNxxx.le.'//cntnam(1:4)
	write(tmpnam(4:6),'(i3.3)')iwdoy
	inquire(file=patheph(1:leneph)//tmpnam(1:14),err=998,
     1		iostat=ios,exist=exsts)
c - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
	wjd = wjdint + wjdf
	if((ndays .GT. 0) .AND. (wjd .GT. vjd)) then
c We're doing number of days, and we have exceeded stop date.
	done = .TRUE.
	loopup = .FALSE.
c - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
	else
c Either doing standard number of predictions for auto mode (ndays.EQ.0),
c or have not exceeded stop date.
c . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
c If file exists, try the next moonrise.  If not, do prediction.
	if(exsts)then
	nrx=nrx+1
	if(nrx .LE. limnum)then
c	  not all prediction times checked; advance to next moonrise
	  delsecs = (hrate*86400.d0) - (2.d0*stepjd)
	  call addtim(wjdint,wsutc,delsecs)
	  loopup = .TRUE.
	else
c	  We've checked all candidates and there are no more
c	  predictions to do; set flag to quit.
	  done = .TRUE.
	  loopup = .FALSE.
	endif
c . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
	else
	pjdint=wjdint
	pjdf=wjdf
	moreprd = nrx .LT. limnum
	loopup = .FALSE.
	endif
c . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
	endif
c - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
	if(loopup) goto 21
	goto 10000
c
998	errmsg=' "getprdtim": Error in inquire; iostat is'
	call error(1,errmsg,ios,ios)

10000	continue
c
c	   the
	   end



	subroutine addtim(rjdint,rsutc,rsecs)
	IMPLICIT DOUBLE PRECISION (A-H,O-Z)
c  Add to input time in integer julian date and utc seconds the third
c  parameter in seconds.  Then adjust so that utc seconds is non-negative
c  and less than one days worth (0 =< rsutc < 86400).
	logical done

	rsutc = rsutc + rsecs
	done=.FALSE.
1	continue
	  if(rsutc .GE. 86400.d0) then
	   rjdint=rjdint+1.d0
	   rsutc=rsutc-86400.d0
	  elseif(rsutc .LT. 0.d0) then
	   rjdint=rjdint-1.d0
	   rsutc=rsutc+86400.d0
	  else
	   done = .TRUE.
	  endif
	if(.NOT.done)goto 1
c
c	   the
	   end


	subroutine calsubl(rjdint,rjdf,rsutc,vlite,ut1c,xpole,ypole,
     1		eosrc,ttltim,lhatrn,azaltrn,lharcv,azalrcv,
     2		t2rise,riset,t2set,sett)
	IMPLICIT DOUBLE PRECISION (A-H,O-Z)

C  This routine calls sublasr for the input time rjdint,rjdf and
c  returns the azaltrn array as well as the riset and sett with
c  corresponding t2rise and t2set.  UTC seconds is also returned.
c 23 Jul 90, 1036 - tlr
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

	DOUBLE PRECISION LHATRN(6),AZALTRN(6),LHARCV(6),AZALRCV(6)
	CHARACTER*1 EOSRC
	common /twopi/twopi,cdr,casr,ctsr

C........... INCLUDE FILES .....................................

	INCLUDE 'CNTINF.INC'

	INCLUDE 'DEFCMN.INC'

	INCLUDE 'DEVCOM.INC'
 
	INCLUDE 'INTRNG.INC'

	INCLUDE 'OBSINF.INC'

	INCLUDE 'PDCMN.INC'

	INCLUDE 'PRDARGCOM.INC'

C...............................................................

	rsutc = dint(rjdf*86400.d0+0.5d0)

cxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
C:LOAD OBSERVATORY COORDINATES
C  THESE PARAMETERS ARE MODIFIED BY POLAR MOTION EACH PASS
      TRANSM(1) = TRAD
      TRANSM(2) = TLONG
      TRANSM(3) = TLAT
      GDLTT     = TGDLAT
      RECEIV(1) = RRAD
      RECEIV(2) = RLONG
      RECEIV(3) = RLAT
      GDLTR     = RGDLAT

C:LOAD TARGET COORDINATES (for the center of the moon, that is)
      REFLCT(1) = 0.D0
      REFLCT(2) = 0.D0
      REFLCT(3) = 0.D0

     0	call subrecr(rjdint,rsutc,vlite,ut1c,xpole,ypole,eosrc,
     1		     ttltim,lhatrn,azaltrn,lharcv,azalrcv)
cxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

C:GET THE PROGRAM QUICKLY TO MOON RISE, TO SAVE TIME
C     COSHA, HRISE REFER TO HOUR ANGLE OF MOON AT RISE
C     RISET & SETT ARE LUNAR RISE & SET TIMES RESPECTIVELY (IN SEC OF JD)
 420  cosha= -dtan(tlat*cdr)*dtan(lhatrn(3)*cdr)
      hrise= -dabs(datan2(dsqrt(1.d0-cosha*cosha),cosha))
      t2rise= HRATE*(HRISE-LHATRN(2)*CDR*15.D0)/CTSR
c check that we're not going to the wrong moonrise (e.g. don't go
c to the previous moonrise when we want the next one.)
	if(azaltrn(3) .LT. minelv) then
c	  Moon is down; make sure we go to future moonrise.
	  if(t2rise .LT. 0.d0) then
	    t2rise = t2rise + hrate*86400.d0
	  endif
	else
c	  Moon is up; make sure we go to past moonrise.
	  if(t2rise .GT. 0.d0) then
	    t2rise = t2rise - hrate*86400.d0
	  endif
	endif

      RISET= RSUTC + t2rise
      t2set = - 2.D0*HRISE*HRATE/CTSR
      SETT= RISET + t2set
c
c	   the
	   end

