	subroutine kustner(day,ex,wy,dif,iflg)
	implicit double precision (a-h,o-z)
      INCLUDE 'LICENSE-BSD3.inc'
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C  Thu Oct 20, 1988 15:49:27
C     1) REPLACE "COMPILER DOUBLE PRECISION" WITH
C        "IMPLICIT DOUBLE PRECISION (A-H,O-Z)".
C     2) CHANGE OPTIONAL COMPILE STATEMENTS TO COMMENTS.
C     3) REMOVE IN-LINE COMMENTS.
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C  Fri Oct 28, 1988 11:07:20
C     1) MOVE DATA STATEMENTS TO "BLKNP.FOR".
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C  Mon Oct 31, 1988 15:06:59
C     1) CHANGE TYPE REAL TO TYPE DOUBLE PRECISION.
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C  Wed Mar 29, 1989 09:36:57
C     1) KLUDGE D.P. RELATIONAL EXPRESSIONS.
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
c  Fri Dec 15 15:37:09 CST 1989
C     1) use path name in open statement.       tlr
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
c  30 Jan 1990
C     1) change "1000 format(2f10.0,i3)" to
c		"1000 format(2d10.1,d4.1)".       tlr
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C 06 Feb 91, 0742: set up for changing 'poldat' file name
C  when using debugger.  Need to make command line option later.
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C 03 Jul 91, 1632: use polar motion file name from /POLCOM/
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C 29 June, 1993: fix adjustment for leap sec, when date is prior to leap sec.
C		 and when the time = time of leap.
C		 rlr.
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C 28 June, 1994 - Change compare for leap sec from GE to GT. rlr
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C  ROUTINE KUSTNER READS AND INTERPOLATES POLAR MOTION AND
C  EARTH ROTATION DATA FOR THE CALLING PROGRAM. 
C  INPUT:
C	DAY 	- JULIAN DATE
C	POLDAT	- POLAR MOTION FILE
C		D	- JULIAN DATE
C		X	- X AXIS P.M. (ARCSEC)
C		Y	- Y AXIS P.M. (ARCSEC)
C		U	- EARTH ROTATION (TIME SEC)
C  OUTPUT:
C	EX	- INTERPOLATED X POLAR MOTION
C	WY	- 	"      Y   "      "
C	DIF	-	"      EARTH ROTATION
C	IFLG	- 0 => OK, 1 => NO INFO
C
C  REVISIONS:
C	11-09-87 - DON'T ASSUME THAT THE CHANGE IN TIME AT A BREAK POINT
C		   SPECIFIED IN /ATMAUX/ IS AN INTEGRAL SECOND.  THIS FIXES
C		   A PROBLEM HANDLING DATA FOR TIMES AROUND 2441317.5 JD
C		   (AND EARLIER JUMPS). 
C		 - ALSO, READ 1 MORE SIGNIFICANT DIGIT IN EACH FIELD USING
C		   BLANK FIELD BETWEEN VALUES.
C		   V01. RLR.
C
c The changes below were copied from the NOVA version; kustner02.
c	25 Oct 89 - Accept fixed time interval (jdstep) in poldat.
c		  - Add one more significant digit to 3 er/pm
c		    values in the table.
c
C$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$ RLR

	character*(*) sccsid
	parameter (sccsid = "@(#)kustner01.f	1.7\t07/01/94")

	INCLUDE 'ATMAUX.INC'

	INCLUDE 'BLD.INC'

	INCLUDE 'DEVCOM.INC'

	INCLUDE 'POLCOM.INC'

	common /kuscmn/ frstim, tt(4), jd1,jd2, jdstep
	logical frstim
	double precision jd1,jd2, jdstep

	double precision adu(4),coef(4)

	iflg= 1
	ex= 0.d0
	wy= 0.d0
	dif= 0.d0

	if (.not.frstim) goto 100
	  frstim= .false.
	  n= 1
	  if(filpol(1:1) .EQ. '/')then
	   open(4,file=filpol(1:lenpol),iostat=ierr)
	  else
	   open(4,file=pathfil(1:lenfil)//filpol(1:lenpol),iostat=ierr)
	  endif
c		no file
	  if (ierr.NE.0) return
	  read(4,1000,err=999,end=999,iostat=ios) jd1,jd2, jdstep
1000	  format(2d10.1,d4.1)
	if(jdstep .EQ. 0.d0) jdstep = 5.d0
	do 10 i10=1,4
10	tt(i10)=(i10-1)*jdstep
	  close(4,iostat=ierr)
c		wrong time

100	if (day.LT.jd1+jdstep .OR. day.GT.jd2-jdstep) return
cccc	write(npsumfil,*) day,n,dd(1),dd(n),jdstep


	if (day.GE.dd(1)+jdstep .AND. day.LE.dd(n)-jdstep) go to 200

c  read the file (again?)
	if(filpol(1:1) .EQ. '/')then
	 open(4,file=filpol(1:lenpol),iostat=ierr)
	else
	 open(4,file=pathfil(1:lenfil)//filpol(1:lenpol),iostat=ierr)
	endif
	read(4,1005)
 1005	format(1x)
	n= 1
 150	read(4,1010,err=999,end=190) dd(n),x(n),y(n),u(n)
 1010	format(f10.0,2f8.0,f9.0)
c		limit record was wrong

	if (n.EQ.1 .AND. day.LT.dd(n)) go to 999

c		too early

	if (dd(n).LT.(day-jdstep)) go to 150

	n= n+ 1
c		copy in another
	if (n.LE.10) go to 150
 190	n= n- 1

c		too few to interp.
	if (n.LT.4) go to 999

c  now set up for interp
 200 	nm1= n- 1
	do 210 i=1,nm1

	  if (day.LT.dd(i) .OR. day.GE.dd(i+1)) go to 210

c		set lower interp limit
	  min= max0(1,i-1)
	  min= min0(min,7)
	  go to 220
 210	continue

c  adjust for leap sec
 220	max= min+ 3
	do 221 i=min,max
 221	adu(i-min+1)= u(i)
	do 225 j=1,ltim,3
cccc	write(npsumfil,*) j, time(j),dd(min),dd(max)
C	if (time(j).LE.dd(min) .OR. time(j).GE.dd(max)) go to 225
	if (time(j).LE.dd(min) .OR. time(j).GT.dd(max)) go to 225
	    call cutabl(time(j)- 1.d-6, datutb)
	    call cutabl(time(j)+ 1.d-6, datuta)
	    timadj= datuta- datutb
cccc	    write(npsumfil,*) day,time(j),ltim,datuta,datutb,timadj
	    if (time(j).GT.day) go to 222

c  early point had jump
	    adu(1)= adu(1)+ timadj
	    if (dd(min+1).LT.time(j)) adu(2)= adu(2)+ timadj
	    if (dd(min+2).LT.time(j)) adu(3)= adu(3)+ timadj
c		only 1 jump possible

cccc	    write(npsumfil,*) "earlier",adu
	    go to 250

c  late point had jump
 222	    adu(4)= adu(4)- timadj
	    if (dd(min+2).GE.time(j)) adu(3)= adu(3)- timadj
	    if (dd(min+1).GE.time(j)) adu(2)= adu(2)- timadj

c		only 1 jump possible
cccc	    write(npsumfil,*) "later",adu
	    go to 250
 225	continue

 250	np= 4
	call entset(tt,day-dd(min),np,coef)

	ex= entrp(coef,x(min),np)
	wy= entrp(coef,y(min),np)
	dif= entrp(coef,adu,np)
	iflg= 0
cccc    write(npsumfil,*) "end",ex,wy,dif,iflg
 999	close(4,iostat = ierr)
c
c	   the
	   end

