      SUBROUTINE DCTPO(IHAND,X,Y,Z,XL,XC)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      INCLUDE 'LICENSE-BSD3.inc'
C
C    CONVERTS DIRECTION COSINES TO POLAR COORDINATES. HANDEDNESS OF THE
C      POLAR COORDINATES IS CHANGED IF IHAND{0. (IHAND.LT.0?)
C    MARK A. POWELL - OCTOBER 28, 1976
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C  Tue Jul 5, 1988 11:24:39
C     1) REPLACE "COMPILER DOUBLE PRECISION" WITH
C        "IMPLICIT DOUBLE PRECISION (A-H,O-Z)".
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C  Fri Mar 31, 1989 16:12:40
C	CHANGE ATAN2 TO DATAN2, AND SQRT TO DSQRT.
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C

	character*(*) sccsid
	parameter (sccsid = "@(#)misclibnp.f	1.15\t01/06/92")

      XL=DATAN2(Y,X)
      R=DSQRT(X*X+Y*Y)
      XC=DATAN2(Z,R)
      IF (IHAND.LT.0) XL=-XL
C
C	   THE
	   END


      SUBROUTINE POTDC(IHAND,XL,XC,X,Y,Z)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C
C    CONVERTS POLAR COORDINATES TO DIRECTION COSINES. HANDEDNESS OF THE
C      COORDINATES SYSTEM IS CHANGED IF IHAND{0.  (ihand<0?)
C    MARK A. POWELL - OCTOBER 28, 1976
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C  Tue Jul 5, 1988 11:24:27
C     1) REPLACE "COMPILER DOUBLE PRECISION" WITH
C        "IMPLICIT DOUBLE PRECISION (A-H,O-Z)".
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C  Fri Mar 31, 1989 16:14:44
C	REPLACE SINE AND COSINE WITH D.P. VERSIONS.
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
      CSXC=DCOS(XC)
      X=CSXC*DCOS(XL)
      Y=CSXC*DSIN(XL)
      Z=DSIN(XC)
      IF (IHAND.LT.0) Y=-Y
C
C	   THE
	   END


	SUBROUTINE READR (IUNIT, RECM1, BUFF, NRECD, LEN, IOS)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C  Mon Aug 15, 1988 16:30:11
C     CONVERT DUMMY "BUFF" TO CHAR*1(500), USE IMPLIED DO TO
C     READ RECORD OF LENGTH "LEN" INTO BUFF.  USE "LEN" TO
C     SET UP LIMITS OF IMPLIED DO.
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C 08 Mar 91, 1707: add parameter to error for calling exit.
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

	INCLUDE 'DEVCOM.INC'

      INTEGER IUNIT,RECM1,NRECD,IOS,RECDNO
	CHARACTER*1 BUFF(1560)

	KF=1
	KL=KF+LEN-1
      DO 10 I=1,NRECD
        RECDNO=RECM1+I
        READ(UNIT=IUNIT,IOSTAT=IOS,ERR=999,REC=RECDNO)
     1		(BUFF(J),J=KF,KL)
	KF=KF+LEN
	KL=KL+LEN
10    CONTINUE
      IOS=IOS+1
	GOTO 10000

ccc 11 Apr 90 - use call to "error"
CCC ADD ERROR MESSAGE, Fri Jul 8, 1988 10:03:05  CCC
999	errmsg = " READR: ERROR IN READ: UNIT, IOSTAT ="
	call error(1,errmsg,iunit,ios)

10000	RETURN
C
C	   THE
	   END 


      SUBROUTINE WRTR ( IUNIT, RECM1, BUFF, NRECD, LEN, IOS )
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C  Mon Aug 15, 1988 16:30:11
C     CONVERT DUMMY "BUFF" TO CHAR*1(500), USE IMPLIED DO TO
C     TO READ RECORD OF LENGTH "LEN" INTO BUFF.  USE "LEN" TO
C     SET UP LIMITS OF IMPLIED DO.  (ACTUALLY USE CHAR*1(1560))
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C 08 Mar 91, 1707: add parameter to error for calling exit.
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

	INCLUDE 'DEVCOM.INC'

      INTEGER IUNIT,RECM1,NRECD,IOS,RECDNO
	CHARACTER*1 BUFF(1560)

	KF=1
	KL=KF+LEN-1
      DO 10 I=1,NRECD
        RECDNO=RECM1+I
        WRITE(UNIT=IUNIT,IOSTAT=IOS,ERR=999,REC=RECDNO)
     1		(BUFF(J),J=KF,KL)
	KF=KF+LEN
	KL=KL+LEN
10    CONTINUE
      IOS=IOS+1
	GOTO 10000

ccc 11 Apr 90 - use call to "error"
CCC ADD ERROR MESSAGE, Fri Jul 8, 1988 10:03:05  CCC
999	errmsg = " WRTR: ERROR IN WRITE: UNIT, IOSTAT ="
	call error(1,errmsg,iunit,ios)

10000	CONTINUE
C
C	   THE
	   END 

	subroutine getargs
c  This routine checks the command line for argument for the
c  normal point program (npt), and sets flags and gets file names
c  if any are found.
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C 08 Mar 91, 1707: add parameter to error for calling exit.
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C 25 Mar 91, 1834: move inquire of filnami to parin, which has the
C  pathnames.  (We might want to put pathnames on the command line.)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C 30 Jul 91,     : add summary file name input. rlr.
C		   add polar motion file input. rlr.
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C 18 May 99,     : add '-c' switch to allow creating old cospar-style
C		   mini-normalpoints rather than ilrs/cstg nps. rlr.
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C 10 Dec 99      : add -f switch and name for control file. rlr.
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

	INCLUDE 'DEVCOM.INC'

	include 'COMARG.INC'

	include 'LIOCMN.INC'

	include 'POLCOM.INC'

	character anarg*255, msg*255
	double precision temp
	logical exists

c  initialize number-of-arguments counter and existance
c  flags for various file names.
	nargs=0
	nflg=.FALSE.
	rflg=.FALSE.
	iflg=.FALSE.
	oflg=.FALSE.
	mflg=.FALSE.
	cflg=.FALSE.
	lflg=.FALSE.
	sflg=.FALSE.
	tflg=.FALSE.
	gflg=.FALSE.
        ntype= -1
	filpol= 'poldat'
	lenpol=itrimlen(filpol)
	ispol= .FALSE.
	sigmult=3.d0

c  top of loop for processing command line arguments
c__________________________________________________________
1	continue
c  get the next argument in the line and exit if there wasn't one.
	call getarg(nargs+1,anarg)
        lena= itrimlen(anarg)

c      be sure there are blanks as far as needed for formats
	anarg(lena+1:255)= '          '

CCC	if(lena .EQ. -1) goto 9 linux fix
	if(lena .LE. 0) goto 9
	if(lena .GE. 255) goto 9
c  there was and argument; increment the argument counter.
	nargs=nargs+1
c  check for what kind of option it is (must be an option flag)
c . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
	if(anarg(1:2) .EQ. '-r')then
c	  process "recalculate" option flag
	  rflg=.TRUE.
c . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
	elseif(anarg(1:2) .EQ. '-n')then
c	  process "normal point" option flag
	  nflg=.TRUE.
c . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
C  Remove COSPAR option with CRD re-write
CC	elseif(anarg(1:2) .EQ. '-c')then
c	  process "cospar normalpoint" option flag
CC	  ntype= 0
c . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
	elseif(anarg(1:2) .EQ. '-i')then
c	  process input file name
CC	  lenb=igetarg(nargs+1,anarg,255)
	  call getarg(nargs+1,anarg)
          lenb= itrimlen(anarg)
	  if(lenb.EQ.-1 .OR. anarg(1:1).EQ.'-') then
	    errmsg=' no input file name after "-i" option.'//char(7)
	    call error(2,errmsg,0,nargs)
	  endif
	  nargs=nargs+1
ccc	  inquire(file=anarg,iostat=inqios,err=101,exist=exists)
ccc101	  errmsg=' Error in inquire of file '//anarg(1:lenb)//char(7)
ccc	  if(inqios .NE. 0)call error(5,errmsg,lenb,inqios)
ccc	  errmsg=' File '//anarg(1:lenb)//' not found'//char(7)
ccc	  if(.NOT.exists)call error(2,errmsg,lenb,inqios)
	  clfilin = anarg(1:lenb)
	  lenfilin=itrimlen(clfilin)
cebug.debug.debug.debug.debug.debug.debug.debug.debug.debug.debug.
ccc	write(ioutdev,'(a)')clfilin
ccc	write(ioutdev,'(2i8)')lenb,lenfilin
ccc	read(indev,'(1x)')
cebug.debug.debug.debug.debug.debug.debug.debug.debug.debug.debug.
	  iflg = .TRUE.
c . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
	elseif(anarg(1:2) .EQ. '-o')then
c	  process output file name
CC	  lenb=igetarg(nargs+1,anarg,255)
	  call getarg(nargs+1,anarg)
          lenb= itrimlen(anarg)
	  if(lenb.EQ.-1 .OR. anarg(1:1).EQ.'-') then
	    errmsg=' no output file name after "-o" option.'//char(7)
	    call error(2,errmsg,1,nargs)
	  endif
	  nargs=nargs+1
	  clfilout = anarg(1:lenb)
	  oflg = .TRUE.
c . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
	elseif(anarg(1:2) .EQ. '-m')then
c	  process normal point file name
CC	  lenb=igetarg(nargs+1,anarg,255)
	  call getarg(nargs+1,anarg)
          lenb= itrimlen(anarg)
	  if(lenb.EQ.-1 .OR. anarg(1:1).EQ.'-') then
	    errmsg=' no mininormal point file name after "-n" option.'
     1							//char(7)
	    call error(2,errmsg,0,nargs)
	  endif
	  nargs=nargs+1
	  clfilmini = anarg(1:lenb)
	  mflg = .TRUE.
c . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
	elseif(anarg(1:2) .EQ. '-c')then
c	  process residual file name
CC	  lenb=igetarg(nargs+1,anarg,255)
	  call getarg(nargs+1,anarg)
          lenb= itrimlen(anarg)
	  if(lenb.EQ.-1 .OR. anarg(1:1).EQ.'-') then
	    errmsg=' no residual file name after "-c" option.'
     1							//char(7)
	    call error(2,errmsg,0,nargs)
	  endif
	  nargs=nargs+1
	  clfilres = anarg(1:lenb)
	  cflg = .TRUE.
c . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
	elseif(anarg(1:2) .EQ. '-l')then
c	  process listing (log) file name
CC	  lenb=igetarg(nargs+1,anarg,255)
	  call getarg(nargs+1,anarg)
          lenb= itrimlen(anarg)
	  if(lenb.EQ.-1 .OR. anarg(1:1).EQ.'-') then
	    errmsg=' no log file name after "-l" option.'//char(7)
	    call error(2,errmsg,2,nargs)
	  endif
	  nargs=nargs+1
	  clfillog = anarg(1:lenb)
	  lflg = .TRUE.
c . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
	elseif(anarg(1:2) .EQ. '-s')then
c	  process summary file name
CC	  lenb=igetarg(nargs+1,anarg,255)
	  call getarg(nargs+1,anarg)
          lenb= itrimlen(anarg)
	  if(lenb.EQ.-1 .OR. anarg(1:1).EQ.'-') then
	    errmsg=' no summary file name after "-s" option.'//char(7)
	    call error(2,errmsg,2,nargs)
	  endif
	  nargs=nargs+1
	  clfilsum = anarg(1:lenb)
	  sflg = .TRUE.
c . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
	elseif (anarg(1:2) .EQ. '-p') then
c	  process polar motion file name
CC	  lenb=igetarg(nargs+1,anarg,255)
	  call getarg(nargs+1,anarg)
          lenb= itrimlen(anarg)
	  if(lenb.EQ.-1 .OR. anarg(1:1).EQ.'-') then
	    errmsg=' no polar motion file name after "-p" option.'
     X		//char(7)
	    call error(2,errmsg,2,nargs)
	  endif
	  nargs=nargs+1
	  filpol=anarg(1:lenb)
	  lenpol=itrimlen(filpol)
	  ispol = .TRUE.
c . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
	elseif (anarg(1:2) .EQ. '-f') then
c	  process control file name
	  call getarg(nargs+1,anarg)
          lenb= itrimlen(anarg)
	  if(lenb.EQ.-1 .OR. anarg(1:1).EQ.'-') then
	    errmsg=' no control file name after "-f" option.'
     X		//char(7)
	    call error(2,errmsg,2,nargs)
	  endif
	  if(lenb.GT.14) then
	    errmsg=' control file name too long.'
     X		//char(7)
	    call error(2,errmsg,2,nargs)
	  endif
	  nargs=nargs+1
	  FNAME(7)=anarg(1:lenb)
c . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
	elseif(anarg(1:2) .EQ. '-t')then
c	  get the normalpoint time bin length in sec
C	  read(anarg(3:lena),'(F6.0)',iostat=ios) nptlen
	  read(anarg(3:8),'(F6.0)',iostat=ios) nptlen
	  if (ios.ne.0) then
	    msg = 'bad character(s) in command line arg '//anarg(1:lena)
	    call error(2,msg,narg,ios)
	  endif
	  if (nptlen .gt. 0) tflg= .TRUE.
c . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
	elseif(anarg(1:2) .EQ. '-g')then
c	  get the sigma multiplier for the least squares fit for normalpoint
C	  read(anarg(3:lena),'(F4.0)',iostat=ios) temp
	  read(anarg(3:6),'(F4.0)',iostat=ios) temp
	  if (ios.ne.0) then
	    msg = 'bad character(s) in command line arg '//anarg(1:lena)
	    call error(2,msg,narg,ios)
	  endif
	  if (temp .ge. 0.) then
	      gflg= .TRUE.
	      sigmult= temp
	  endif
c . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
	else
	  errmsg=' Bad option flag: "'//anarg(1:lena)//'"'//char(7)
	  call error(2,errmsg,lena,nargs)
	endif
c . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
	goto 1
c__________________________________________________________
c  bottom of loop.

c  exit from loop.
9	continue
			if(nargs .GT. 0) then
c  Check for certain errors in command line input.
	if(.NOT.iflg)then
	  errmsg=' Error: must have input file in command line'
     1							//char(7)
	  call error(2,errmsg,0,nargs)
	elseif(lenfilin .LE. 3) then
	  errmsg=' Error: input file name too short.'//char(7)
	  call error(2,errmsg,0,lenfilin)
	elseif(rflg .AND. nflg) then
	  errmsg=' Both -n and -r flags set'//char(7)
	  call error(2,errmsg,0,nargs)
	endif

c  Default is flag for normal point.
	if(.NOT.nflg .AND. .NOT.rflg .AND. iflg) then
	  nflg = .TRUE.
	endif

c  Create default file names if necessary.
	if(.NOT.oflg)then
	  clfilout=clfilin(1:lenfilin-2)//'lo'
	endif
	if(.NOT.lflg)then
	  clfillog=clfilin(1:lenfilin-2)//'ll'
	endif
	if(.NOT.mflg .AND. nflg)then
	  clfilmini=clfilin(1:lenfilin-2)//'lm'
	endif
	if(.NOT.cflg) then
	  clfilmini=clfilin(1:lenfilin-2)//'rsc'
	endif
			endif
c
c	   the
	   end


	function iprdargs(autoprd, ndays)
c  This function by "euler", the main program for "prd", the prediction
c  program.  If there are no command line arguments, then autoprd is
c  .FALSE. and the function value is zero.  If there are arguments, then
c  the number of them is returned as the value of the function, and
c  autoprd is .TRUE. if the first argument is "-a".  If the
c  first argument is not "-a", the program prints an error message
c  and stops.
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c 14 Jun 90, 1122: This function added to misclibnp for prd (euler).
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
C 08 Mar 91, 1707: add parameter to error for calling exit.
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
c 31 May 91, 0931: Let user specify on command line the number of days
c	from current time to do automatic predictions.  The argument
c	"ndays" is the number of days to do predictions.  If ndays is
c	zero, the predictions are done as before.  If ndays is greater
c	than zero, ndays is added to the julian date, and predictions
c	are done as long as the time of moon rise is less than this
c	new date.
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C 05 Jul 91: provide for polar motion file name on command line.
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C 10 Dec 01: Provide -f switch and control file name on command line. rlr.
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
	INCLUDE 'POLCOM.INC'

	INCLUDE 'LIOCMN.INC'

	character anarg*255, msg*80, d1*1,d2*1
	character ctrlfile*14
	logical autoprd

	continue
5001	format(i1)

	autoprd=.FALSE.
	numargs=iargc()
	iprdargs=numargs
	narg=0
	filpol = 'poldat'
	ctrlfile = 'default'
	lenpol=itrimlen(filpol)
	ispol = .FALSE.
c===================================================================
1	continue
	if(narg .GE. numargs) goto 2
	  narg=narg+1
CC	  lena=igetarg(narg,anarg,255)
	  call getarg(narg,anarg)
          lena= itrimlen(anarg)
c + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + +
	  if(anarg(1:2) .EQ. '-a') then
	    autoprd=.TRUE.
c see if number of days is specified.
	ndays = 0
c------------------------------------------------------------------
	if(lena .GT. 2) then
	  if(lena .GT. 6) then
	    msg=' command line argument too long: '//anarg(1:lena)
	    call error(2,msg,narg,lena)
	  endif
	  d1 = anarg(3:3)
	  if(llt(d1,'0') .OR. lgt(d1,'9')) then
	    msg=' Bad char in command line argument: '//d1
	    call error(2,msg,narg,lena)
	  endif
	  read(d1,5001)idig1
c	  check for 1 or 2 digits:
c..................................................................
	  if(lena .EQ. 3) then
c	    one digit
	    ndays = idig1
	  else
c..................................................................
c	    two digits
	    d2 = anarg(4:4)
	    if(llt(d2,'0') .OR. lgt(d2,'9')) then
	      msg=' Bad char in command line argument: '//d2
	      call error(2,msg,narg,lena)
	    endif
	    read(d2,5001)idig2
	    ndays = 10*idig1 + idig2
	  endif
c..................................................................
	endif
c------------------------------------------------------------------
c + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + +
	  elseif (anarg(1:2) .EQ. '-p') then
	    narg=narg+1
	    if(narg.GT.numargs) then
	      msg=' No file name follows "-p" on command line.'
	      call error(2,msg,narg,numargs)
	    endif
CC	    lena=igetarg(narg,anarg,255)
	    call getarg(narg,anarg)
            lena= itrimlen(anarg)
	    if(anarg(1:2).EQ.'-a' .OR. anarg(1:2).EQ.'-f') then
	      msg=' No file name follows "-p" on command line.'
	      call error(2,msg,narg,numargs)
	    endif
	    filpol=anarg
	    lenpol=itrimlen(filpol)
	    ispol = .TRUE.
c + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + +
c + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + +
CC Prediction control file name is overridden...
	  elseif (anarg(1:2) .EQ. '-f') then
	    narg=narg+1
	    if(narg.GT.numargs) then
	      msg=' No file name follows "-f" on command line.'
	      call error(2,msg,narg,numargs)
	    endif
	    call getarg(narg,anarg)
            lena= itrimlen(anarg)
	    if(anarg(1:2).EQ.'-a' .OR. anarg(1:2).EQ.'-p') then
	      msg=' No file name follows "-p" on command line.'
	      call error(2,msg,narg,numargs)
	    endif
	    if(lenb.GT.14) then
	      msg=' Control file name too long.'
     X		//char(7)
	      call error(2,msg,2,nargs)
	    endif
	    FNAME(6)=anarg
c + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + +
	  else
	    msg=' Bad command line argument: "'//anarg(1:lena)//'"'
	    call error(2,msg,narg,lena)
	  endif
c + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + +
	goto 1
2	continue
c===================================================================
cccccc===================================================================
ccccc	if(numargs .GT. 0) then
ccccc	  narg=narg+1
ccccc	  lena=igetarg(narg,anarg,255)
ccccc	  if(anarg(1:2) .EQ. '-a') then
ccccc	    autoprd=.TRUE.
ccccc	  else
ccccc	    msg=' Bad command line argument: "'//anarg(1:lena)//'"'
ccccc	    call error(2,msg,narg,lena)
ccccc	  endif
cccccc
cccccc see if number of days is specified.
ccccc	ndays = 0
cccccc------------------------------------------------------------------
ccccc	if(lena .GT. 2) then
ccccc	  if(lena .GT. 4) then
ccccc	    msg=' command line argument too long: '//anarg(1:lena)
ccccc	    call error(2,msg,narg,lena)
ccccc	  endif
ccccc	  d1 = anarg(3:3)
ccccc	  if(llt(d1,'0') .OR. lgt(d1,'9')) then
ccccc	    msg=' Bad char in command line argument: '//d1
ccccc	    call error(2,msg,narg,lena)
ccccc	  endif
ccccc	  read(d1,5001)idig1
cccccc	  check for 1 or 2 digits:
cccccc..................................................................
ccccc	  if(lena .EQ. 3) then
cccccc	    one digit
ccccc	    ndays = idig1
ccccc	  else
c..................................................................
c	    two digits
ccccc	    d2 = anarg(4:4)
ccccc	    if(llt(d2,'0') .OR. lgt(d2,'9')) then
ccccc	      msg=' Bad char in command line argument: '//d2
ccccc	      call error(2,msg,narg,lena)
ccccc	    endif
ccccc	    read(d2,5001)idig2
ccccc	    ndays = 10*idig1 + idig2
ccccc	  endif
c..................................................................
ccccc	endif
c------------------------------------------------------------------
ccccc	endif
cccccc===================================================================
c
c	   the
	   end


	subroutine findext(filnam,lenth,lenpre)
c  This subroutine finds the extension of the file name "filnam".
c  The subroutine accomplishes this by returning the length of the
c  prefix (up to the last ".").  If filnam doesn't contain a ".",
c  the length returned is -1.  A prefix length of zero means the file
c  name begins with a "." and has no other "." in it.  The entire
c  length of filnam is returned as "lenth".
	character*(*) filnam

c  get the length of filnam up to the first blank or null.
	lenth=itrimlen(filnam)
c - - - - - - - - - - - - - - - - - - - - - - - - - - -
	if(lenth .GT. 0) then
	  lenpre=-1
	  do 1 i1=lenth,1,-1
	  if (filnam(i1:i1) .EQ. '.') then
	    lenpre = i1-1
	    goto 2
	  endif
1	  continue
2	  continue
	else
c - - - - - - - - - - - - - - - - - - - - - - - - - - -
c  here the length of the file name is zero (or less!).
	  lenpre=-9
	endif
c - - - - - - - - - - - - - - - - - - - - - - - - - - -
c
c	   the
	   end


	subroutine gotoeof(idevic,nrec)
c  This routine sets the file pointer for unit "idevic" after
c  the end of the last record.  The file connected to
c  idevic will then be ready for appending.
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C 08 Mar 91, 1707: add parameter to error for calling exit.
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

	INCLUDE 'DEVCOM.INC'

	continue
5001	format(1x)
	nrec=0
1	continue
	read(idevic,5001,err=9,end=2,iostat=ios)
	nrec=nrec+1
	goto 1
2	continue
	goto 10
9	continue
	errmsg=' error when trying to set device for appending'
	call error(1,errmsg,idevic,ios)
10	continue
c
c	  the
	  end

