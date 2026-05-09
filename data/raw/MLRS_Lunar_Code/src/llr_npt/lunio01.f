      SUBROUTINE LUNIO(MODE,INFTYP,NRECD,IRECD)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      INCLUDE 'LICENSE-BSD3.inc'
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C  Tue May 17, 1988 15:33:45
C     1) REPLACE "COMPILER DOUBLE PRECISION" WITH
C        "IMPLICIT DOUBLE PRECISION (A-H,O-Z)".
C     2) MAKE SURE STATEMENTS DON'T EXTEND BEYOND COLUMN 72.
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C  Wed May 18, 1988 14:32:05
C     1) MAKE 'FNAME(6)' CHARACTER*14 INSTEAD OF 'INTEGER FNAME(7,6)'
C     2) REMOVE DATA STATEMENTS OF COMMON ITEMS TO BLOCK DATA
C        SUBPROGRAM YCLEPT "BLKGPRED".
C     3) IN OPEN, USE FNAME(INFTYP)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C  Wed Aug 3, 1988 13:29:32
C     MAKE "IRECD" A CHARACTER VARIABLE.
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C  Mon Aug 15, 1988 16:26:30
C     ADD RECORD LENGTH TO CALLS TO "READR" AND "WRTR".
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
c  Tue Dec 19 07:54:32 CST 1989; use pathname.   tlr
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
cc 20 Apr 90 - tlr
c  change comments to indicate a 7th file type
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
cc 23 Apr 90 - tlr:  For file types 6 & 7, only read from or
c  write to the first record.
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C 26 Nov 90, 0840 - add variable idbx to facilitate debugging.
C    Feb 91 - further changes to debugging code.
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C 08 Mar 91, 1713: add parameter to error for calling exit.
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C.....
C     LUNIO HANDLES I/O FOR GENERALIZED LUNAR PARAMETER FILES
C     INPUT: MODE    1 = READ FROM DESIGNATED FILE
C                    2 = WRITE TO DESIGNATED FILE
C            INFTYP  INFORMATION TYPE; EACH IMPLIES A SPECIFIC FILE NAME
C                    AND RECORD LENGTH
C                    1 = SITE INFO
C                    2 = REFLECTOR INFO
C                    3 = LUNAR SURFACE INFO
C                    4 = LUNAR MODEL PARAMETER INFO
C                    5 = LUNAR PREDICTION CONTROL PARAMETER INFO
C                    6 = PREDICTION CONTROL SET DESTINED FOR "prd"
C                    7 = NORMAL POINT/RECALC CONTROL SET FOR "npt"
C            NRECD   RECORD NUMBER 0 TO N (fortran records 1 to (n+1))
C                    0 WILL PRESUMABLY CONTAIN HEADER INFO AND
C                    1 TO N WILL EACH CONTAIN A COMPLETE SET OF INFO
C                    AN ARRAY TO BE WRITTEN
C     OUTPUT:IRECD   AN ARRAY TO BE READ FROM FILE
C
C  REVISIONS:
C	07-20-87 - EXPAND SIZE OF LREFL.PC & GIVE IT A NEW NAME FOR
C		   CHANGE TO ALLOW MULTIPLE REFLECTOR COORDINATE SETS. V01. RLR.
C.....

	character sccsid*(*)
	parameter (sccsid = "@(#)lunio01.f	1.8\t06/18/91")

	INCLUDE 'LIOCMN.INC'

	INCLUDE 'DEVCOM.INC'

	CHARACTER IRECD*(*)
	
	CHARACTER*18 pathtmp, pathtmps, efnam*14
	save pathtmps
	
C...THIS IS WHAT IS ACTUALLY NOW IN BLOCK DATA SUBPROGRAM.
cc 20 Apr 90
C      DATA RCDLEN/200,152,130,500,252, 252,252/
C      DATA FNAME/"LSITES.PC    ","LREFL1.PC    ",
C     $           "LFEAT.PC     ","LMODEL.PC    ",
C     $           "LCONTOL.PC   ","PRDCTRL.PC   ",
C     $                           "NPTCTRL.PC   "/
C      DATA OFFSET/0,1,11,0,0, -1,-1/
c	(The last two offsets don't matter.)
C.....

C  GET THE RECORD LENGTH FOR THIS TYPE OF RECORD
	LEN = RCDLEN(INFTYP)
C.....
C.....HEADER, ALTHOUGH RECD 0, CAN BE LONGER THAN ONE RECD. HANDLE OFFSETS
      IF (NRECD.GT.0) GO TO 10
C  HEADER (ONE OR MORE PHYSICAL RECORDS PER HEADER)
      MRECD = NRECD
      NOREC = OFFSET(INFTYP)+1
	write(*,*) "no 1: ",inftyp,offset(inftyp),mrecd,nrecd,norec
      GO TO 11
C  DATA RECORDS (ONE PHYSICAL RECORD PER DATA RECORD)
   10 MRECD = NRECD+OFFSET(INFTYP)
      NOREC = 1
	write(*,*) "no 2: ",inftyp,offset(inftyp),mrecd,nrecd,norec
c for control files, force data into first record.
11	if(inftyp.EQ.6 .OR. inftyp.EQ.7)then
	  mrecd = 0
	  norec = 1
	write(*,*) "no 3: ",inftyp,offset(inftyp),mrecd,norec
	endif
C.....
C.....OPEN THE APPROPRIATE FILE
C.....(OPEN STATMENT CHANGED TO FORTRAN77 (JULY 88?))

c debug debug debug debug debug debug debug debug debug debug
c set the value of idbx to 1 with the debugger to use alternate
c prediction or normal point control file.
		idbx = 0
		if(idbx .NE. 0) then
		  pathtmp = pathfil
		  pathtmps= pathfil
		  if(idbx .EQ. 1)pathtmp = '/users/tlr/np.dr/'
		  pathfil = pathtmp
		  lenfil  = itrimlen(pathfil)
		  efnam=fname(inftyp)
		  if(idbx .EQ. 2)efnam='TLRNPTCTRL.PC'
		  if(idbx .EQ. 3)efnam='jTLRNPTCTRL.PC'
		  fname(inftyp)=efnam
		endif
c debug debug debug debug debug debug debug debug debug debug

	IERR=0
     0	errmsg=' ERROR IN SUBROUTINE LUNIO OPENING FILE "'//
     1		FNAME(INFTYP)//'"'
     0	OPEN(13,IOSTAT=IOS,ERR=9911,
     1		FILE=pathfil(1:lenfil)//FNAME(INFTYP),STATUS='OLD',
     2		ACCESS='DIRECT',FORM='UNFORMATTED',RECL=LEN)
c debug debug debug debug debug debug debug debug debug debug
		if(idbx .NE. 0) then
		  pathfil = pathtmps
		  lenfil  = itrimlen(pathfil)
		  idbx = 0
		endif
c debug debug debug debug debug debug debug debug debug debug
	IERR=1
9911  IF (IERR.NE.1) CALL ERROR(1,errmsg,IERR,IOS)
C.....
C.....READ
      IF (MODE.NE.1) GO TO 20
	errmsg=' ERROR IN SUBROUTINE LUNIO; "READR" ERROR NUMBER BELOW.'
	write(*,*) "ierr= ",ierr
      CALL READR(13,MRECD,IRECD,NOREC,LEN,IERR)
	write(*,*) MRECD, IRECD, NOREC, LEN, IERR,pathfil(1:lenfil)//FNAME(INFTYP)
c  (readr should have called error)
      IF (IERR.NE.1) CALL ERROR (1,errmsg,IERR,IERR)
      GO TO 50
C.....
C.....WRITE
   20 IF (MODE.NE.2) GO TO 50
	errmsg=' ERROR IN SUBROUTINE LUNIO; "WRTR" ERROR NUMBER BELOW.'
      CALL WRTR(13,MRECD,IRECD,NOREC,LEN,IERR)
      IF (IERR.NE.1) CALL ERROR (1,errmsg,IERR,IERR)
C.....
C.....CLOSE
50	IERR=0
	errmsg=' ERROR IN SUBROUTINE LUNIO CLOSING FILE "'//
     1		FNAME(INFTYP)//'"'
	CLOSE (13)
	IERR=1
9950	CONTINUE
c  (this is supposed to be error in close)
	IF (IERR.NE.1) CALL ERROR (1,errmsg,IERR,IOS)
      RETURN

C	   THE
	   END

