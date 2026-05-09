	subroutine parin(vlite,obcdsv,vers,modeop,filnami,filtmp)
	implicit double precision (a-h,o-z)
      INCLUDE 'LICENSE-BSD3.inc'
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C  Mon Aug 29, 1988 10:45:30
C     1) REPLACE "COMPILER DOUBLE PRECISION" WITH
C        "IMPLICIT DOUBLE PRECISION (A-H,O-Z)".
C     2) MAKE SURE STATEMENTS DON'T EXTEND BEYOND COLUMN 72.
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C  Wed Oct 19, 1988 14:32:48
C     1) REMOVE IN-LINE COMMENTS.
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C  Thu Oct 27, 1988 09:20:59
C     1) MOVE DATA STATEMENTS OF VARIABLES IN COMMON TO BLKNP.FOR
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C  Mon Oct 31, 1988 15:59:07
C     1) CHANGE TYPE REAL TO TYPE DOUBLE PRECISION.
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
ccc	Wed Dec  6 14:19:46 CST 1989 (start)  -  tlr
c  provide for command line input
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c  Fri Dec 15 15:52:04 CST 1989, use pathnames in opens.
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cc 20 Apr 90 - tlr
c  take control file information from file type 7: NPTCTRL.PC
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c 25 Mar 91, 1736:
c   o  use pathnames in all open statements.
c   o  put inquire of filnami here instead of in "getargs".
c   o  add debugging aids.
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c 01 May 91, 1756: in command line case:
c   o  use length of output file name in inquire statement;
c   o  open filtmp in patheph directory;
c ?9 May 91? - accommodate file names beginning with "/".
c 03 Dec 91  - differentiate header for normalpoint and recalc. rlr.
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

C  READ CONTROL PARAMETERS,
C       MODEL PARAMETERS, AND
C       REFLECTOR COORDINATES
C  FROM DISK;
C  TRANSFER PARAMETERS FROM INPUT TO TARGET COMMONS;
C  GET I/P & O/P FILE NAMES FROM OPERATOR AND OPEN FILES
C
C  OUTPUT:
C	MODEOP - 0 FOR NORMAL POINTING RUN
C		 1 FOR RESIDUAL RECALULATION ("NEWRES")
C
C  REVISIONS:
C	12-09-86 - NEW LUNAR '86 DATA FORMAT. 
C		   COMBINE FUNCTION OF NORMAL PT & NEWRES INTO ONE PGM.
C		   OUTPUT IS WRITTEN TO A FILE RENAMED TO THE OLD FILE NAME.
C		   OUTPUT LISTING FILE IS CREATED. V10. RLR.
C	07-20-87 - READ PROPER REFLECTOR COORDINATE SET FROM FILE LREFL1.PC
C		   V 11. RLR.
c
c	30 Oct 89 - Use lunar prediction files (replace call to
c		    "datio(1,5,8)" with "datio(1,6,1)".)
c		    Provide for variable jdstep.
c		    (this is before converting to pparin)
c							  tlr
C
C$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$ RLR 10/85

	character*(*) sccsid
	parameter (sccsid = "@(#)parin11.f	1.12\t01/06/92")

C********************************************

C   NO REFERENCES TO VARIABLES IN /BLD/?
	INCLUDE 'BLD.INC'

	INCLUDE 'CETB1.INC'

	INCLUDE 'CNTINF.INC'

	include 'COMARG.INC'

	INCLUDE 'DEFCMN.INC'

	INCLUDE 'DEVCOM.INC'

	INCLUDE 'MODINF.INC'

	INCLUDE 'OBSINF.INC'

	INCLUDE 'PASSIT.INC'

	INCLUDE 'PDCMN.INC'

	INCLUDE 'U1UC.INC'

	INCLUDE 'USNOC.INC'

C********************************************

	character*9 vers
	character*1 ichr
      character*255 filnami,filnamo,filtmp
      character*255 pathop
c debug debug debug debug debug debug debug debug debug debug debug
      character*255 pathtmp, filetmp*28
c debug debug debug debug debug debug debug debug debug debug debug
      double precision obcdsv(10)
      double precision pcs(14)
	logical exists

C...FOR KUSTNER TO SAVE INFO
	common/kuscmn/ frstim,tt(4),jd1,jd2,jdstep
	logical frstim
	double precision jd1,jd2,jdstep

      equivalence (pcs(1),k1)
      
C----------------------------------------------------------
C:READ CONTROL, MODEL, & REFLECTOR FILES
C	CONTROL FILE(MIT)
cc Fri Apr 20 13:41:28 CDT 1990 - tlr
c  use file type 7 (NPTCTRL.PC) for recalc/normal point.
	call datio(1,7,1)
CC	write(*,*) "datio ctrl: nmodel=",nmodel
CC	write(*,*) CNTNAM,OBSPIK,RFLPIK,NMODEL,RJD1,RJD2,
CC     1          STEPJD,MINELV,DEBUG,RSETPK,IEPH,FILEPH
CC	write(*,1111) nmodel, rsetpk
CC 1111	format("octal: ",o8,o8)
c.debug.debug.debug.debug.debug.debug.debug.debug.debug.debug.debug.
	idbx = 0
	if(idbx .NE. 0) then
c  With new FORTRAN, xdb can't see names in common; use intermediate
c  variables to examine/change common variables.
	  filetmp = fileph
	  fileph  = filetmp
	  lentmp  = itrimlen(fileph)
	  leneph  = lentmp
	  pathtmp = patheph
	  patheph = pathtmp
	  lentmp  = itrimlen(patheph)
	  leneph  = lentmp
	endif
c.debug.debug.debug.debug.debug.debug.debug.debug.debug.debug.debug.

C	PARAMETER FILE
      call datio(1,4,nmodel)

C	REFLECTOR COORDINATES
      call datio(1,2,rsetpk)

C:INITIALIZE EPHEMERIS & PARAMETERS
C  COPY FROM MODEL INFO COMMON INTO LOCAL COMMONS
      emrat = emrati
      k2love = k2lovi
      uttide = ut1tid
      do 10 ii=1,14
10    pcs(ii) = cs(ii)

      vlite=vlight*86400.d0/aukm
      k2mc3=(kk+kk)/(vlite*vlite*vlite)

C  RECORD WHICH COORDINATE SETS ARE AVAILABLE
      do 50 nobs=1,10
        obcdsv(nobs)= 0.
        if (obspik(nobs) .neqv. .true.) go to 50
        call datio(1,1,nobs)
        obcdsv(nobs)= obscode
   50 continue

c======== case of interactive or command line input. ============
			if(nargs .NE. 0) then
C..set mode of op
	if(nflg)then
	  modeop= 0
	  filnami=clfilin
	  filtmp=clfilout
	  filnamo=clfilmini
	else
	  modeop= 1
	  filnami=clfilin
	  filtmp=clfilout
	endif

C..GET INPUT & OUTPUT FILE
	  lenami = itrimlen(filnami)
	if(filnami(1:1) .EQ. '/') then
	  pathop = filnami
	else
	  pathop = patheph(1:leneph)//filnami
	endif
	  inquire(file=pathop,iostat=inqios,err=101,exist=exists)
101	  errmsg=' Error in inquire of file '//pathop//char(7)
	  if(inqios .NE. 0)call error(1,errmsg,lenb,inqios)
	  errmsg=' File '//pathop//' not found'//char(7)
	  if(.NOT.exists)call error(2,errmsg,lenb,inqios)

	open(21,file=pathop,iostat=ierr)

	if(filtmp(1:1) .EQ. '/') then
	  pathop = filtmp
	else
	  pathop = patheph(1:leneph)//filtmp
	endif
	open(22,file=pathop, iostat=ierr)

C  GET & OPEN THE MINI-NORMALPOINT FILE (if this is not recalc)
	if (modeop .EQ. 0) then
	  if(filnamo(1:1) .EQ. '/') then
	    pathop = filnamo
	  else
	    pathop = patheph(1:leneph)//filnamo
	  endif
	  open(23,file=pathop,iostat=ierr)
	endif
c================================================================
			else
C..SEND GREETING TO CRT
      if (ioutdev .gt. 0) write(ioutdev,1000) vers

C..GET MODE OF OP
      modeop= 1
	ichr=' '
CC      type *, 'Create normal points? (Y OR N) <Y>:'
      write(*,*) 'Create normal points? (Y OR N) <Y>:'
      read(indev,1026) ichr
      if (ichr.EQ.'Y' .OR. ichr.EQ.'y' .OR. ichr.EQ.' ') modeop= 0

C..GET INPUT & OUTPUT FILE NAMES
CC      type *, "Input/Output lunar data file: "
      write( *,* ) "Input/Output lunar data file: "
      read(indev,1025) filnami
	open(21,file=patheph(1:leneph)//filnami,iostat=ierr,
     1				status='old')

C  OPEN THE TEMPORARY OUTPUT FILE.
	filtmp = 'npout'
	open(22,file=patheph(1:leneph)//filtmp,iostat=ierr)

C  GET & OPEN THE MINI-NORMALPOINT FILE
	if (modeop .EQ. 0) then
CC	  type *, "Output mini normal point file:"
CC	  write( *,* ) "Output mini normal point file:"
	  write( *,* ) "Output CRD normal point file:"
	  read(indev,1025) filnamo
1025	  format(a)
1026	  format(a1)
	  open(23,file=patheph(1:leneph)//filnamo,iostat=ierr)
	endif
			endif
c================================================================

C..WRITE OUT GREETINGS TO LISTING FILE
	if (modeop.eq.0) then
	    write(nplogfil,1000) vers
 1000       format(10x,
     x		'Lunar Normalpoint / Residual Calculation Program',
     x		5x,a9//)
	else
	    write(nplogfil,1001) vers
 1001       format(10x,
     x		'Lunar Residual Calculation / Normalpoint Program',
     x		5x,a9//)
	endif
      	if (ioutdev .gt. 0) 
     x    write(ioutdev,1010) filnami(1:itrimlen(filnami))
        write(nplogfil,1010) filnami(1:itrimlen(filnami))
 1010   format(1x,'Input/Output Lunar Data File:',/,a,/)
	if (modeop.eq.0) then
            if (ioutdev .gt. 0)
     x	    write(ioutdev,1011) filnamo(1:itrimlen(filnamo))
	    write(nplogfil,1011) filnamo(1:itrimlen(filnamo))
	endif
 1011   format(1x,'Output Mini-Normalpoint File: ',/,a,/)
C
c	   the
	   end

