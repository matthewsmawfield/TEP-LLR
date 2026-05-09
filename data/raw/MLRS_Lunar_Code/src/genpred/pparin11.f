	subroutine pparin(vlite)
ccc	subroutine pparin(vlite,obcdsv,vers,filnami,filtmp)
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
c  7 Nov 89 - tlr
c	convert to NOVA-style pparin.
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
c 06 Nov 90, 1600: add "ichngdeb" and code around 127 to make it
c	easier to set the debug vector from xdb.
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

C  READ CONTROL PARAMETERS,
C       MODEL PARAMETES, AND
C       REFLECTOR COORDINATES
C  FROM DISK;
C  TRANSFER PARAMETERS FROM INPUT TO TARGET COMMONS;
C  GET I/P & O/P FILE NAMES FROM OPERATOR AND OPEN FILES
C
C  OUTPUT:
C	MODEOP:  0 FOR NORMAL POINTING RUN
C		 1 FOR RESIDUAL RECALULATION ("NEWRES")
c		-1 for prediction program
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
	parameter (sccsid = "@(#)pparin11.f	1.4\t02/21/91")

C********************************************

C   NO REFERENCES TO VARIABLES IN /BLD/?
	INCLUDE 'BLD.INC'

	INCLUDE 'CNTINF.INC'

	INCLUDE 'CETB1.INC'

	INCLUDE 'DEFCMN.INC'

	INCLUDE 'DEVCOM.INC'

	INCLUDE 'EPHREC.INC'

	INCLUDE 'MODEOP.INC'

	INCLUDE 'MODINF.INC'

	INCLUDE 'OBSINF.INC'

	INCLUDE 'PASSIT.INC'

	INCLUDE 'PDCMN.INC'

	INCLUDE 'U1UC.INC'

	INCLUDE 'USNOC.INC'

C********************************************

ccc	character*1 ichr
ccc	character*56 locfil
ccc	character*28 filnami,filnamo,filnaml,filtmp
      double precision pcs(14),cns(18),cut(5)

      common/cnnsnn/c30,c31,c32,c33,c40,c41,c42,c43,c44,
     1              s30,s31,s32,s33,s40,s41,s42,s43,s44

C...FOR KUSTNER TO SAVE INFO
	common/kuscmn/ frstim,tt(4),jd1,jd2,jdstep
	logical frstim
	double precision jd1,jd2,jdstep

      equivalence (pcs(1),c1), (cns(1),c30), (cut(1),utdot)
      
C----------------------------------------------------------
C:READ CONTROL, MODEL, & REFLECTOR FILES
C	CONTROL FILE(MIT)
ccc      call datio(1,5,8)
	call datio(1,6,1)

C  PARAMETER FILE
c  read into /modinf/ starting at "nammod"
      call datio(1,4,nmodel)

C  REFLECTOR COORDINATES
c  read into /rflinf/ starting at "rfset"
      call datio(1,2,rsetpk)

C:INITIALIZE EPHEMERIS & PARAMETERS
C  COPY FROM MODEL INFO COMMON INTO LOCAL COMMONS
      emrat = emrati
      k2love = k2lovi
      uttide = ut1tid
	ephsrc = ieph

      do 10 ii=1,14
10    pcs(ii) = cs(ii)

cdbdbdbdbdbdbdbdbdbdbdbdbdbdbdbdbdbdbdbdbdbdbdbdbdbdbdbdbdbdbdbdb
	ichngdeb = 0
	if(ichngdeb .NE. 0)then
	  do 127 i127=1,20
127	  debug(i127)=.TRUE.
	endif
cdbdbdbdbdbdbdbdbdbdbdbdbdbdbdbdbdbdbdbdbdbdbdbdbdbdbdbdbdbdbdbdb
	do 11 i11=1,20
11	deb(i11)=debug(i11)

      vlite=vlight*86400.d0/aukm
      k2mc3=(kk+kk)/(vlite*vlite*vlite)

C  RECORD WHICH COORDINATE SETS ARE AVAILABLE
	do 50 nobs=1,10
	  if (obspik(nobs)) then
            call datio(1,1,nobs)
	  endif
50	continue

c  set up receiver flag
	if(trad.NE.rrad .OR. tlat.NE.rlat .OR. tlong.NE.rlong) irecvr=2

cccccccccccc normal point code below cccccccccccccccccccccccccccccccccyy
ccc possible future task: have npt call pparin, not parin. cccccccccccyy
cC..GET MODE OF OP
c      modeop= 1
c	ichr=' '
c      type *, 'Create normal points? (Y OR N) <Y>:'
c      read(indev,1026) ichr
c      if (ichr.EQ.'Y' .OR. ichr.EQ.'y' .OR. ichr.EQ.' ') modeop= 0
c
cC..GET INPUT & OUTPUT FILE NAMES
c      type *, "Input/Output lunar data file: "
c      read(indev,1025) filnami
c	open(21,file=filnami,iostat=ierr)
c
cC  OPEN THE TEMPORARY OUTPUT FILE.
c      do 20 i=1,56
c        k = 56 - i + 1
c        if (locfil(k:k).EQ.':') go to 25
c 20   continue
c      k= 0
c25	continue
c	filtmp = 'npout'
c	open(22,file=filtmp,iostat=ierr)
c
cC  GET & OPEN THE MINI-NORMALPOINT FILE
c      if (modeop.ne.0) go to 30
c      type *, "Output mini normal point file:"
c      read(indev,1025) filnamo
c 1025 format(a28)
c 1026 format(a1)
c	open(23,file=filnamo,iostat=ierr)
c
cCCCC..GET LISTING FILE NAME.
c30	continue
c
cC..WRITE OUT GREETINGS TO LISTING FILE
c	write(nplogfil,1000) vers
c 1000 format(10x,'Lunar Normalpoint / Residual Calculation Program',
c     x		5x,a9//)
c      write(ioutdev,1010) filnami
c      write(nplogfil,1010) filnami
c 1010 format(1x,'Input/Output Lunar Data File: ',a28/)
c	if (modeop.eq.0) then
c	  write(ioutdev,1011) filnamo
c	  write(nplogfil,1011) filnamo
c	endif
c 1011 format(1x,'Output Mini-Normalpoint File: ',a28/)
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccyy

c	   the
	   end

