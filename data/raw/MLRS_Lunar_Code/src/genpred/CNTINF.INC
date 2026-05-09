C  INCLUDE FILE MADE Tue Jul 12, 1988 15:15:35 FOR
C  SUBPROGRAMS "PRDCT", "CNTUPD", AND "DATIO".
cc Apr 24 1518 CDT 1990, make ephemeris file name (FILEPH) 255 bytes.
c (That means control files PRDCTRL and NPTCTRL must be bigger (500 bytes))
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cc 30 May 90, 0908 - make FILEPH 28 bytes again; PRDCTRL etc 250 bytes.
C  2/97 - expand to 252 (multiple of 4 char) for Linux. rlr.
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
	CHARACTER*18 CNTNAM
	CHARACTER*30 FILEPH
        CHARACTER*18 IZZZFL
	LOGICAL OBSPIK,RFLPIK
	DOUBLE PRECISION MINELV
	LOGICAL DEBUG
	INTEGER RSETPK
	COMMON/CNTINF/CNTNAM,OBSPIK(10),RFLPIK(5),NMODEL,RJD1,RJD2,
     1		STEPJD,MINELV,DEBUG(20),RSETPK,IEPH,FILEPH,
     2		IZZZFL
