	SUBROUTINE WRLUN(OMC,AZ,ALT)
	IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      INCLUDE 'LICENSE-BSD3.inc'
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C  Wed Oct 19, 1988 15:39:19
C     1) REPLACE "COMPILER DOUBLE PRECISION" WITH
C        "IMPLICIT DOUBLE PRECISION (A-H,O-Z)".
C     2) MAKE SURE STATEMENTS DON'T EXTEND BEYOND COLUMN 72.
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C  Thu Nov 3, 1988 15:28:51
C     1) COMMENT OUT CALLS TO "FLA1".
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C  Mon Jan 30, 1989 09:44:28
C     1) REPLACE OLD CALLS TO "FLA1" W/ INTERNAL WRITES OF SCALED INTEGERS
C     2) IP AND BH IN /PASSIT/ ARE NOW CHAR*130 SIMPLE VARIABLES.
c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c
c  22 Nov 89
c     1) Fix write of xaz,alt,omc to internal file.
c	 (was writing to single character.)
c     2) Fix the format: (3i10) <- (2i10,i9).
c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c
c  Tue Dec  5 08:22:03 CST 1989
c	write xaz to internal file w/ f12.0 format, then move 
c	digits sans radix to ip.
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

C  WRITE LUNAR DATA TO OUTPUT FILE
C  REVISIONS:
C       09/29/85 - SET UP TO WRITE ALL TYPES OF LUNAR RECDS - RLR
C       12/03/86 - LUNAR '86 FORMAT. V10. RLR.
C	02/09/87 - FIX FOR NEGATIVE AZ. V11. RLR.
C	05/01/87 - Lunar format '86, mod 1. (mini-normalpoint). V12. rlr.
C$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$

	character*(*) sccsid
	parameter (sccsid = "@(#)wrlun12.f	1.6\t04/11/90")

	INCLUDE 'MODEOP.INC'

	INCLUDE 'PASSIT.INC'

	character temp*12

C  HANDLE NEGATIVE AZ:
	XAZ= AZ
	IF (XAZ.LT.0.D0) XAZ= 360.D0+ XAZ

C  HEADER RECORDS. IF RUN HEADER, WRITE RUN & BURST HEADERS OUT.
        IF (IP(1:1).NE.'1') GO TO 10
        WRITE(22,1010) IP
C...BUT ONLY IF IT'S NOT A NORMALPOINT HEADER
        IF (BH(68:68).NE.'N') THEN
	   WRITE(22,1010) BH
	ENDIF
        GO TO 40

C  DETAIL OR ARCHIVAL NORMALPOINT RECORDS
 10     IF (IP(1:1).NE.'3' .AND. IP(1:1).NE.'4') GO TO 20
C...DON'T UPDATE ON NORMALPOINT RUN
	IF (MODEOP.EQ.0) GO TO 30
C...WRITE BURST HEADER FOR NORMALPOINTS HERE...
	IF (IP(1:1).EQ.'4') THEN
	   WRITE(22,1010) BH
	ENDIF
cc........................................................
cccc	IXAZ=XAZ*1.d7
c	IXAZ=XAZ*1.d6
c	jxaz=(xaz*1.d7)-(10.d0*ixaz)
cc.. Don't do above 2-integer stunt.  Instead, do ........
cc.. internal write and then move characters. ............
	zxaz=xaz*1.d7
cc........................................................
	IALT=ALT*1.d7
	IOMC=OMC*1.d3
cc........................................................
c	WRITE(IP(91:119),5701)IXAZ,jxaz,IALT,IOMC
c5701	FORMAT(i9,i1,I10,i9)
cc.. Two-integer code above replaced by code below. ......
	write(temp,5702)zxaz
	ip(91:100)=temp(2:11)
	WRITE(IP(101:119),5701)IALT,IOMC
5701	FORMAT(I10,i9)
5702	FORMAT(f12.0)
cc........................................................
        GO TO 30

C  MINI-NORMALPOINT RECORD
 20     IF (IP(1:1).NE.'5') GO TO 30
cc........................................................
cccc	IXAZ=XAZ*1.d7
c	IXAZ=XAZ*1.d6
c	jxaz=(xaz*1.d7)-(10.d0*ixaz)
cc........................................................
cc.. Don't do above 2-integer stunt.  Instead, do ........
cc.. internal write and then move characters. ............
	zxaz=xaz*1.d7
cc........................................................
	IALT=ALT*1.d7
	IOMC=OMC*1.d3
cc........................................................
cc	WRITE(IP(81:109),5701)IXAZ,jxaz,IALT,IOMC
cc.. Two-integer code above replaced by code below. ......
	write(temp,5702)zxaz
	ip(81:90)=temp(2:11)
	WRITE(IP(91:109),5701)IALT,IOMC
cc........................................................

C  'Z', CAL, AND COMMENT RECORDS
30	CONTINUE
	WRITE(22,1010) IP
1010	FORMAT(A130)

40     CONTINUE
C
C	   THE
	   END

