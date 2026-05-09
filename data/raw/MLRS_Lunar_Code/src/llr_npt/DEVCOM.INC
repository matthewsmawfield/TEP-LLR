C  INCLUDE FILE MADE Fri Aug 5, 1988 15:52:12 TO PROVIDE
C  SINGLE POINT DEFINITION OF CONSTANTS SUCH AS CONSOLE
C  INPUT AND OUTPUT DEVICES.
c  Thu Dec 14 19:04:59 CST 1989, add pathname parameters for opening files.
c  Apr 1990: add ierrout (set to 7 in common block.); add errmsg.
c  02 Oct 90, 1341: include log pathname
c  30 Jul 91      : include summary file name
	INTEGER EULOG
cccccc	character*255 patheph, pathfil, errmsg*80
	character*255 patheph, pathfil, pathlog, errmsg*80
	COMMON /DEVCOM/ INDEV, IOUTDEV, ierrout, NPLOGFIL, NPSUMFIL, 
     1	EULOG, patheph, pathfil, pathlog, leneph, lenfil, lenlog, errmsg 
ccccc1		patheph, pathfil, leneph, lenfil,  errmsg

