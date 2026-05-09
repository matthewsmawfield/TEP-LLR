      SUBROUTINE LUNIO(MODE,INFTYP,NRECD,IRECD)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C.....
C     LUNIO HANDLES I/O FOR GENERALIZED LUNAR PARAMETER FILES
C     INPUT: MODE    1 = READ FROM DESIGNATED FILE
C                    2 = WRITE TO DESIGNATED FILE
C            INFTYP  INFORMATION TYPE; EACH IMPLIES A SPECIFIC FILE NAME
C                    AND RECORD LENGTH
C                    1 = SITE INFO
C                    2 = REFLECTOR INFO
C                    3 = LUNAR SURFACE INFO (Not functional here)
C                    4 = LUNAR MODEL PARAMETER INFO
C                    5 = LUNAR PREDICTION CONTROL PARAMETER INFO
C                    6 = PREDICTION CONTROL SET DESTINED FOR "prd"
C                    7 = NORMAL POINT/RECALC CONTROL SET FOR "npt"
C            NRECD   RECORD NUMBER 0 TO N (fortran records 1 to (n+1))
C                    0 WILL PRESUMABLY CONTAIN HEADER INFO AND
C                    1 TO N WILL EACH CONTAIN A COMPLETE SET OF INFO
C                    AN ARRAY TO BE WRITTEN
C     OUTPUT:IRECD   AN ARRAY TO BE READ FROM FILE
C                    Actually, the data is copied to the proper common.
C
C  REVISIONS:
C	07-20-87 - EXPAND SIZE OF LREFL.PC & GIVE IT A NEW NAME FOR
C       03-29-01 - Rewrite to use ascii files. rlr.
C.....

	character sccsid*(*)
	parameter (sccsid = "@(#)lunio01.f	1.8\t06/18/91")

	INCLUDE 'LIOCMN.INC'
	INCLUDE 'DEVCOM.INC'
	INCLUDE 'RFLINF.INC'
	INCLUDE 'OBSINF.INC'
	INCLUDE 'CNTINF.INC'
	INCLUDE 'MODINF.INC'
CC	INCLUDE 'LABELS.INC'

	CHARACTER IRECD(100)*18
	
	CHARACTER*18 pathtmp, pathtmps, efnam*14
        CHARACTER*255 filea, fileb
        INTEGER recdno
        CHARACTER*78 recdval
	save pathtmps
	
C...THIS IS WHAT IS ACTUALLY NOW IN BLOCK DATA SUBPROGRAM.
cc 20 Apr 90
C      DATA FNAME/"LSITES.PCA   ","LREFL1.PCA   ",
C     $           "LFEAT.PCA    ","LMODEL.PCA   ",
C     $           "LCONTOL.PCA  ","PRDCTRL.PCA  ",
C     $                           "NPTCTRL.PCA  "/
C.....

C.....
C.....OPEN THE APPROPRIATE FILE
	IERR=0
     	errmsg=' ERROR IN SUBROUTINE LUNIO OPENING FILE "'//
     1		FNAME(INFTYP)//'"'
     	OPEN(13,IOSTAT=IOS,ERR=9911,
     1		FILE=pathfil(1:lenfil)//FNAME(INFTYP),STATUS='OLD')
    
	IERR=1
9911  IF (IERR.NE.1) CALL ERROR(1,errmsg,IERR,IOS)
C.....
C..... Position to proper record
      IF (MODE.NE.1) GO TO 20
 100    read (13,1000,END=1) recdno, recdval
 1000   format (i2,a78)
        if (recdno .ne. nrecd) go to 100
C.....READ
CC        GO TO (1,2,3,4,5,6,7) INFTYP
        GO TO (1,2,3,4,5,6,6) INFTYP
C.....Telescope coordinates - transmit and receive
 1        if (nrecd .eq. 0) then
            read (recdval, 1010, END=101) irecd(1)
            do i=2,10
              read (13, 1015) recdno,irecd(i)
            enddo
 1010       format (1x,a18)
 1015       format (i2,1x,a18)
            go to 199
          endif
 101      read (recdval, 1100) obscode, obsite
          read (13, 1105) recdno, trad, tlat, tlong, tgdlat
          read (13, 1105) recdno, rrad, rlat, rlong, rgdlat
 1100     format (1x,f8.0,1x,a18)
 1105     format (i2,4(1x,f15.8))
          go to 199
C.....Reflector coordinates
 2        if (nrecd .eq. 0) then
            read (recdval, 1010, END=201) irecd(1)
            do i=2,10
              read (13, 1015) recdno,irecd(i)
            enddo
            go to 199
          endif
 201      read (recdval, 1200) rfset
          read (13, 1205) (recdno,serad(I),selat(I),selon(I), I=1,5)
 1200     format (1x,a18)
 1205     format (4(i2,1x,f15.9,2(1x,f15.8),/),i2,1x,f15.9,2(1x,f15.8))
          go to 199
C.....Surface feature coordinates
 3        read (recdval, 1300) srfnam
          read (13, 1305) recdno, sfrd,dcs1,dcs2
 1300     format (1x,a18)
 1305     format (i2,3(1x,f15.8))
          go to 199
C.....Model parameters
 4        if (nrecd .eq. 0) then
            read (recdval, 1010, END=401) irecd(1)
            do i=2,10
              read (13, 1015) recdno,irecd(i)
            enddo
            go to 199
          endif
 401      read (recdval, 1400) nammod
          read (13, 1405) recdno, twoway,polmotn, TOLERR,ITERMAX,
     1          recdno, VLIGHT, AUKM,EMRATI,
     2          recdno, K2LOVI,PMTITL,UT1TID
 1400     format (1x,a18)
 1405     format (i2,1x,2l2,1x,d10.3,1x,i5,/,
     1           i2,3(1x,d20.12),/,
     2           i2,1x,d10.3,1x,a18,1x,l2)
          go to 199
C.....Control parameters
 5        if (nrecd .eq. 0) then
            read (recdval, 1010, END=6) irecd(1)
            do i=2,10
              read (13, 1015) recdno,irecd(i)
            enddo
            go to 199
          endif
 6        read (recdval, 1500) cntnam
          read (13, 1505) recdno, (obspik(I),I=1,10),(rflpik(I),I=1,5),
     1          recdno, nmodel,rjd1,rjd2,stepjd,minelv,
     2          recdno, (debug(I),I=1,20),
     3          recdno, rsetpk,ieph,fileph
 1500     format (1x,a18)
 1505     format (i2,1x,10l2,5x,5l2,/,i2,1x,i8,4(1x,D15.8),/,
     1            i2,1x,20l2,/,i2,1x,2(1x,i8),1x,a18)
 199    CONTINUE
      IF (IERR.NE.1) CALL ERROR (1,errmsg,IERR,IERR)
      CLOSE (13)
      GO TO 50
C.....
C.....WRITE
   20 IF (MODE.NE.2) GO TO 50
C.....OPEN THE OUTPUT TEMPORARY FILE
	IERR=0
     	errmsg=' ERROR IN SUBROUTINE LUNIO OPENING FILE "'//
     1		"tmp_"//FNAME(INFTYP)//'"'
     	OPEN(14,IOSTAT=IOS,ERR=9912,
     1		FILE=pathfil(1:lenfil)//"tmp_"//FNAME(INFTYP))
    
	IERR=1
9912  IF (IERR.NE.1) CALL ERROR(1,errmsg,IERR,IOS)
C.....Copy unchanged data to output file.
 200    read (13,1020,END=205) recdno, recdval
 1020   format (i2,a78)
        if (recdno .ne. nrecd) then
          write (14,1020) recdno, recdval
          go to 200
        endif
CC        GO TO (11,12,13,14,15,16,17) INFTYP
 205    GO TO (11,12,13,14,15,16,16) INFTYP
C.....Telescope coordinates - transmit and receive
C 11       if (recdno .eq. 0) then
 11       if (nrecd .eq. 0) then
            do i=1,10
              write (14, 1015) nrecd,irecd(i)
            enddo
            nout= 9
            go to 299
          endif
          write (14, 1101) nrecd, obscode, obsite
          write (14, 1105) nrecd, trad, tlat, tlong, tgdlat
          write (14, 1105) nrecd, rrad, rlat, rlong, rgdlat
 1101     format (i2,1x,f8.0,1x,a18)
          nout= 3
          go to 299
C.....Reflector coordinates
 12       if (nrecd .eq. 0) then
            do i=1,10
              write (14, 1015) nrecd,irecd(i)
            enddo
            nout= 9
            go to 299
          endif
          write (14, 1015) nrecd, rfset
          write (14, 1205) (nrecd,serad(I),selat(I),selon(I), I=1,5)
          nout= 6
          go to 299
C.....Surface features
 13       write (14, 1015) nrecd, srfnam
          write (14, 1305) nrecd, sfrd,dcs1,dcs2
          nout= 2
          go to 299
C.....Model parameters
 14       if (nrecd .eq. 0) then
            do i=1,10
              write (14, 1015) nrecd,irecd(i)
            enddo
            nout= 9
            go to 299
          endif
          write (14, 1015) nrecd, nammod
          write (14, 1405) nrecd, twoway,polmotn, TOLERR,ITERMAX,
     1          nrecd, VLIGHT, AUKM,EMRATI,
     2          nrecd, K2LOVI,PMTITL,UT1TID
          nout= 4
          go to 299
C.....Control parameters
 15       if (nrecd .eq. 0) then
            do i=1,10
              write (14, 1015) nrecd,irecd(i)
            enddo
            nout= 9
            go to 299
          endif
 16       write (14, 1015) nrecd, cntnam
          write (14, 1505) nrecd, (obspik(I),I=1,10),(rflpik(I),I=1,5),
     1          nrecd, nmodel,rjd1,rjd2,stepjd,minelv,
     2          nrecd, (debug(I),I=1,20),
     3          nrecd, rsetpk,ieph,fileph
          nout= 5
 299  continue
      errmsg=' ERROR IN SUBROUTINE LUNIO; "WRTR" ERROR NUMBER BELOW.'
      IF (IERR.NE.1) CALL ERROR (1,errmsg,IERR,IERR)
C.....Bypass data we just wrote
        do i=1,nout
          read (13,1020,END=290) recdno, recdval
          if (recdno .ne. nrecd) write (14,1020) recdno, recdval
        enddo
C.....Copy unchanged data to output file.
 210    read (13,1020,END=290) recdno, recdval
        write (14,1020) recdno, recdval
        go to 210
 290    CLOSE (13)
        CLOSE (14)
        filea = pathfil(1:lenfil)//"tmp_"//FNAME(INFTYP)
        istat= nullterm (filea,255)
        fileb = pathfil(1:lenfil)//FNAME(INFTYP)
        ierr= nullterm (fileb,255)
        ierr= rename (filea, fileb)
        if (ierr .LT. 0) then
           call perror(' lunio - call rename') 
        endif
C.....
C.....End
50      CONTINUE
      RETURN

C	   THE
	   END

