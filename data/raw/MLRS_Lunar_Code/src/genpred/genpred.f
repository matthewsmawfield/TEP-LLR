C        PROGRAM EULER(INPUT,OUTPUT,TTY,UTPOLE,POLDAT,MITEPHB,
C     X TAPE5=INPUT,TAPE6=OUTPUT,TAPE10=TTY) 
        IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C *****************************************************************************
c Copyright (c) 2017, The University of Texas at Austin
c All rights reserved.
c
c Redistribution and use in source and binary forms, with or without
c modification, are permitted provided that the following conditions are met:
c
c 1. Redistributions of source code must retain the above copyright notice,
c this list of conditions and the following disclaimer.
c
c 2. Redistributions in binary form must reproduce the above copyright notice,
c this list of conditions and the following disclaimer in the documentation
c and/or other materials provided with the distribution.
c
c 3. Neither the name of the copyright holder nor the names of its contributors
c may be used to endorse or promote products derived from this software without
c specific prior written permission.
c
c THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
c AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
c IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
c ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE
c LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
c CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
c SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
c INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
c CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
c ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF
c THE POSSIBILITY OF SUCH DAMAGE.
c
C *****************************************************************************
C:%TITLE PREDICT -- LUNAR PREDICTION PACKAGE
C:%NOSOURCE
C:%G LUNAR RANGE AND POSITION PREDICTION PROGRAM
C:%T PURPOSE
C:THIS PROGRAM IS DESIGNED TO PRODUCE LUNAR RANGE AND POSITION PREDICTIONS
C:UNDER A GENERALIZED SET OF CONDITIONS.
C:FILES DETAILING MODEL PARAMETERS, REFLECTOR, TRANSMITTER, & RECEIVER
C:COORDINATES, AND CONTROL PARAMETERS ARE PRODUCED BY PROGRAM 'PREDCNTRL'
C:WHICH INITIATES THIS PROGRAM.
C
C$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$ RA & RLR
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C  Wed Mar 8, 1989 09:23:56
C     1) Add CHAR*1 EOSRC to call to INITRNG
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C  Fri Mar 10, 1989 07:57:10
C     1) MAKE SURE ALL CONSTANTS ARE DOUBLE PRECISION.
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C  Tue Mar 28, 1989 10:57:00
C     1) CONFORM TO D.P. SECONDS IN JDTGR ARGUEMENT LIST.
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C  Tue Oct 23, 1989
C     1) RECORD RA AND DEC FOR STAR OFFSETS.
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C  Oct 89 - tlr
C     1) initial & final lunar RA and DEC in ephemeris header
c     2) set modeop = -1
c     3) fix "14:29:60.0" problem
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C  3-7 Nov 89 - tlr
C     1) change to call sublasr03
C     2) change to call pmodrep
C     3) change to call pparin
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
c 30 May 90, 0953 - tlr
c  o Test RFLPIK for TRUE or FALSE, not "Y" or "N".
c  o Compare "secs1" from jdtgr with double precision 59.95d0,
c    not single precision.
c  o Don't let pparin change values of rjd1, rjd2.
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C 14 Jun 90, 1432 - tlr
c  o When program is done, append line w/ ephemeris span and
c    current time to "Lunar_Pred_Recd".
c  o Form the log file name from the year and day of year
c    catenated with "lun_pred_log".
c  o If there is a command line arguement '-a', do automatic prediction.
c c c c c c c c c c c c c c c c c c c c c c c c c c c c c c
c 06 Jul 90 through approx. 24 July 90 - tlr
c  o Changes to automatic predictions:
c     - Do predictions based on UTC calling gmtime through c routine.
c     - Use day-of-year in file names.
c  o Check for moon rise/set only for center of moon, not each reflector.
c  o Put Gregorian dates in prediction file.
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C (This is when we started making the program capable of using either
C  mit or jpl ephemerides.)
C 29 Aug 90, 1052 -
C        o  Remove "/us1" from pathnames.
C        o  Put ephinit after pparin so that the latter can
C           tell the former which ephemeris to open.
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C 02 Oct 90, 1130 -
C  o Automatic predictions: in log file and output ephemeris heading,
C    make ending time exact, not initial estimate.
C  x (don't)Call "chmod" for temporary files; give everyone r,w privileges.
C  o Activate pathname for log file in DEVCOM.INC.
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C 28 Nov 90, 1546: chmod problem already taken care of with scratch
C                   temporary log file.
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C 07 Feb 91, 1028: add eosrc to /mdrep/ so pmodrep can print it.
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C 07 Feb 91, 1521: for manual as well as auto, do pmodrep at moonrise.
C                   (so eosrc will be defined.)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C 18 Feb 91, 1733: open final log file (77) w/ access='append"
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C 21 Feb 91, 1624: convert to new pathnames; use pathname in rename.
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C 08 Mar 91, 1048: add parameter to error for calling HP-UX exit.
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C 31 May 91, 1354: Add "ndays" parameter to call to "iprdargs" to allow
C  doing predictions for given number of days after program start time.
C                                 finished: 13 Jun 91, 1331
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C 27 Jun 91, 1031: LRADEC(1-4) - make ra in hours, dec in degrees.
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C 03 Jul 91, 1459: add code for getting polar motion file from user
C        or command line.
C 14 Dec 91         : Close the temp file (unit 24) on exit.  This should
C        keep the tmp.F file from staying around.
C 10 Dec 01      : Allow input of prediction file name from command line. rlr.
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
        character*(*) sccsid
        parameter (sccsid = "%W%\t%G%")
        DOUBLE PRECISION LHATRN(6),AZALTRN(6),LHARCV(6),AZALRCV(6)
        double precision GED2GEC
        double precision rnjdint, rnjdf
        double precision lastgast
        LOGICAL RISEN, moreprd, done
        INTEGER D1(5), D2(5), idf(6),idl(6)
        integer velout
        CHARACTER*1 EOSRC,nul,ephfp*30, ephfpx*30, logstring*26
        character*20 tempeph, templog
C        character*81 line
        character*101 line
        character*256 ephname
        character*10 ephlabel
chmod chmod chmod chmod chmod chmod chmod chmod chmod chmod chmod
ccc        integer chmod
ccc$ALIAS  chmod = 'chmod' (%VAL,%VAL)
chmod chmod chmod chmod chmod chmod chmod chmod chmod chmod chmod

C........... INCLUDE FILES .....................................

        INCLUDE 'BLD.INC'

        INCLUDE 'BOUNCE.INC'

        INCLUDE 'CETB1.INC'

        INCLUDE 'CETB4.INC'

        INCLUDE 'CNTINF.INC'

        INCLUDE 'DDMU1.INC'

        INCLUDE 'DEFCMN.INC'
 
        INCLUDE 'DEPHREC.INC'

        INCLUDE 'DEVCOM.INC'

        INCLUDE 'DIRCTN.INC'

        INCLUDE 'EPHREC.INC'

        INCLUDE 'INTRNG.INC'

        INCLUDE 'LIB.INC'

        INCLUDE 'MODEOP.INC'

        INCLUDE 'MODINF.INC'

        INCLUDE 'OBSINF.INC'

        INCLUDE 'PDCMN.INC'

        INCLUDE 'POLCOM.INC'

        INCLUDE 'PRDARGCOM.INC'

        INCLUDE 'RFLINF.INC'

        INCLUDE 'U1UC.INC'

        INCLUDE 'USNOC.INC'

        COMMON /EOP/ XPOLE, YPOLE, UT1C
                                                                                
C...............................................................

        DOUBLE PRECISION PCS(14)
        EQUIVALENCE (PCS(1),C1)

        COMMON/TWOPI/TWOPI,CDR,CASR,CTSR

        COMMON /FALOOK/ FA(26)

COOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOO
C  The following named common blocks were included in main program
C  EULER so that they would remain defined during execution of the
C  program.  In ABSoft (MicroSoft) FORTRAN, named common blocks may
C  become undefined if the highest level subprogram in which they
C  are defined executes a RETURN or END.
C...............................................................

        INCLUDE 'PVHLD.INC'

        INCLUDE 'ATMAUX.INC'

        INCLUDE 'CETB2.INC'

        INCLUDE 'DATA2.INC'

        INCLUDE 'PRMAT.INC'

        INCLUDE 'PNT.INC'

        INCLUDE 'FUNBLK.INC'

C...FOR KUSTNER TO SAVE INFO
        common/kuscmn/ frstim,tt(4),jd1,jd2,jdstep
        logical frstim
        double precision jd1,jd2,jdstep

C...............................................................
        INTEGER RFLSET(5)
c put these in common so pmodrep can use them.
        DOUBLE PRECISION JDINT1,JDF1,JDINT2,JDF2
        common /mdrep/ rflset,jdint1,jdf1,jdint2,jdf2, eosrc

        integer gunit(6)
        character*2 gunitasc
        character*17 mvern
        logical first(6), firstopen, recd31
        character*10 refname(6)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

c  Only open debugging file if at least one debug flag is true.
        logical anydeb

C       DATA mvern /"v1.00 16 Nov 2006"/
        DATA mvern /"v1.02 12 May 2010"/
C       V 1.00 - Nov 2006
C       v 1.01 - Unofficial. DE-421 with eph version 2. Discarded.
C       v 1.02 - Take CPF ephemeris name from nammod, the model name.
C
C  MOON HAS NOT YET RISEN:
        DATA RISEN /.FALSE./
        DATA RFLSET /0,-1,2,3,4/
        DATA first /.true.,.true., .true., .true., .true., .true./
        DATA recd31 /.false./
        DATA lastgast /-100.d0/
        DATA refname /"luncenter ","apollo11  ","luna17    ",
     x                "apollo14  ","apollo15  ","luna21    "/
C-----------------------------------------------------------------

        nul=char(0)
        done=.FALSE.
        firstopen= .true.
        ieph_sub_seq= 1
	iterdir= 1
CC        ephname= 'LUNxxx.le.'

cXXXXXXXXXXXXXXXXXXXXXXXXXXXXX << 1 >> XXXXXXXXXXXXXXXXXXXXXXXXXXX
c get date, time, day-of-year (these are UTC values!)
        nsecs=icgetime(nowdoy,nowyr,nowmon,nowda,nowhr,nowmin,nowsec)
        errmsg = ' Error in icgetime calling time'
        if(nsecs .EQ. -1) call error(1,errmsg,nsecs,nsecs)
        rnowsec = nowsec
        call grtjd(nowyr,nowmon,nowda,
     1                nowhr,nowmin,rnowsec,rnjdint,rnjdf)
        call grtdoy(nowyr,nowmon,nowda,nowdoy)
C genpred changes:
         rnjdf= 0.d0
          pjdint = rnjdint- 1.d0
          pjdf = rnjdf

cXXXXXXXXXXXXXXXXXXXXXX << 2 >> XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
c check for command line argument of '-a' (or '-a5', -a12', etc.), and
c for "-p <filename>" for polar motion data file, and "-f filename" for
c prediction control file.
        nargs=iprdargs(autoprd, ndays, ephname, iterdir, recd31)

        if(autoprd) then
c          if so, divert console output to scratch file.
          ioutdev=-1
        endif

CXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
2        continue
CXXXXXXXXXXXXXXXXXXXXXXXXXX < 3 > XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
CC        pathfil='data/lib/'
CC        patheph='data/llr/pred/'
CC        pathlog='data/llr/pred/'
        pathfil='data/lib/'
        patheph='data/pred/'
        pathlog='data/pred/'
        lenfil=itrimlen(pathfil)
        leneph=itrimlen(patheph)
        lenlog=itrimlen(pathlog)

CXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX < 3 1/2 > XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
C  OPEN LOG FILE NAMED "EULER-LOGFILE".  THIS FILE MAY BE PRINTED AFTER
C  THE USER EXITS THE EULER PROGRAM.  IN THE FUTURE, WE MAY WANT TO GIVE
C  THE USER THE OPTION OF CHANGING THE OUTPUT DEVICE TO 6: LINE PRINTER.
C  (P.S. THE DEFAULT VALUE OF "EULOG" IS 11.)
c   (Now consider this.  The program "EULER" references subprograms in
c   the normal point directory.  Those subprograms are conditioned to
c   write to unit "nplogfil" which is equal to 10.  Since everyone
c   should write to the same unit, we redefine "nplogfil" to be eulog
c   before any calls to normal point routines are made.)
        nplogfil = eulog
c Open a temporary log file. At the end of the program, the log file
c will be renamed with a name involving the year and day of year.
        templog(1:15)='xxxdel_log_temp'
        templog(16:16)=char(0)
        open(eulog,status='scratch')

CXXXXXXXXXXXXXXXXXXXXXXXXX < 4 > XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
        MODEOP = -1
        SITENO = 999
        EPHSRC = 888
        MINELV = 0.D0

CXXXXXXXXXXXXXXXXXXXXXXXX < 5 > XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
c:read polar motion & ut1-utc polynomial coefficients
        call polin(pmtitl,pcs)

CXXXXXXXXXXXXXXXXXXXXXX < 6 > XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
c:write greeting
        if (ioutdev .gt. 0) write(ioutdev,5009) mvern
        write(eulog,5009) mvern
5009    format(18x,'Generalized Solar System Prediction Package',1x,A17,
     1                /18x,43('*')
     2                //,' Enter output ephemeris file name: ')
5010    format(1x,'Enter polar motion file name ',
     1                                '(just return for "poldat")')
5011    format(' POLAR MOTION FILE NAME')

c:get output ephemeris file name (manual or auto mode).  in manual mode,
c:also get polar motion file name (default is "poldat")
        if(autoprd) then
          ephfp  ='LUNxxx_del.le'
          tempeph='LUNxxx_del.le'
          lenteph=itrimlen(tempeph)
        else
          read(indev,'(a)') ephfp
          write(eulog,'(1x,a20)') ephfp
          ephname= ephfp
          if(.NOT.ispol) then
            if (ioutdev .gt. 0) write(ioutdev,5010)
            read(indev,'(a)')filpol
            if(filpol(1:1) .EQ. ' ')filpol='poldat'
            lenpol=itrimlen(filpol)
          endif
        endif

CXXXXXXXXXXXXXXXXXXXXXXXXX < 7 > XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
        call pparin(vlite)
        if (iterdir .eq. -1) write (eulog,5100) 
 5100   format(/" Predicting backwards in time"/)
ccc Set ephemeris source (EPHSRC) and site number (SITENO) if necessary. ccc

CXXXXXXXXXXXXXXXXXXXXXXXXX < 7 1/2 > XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
c put ephinit here so pparin can tell it which ephemeris to open.
      REM=AUKM
      CALL EPHINIT(AUKM,REM)

CXXXXXXXXXXXXXXXXXXXXXXXX < 8 > XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
c  Only open debugging file if at least one debug flag is true.
        anydeb=.FALSE.
        do 1 i1=1,20
        if(deb(i1))anydeb=.TRUE.
1        continue
c  change the definition of dunit and open debug file.
        dunit = 33
        if(anydeb)open(dunit,file='Prd-debug')

CXXXXXXXXXXXXXXXXXXXXXXXXX < 9 > XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
c get timespan of prediction from user or from current date, depending
c on whether we are in auto mode or not.  Subroutine "getprdtim" expects
c (pjdint,pjdf) to be the program execution time with each call.
        if(autoprd)then
            pjdf= 0.5d0
            pjdint= pjdint+ 1.0d0
          if (pjdint+pjdf .gt. rnjdint+rnjdf + ndays) then
            moreprd= .false.
            done= .true.
          else
            moreprd= .true.
            done= .false.
          endif  

c        ***************************************************
c        Quit if there are no more predictions to be done.
          if ( done ) goto 10001
c        ***************************************************
c          .... make sure initial time is integer multiple of stepjd
          rjdfstp=pjdf*86400.d0/stepjd
          ijdfstp=dint(rjdfstp)
          if((rjdfstp-ijdfstp) .GT. 0.01d0) ijdfstp=ijdfstp+1
          pjdf=ijdfstp*stepjd/86400.d0
c          .....................................................
          jdint1=pjdint
          jdf1=pjdf
          jdint2=jdint1
          putc2=dint(jdf1*86400.d0 + 0.5d0)
c          final time not actually used for automatic predictions.
          call addtim(jdint2,putc2,44700.d0)
        else
          moreprd = .FALSE.
          call dateqry(1,jdint1,jdf1)
          call dateqry(2,jdint2,jdf2)
          pjdint = jdint1
          pjdf = jdf1
          ndays= jdint2-jdint1
        endif
        rjd1 = jdint1 + jdf1
        rjd2 = jdint2 + jdf2
          call jdtgr(pjdint,pjdf,ipyr,ipmon,ipda,iphr,ipmin,psec)
          if(psec .GE. 59.95d0) call jdtgr(pjdint,pjdf+1.d-9,
     1                                ipyr,ipmon,ipda,iphr,ipmin,psec)
          call grtdoy(ipyr,ipmon,ipda,ipdoy)

CXXXXXXXXXXXXXXXXXXXXXXXXXXX < 10 > XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
c:initialize output ephemeris file
        recd(1) = jdint1
        recd(2) = jdf1
        istep= idint(stepjd)

CXXXXXXXXXXXXXXXXXXXXXXXXXXXXX et cetera XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX

C  INSURE TIME STARTS ON EVEN SECOND.
      TJDINT = JDINT1
      SUTCT = DINT((JDF1)*86400.D0+.5D0)
      NRCD = 1
      NP = 0

C We want xmit and xcv angles, etc.
      IRECVR= 2

      if (firstopen) then
      firstopen= .false.
       do i= 0,5
        gunit(i+1) = 15+i

C  Define name based on ILRS standards
        ilen= itrimlen(refname(i+1))
        iyr= mod(ipyr,100)
        write(ephfpx,1001) refname(i+1)(1:ilen),iyr,ipmon,ipda,
     &                     (ipdoy+500)*10+ ieph_sub_seq
 1001   FORMAT (a,'_cpf_',i2,i2,i2,'_',i4,'.utx')

C  provide option to create user-defined id and station|geocenter to name
        ilen= itrimlen(ephname)
        if (ilen .gt. 0 .and. ilen .lt. 256) then
          ephfpx= 'Lx-'//CNTNAM(1:4)//'-'//ephname
          if (i .eq. 0) then
            write(ephfpx(2:2),'(a1)') "C"
          else
            write(ephfpx(2:2),'(i1)') (i-1)
          endif
        endif
        ilen= itrimlen(ephfpx)
        do j=1,ilen
          if (ephfpx(j:j) .eq. ' ') ephfpx(j:j)= '0'
        enddo
        open(gunit(i+1),file=ephfpx)
       enddo
      endif


C  TIME LOOP
C:TIME = BEGIN TIME
C:DO WHILE TIME < ENDTIME
100        NP=NP+1
        FJDUTC = SUTCT/86400.D0
        TJDUTC = TJDINT + FJDUTC
      CALL JDTGR(TJDINT,FJDUTC,D1(3),D1(1),D1(2),D1(4),D1(5),SECS1)
        if(secs1 .GE. 59.95d0)then
          call jdtgr(tjdint,fjdutc+1.d-9,
     2                        d1(3),d1(1),d1(2),d1(4),d1(5),secs1)
        endif
        d1(3)= d1(3)+ 1900
      RECD(1) = TJDINT
      RECD(2) = SUTCT
      DRECD(1) = TJDINT
      DRECD(2) = SUTCT
c Save the ephemeris time in case it is the first
c for output to "Lunar_Pred_Recd", and for constructing file name.
        if(nrcd .EQ. 1)then
          rjdff = tjdutc
          rjdist = tjdint
          rsutcst = sutct
           idf(1)  = d1(3)
           idf(2)  = d1(1)
           idf(3)  = d1(2)
           idf(4)  = d1(4)
           idf(5)  = d1(5)
           idf(6)  = secs1 + 0.5d0
        endif

C  REFLECTOR LOOP
C:DO FOR CENTER OF MOON & EACH REFLECTOR
      DO 400 NREFL=-1,4

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
C Get the irritating geocentric latitude for time SDE calculations, etc.
      IF (TLAT .EQ. 0.d0 .OR. TLAT == TGDLAT) TLAT= GED2GEC(TGDLAT)
      IF (RLAT .EQ. 0.d0 .OR. RLAT == RGDLAT) RLAT= GED2GEC(RGDLAT)

C:LOAD TARGET COORDINATES
      IF (NREFL.NE.-1) GO TO 200
      REFLCT(1) = 0.D0
      REFLCT(2) = 0.D0
      REFLCT(3) = 0.D0
      GO TO 210

  200 NRFP1 = NREFL+1
      IF (.NOT. RFLPIK(NRFP1)) GO TO 400
      REFLCT(1) = SERAD(NRFP1)
      REFLCT(2) = SELON(NRFP1)
      REFLCT(3) = SELAT(NRFP1)
210   continue

      call sublasr(tjdint,sutct,vlite,ut1c,xpole,ypole,eosrc,
     1                     ttltim,lhatrn,azaltrn,lharcv,azalrcv,
     2                     corr12,corr23,dtadt)

C  RECORD THE RANGE 
375        continue
        IF (NREFL.GE.0) RECD(31+NREFL) = TTLTIM

        RISEN= .TRUE.
        if(nrcd.EQ.1 .AND. nrefl.EQ.-1) then
c        1st prediction, and for center of moon
          if(autoprd) then
c            is also automatic prediction.
            call grtdoy(d1(3)-1900,d1(1),d1(2),idoyprd)
            ephfp='Lxxx'//ephname//CNTNAM(1:4)
            write(ephfp(2:4),'(i3.3)')idoyprd
            irx = rename(patheph(1:leneph)//tempeph(1:lenteph)//nul,
     1                   patheph(1:leneph)//ephfp(1:14)        //nul)
            write(eulog,'('' Auto mode: '',a20)')ephfp
            write(eulog,'('' Polar motion file: '',a)') filpol
          endif
c  we do all (auto & manual) pmodreps here so they know what eosrc is.
          call pmodrep
          if (eosrc .NE. 'T') call system("XClunprederr")
        endif

C:PRINT REFLECTOR NUMBER, AZ, ALT, & RANGE ONTO TTY & LPT
      IF (NREFL.NE.-1) GO TO 380
        WRITE(EULOG,1035) NRCD,TJDUTC,D1,SECS1
 1035 FORMAT(/1X,'RECORD:',I3,3X,'TIME:',F15.6, 3X,2(i2.2,'/'),i4.4,
     $                2X,2(i2.2,':'),F4.1/
     2        3X,'REFL',3X,'AZ',9X,'EL',9X,'LHA',8X,
     3                                'RA',9X,'DEC',6X,'RANGE')

C SAVE STARTING AND ENDING RA & DEC FOR STAR OFFSET.
c        (save as hours and degrees)
        IF(NRCD .EQ. 1) THEN
          LRADEC(1) = LRA
          if(lradec(1) .LT. 0.d0) lradec(1)=lradec(1)+24.d0
          IF(LRADEC(1) .GE. 24.d0) LRADEC(1) = lradec(1)-24.d0
          LRADEC(2)=LHATRN(3)
        ELSE
          LRADEC(3)=LRA
          if(lradec(3) .LT. 0.d0) lradec(3)=lradec(3)+24.d0
          IF(LRADEC(3) .GE. 24.d0) LRADEC(3) = lradec(3)-24.d0
          LRADEC(4)=LHATRN(3)
        ENDIF

380        continue

        if(lra .LT. 0.d0) lra = lra+24.d0
        IF(LRA .GE. 24.d0) LRA = lra-24.d0
        WRITE(EULOG,1050) NREFL,AZALTRN(2),AZALTRN(3),
     2                LHATRN(2),LRA,LHATRN(3),TTLTIM
C     2                LHATRN(2),LRA,LHATRN(3),TTLTIM, range
        WRITE(EULOG,1050) NREFL,AZALRCV(2),AZALRCV(3),
     2                LHARCV(2),LRA,LHARCV(3),TTLTIM
 1050 FORMAT(5X,I1,2F11.5,2F12.6,F11.5,2F16.11)

C:IF DESIRED, PRINT LHA & DEC & AXIS RATES
      IF (DEB(20)) CALL LRATE(NRCD,NREFL,LHATRN(2),LHATRN(3))

C:ENDDO
400        CONTINUE

C:WRITE RECORD ON OUTPUT EPHEMERIS FILE, save the jd w/ utc in case it is
c the last, and increment the record counter.
        rjdll=tjdutc
        if(autoprd) then
          jdint2=tjdint
          jdf2=fjdutc
          rjd2=tjdint+jdf2
        endif
           idl(1)  = d1(3)
           idl(2)  = d1(1)
           idl(3)  = d1(2)
           idl(4)  = d1(4)
           idl(5)  = d1(5)
           idl(6)  = secs1 + 0.5d0
      do i=0, 5
        if (first(i+1)) then
          first(i+1)= .false.
          velout= 0; ! usually no velocity output
C        Nedd to get the ivers and eph_seq from somewhere...
          ivers= 1
          ieph_seq= (ipdoy+500)*10+ ieph_sub_seq
          icompat= 0
          ittype= 2
          iref_frame= 0
          iratype= 0
          icofmapp= 0
          if (i .EQ. 0) iratype= 1
CC          if (ieph .eq. 0) then
CC            ephlabel= "mitpep    "
CC          else if (ieph .eq. 1) then
CC            ephlabel= "jpl_de-421"
CC          else
CC            ephlabel= "unknown   "
CC          endif
          ephlabel= nammod(1:10)

          call jdtgr(drecd(1),drecd(2)/86400.d0,
     x                ipyr,ipmon,ipda,iphr,ipmin,psec)
          call jdtgr(drecd(1)+ndays,(drecd(2)-stepjd+0.001)/86400.d0,
     x                ipeyr,ipemon,ipeda,ipehr,ipemin,pesec)
        write (gunit(i+1),1200) ivers,nowyr+1900,nowmon,nowda,nowhr,
     x      ieph_seq, refname(i+1), ephlabel
 1200   format ("H1 CPF",1x,i2,1x," UTX",1x,i4,3(1x,i2),
     x      2x,i4,1x,A10,1x,A10)
        write (gunit(i+1),1201) 100+i-1, 100+i-1, 0,
     x                ipyr+1900,ipmon,ipda,iphr,ipmin,int(psec),
     x                ipeyr+1900,ipemon,ipeda,ipehr,ipemin,int(pesec),
     x                int(stepjd),icompat,ittype, iref_frame, iratype,
     &                icofmapp
 1201  format ("H2 ", I8, 1x, i4, 1x, i8, 1x, 2(i4, 5(1x,i2),1x), 
     x          i5, 1x,i1,1x,i1,1x,i2,1x,i1,1x,i1)
        write (gunit(i+1),1202) 
 1202        format ("H9")
        endif
        
        dmjd= (drecd(1)-2400000.d0) + (drecd(2)/86400.d0-0.5d0)
        mjd= int(dmjd)
        dsod= (dmjd- mjd)*86400.d0

        write (gunit(i+1),1210) mjd,dsod,(drecd(k),k=3+(i*12),5+(i*12))
 1210        format ("10 1 ",i5, 1x,f7.1,1x,"0",3(1x,f18.3))
        if (velout .eq. 1) then
          write (gunit(i+1),1212) (drecd(k),k=6+(i*12),8+(i*12))
 1212          format ("20 1 ",3(1x,f18.6))
        endif

        write (gunit(i+1),1211) mjd,dsod,(drecd(k),k=9+(i*12),11+(i*12))
 1211        format ("10 2 ",i5, 1x,f7.1,1x,"0",3(1x,f18.3))
        if (velout .eq. 1) then
          write (gunit(i+1),1213) (drecd(k),k=12+(i*12),14+(i*12))
 1213          format ("20 2 ",3(1x,f18.6))
        endif
        write (gunit(i+1),1214) (drecd(74+k+3*i), k=1,3), 
     x                corr12*86400.d9
 1214        format ("30 1 ",3(1x,f7.0),1x,f5.1,1x)
        if (recd31) then
           write (gunit(i+1),1215) (drecd(96+k), k=1,3)
 1215      format ("31 1 ",3(1x,f20.10))
        endif
        if (i .eq. 0 .AND. ieph .ne. 0) then
        write (gunit(i+1),1216) mjd, dsod, (drecd(92+k), k=1,3),
     x                drecd(96)/CDR/15.d0
 1216        format ("60 ",i5, 1x,f7.1,3(1x,f18.12), f18.12)
        lastgast= drecd(93)
        endif

      enddo
      NRCD= NRCD+ 1

        if (tjdutc.GE.rjd2 .AND. .NOT.autoprd) go to 500
c (if this is autoprd, keep going until moon set.)

      SUTCT = DINT(SUTCT+STEPJD+0.1D0)
C  genpred: create new file at UTC day boundary. 
      if (dabs(SUTCT- 129600.d0) .lt. 14.d0) go to 500
      GO TO 100

C:COMPLETE OUTPUT EPHEMERIS PROCESSING
500        CONTINUE
        RECD(1)=JDINT2
        RECD(2)=JDF2
        RECD(3)=obscode
        do 501 i501=1,6
        i=i501+3
        j=i501+3+6
        recd(i)=idf(i501)
501        recd(j)=idl(i501)

c leave a trace of this prediction run.
        open(77,access='append',
     1                file=pathlog(1:lenlog)//'Lunar_Pred_Recd')
        write(77,5077)nowyr+1900,nowmon,nowda,nowhr,nowmin,nowsec,
     1                rjdff,rjdll
5077    format(i4.2,2('/',i2.2),', ',i2.2,':',i2.2,':',i2.2,
     1         ' ephemeris spanning ',f13.5,'-',f13.5,' created')
        close(77)

c   rename the log file.
        call addtim(rjdist,rsutcst,0.d0)
        rjdfst=rsutcst/86400.d0
        call jdtgr(rjdist,rjdfst,ipyr,ipmon,ipda,iphr,ipmin,psec)
        if(psec .GE. 59.95d0) call jdtgr(rjdist,rjdfst+1.d-9,
     1                                ipyr,ipmon,ipda,iphr,ipmin,psec)
        call grtdoy(ipyr,ipmon,ipda,ipdoy)
c        (We are prepared for the twenty first century!)
        if(ipyr .GT. 60) then
          ipyr=ipyr+1900
        else
          ipyr=ipyr+2000
        endif
c        construct permanent log file name; rename file.
        logstring(26:26)=char(0)
        logstring(5:5)='-'
        logstring(13:25)='.lun_pred_log'
        write(logstring(1:4),'(i4)')ipyr
        write(logstring(6:8),'(i3.3)')ipdoy
        write(logstring(9:12),'(a4)')CNTNAM(1:4)
c..................................................................
c Move the log information to the new log file and make sure the
c FINAL TIME is exact, not estimate.
        open(77,file = pathlog(1:lenlog)//logstring)
        rewind (eulog)
600        continue
        read(eulog,'(a)')line
        if(line(1:11) .EQ. ' FINAL TIME')goto 601
        lenlin=nullterm(line,101)
        if(lenlin .EQ.0)then
          write(77,'(1x)')
        else
          write(77,'(a)')line(1:lenlin)
        endif
        goto 600

601        continue
        call jdtgr(jdint2,jdf2,d2(3),d2(1),d2(2),d2(4),d2(5),secs2)
          if(secs2.GE.59.95d0) call jdtgr(jdint2,jdf2+1.d-9,
     1                        d2(3),d2(1),d2(2),d2(4),d2(5),secs2)
C Y2k 4-digit year
        d2(3)= d2(3)+ 1900
        write(77,1000) rjd2,d2,secs2
1000        format(1x,'FINAL TIME: ',  f15.4, 3x,2(i2.2,'/'),i4.4,2x,
     4                                             2(i2.2,':'),f4.1)

602        continue
        read(eulog,'(a)',end=603)line
        lenlin=nullterm(line,101)
        if(lenlin .EQ.0)then
          write(77,'(1x)')
        else
          write(77,'(a)')line(1:lenlin)
        endif
        goto 602
603        continue
c..................................................................
        close (eulog)
        close (77)
        close (24)

ccccccccccccccccccc check for more predictions (automatic case only)
        if(moreprd) goto 2
c
10001 continue
      do i=0, 5
        write (gunit(i+1),1220) 
 1220   format ("99")
      enddo
c  Terminate the program through HP-UX exit.
        call ceror(0)
c
c           the
           end

        double precision function ged2gec (gdlat)
        double precision gdlat
        double precision if, f,f2,f3, gl
        data if /298.257D+00/

        COMMON/TWOPI/TWOPI,CDR,CASR,CTSR
        double precision TWOPI,CDR,CASR,CTSR

        f= 1.d0/if
        f2= f*f
        f3= f*f2
        gl= gdlat*CDR
        dlat = (f+f2/2.d0)*dsin(2.d0*gl)- 
     x                ((f2+f3)/2.0d0)*dsin(4.d0*gl)+ 
     x                f3/3.d0*dsin(6.d0*gl)
        ged2gec= (gl- dlat)/CDR
CC        write (*,*) "geod, geoc, dlat", GDLAT,GED2GEC,dlat/CDR
        return
        end
