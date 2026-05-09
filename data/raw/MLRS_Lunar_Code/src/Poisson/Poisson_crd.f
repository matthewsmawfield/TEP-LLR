C
C  Poisson Analysis System.
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
C
C  Author:  Randall Ricklefs, McDonald Observatory, Univ of Texas
C  Date: 1988
C
C  This routine is driver for the poisson data filtering subroutine suite.
C  Essentially, this program reads the time and range residuals from
C  a range data file and then reads the filtering parameters for that
C  satellite from one of several parameter files.  The data is then subjected
C  to the Poisson analysis routines, which produces the slope and bin number
C  of the most anomalous bin (if there is one).  If there is a substantially
C  anomalous bin, the data from that bin pair is tagged for further processing.
C  The tags are written to a copy of the input data file called name.out.
C
C  Program prompts for input data file name.  To quit program, type 'qq'
C  for file name or type '^C'.
C
C  Parameters:
C      MAXOBS  - maximum number of observations expected per pass.
C
C  Revisions:
C      11/29/89 - specialize for lunar ranging. Poisson. rlr.
C v1.3 08/09/90 - 1) Allow user to specify alternate command parameter file
C      with '-c comname' switch;
C         2) Specify through the command file the minimum number
C      of points per bin that is acceptable.
C v1.4 08/10/90 - 1) Use mint and maxt to limit data used.
C v1.5  08/15/90 - 1) Flag cases where the highest number of points in a
C      time bin is available in more than one bin.  In such
C      a case, do not select any data.
C         2) warn and create npt break records for multiple time bins.
C         3) set up problem warning (for batch processing).
C v1.6 08/20/90 - 1) Recognize multi-slope ties describing identical or
C      nearly-identical data points.
C         2) Recognize when wider residual bin-width should be used,
C      and use it.
C         3) Display date/time of run on listing
C v1.7  09/12/90 - 1) Make slope tie a comment, not a warning.
C         2) clean up lunar time breaks.
C v1.8  10/12/90 - 1) expand t & r in avg_data to "MAXOBS";
C         2) checked for too many points in readres;
C         3) changed error handling in writeres;
C         4) changed temp file length in filnme to 256.
C v1.8a 09/18/96 - 1) exoand MAxOBS from 5000 to 10000. rlr.
C                  2) put limits on ntimebins to prevent breakt from
C                     overflowing.
C                  3) Put in sPoisson fix for day wrap-around.
C v1.9  01/02/97 - Drop records with residuals > maxres & < minres, rlr.
C v1.9a 06/20/97 - Drop record at read rather than write, so that stats
C         are reasonable. Required reading parameters before data. rlr.
C v2.0  07/xx/07 - CRD version. Combines slr and llr versions. rlr.
C v2.1  05/20/08 - Fix nu,ber of input and output points reported;
C                - H4 output to reports are from standard variables.
C                - cleaned up source with fpret. rlr
C v2.2  07/02/08 - CRD v 1.00. rlr.
C v2.2a 01/07/10 - Add version # to output file. rlr.
C
      IMPLICIT none
      INTEGER*2 maxobs
      PARAMETER (maxobs=10000)
C
      CHARACTER*256 comname,inname,listname,sumname,outname
C
      REAL*8 arr x(maxobs),arr y(maxobs),breakt(5),binlength,binwidth,
     +    calsdavg,maxslope,minslope,mint,maxt,minyd,maxyd,minres,maxres
C
      INTEGER*2 minpts,ntot,psens,targid
C
      INTEGER*4 nsignal,ioerr
C
      INTEGER len0,len,trimlen
C
      LOGICAL selected(maxobs),warn
C
      COMMON /tbcmn/ debug
      LOGICAL debug
C
      CHARACTER*80 sccs_poisson
C
      COMMON /versions/ version
      CHARACTER*20 version
C
      DATA comname /'data/lib/filter_crd.pc'/
      DATA sumname /'data/tmp/trash'/
      DATA warn /.false./
      DATA debug /.false./
      DATA breakt /5*1.e9/
      DATA version /"2.2a"/
C
      DATA sccs_poisson /'%W%\t%G%'/
C
C get the needed file names
C
      CALL getnames(comname,inname,outname,listname,sumname)
C
C open listing file named input file name+'.lst'
C
      OPEN(2,file=listname,status='unknown',access='append',
     +    iostat=ioerr,err=900)
      OPEN(5,file=sumname,status='unknown',access='append',iostat=ioerr,
     +    err=901)
      len0= trimlen(version)
      len= trimlen(inname)
      WRITE(2,1000) version(1:len0),inname(1:len)
 1000 FORMAT(25x,
     +    'Poisson Filtering (CRD)   v',a,//'Input data file name: ',a/)
C
C get the time/residual set to work with
C
      minres= -1.d+100
      maxres= 1.d+100
      CALL readres(inname,targid,arrx,arry,selected,ntot,calsdavg,
     +    minres,maxres)
C
C default control parameter names.  See timebin for description of parameters
C
      binlength= 180.
      binwidth= 0.5
      minslope= -0.02
      maxslope= +0.02
      minres= -10000
      maxres= 10000
C
C now, get the parameters from file, if available
C
      CALL readparam(comname,targid,binlength,binwidth,minslope,
     +    maxslope,minyd,maxyd,psens,minpts,mint,maxt,minres,maxres)
C
C C     write(*,*) "calsd",calsdavg
C
C   toss out-of-range points (esp. needed for llr)
C C     write(*,*) "pre-oor",ntot
C       i= 1
C       do ii=1,ntot
C         if ((arry(i).gt.maxres .or. arry(i).lt.minres)) then
C           ntot= ntot- 1
C           do j= i, ntot
C             arrx(j)= arrx(j+1)
C             arry(j)= arry(j+1)
C           enddo
C         else
C           i= i+ 1
C         endif
C       enddo
C C     write(*,*) "post-oor",ntot
C
C   main routine: identify optimal slope and bin
C
      CALL timebin(arrx,arry,ntot,binlength,binwidth,minslope,maxslope,
     +    minyd,maxyd,psens,minpts,mint,maxt,minres,maxres,calsdavg,
     +    selected,warn,breakt)
C
C record the results on file input file name+'.out'
C
      CALL writeres(outname,targid,arrx,arry,selected,ntot,breakt,
     +    nsignal,minres,maxres)
C
C close the listing file -- all the others HAVE been closed
C
      CLOSE (2)
C
C Stop the boat..  Get off at correct dock ( :-) ).
C
      IF (nsignal.gt.0.and..not.warn) THEN
C
C ..normal return if data marked
C
        CALL exit(0)
      ELSE IF (nsignal.gt.0.and.warn) THEN
C
C ..it's ok, but should be looked at
C
        CALL exit(1)
      ELSE
C
C ..abnormal return if no data found.A
C
        CALL exit(2)
      ENDIF
C
  900 CALL error('listing file open error',ioerr)
  901 CALL error('summary file open error',ioerr)
C
      END
C
      SUBROUTINE getnames(comname,inname,outname,listname,sumname)
C
C get the input, output, and listing file names from the command line
C ot the user.
C
      CHARACTER*256 argv(10),comname,inname,listname,sumname,outname
C
      INTEGER argc,argl(12),trimlen
C
      LOGICAL exists
C
      argc= iargc()
      IF (argc.ge.6.and.argc.le.12) THEN
        DO 10 i=1,argc
          CALL getarg(i,argv(i))
   10   argl(i)= trimlen(argv(i))
C
C now get the arguments we want!
C
        i= 0
   15   i= i+ 1
        IF (argv(i)(1:2).eq.'-c'.and.argl(i).eq.2) THEN
          i= i+ 1
          comname= argv(i)(1:argl(i))
          GO TO 15
        ENDIF
        IF (argv(i)(1:2).eq.'-i'.and.argl(i).eq.2) THEN
          i= i+ 1
          inname= argv(i)(1:argl(i))
          GO TO 15
        ENDIF
        IF (argv(i)(1:2).eq.'-o'.and.argl(i).eq.2) THEN
          i= i+ 1
          outname= argv(i)(1:argl(i))
          GO TO 15
        ENDIF
        IF (argv(i)(1:2).eq.'-l'.and.argl(i).eq.2) THEN
          i= i+ 1
          listname= argv(i)(1:argl(i))
          GO TO 15
        ENDIF
        IF (argv(i)(1:2).eq.'-s'.and.argl(i).eq.2) THEN
          i= i+ 1
          sumname= argv(i)(1:argl(i))
          GO TO 15
        ENDIF
        IF (i.le.10) GO TO 15
      ELSE
   30   CALL filnme
     +      ('Command parameter file name: ','.pc',comname,exists)
        IF (.not.exists) THEN
          WRITE(*,'(''    File does not exist.  Try again.'')')
          GO TO 30
        ENDIF
        inname= 'test.lu'
   40   CALL filnme('Input data file name: ','.lu',inname,exists)
        IF (.not.exists) THEN
          WRITE(*,'(''    File does not exist.  Try again.'')')
          GO TO 40
        ENDIF
        len= trimlen(inname)
        outname= inname(1:len-2)//'l8'
        CALL filnme('Output data file name: ','.lf',outname,exists)
        listname= inname(1:len-2)//'l9'
        CALL filnme('Listing file name: ','.ll',listname,exists)
        sumname= inname(1:len-2)//'l0'
        CALL filnme('Summary file name: ','.ll',sumname,exists)
      ENDIF
C
      RETURN
      END
      SUBROUTINE readres(inname,targid,arrx,arry,selected,ntot,calsdavg,
     +    minres,maxres)
C
C READRES
C Author: R. L. Ricklefs, Univ of Texas McDonald Observatory
C read time / residual data from file.
C input
C     none
C output
C     inname  - input data file name
C     targid  - SIC or 100+lunar reflector number (0-4)
C     arrx    - time array (sec of day)
C     arry    - residual array (nsec)
C     selected  - status array: true => ok data, false => noise
C     ntot    - total number of data points
C revisions
C
      IMPLICIT none
      INTEGER*2 maxobs
      PARAMETER (maxobs=10000)
C
      CHARACTER*256 inname
C
      CHARACTER*512 line
C
      INTEGER*2 ntot,targid
C
      INTEGER*4 day,fsec,hour,iburst,ioerr,min,mon,sec,year,qual
C
      REAL*8 arrx(ntot),arry(ntot),calsd,calsdavg,calsdsum,dx,minres,
     +    maxres
C
      LOGICAL selected(ntot)
C
      INCLUDE '../include/crd.inc'
      INCLUDE '../include/crd_MLRS.inc'
C
      COMMON /tbcmn/ debug
      LOGICAL debug
C
C open file
C
      OPEN(1,file=inname,status='old',err=98,iostat=ioerr)
C
      ntot= 1
      iburst= 0
      calsdsum= 0.d0
   20 READ(1,'(a)',end=999,err=99,iostat=ioerr) line
      IF (line(1:2).eq.'h3') THEN
        CALL read_h3 (line)
        targid= sic
C
C C     write(*,*) "read",sic,targid
C
      ENDIF
      IF (line(1:2).eq.'h4') then
        CALL read_h4 (line)
        WRITE(2,'(''Data for '',i4,2(''/'',i2),
     +1x,i2,'':'',i2/)
     +') start_year, start_mon, start_day, start_hour, start_min
      ENDIF
      IF (line(1:2).eq.'40') THEN
        iburst= iburst+ 1
        CALL read_40 (line)
        calsdsum= calsdsum+ cal_rms/1000.d0
      ENDIF
      IF (line(1:2).eq.'93') THEN
        CALL read_93 (line)
        arrx(ntot)= d93_sec_of_day
        IF (target_type.eq.1) THEN
          arry(ntot)= range_omc
        ELSE IF (target_type.eq.2) THEN
C
C Don't copy in k'only record's dummy O-C_post. Use out of ranger val
C
          IF (dabs(range_omc_post-(-1.d0))<1.d-6) THEN
            arry(ntot)= range_omc
          ELSE
            arry(ntot)= range_omc_post
          ENDIF
        ELSE
          arry(ntot)= range_omc ! for now
        ENDIF
C
C   handle wrap-around to next day
C
        IF (ntot.gt.1) THEN
          IF (arrx(ntot).lt.arrx(1)) arrx(ntot)= arrx(ntot)+ 86400.d0
        ENDIF
C
C assume it's not data:
C
        selected(ntot)= .false.
        ntot= ntot+ 1
        IF (ntot.gt.maxobs) CALL error("Too many observations!",1)
      ENDIF
C
      GO TO 20
C
  999 ntot= ntot- 1
C
C C     write(*,*) "Po ntot= ",ntot
C
      IF (iburst.gt.0) THEN
        calsdavg= calsdsum/iburst
      ELSE
        calsdavg= 0.d0
      ENDIF
C
C  if (debug)
C    write(2,'(''Cal average= '',f6.3,'' over '',i3,'' bursts''/)')
C .calsdavg,iburst
C
      RETURN
C
C error calls
C
   97 CALL error('decoding input data file record',ioerr)
   98 CALL error('opening input data file',ioerr)
   99 CALL error('reading input data file',ioerr)
C
      END
C
      SUBROUTINE readparam(comname,insic,binlength,binwidth,minslope,
     +    maxslope,minyd,maxyd,psens,minpts,mint,maxt,minres,maxres)
C
C READPARAM
C Author: R. L. Ricklefs, Univ of Texas McDonald Observatory
C
C Read control parameters for Poisson filtering run.
C
C input
C     comname - name of command parameter file
C output
C     binlength       - length of time bin in sec
C     binwidth- width of residual semi-bin in nsec
C     minslope- maximum negative slope for data
C     maxslope- maximum positive slope for data
C     minyd   - minimum acceptable range residual
C     maxyd   - maximum acceptable range residual
C     psens   - fudge factor to desensitize Poisson filter
C     minpts  - the minimum number of points in a time bin for
C       the bin to be accepted as data
C     mint    - minimum time of data used in filtering
C     maxt    - maximum time of data used in filtering
C     minres  - minimum o-c range residual used in Poisson analysis
C     maxres  - maximum o-c range residual used in Poisson analysis
C revisions
C
      IMPLICIT none
      CHARACTER*256 comname
C
      INTEGER*2 minpts,psens,insic,sic
C
      INTEGER*4 ioerr
C
      REAL*8 binlength,binwidth,maxres,maxslope,maxt,maxyd,minres,
     +    minslope,mint,minyd
C
C open file, read parameters, and close
C
      OPEN(3,file=comname,iostat=ioerr,err=900)
   10 READ(3,
     +    '(i4,6(1x,f10.0),2(1x,i4),4(1x,f10.0))',iostat=ioerr,err=905,
     +    END=20) sic,binlength,binwidth,minslope,maxslope,minyd,maxyd,
     +    psens,minpts,mint,maxt,minres,maxres
      IF (sic.eq.insic) GO TO 20
      GO TO 10
C
   20 CLOSE(3)
      WRITE(2,1000)insic,binlength,binwidth,minslope,maxslope,
     +    minyd,maxyd,psens,minpts,mint,maxt,minres,maxres
 1000 FORMAT('Satellite ID: ',i4/'Bin Length (seconds): ',f10.3/
     +    'Bin Width (nsec):     ',f10.3/
     +    'Minimum Slope (nsec/hour): ',f12.3/
     +    'Maximum Slope (nsec/hour): ',f12.3/
     +    'Minimum Acceptable Residual (nsec): ',f12.3/
     +    'Maximum Acceptable Residual (nsec): ',f12.3/
     +    'Poisson Sensitivity Reduction (points): ',i3/
     +    'Minimum Points in Time Bin (points):    ',i3/
     +    'Minimum Time Used (hours): ',f8.4/
     +    'Maximum Time Used (hours): ',f8.4/
     +    'Minimum O-C Residual (nsec): ',f12.3/
     +    'Maximum O-C Residual (nsec): ',f12.3/)
C
C change from nsec/hour to nsec/sec
C
      minslope= minslope/3600.d0
      maxslope= maxslope/3600.d0
C
C ...and from hours to seconds
C
      mint= mint*3600.d0
      maxt= maxt*3600.d0
C
      RETURN
C
C error handling
C
  900 CALL error('parameter file error on open',ioerr)
  905 CALL error('reading parameter file',ioerr)
      END
C
      SUBROUTINE writeres(outname,targid,arrx,arry,selected,ntot,breakt,
     +    nsignal,minres,maxres)
C
C  WRITERES
C  Author: R. L. Ricklefs, Univ of Texas, McDonald Observatory
C
C  Copy input data file to output adding the new status flag, indicating
C  whether the observation appears to be data or noise.
C
C  NOTE: this routine has commented-out code to write different types
C  of files used by microcomputer program to plot this data.
C
C  input
C      inname  - input data file name
C      targid  - SIC or 100+lunar reflector number
C      arrx    - observation time
C      arry    - observation residual
C      selected- observation status
C      ntot    - number of observations
C
C  output
C      nsignal - number of signal photons marked.
C
C History
C      12/xx/97 - changed to filter output data based on a min and max
C                  residual. rlr.
C
      IMPLICIT none
      CHARACTER*256 outname
C
      CHARACTER*512 hold_line,line,type1recd,type2recd
C
      INTEGER*2 hour,min,ntot,sec,targid
C
      INTEGER*4 nline,nsignal,ioerr,nbreak
C
      REAL*8 arrx(ntot),arry(ntot),breakt(5),minres,maxres,fsec,res,
     +    time3
C
      LOGICAL selected(ntot),hold
C
      INTEGER lengtht,trimlen
C
      COMMON /versions/ version
      CHARACTER*20 version
C
      INCLUDE '../include/crd.inc'
      INCLUDE '../include/crd_MLRS.inc'
C
C open the output file
C
      OPEN(4,file=outname,status='unknown',err=98,iostat=ioerr)
      REWIND (1)
C
C copy from input file to output file, replacing the 'quality' (filter) flag
C
CC      nline= 1
      nline= 0
      nbreak= 1
      nsignal= 0
      hold= .false.
   20 READ(1,'(a)',end=999,err=99,iostat=ioerr) line
C
C   range record -- assumes each is paired with a '93' filter record.
C
      IF (line(1:2).eq.'10') THEN
        CALL read_10 (line)
C
C          This assumed a '10' is always followed by a '93'
C
        hold= .true.
      ENDIF
      IF (line(1:2).eq.'93') THEN
        CALL read_93 (line)
        IF (target_type.eq.1) THEN
          res= range_omc
        ELSE IF (target_type.eq.2) THEN
C
C Don't copy in k'only record's dummy O-C_post. Use out of ranger val
C
          IF (dabs(range_omc_post-(-1.d0))<1.d-6) THEN
            res= range_omc
          ELSE
            res= range_omc_post
          ENDIF
        ELSE
          res= range_omc ! for now
        ENDIF
C
C ...Tag residuals that are out of bounds (same done in readres)
C         Don't assume '93' has accompanying '10'.
C C       if (dabs(res) >= 1.e-30) then
C         NOW write the '10' record.
C
        IF (hold) THEN
C
C C     write(*,*) nline, selected(nline),res
C
          IF (selected(nline+1)) THEN
            filter_flag= 2
            nsignal= nsignal+ 1
          ELSE
            filter_flag= 1
          ENDIF
          IF ((res.gt.maxres.or.res.lt.minres)) filter_flag= 3
          CALL write_10 (hold_line)
          lengtht= trimlen(hold_line)
          WRITE(4,'(a)') hold_line(1:lengtht)
          hold= .false.
C
C Keep nline in sync. Out of bound records are not counted in nline
C if (filter_flag .lt. 3) nline= nline+ 1
C
        ENDIF
        nline= nline+ 1
C
C ....see if there needs to be a time break here
C C       time3= hour*3600.d0+ min*60.d0+ sec+ fsec/1.d7
C
C ....if so, create marked run header & subheader here for normalpoint
C C       if (time3 .gt. breakt(nbreak)) then
C C         nbreak= nbreak+ 1
C C         type1recd(3:16)= line(3:16)
C C         type1recd(68:68)= 'F'
C C         write(4,'(a)') type1recd
C C         type2recd(3:16)= line(3:16)
C C         type2recd(68:68)= 'F'
C C         write(4,'(a)') type2recd
C C       endif
C
      ELSE IF (line(1:2).eq.'95') THEN
              CALL read_95 (line)
              psnvers= version
              CALL write_95 (line)
              lengtht= trimlen(line)
              WRITE(95,'(a)') line(1:lengtht)
      ENDIF
C
      IF (.not.hold) THEN
        lengtht= trimlen(line)
        WRITE(4,'(a)') line(1:lengtht)
      ENDIF
      GO TO 20
C
C finish up
C
  999 CLOSE(1)
      CLOSE(4)
      WRITE(2,1000) nline, nsignal
      WRITE(5,1000) nline, nsignal
 1000 FORMAT(/"Poisson Filtering: Total obs=",i6,",  Selected obs=",i6)
      RETURN
C
C error calls
C
   98 CALL error('opening output file',ioerr)
   99 CALL error('writing output file (re-reading input)',ioerr)
C
      END
C
      SUBROUTINE timebin(arrx,arry,ntot,binlength,binwidth,minslope,
     +    maxslope,minyd,maxyd,psens,minpts,mint,maxt,minres,maxres,
     +    calsdavg,selected,warn,breakt)
C
C TIME BIN
C Author: R. Ricklefs
C
C Poisson analysis subsystem
C Purpose:  Do Poisson statistics on every time bin of the proper length
C
C input:
C     arrx    - observation time in sec of day (usually)
C     arry    - observation residual in nsec (usually)
C     ntot    - number of observations
C     binlength       - length of time bin in units of arrx
C     binwidth- width of residual semi-bin in units of arry
C     minslope- max negative slope for data in units of arry/arrx
C     maxslope- max positive slope for data in units of arry/arrx
C     minyd   - minimum acceptable range residual in units of arry
C     maxyd   - maximum acceptable range residual in units of arry
C     psens   - fudge factor to desensitize Poisson filter in obs.
C     minpts  - minimum points in a time bin
C     mint    - minimum time used in filtering
C     maxt    - maximum time used in filtering
C     calsdavg- average cal standard deviation (nsec)
C
C local:
C     nintbin - number of points in current time bin
C
C output:
C     selected- observation status (true=> data, false=>noise)
C     warn    - warn the user of a potential problem
C     breakt  - time breaks to be put into lunar data file
C
C Revisions:
C     08/10/90- add mint and maxt limits on data used;
C       handle case of multiple max bin with same
C       number of points.
C
      IMPLICIT none
      INTEGER*2 anompts,i,maxexp,minpts,nintbin,nlow,ntot,optbin,psens
C .,       ntimebin
      INTEGER*4 ntimebin
C
      REAL*8 arrx(ntot),arry(ntot),binlength,binwidth,breakt(5),
     +    calsdavg,eqvres,hourhigh,hourlow,maxslope,maxt,maxx,maxy,
     +    maxyd,minslope,mint,minx,miny,minyd,optslope,refbinwidth,sd,
     +    snratio,xlow,minres,maxres
C
      LOGICAL selected(ntot),tieres,tieslp,warn,widenbin
C
      COMMON /tbcmn/ debug
      LOGICAL debug
C
      refbinwidth= binwidth
C
C Find first point to use
C
      nlow= ntot+ 1     ! 052008
    1 DO 5 i=1,ntot
        IF (arrx(i).gt.mint) THEN
          nlow= i
          GO TO 6
        ENDIF
    5 CONTINUE
    6 nintbin= 0
      xlow= arrx(nlow)- binlength
      ntimebin= 0
C
C find the point at which to start the next time span
C
   10 nlow= nlow+ nintbin
C
CC      write(*,*) "nlow,ntot,arrx,maxt",nlow,ntot,arrx(nlow),maxt
C   none: quit
C
      IF (nlow.ge.ntot.or.arrx(nlow).gt.maxt) GO TO 999
      ntimebin= ntimebin+ 1
C
C C?    if (ntimebin .gt. 6) then
C C            call error('Too many time bins (>6)',ntimebin)
C C     endif
C
C   set start value of x array (time) bin
C
      xlow= xlow+ binlength
      IF (ntimebin.gt.1) breakt(ntimebin-1)= xlow
C
C get all points within the time bin (xlow through xlow+ bin length)
C
      DO 50 i=nlow,ntot
        IF (arrx(i)-xlow.gt.binlength.or.arrx(i).gt.maxt) THEN
          nintbin= i- nlow
          GO TO 60
        ENDIF
   50 CONTINUE
C
C C     write (*,*) xlow,ntot,nlow,nintbin
C   here for last bin (end of data)
C   calculate number of points in bin
C C     nintbin= ntot- nlow+ 1
C
      nintbin= ntot- nlow+ 1
C
C C     write(*,*) nintbin
C
C   nothing in this group? [increment time bin count?]
C
   60 IF (nintbin.le.1) GO TO 10
C
C test output
C
      IF (debug) THEN
        hourlow= arrx(nlow)/3600.d0
        hourhigh= arrx(nlow+nintbin-1)/3600.d0
        WRITE(2,1000) hourlow, hourhigh, nintbin
 1000 FORMAT(/'==> Time bin: ',f8.3,'UTC to ',f8.3,'UTC with ',i4,
     +      ' points.')
      ENDIF
C
C iterate through range of slopes minslope to maxslope
C
      CALL chgslope(arrx(nlow),arry(nlow),nintbin,binwidth,minslope,
     +    maxslope,minx,maxx,miny,maxy,minyd,maxyd,psens,minres,maxres,
     +    maxexp,anompts,optbin,optslope,eqvres,sd,tieslp,tieres,
     +    widenbin)
C
C C        write (*,*) "from changelsope:",
C C     .          nlow,arrx(nlow),arry(nlow),nintbin,binwidth,
C C     .          minslope,maxslope,minx,maxx,miny,maxy,
C C     .          minyd,maxyd,psens,maxexp,
C C     .          anompts,optbin,optslope,eqvres,
C C     .          sd,tieslp,tieres,widenbin
C
C   now what do we do with the results???
C   Tag data in selected bin.
C
      IF (anompts.gt.maxexp.and.anompts.ge.minpts) THEN
        CALL tagdata(arrx(nlow),arry(nlow),nintbin,minx,maxx,miny,maxy,
     +      minres,maxres,optslope,binwidth,optbin,selected(nlow),
     +      snratio)
        WRITE(2,1010) ntimebin,nintbin,maxexp,anompts,calsdavg,
     +      binwidth,optslope*3600.,eqvres,sd,snratio
        WRITE(5,1010) ntimebin,nintbin,maxexp,anompts,calsdavg,
     +      binwidth,optslope*3600.,eqvres,sd,snratio
 1010 FORMAT(/'Poisson: Time bin: ',i5,', Points: total ',i4,', ',
     +      'expected: ',i3,', selected: ',i3/'Poisson: Cal SD: ',f5.3,
     +      ' nsec, Bin width: ',f5.2,' nsec, '/'Poisson: Slope:',f9.2,
     +      ' nsec/hour, ','O-C: ',f7.2,' nsec, SD: ',f5.3,
     +      ' nsec, S:N: ',f4.1)
        IF ((anompts-maxexp).le.2) THEN
          WRITE(2,1015)
          WRITE(5,1015)
 1015 FORMAT('Poisson: ==> Data is weak. Number of points is ',
     +        'near noise level. Check results.')
          warn= .true.
        ENDIF
        IF (snratio.lt.1.0) THEN
          WRITE(2,1020)
          WRITE(5,1020)
 1020 FORMAT('Poisson: ==> Data is weak. Signal to noise ratio ',
     +        'is quite low. Check results')
          warn= .true.
        ENDIF
      ELSE
        WRITE(2,1025) ntimebin,nintbin,maxexp
        WRITE(5,1025) ntimebin,nintbin,maxexp
 1025 FORMAT("Poisson: Time bin ",i5,", Total points: ",i4,", ",
     +      "Selected points: NONE"/"Poisson: Expected points: ",i5)
        IF (anompts.gt.maxexp.and.anompts.lt.minpts) THEN
          WRITE(2,1030)
          WRITE(5,1030)
 1030 FORMAT('Poisson: ==> Too few points in significant bin')
        ENDIF
      ENDIF
C
      IF (tieres.and..not.tieslp) THEN
        WRITE(2,1035)
        WRITE(5,1035)
 1035 FORMAT('Poisson: ==> Tie in significant residuals - ',
     +      'please check the results.')
        warn= .true.
      ELSEIF (tieslp.and..not.tieres) THEN
        WRITE(2,1040)
        WRITE(5,1040)
 1040 FORMAT('Poisson: ==> Tie in significant slopes - ',
     +      'please check the results.')
        warn= .true.
      ELSEIF (tieres.and.tieslp) THEN
        WRITE(2,1045)
        WRITE(5,1045)
 1045 FORMAT('Poisson: ==> Tie in significant residuals and slopes ',
     +      '- please check the results.')
        warn= .true.
      ENDIF
C
C  If there is a need to widen the bin, but only once...
C      if (widenbin .and. refbinwidth.eq.binwidth) then
C        binwidth= binwidth*1.5d0
C        write(2,1022) binwidth
C 1022   format(/'Bin width too narrow - widening to ',f5.2,' nsec'/)
C        go to 1
C      endif
C
C  Warn of need to widen bin
C
      IF (widenbin) THEN
        WRITE(2,1050)
C
C C       write(5,1050)
C
 1050 FORMAT('Poisson: ==> Bin width appears to be too small - ',
     +      'please check the results.')
        warn= .true.
      ENDIF
C
C repeat for next time bin
C
      GO TO 10
C
  999 IF (ntimebin.gt.1) THEN
        WRITE(2,1060) ntimebin
        WRITE(5,1060) ntimebin
 1060 FORMAT('Poisson: ==> More than one time bin (',i3,
     +      '). Is the break optimal?')
        warn= .true.
      ENDIF
      WRITE(2,'(/)')
      WRITE(5,'(/)')
C
      RETURN
      END
C
      SUBROUTINE chgslope(arrx,arry,n,binwidth,minslope,maxslope,minx,
     +    maxx,miny,maxy,minyd,maxyd,psens,minres,maxres,maxexp,anompts,
     +    optbin,optslope,optres,optsd,tieslp,tieres,widenbin)
C
C CHG-SLOPE
C Author: R. Ricklefs, Univ of Texas McDonald Observatory
C
C Poisson analysis subsystem
C Purpose:  Vary the slope of the data to find the best data bin.
C
C input:
C     arrx    - observation time in sec of day (usually)
C     arry    - observation residual in nsec (usually)
C     n       - number of observations
C     binwidth- width of residual semi-bin in units of arry
C     minslope- maximum negative slope for data in units of arry/arrx
C     maxslope- maximum positive slope for data in units of arry/arrx
C     minyd   - minimum acceptable range residual in units of arry
C     maxyd   - maximum acceptable range residual in units of arry
C     psens   - fudge factor to desensitize Poisson filter in obs.
C
C output:
C     minx    - minimum observation time
C     maxx    - maximum observation time
C     miny    - minimum observation residual
C     maxy    - maximium observation residual
C     maxexp  - maximum number of observations expected in a bin
C     anompts - maximum number of points for any bin exceeding
C       Poisson limits.
C     optbin  - bin exceeding Poisson limit by greatest amount
C     optslope- slope associated with optbin
C     optres  - equivalent O-C residual of data (i.e., at center of
C       time bin)
C     optsd   - standard deviation of the chosen data in nsec.
C     tieslp  - .true. if there is a tie in slope
C     tieres  - .true. if there is a tie in residual bins
C     widenbin- 3 or more significant adjacent bins found at best
C       slope.  Widen the binwidth variable.
C
C Revisions:
C
      IMPLICIT none
C
      INTEGER*2 maxcan
      PARAMETER (maxcan=500)
C
      INTEGER nrej
C
      INTEGER*2 anompts,cananom(maxcan),canbin(50,maxcan),
     +    canexp(maxcan),histarr(10000),i,j,k,m,mo2,maxanom,maxbin(50),
     +    maxexp,maxpts,n,nbest,nbins,ncanbin(maxcan),
     +    nconsbin,nconscan(maxcan),nmax,npickres,npicksl,nslopes,
     +    optbin,psens,bestsl(100),np
C
      REAL*8 arrx(n),arry(n),avgr,avgt,
     +    binwidth,canslope(maxcan),canres(50,maxcan),eqvres,maxslope,
     +    maxx,maxy,maxyd,minslope,minx,miny,miny0,minyd,optres,
     +    optslope,savr,optsd,savt,sd,slope,slopeincr,minres,maxres,
     +    totslopes,tmpx
C
      LOGICAL tieres,tieslp,widenbin
C
      COMMON /tbcmn/ debug
      LOGICAL debug
C
C find min and max in resid & time
C
      minx= 1.d30
      maxx= -1.d30
      miny= 1.d30
      maxy= -1.d30
      np= 0
      anompts= 0
      maxexp= 0
      DO 10 i=1,n
        IF (arry(i).gt.minres.and.arry(i).lt.maxres) THEN
          IF (arrx(i).lt.minx) minx= arrx(i)
          IF (arrx(i).gt.maxx) maxx= arrx(i)
C
C C         if (arry(i).lt.miny) miny= arry(i)
C C         if (arry(i).gt.maxy) maxy= arry(i)
C   Take as minimum if it is not out of bounds (i.e. beyond +/- 1.d5 nsec)
C
          IF (arry(i).lt.miny.and.arry(i).gt.-1.d5) miny= arry(i)
          IF (arry(i).gt.maxy.and.arry(i).lt.1.d5) maxy= arry(i)
          np= np+ 1
        ENDIF
   10 CONTINUE
C
C make sure we have meaningful number of bins even in absence of noise
C
CC        write(*,*) "prelim miny, maxy, n",miny,maxy,n
      IF ((maxy-miny).lt.10*binwidth) THEN
        miny= dmin1(miny,maxy-5*binwidth)
        maxy= dmax1(maxy,miny+5*binwidth)
      ENDIF
CC        write(*,*) "final miny, maxy",miny,maxy
C
C initialize increment for slope
C
      if (maxx .lt. minx) THEN
        tmpx= minx
        minx= maxx
        maxx= tmpx
      endif
      totslopes= (maxx-minx)/(2*binwidth)
CC	write(*,*) "nts: ",totslopes
      IF (totslopes .gt. 10000) THEN	! 06/09/08
        WRITE(2,*) 'Too many slopes (>10000)',totslopes
        RETURN
CC        CALL error('Too many slopes (>10000)',totslopes)
      ENDIF
      IF (totslopes .lt. 1) THEN	! 06/09/08
      WRITE(2,*) 'Too few slopes (<1)',totslopes
      RETURN
CC        CALL error('Too few slopes (<1)',totslopes)
      ENDIF
      slopeincr= 2*binwidth/(maxx-minx)
CC        write(*,*) "Maxx, minx, slopeincr",maxx,minx,slopeincr
      nslopes= 1
      DO 15 i=1,maxcan
   15 ncanbin(i)= 0
      tieslp= .false.
      tieres= .false.
      widenbin= .false.
C
C increment slope by one bin width at each end of time swath
C  Be sure hist arr doesn't overflow!!
C
      DO 50 slope= minslope, maxslope, slopeincr
        IF (debug) WRITE(2,'(''Slope='',f10.5)') slope
C
C bin the data into bins of "bin width" in residual space (y)
C
        CALL bindata(arrx,arry,n,minx,maxx,miny,maxy,slope,histarr,
     +      binwidth,miny0,m)
C
C if there are enough bins, compare number per bin against maximum expected
C by Poisson statistics
C
        IF (m.ge.3) THEN
          mo2= m/2
CC          nmax= maxpts(n,mo2)+ psens
          nmax= maxpts(np,mo2)+ psens
          IF (debug) WRITE(2,1010) mo2,nmax
 1010 FORMAT('For ',i4,' bins, max expected in any bin=',i6)
C
C find most outrageous bins for this slope
C
          CALL findanom(histarr,m,nmax,maxbin,nbins,maxanom,nconsbin)
C
C copy to storage every bin of this slope which has the highest known
C number of anomalies
C
          IF (maxanom.ge.anompts-1) THEN
C
C if (maxanom .ge. anompts) then
C
            DO 20 i=1,nbins
C
C find the residual at 1/2 the bin time span
C
              eqvres= miny0+ (maxbin(i)+1)*binwidth+ slope*(maxx-minx)/2
              IF (eqvres.ge.minyd.and.eqvres.le.maxyd) THEN
                anompts= maxanom
                cananom(nslopes)= maxanom
                canslope(nslopes)= slope
                canexp(nslopes)= nmax
                ncanbin(nslopes)= ncanbin(nslopes)+ 1
                canbin(ncanbin(nslopes),nslopes)= maxbin(i)
                canres(ncanbin(nslopes),nslopes)= eqvres
                nconscan(nslopes)= nconsbin
              ENDIF
   20       CONTINUE
            IF (ncanbin(nslopes).gt.0) THEN
              nslopes= nslopes+ 1
              IF (nslopes.gt.maxcan) CALL
     +            error('Too many slopes (>500)',nslopes)
            ENDIF
          ENDIF
        ENDIF
   50 CONTINUE
C
C Now find the best bin of the best slope
C
      nslopes= nslopes- 1
      nbest= 0
      bestsl(1)= -1     ! Just in case... 052008
C
C (do loop will not run if no anonamlies found)
C
      DO 100 i=1,nslopes
        IF (cananom(i).ge.(anompts-1)) THEN
C
C        if (cananom(i).eq.anompts) then
C          if (debug) write(2,1015) i, canslope(i),
C     .        (canbin(j,i),canres(j,i),j=1,ncanbin(i))
C 1015     format(i3, f9.6, 3x, 5(:,i4,f9.4))
C
          nbest= nbest+ 1
          bestsl(nbest)= i
        ENDIF
  100 CONTINUE
C
C Is there an unambiguos winner?
C
      IF (nbest.eq.1.and.ncanbin(bestsl(1)).eq.1) THEN
C
C Yes...
C
        k= bestsl(1)
        maxexp= canexp(k)
        optbin= canbin(1,k)
        optslope= canslope(k)
        optres= canres(1,k)
        IF (nconscan(k).gt.2) widenbin= .true.
        CALL avgdata(arrx,arry,n,minx,maxx,miny,maxy,canslope(k),
     +      binwidth,canbin(1,k),avgt,avgr,sd,nrej)
        optsd= sd
        IF (debug) WRITE(2,1030) avgt,avgr,sd,nrej
 1030 FORMAT('Average time, residual, std dev, and nrej:',2f10.2,f7.4,
     +      i3)
      ELSE
C
C No... Don't pick anything except a representative 'max expected'
C
        anompts= 0

C Protect yourself in case there are no slopes.
        IF (bestsl(1) .ge. 1 .and. bestsl(1) .lt. maxcan) THEN
          maxexp= canexp(bestsl(1))
         ELSE
          maxexp= 0
         ENDIF
      ENDIF
C
C warn user about ties
C
      IF (nbest.gt.1) tieslp= .true.
      npicksl= 1
      npickres= 1
      DO 110 i=1,nbest
        IF (ncanbin(bestsl(i)).gt.1) tieres= .true.
  110 CONTINUE
C
C  try to resolve tie on multiple slopes
C      if (tieslp .and. .not.tieres) then
C        npicksl= 1
C        npickres= 1
C        do 150 i=1,nbest
C          k= bestsl(i)
C          call avgdata(arrx,arry,n,
C      minx,maxx,miny,maxy,canslope(k),
C     .binwidth,canbin(1,k),avgt,avgr,sd,nrej)
C          if (debug)
C     .      write(2,1030) avgt,avgr,sd,nrej
C          if (i.eq.1) then
C            savt= avgt
C            savr= avgr
C            optsd= sd
C          else
C  ...can't deal with great data selection differences at different slopes
C            if (dabs(savt-avgt).gt.10 .or.
C     .dabs(savr-avgr).gt.binwidth/3.d0) go to 999
C            if (dabs(savt-avgt).le.10 .and.
C     .dabs(savr-avgr).le.binwidth/3.d0) go to 150
C            if (sd.lt.optsd .and.
C     .dabs(savr-avgr).le.binwidth/3.d0) then
C      npicksl= i
C      optsd= sd
C            endif
C          endif
C 150    continue
C      elseif (tieslp .or. tieres) then
C
      IF (tieslp.or.tieres) THEN
        npicksl= 1
        DO 200 i=1,nbest
          k= bestsl(i)
          npickres= 1
          DO 190 j=1,ncanbin(k)
            CALL avgdata(arrx,arry,n,minx,maxx,miny,maxy,canslope(k),
     +          binwidth,canbin(j,k),avgt,avgr,sd,nrej)
            IF (debug) WRITE(2,1080) canslope(k),canbin(j,k),
     +          avgt,avgr,sd,nrej
 1080 FORMAT('Slope, bin, avg time, resid, SD & nrej:',f7.4,i4,2f9.2,
     +          f7.4,i3)
            IF (i.eq.1.and.j.eq.1) THEN
              optsd= sd
            ELSE
C
C ...can't deal with great data selection differences at different slopes
C             if (dabs(savt-avgt).le.10 .and.
C    .  dabs(savr-avgr).le.binwidth/3.d0) go to 150
C             if (sd.lt.optsd .and.
C    .  dabs(savr-avgr).le.binwidth/3.d0) then
C
              IF (sd.lt.optsd.and.nrej.le.1) THEN
                npicksl= i
                npickres= j
                optsd= sd
              ENDIF
            ENDIF
  190     CONTINUE
  200   CONTINUE
      ENDIF
C
C ...they were are essentially the same data.  Take one with smallest sd.
C
      IF (tieslp.or.tieres) THEN
        IF (debug) WRITE(2,'(''Candidate picked '',2i3)
     +') npicksl,npickres
        k= bestsl(npicksl)
        anompts= cananom(k)
        maxexp= canexp(k)
        optbin= canbin(npickres,k)
        optslope= canslope(k)
        optres= canres(npickres,k)
        IF (nconscan(k).gt.2) widenbin= .true.
      ENDIF
C
  999 RETURN
      END
C
      SUBROUTINE bindata(arrx,arry,n,minx,maxx,miny,maxy,slope,histarr,
     +    binwidth,miny0,m)
C
C BINDATA
C Author: R. L. Ricklefs, Univ of Texas, McDonald Observatory
C
C Bin the data by incrementing histogram cells corresponding to the
C residual of each observation residual.  This residual is not necessarily
C the residual in arry, but a residual based on tilting the data by the
C slope argument.  Note that tilting a 'square' residual space by a slope
C will result in bins on the top and bottom of the residual space being
C triangular and short.  These cannot be included in the analysis as they
C will bias the Poisson analysis toward accepting data peaks that are really
C not meaningful.
C Input
C     arrx    - observation time in sec of day (usually)
C     arry    - observation residual in nsec (usually)
C     n       - number of observations
C     minx    - minimum observation time
C     maxx    - maximum observation time
C     miny    - minimum observation residual
C     maxy    - maximium observation residual
C     slope   - slope of data in units of arry/arrx
C     binwidth- width of residual semi-bin in units of arry
C
C Output
C     histarr - histogram array for residual space
C     miny0   - minimum range residual corrected for tilt,
C       in units of arry
C     m       - numbers of bins used corrected for tilt
C NOTE: m here is based on bins of size bin width.  This gives m twice
C the size of m in the maxpts routine.  The reason for this is that we
C look at the data 2 bins at a time.  In other workds, the bins are of
C size 2*binwidth.  This allows the bins to overlap, removing edge effects
C from data straddling 2 bins.  (If data straddles more than 2 bins,
C the iteration through many slopes should still allow all the data to be
C recognized.
C
C Revisions:
C
      IMPLICIT none
      INTEGER*2 binno,histarr(10000),i,m,n
C
      REAL*8 binwidth,slope,arrx(n),arry(n),maxx,maxy,maxy0,minx,miny,
     +    minyx,miny0
C
C When we have non-zero slope, the data in the corners needs to be ignored
C or it will skew the results.  So calculate the min & max residuals for
C the most extreme full-size bins.
C
      IF (slope.gt.0.d0) THEN
        maxy0= maxy- (maxx-minx)*slope
        miny0= miny
      ELSE
        maxy0= maxy
        miny0= miny- (maxx-minx)*slope
      ENDIF
C
C calculate m, the number of bins
C
      m= idint((maxy0-miny0)/binwidth)+ 1
C
C are there too few or to many?
C
      IF (m.lt.3) RETURN
C      IF (m.gt.2500) CALL error('More than 10000 bins!',m) fix 051908
      IF (m.gt.10000) CALL error('More than 10000 bins!',m)
C
C init histogram array
C
      DO 10 i=1,m
   10 histarr(i)= 0
C
C put each observation into proper histogram bin
C
      DO 100 i=1,n
        minyx= miny0+ (arrx(i)-minx)*slope
        binno= idint((arry(i)-minyx)/binwidth)+ 1
        IF (binno.ge.1.and.binno.le.m) histarr(binno)= histarr(binno)+ 1
  100 CONTINUE
 1111 FORMAT(10i6)
C
      RETURN
      END
C
      INTEGER*2 function maxpts(n,m)
C
C MAXPTS
C Author: R. L. Ricklefs      09/29/88
C
C Poisson analysis subsystem
C Find the largest number of points that is expected in any one bin given that
C there are n points divided among m bins.
C
C input:
C     n       - number of observations
C     m       - number of bins into which observations can be divided
C
C output:
C     maxpts  - maximum number of points expected in any one of those bins
C
C References:
C     Numerical Recipes - The Art of Scientific Computing by
C     W. H. Press, et. al.
C       pp 206-208
C Revisions:
C
      IMPLICIT none
      INTEGER*2 m,n,nbins,npts
C
      REAL oprob,prob,rlambda
C
      rlambda= float(n)/float(m)
      prob= exp(-rlambda)
      oprob= -1.e6
      npts= 0
C
   10 nbins= int(prob*m+0.5)
      IF (nbins.gt.0.or.prob.gt.oprob) THEN
        npts= npts+ 1
        oprob= prob
        prob= prob*rlambda/npts
      ELSE
        maxpts= npts- 1
        RETURN
      ENDIF
      GO TO 10
C
      END
C
      SUBROUTINE findanom(histarr,m,nmax,maxbin,nbins,maxanom,nconsbin)
C
C  FIND ANOM
C Author: R. Ricklefs
C
C  Poisson analysis subsystem
C  Purpose:  To find data bins, taken 2 at a time, which
C  exceed the expected maximum number per bin computed from Poisson
C  statistics.
C  input:
C      hist arr- array of observation histogram
C      m       - number of bins into which observations can be divided
C      n max   - maximum number of points expected in any one of those bins
C
C  output:
C      maxbin  - bins have greatest excursion from Poisson statistics
C      nbins   - number of bins used in maxbin
C      maxanom - number of observations in optbin.
C      nconsbin- number of consecutive significant bins
C
C  internal:
C      anom bins- repository for anomalies: col 1 contains bin #;
C         col 2 the bin total
C      n anom  - number of anomalies found
C
C  Revisions:
C
      IMPLICIT none
      INTEGER*2 anombin(2500,2),histarr(10000),i,m,mpairs,maxanom,nanom,
     +    nbins,nconsbin,nconscan,nmax,nsum,maxbin(50)
C
      COMMON /tbcmn/ debug
      LOGICAL debug
C
      nanom= 0
      nbins= 0
      mpairs= m- 1
C
C taking bins 2 at a time, save a list of those having more obs that expected
C from Poisson statistics.
C
      DO 100 i= 1,mpairs
        nsum= histarr(i)+ histarr(i+1)
        IF (nsum.gt.nmax) THEN
          nanom= nanom+ 1
          IF (nanom.gt.2500) CALL
     +        error('Findanom - too many possible data bins (>2500)',i)
          anombin(nanom,1)= i
          anombin(nanom,2)= nsum
        ENDIF
  100 CONTINUE
C
C pick the max bins.
C
      maxanom= 0
      DO 150 i=1,nanom
        IF (anombin(i,2).gt.maxanom) THEN
          maxanom= anombin(i,2)
C
C maxbin= anombin(i,1)
C
        ENDIF
  150 CONTINUE
C
C Pick all bins with the max number of points
C
      DO 160 i=1,nanom
        IF (anombin(i,2).eq.maxanom) THEN
          nbins= nbins+ 1
          IF (nbins.gt.50) CALL
     +        error('Findanom - too many possible anomaly bins (>50)',i)
          maxbin(nbins)= anombin(i,1)
        ENDIF
  160 CONTINUE
C
C find the number of consecutive significant bins
C
      nconsbin= 1
      nconscan= 1
      DO 170 i=2,nanom
        IF (anombin(i,1).eq.(anombin(i-1,1)+1)) THEN
          nconscan= nconscan+ 1
        ELSEIF (nconscan.gt.nconsbin) THEN
          nconsbin= nconscan
          nconscan= 1
        ELSE
          nconscan= 1
        ENDIF
C
C ...inelegant way to handle last bin
C
        IF (i.eq.nanom.and.nconscan.gt.nconsbin) nconsbin= nconscan
  170 CONTINUE
C
C print the results
C
      IF (nanom.gt.0.and.debug) THEN
        WRITE(2,1000) (anombin(i,1),i=1,nanom)
 1000 FORMAT('Bins with possible data:'/5(:,'  Bin no: ',10i5/))
        WRITE(2,1010) (anombin(i,2),i=1,nanom)
 1010 FORMAT(5(:,'  Points: ',10i5/))
        WRITE(2,1020) nconsbin
 1020 FORMAT('Number of consecutive significant bins: ',i4)
        WRITE(2,1025) maxanom,(maxbin(i),i=1,nbins)
 1025 FORMAT('==> Most likely data bins (contain ',i5,
     +      ' points): ',5(:,10i5/))
      ENDIF
C
      RETURN
      END
C
      SUBROUTINE tagdata(arrx,arry,n,minx,maxx,miny,maxy,minres,maxres,
     +    slope,binwidth,optbin,select,snratio)
C
C TAGDATA
C Author: R. L. Ricklefs, Univ of Texas McDonald Observatory
C Use everything learned to tag observations as either 'data' or 'noise'.
C Then form the signal:noise ratio for this group of data.
C Routine is very similar to BINDATA.
C Input
C     arrx    - observation time in sec of day (usually)
C     arry    - observation residual in nsec (usually)
C     n       - number of observations
C     minx    - minimum observation time
C     maxx    - maximum observation time
C     miny    - minimum observation residual
C     maxy    - maximium observation residual
C     slope   - slope of data in units of arry/arrx
C     binwidth- width of residual semi-bin in units of arry
C     optbin  - bin have greatest excursion from Poisson statistics
C     select  - observation status (true=> data, false=>noise)
C Output
C     select  - observation status (true=> data, false=>noise)
C     snratio - signal:noise ratio
C
C Revisions:
C
      IMPLICIT none
      INTEGER*2 binno,i,n,nnoise,nsignal,optbin
C
      REAL*8 binwidth,slope,arrx(n),arry(n),maxx,maxy,maxy0,minx,miny,
     +    minyx,miny0,snratio,minres,maxres
C
      LOGICAL select(n)
C
      IF (slope.gt.0.d0) THEN
        maxy0= maxy- (maxx-minx)*slope
        miny0= miny
      ELSE
        maxy0= maxy
        miny0= miny- (maxx-minx)*slope
      ENDIF
      nnoise= 0
      nsignal= 0
C
      DO 100 i=1,n
        IF (arry(i).gt.minres.and.arry(i).lt.maxres) THEN
          minyx= miny0+ (arrx(i)-minx)*slope
          binno= idint((arry(i)-minyx)/binwidth)+ 1
          IF (binno.eq.optbin.or.binno.eq.(optbin+1)) THEN
            select(i)= .true.
            nsignal= nsignal+ 1
          ELSE
            select(i)= .false.
            nnoise= nnoise+ 1
          ENDIF
        ENDIF
  100 CONTINUE
C
      IF (nnoise.gt.0.and.(maxy0-miny0-2*binwidth).ne.0) THEN
        snratio=
     +      dmax1((nsignal/2*binwidth)/(nnoise/(maxy0-miny0-2*binwidth))
     +      ,0.d0)
        snratio= dmin1(snratio,99.9d0)
      ELSE
        snratio= 99.9d0
      ENDIF
C
      RETURN
      END
C
      SUBROUTINE avgdata(arrx,arry,n,minx,maxx,miny,maxy,slope,binwidth,
     +    optbin,avgt,avgr,sd,nrej)
C
C AVGDATA
C Author: R. L. Ricklefs, Univ of Texas McDonald Observatory
C Find the true average of observations in this bin and slope.
C Routine is very similar to BINDATA.
C Input
C     arrx    - observation time in sec of day (usually)
C     arry    - observation residual in nsec (usually)
C     n       - number of observations
C     minx    - minimum observation time
C     maxx    - maximum observation time
C     miny    - minimum observation residual
C     maxy    - maximium observation residual
C     slope   - slope of data in units of arry/arrx
C     binwidth- width of residual semi-bin in units of arry
C     optbin  - bin have greatest excursion from Poisson statistics
C
C Output
C     avgt    - average of observation times
C     avgr    - average of residuals
C
C Revisions:
C
      IMPLICIT none
      INTEGER*2 maxobs
      PARAMETER (maxobs=10000)
C
      INTEGER*2 binno,i,n,optbin
C
      INTEGER nsignal,nrej
C
      REAL*8 a,b,c,avgr,avgra,avgt,binwidth,slope,arrx(n),arry(n),
     +    maxx,maxy,maxy0,minx,miny,minyx,miny0,r(maxobs),
     +    sd,se,sumr,sumt,t(maxobs)
C
      IF (slope.gt.0.d0) THEN
        maxy0= maxy- (maxx-minx)*slope
        miny0= miny
      ELSE
        maxy0= maxy
        miny0= miny- (maxx-minx)*slope
      ENDIF
      nsignal= 0
      sumt= 0.d0
      sumr= 0.d0
C
      DO 100 i=1,n
        minyx= miny0+ (arrx(i)-minx)*slope
        binno= idint((arry(i)-minyx)/binwidth)+ 1
        IF (binno.eq.optbin.or.binno.eq.(optbin+1)) THEN
          sumt= sumt+ arrx(i)
          sumr= sumr+ arry(i)
          nsignal= nsignal+ 1
          t(nsignal)= arrx(i)
          r(nsignal)= arry(i)
        ENDIF
  100 CONTINUE
C
C get the averages
C
      avgt= sumt/nsignal
      avgra= sumr/nsignal
      CALL lsq(t,r,nsignal,2,a,b,c,sd,se,nrej)
      avgr= a+ avgt*b
C
C write(2,*) a,b,c,sd,se,nrej,avgt,avgra,avgr
C
      RETURN
      END
C
      SUBROUTINE error(msg,ierr)
C
C   ERROR
C   DISPLAYS ERROR MESSAGES FROM CALLING ROUTINE AND HALTS PROGRAM
C
C $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$ rlr
C
      IMPLICIT none
      CHARACTER msg*(*)
      INTEGER*4 ierr
      INTEGER*2 ilen
      INTEGER trimlen
C
      ilen= trimlen(msg)
C
      WRITE(2,'(''Error!  '',i6,'' - '',A)') ierr,msg(1:ilen)
      WRITE(5,'(''Error!  '',i6,'' - '',A)') ierr,msg(1:ilen)
      CLOSE(1)
      CLOSE(2)
      CLOSE(4)
      CLOSE(5)
      CALL exit(3)
      END
C
      SUBROUTINE filnme(quest,ext,filnam,exists)
C
C   THIS ROUTINE PROMPTS THE OPERATOR FOR A FILE NAME, GIVING A DEFAULT NAME.
C   THIS NAME IS READ, SUFFIXED WITH AN EXTENSION IF AVAILABLE AND PACKED
C   TWO CHARACTERS/WORD FOR OUTPUT.
C   IF THE FILENAME READ IN ALREADY CARRIES AN EXTENSION, IT IS USED IN
C   PREFERENCE TO THAT SPECIFIED IN THE INPUT PARAMETERS.
C
C   INPUT:
C        QUEST  -- THE PROMPT QUESTION
C        EXT    -- CHARACTER STRING EXTENSION TO BE ADDED TO NAME.
C                  ZERO LENGTH STRING INDICATES NO EXT.
C        FILNAM -- DEFAULT FILE NAME
C
C   OUTPUT:
C        FILNAM -- NEW PACKED FILENAME.
C       EXISTS  -- STATUS OF FILE: .T.=> EXISTS
C $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$ RLR 8/81
C
      IMPLICIT none
      CHARACTER quest*(*),ext*3,filnam*(*),filtmp*256
      INTEGER lengtht,lengthe,lengthq,trimlen
      INTEGER*4 ierr
      LOGICAL exists
C
C WRITE QUESTION IF ONE WAS PASSED
C
   10 lengthq= len(quest)
      lengtht= trimlen(filnam)
      IF (lengthq.gt.0) THEN
        WRITE(*,1000) quest,filnam(1:lengtht)
 1000 FORMAT(1x,a,'<',a,'>: ')
      ENDIF
C
C Obtain, check and return operator response
C
      READ(*,1010) filtmp
 1010 FORMAT(a)
      IF (filtmp(1:2).eq.'qq') STOP
C
C ..trim off blanks
C
      lengtht= trimlen(filtmp)
      IF (lengtht.eq.0) GO TO 9999
      filtmp= filtmp(1:lengtht)
C
C ..add extension if it exists
C
      lengthe= trimlen(ext)
      IF (lengthe.gt.0) filtmp= filtmp//'.'//ext(1:lengthe)
C
C ..see if file exists and name is valid
C
      inquire(file=filtmp,iostat=ierr,exist=exists)
      IF (ierr.ne.0) THEN
        WRITE(*,1015)
 1015 FORMAT(8x,'Inappropriate')
        GO TO 10
      ENDIF
      filnam= filtmp
C
 9999 RETURN
      END
C
      INTEGER function trimlen(array)
C
C find the length of an array to the first blank (or end of array)
C
      CHARACTER array*(*)
C
      length= len(array)
      trimlen= length
      DO 10 i=1,length
        lenbk= length+1-i
        IF (array(lenbk:lenbk).eq.' ') GO TO 10
        trimlen= lenbk
        GO TO 90
   10 CONTINUE
      trimlen= 0
   90 RETURN
      END
