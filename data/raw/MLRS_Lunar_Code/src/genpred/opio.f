CC	BLOCKDATA
CC        DATA INDEV /5/, OUTDEV /6/, ierrout /0/
CC	END

	SUBROUTINE ALPHAIN(QUEST,LEGALS,INDX)
      INCLUDE 'LICENSE-BSD3.inc'

C  ALPHAIN PROMPTS THE CONSOLE OPERATOR FOR AN ALPHA-NUMERIC RESPONSE,
C  AND THEN CHECKS AND RETURNS THE OPERATOR'S RESPONSE.
C
C  ON INPUT --
C	QUEST  - (CHARACTER STRING) PROMPT QUESTION,
C	         Null string => no prompt question.
C	LEGALS - (character STRING) all valid responses;
C	INDX   - (integer) INDEX OF DEFAULT RESPONSE.
C
C  ON OUTPUT --
C	QUEST,LEGALS - UNCHANGED
C	INDX   - index of chosen response
C
C	CALLS -- SCANA
C                                                                   RWH 6/82
C$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C  Wed Aug 3, 1988 10:31:05
c  provide for null character.
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

CCC    Mon Aug 8, 1988 10:16:13,  ADD INCLUDE FOR CONSOLE I/O
	INCLUDE 'DEVCOM.INC'

	character*(*) sccsid
	parameter (sccsid = "@(#)opio.f	1.8\t04/29/91")

	character*1 null
C	parameter (null=char(0))

	character QUEST*(*),legals*(*),line*80,def*1
	character scana*1,temp*1

C  WRITE QUESTION IF ONE WAS PASSED
	null= char(0)
 10	length= len(quest)
 	if (quest(1:1) .NE. null) then
	  call questout (IOUTDEV,quest)
 	  def= legals(indx:indx)
 	  write(IOUTDEV,1000) legals,def
 1000	  format(1x,'(',a,') <',a,'>: ',$)
CC 	  write(IOUTDEV,1000) quest,legals,def
CC 1000	  format(1x,A,'(',a,') <',a,'>: ',$)
	endif

C  Obtain, check and return operator response
CC	read(INDEV,1010) line
	read(*,1010) line
 1010	format(a80)
CC 	if (line(1:2).eq.'qq') call ceror(0)
 	temp= scana(line,1,80)
	if (temp.eq.' ') go to 30
 	length= len(legals)
 	do 20 i=1,length
 	  if (temp.ne.legals(i:i)) go to 20
 	    indx= i
 	    go to 30
 20	continue
 	write(ierrout,1015)
 1015	format(1x,'Please choose one of the options presented.')
 	go to 10

 30	continue
C
C	   THE
	   END




        SUBROUTINE FILENAME(QUEST,EXT,FILNAM,EXISTS)

C  THIS ROUTINE PROMPTS THE OPERATOR FOR A FILE NAME, GIVING A DEFAULT NAME.
C  THIS NAME IS READ, SUFFIXED WITH AN EXTENSION IF AVAILABLE AND PACKED
C  TWO CHARACTERS/WORD FOR OUTPUT.
C  IF THE FILENAME READ IN ALREADY CARRIES AN EXTENSION, IT IS USED IN
C  PREFERENCE TO THAT SPECIFIED IN THE INPUT PARAMETERS.
C
C  INPUT:
C       QUEST  -- THE PROMPT QUESTION
C       EXT    -- CHARACTER STRING EXTENSION TO BE ADDED TO NAME.
C                 ZERO LENGTH STRING INDICATES NO EXT.
C       FILNAM -- DEFAULT FILE NAME 
C
C  OUTPUT:
C       FILNAM -- NEW PACKED FILENAME.
C	EXISTS  -- STATUS OF FILE: .T.=> EXISTS
C$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$ RLR 8/81
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C  Wed Aug 3, 1988 09:34:27
c  provide for null character.
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
CCC    Mon Aug 8, 1988 10:26:42,  ADD INCLUDE FOR CONSOLE I/O
	INCLUDE 'DEVCOM.INC'
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

	character*1 null
C	parameter (null=char(0))

	character quest*(*),ext*3,filnam*(*),filtmp*255
	character belchar*1, msg*80
	integer ierr,lengtht,lengthq,itrimlen
	logical exists

	null= char(0)
	belchar=char(7)
C  WRITE QUESTION IF ONE WAS PASSED
 10	lengthq= len(quest)
 	lengtht= itrimlen(filnam)
 	if (quest(1:1) .NE. null) then
	  write(IOUTDEV,1000) quest,filnam(1:lengtht)
	endif
 1000	  format(1x,A,'<',A,'>: ',$)

C  Obtain, check and return operator response
C	read(INDEV,1010) filtmp
	read(*,1010) filtmp
 1010	format(a)
CC 	if (filtmp(1:2).eq.'qq') call ceror(0)

C  ..trim off blanks
 	lengtht= itrimlen(filtmp)
CC 	if (lengtht.eq.0) go to 9999
C 	if (lengtht.eq.0) then
 	if (lengtht.eq.0 .or. lengtht.eq.255) then
           filtmp= filnam
 	   lengtht= itrimlen(filtmp)
	end if
 	filtmp= filtmp(1:lengtht)

c ...add the extension, if there is one.
	lengthq=len(ext)
	if(lengthq .NE. 0)filtmp=filtmp//'.'//ext(1:lengthq)

C  ..see if file exists and name is valid
	inquire(file=filtmp,iostat=ierr,exist=exists)
	if(ierr .NE. 0) then
	  lengthq=len(filtmp)
	  msg=' Error when using "inquire" on file '//filtmp(1:lengthq)
	  write (*,*) msg
CC	  call error(1,msg,0,ierr)
	endif
CC	if (.NOT.exists) then
CC	  write(ierrout,1015) belchar
CC 1015	  format(1x,'Named file not found.',a1)
CC 	endif
CC	if(.NOT.exists) go to 10
 	filnam= filtmp

9999	continue
C
C	   THE
	   END


	DOUBLE PRECISION FUNCTION CSCANF(ARRAY,LOCa,LOCz,IERR)
c
c  CSCANF is a function that scans an input character array (a1) to
c  find the numeric value it contains.  The value is returned as 
c  a double precision real number.  The value in ARRAY can be integer,
c  decimal, or even exponential.
c
c$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$ rlr
CCC    Mon Aug 8, 1988 10:46:27,  ADD INCLUDE FOR CONSOLE I/O
ccc	INCLUDE 'DEVCONST.INC'
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

	character*1 array(*)
	integer loc1,loc2,ierr
	double precision val,tens
	integer exval
	INTEGER MSIGN,ESIGN
	LOGICAL FRACT,EXPO,aftern,aftere

	cscanf = -999.d0
	loc1 = loca
	loc2 = locz
	ierr= 1
	aftere= .false.
	aftern= .false.
	MSIGN= 1
	VAL= 0.
	FRACT= .FALSE.
	EXPO= .FALSE.
	ESIGN= 1
	EXVAL= 0
	TENS= 10.

c  if loc2<loc1, take error exit.
	if(loc2 .LT. loc1) goto 900
c  skip over leading blanks. (The "j9=loc2+1" takes care of
c  lines that are only blanks.)
	j9=loc2+1
	do 9 i9=loc1,loc2
	  if(array(i9) .EQ. ' ')goto 9
	  j9=i9
	  goto 19
9	continue
19	loc1=j9

c  if loc2<loc1, take error exit.
	if(loc2 .LT. loc1) goto 900
	DO 10 I=LOC1,LOC2
	  IF (EXPO) go to 8
	  IF ((ARRAY(I).EQ.'D'.OR.ARRAY(I).EQ.'E').and.aftern) GO TO 8
	  IF ((ARRAY(I).EQ.'d'.OR.ARRAY(I).EQ.'e').and.aftern) GO TO 8
	  IF (ARRAY(I).EQ.'+'.and..not.aftern) then
	 	MSIGN= +1
	 	go to 10
	  elseif (ARRAY(I).EQ.'-'.and..not.aftern) then
	  	MSIGN= -1
	  	go to 10
	  elseif (ARRAY(I).EQ.'.'.and..not.fract) then
	  	FRACT= .TRUE.
	  	go to 10
	  endif
	  if (array(i).eq.' '.and.aftern) go to 15
	  IF (ARRAY(I).LT.'0'.OR.ARRAY(I).GT.'9') GO TO 900
	  IF (FRACT) GO TO 5
	  VAL= VAL*10+ ICHAR(ARRAY(I))-48
	  aftern= .true.
	  GO TO 10

 5	  VAL= VAL+ (ICHAR(ARRAY(I))-48)/TENS
	  TENS= TENS*10.
	  GO TO 10

 8	  EXPO= .TRUE.
	  IF (ARRAY(I).EQ.'-'.and..not.aftere) then
	  	ESIGN= -1
	  	go to 10
	  elseif (array(i).eq.'+'.and..not.aftere)then
	  	esign= +1
	  	go to 10
	  endif
	  if (array(i).eq.' '.and.aftere) go to 15
	  if (ARRAY(I).LT.'0'.OR.ARRAY(I).GT.'9') GO TO 10
	  EXVAL= EXVAL*10+ ICHAR(ARRAY(I))-48
	  aftere= .true.
 10	CONTINUE

 15	CSCANF= MSIGN*VAL* 10.**(ESIGN*EXVAL)
 	go to 9999
 900	ierr= -1
9999  con  t  inue
C	   h
	   end


	SUBROUTINE INGRIN(QUEST,RESPON,MIN,MAX)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C  10 Nov 89 -  Completely replace ingrin with better
c		code adapted to HP.   tlr
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
c 31 May 90, 1653 - convert back to using improved "cscanf"
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

C  INGRIN PROMPTS THE CONSOLE OPERATOR
C  FOR A DOUBLE-PRECISION REAL NUMBER, (THAT IS TO SAY, AN INTEGER),
C  AND THEN CHECKS AND RETURNS THE OPERATOR'S RESPONSE.
C
C  ON INPUT --
C	QUEST  - (CHARACTER STRING) PROMPT QUESTION,
C	         Null string => no prompt question.
C	RESPON - (integer) DEFAULT RESPONSE
C	MIN    - (integer) LOWEST ACCEPTABLE VALUE
C	MAX    - (integer) HIGHEST ACCEPTABLE VALUE
C
C  ON OUTPUT --
C	QUEST,MIN,MAX - UNCHANGED
C	RESPON - UPDATE W/ OPERATOR RESPONSE IF DEFAULT NOT TAKEN
C
C					tlr 10 Nov 89
C$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C  Wed Aug 3, 1988 09:34:27
c  provide for null character.
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
CCC    Mon Aug 8, 1988 10:30:13,  ADD INCLUDE FOR CONSOLE I/O
	INCLUDE 'DEVCOM.INC'
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

	character*1 null
C	parameter (null=char(0))
	character QUEST*(*),line*80
	integer RESPON,MIN,MAX,TEMP,length
	integer itrimlen
	double precision CSCANF

	null= char(0)
C  WRITE QUESTION IF ONE WAS PASSED
 10	length= len(quest)
 1000	  format(1x,A,'<',i8,'>: ',$)
 	if (quest(1:1) .NE. null) write(IOUTDEV,1000) quest,respon

C  Obtain, check and return operator response
CC	read(INDEV,1010) line
	read(*,1010) line
 1010	format(a80)
CC 	if (line(1:2).eq.'qq') call ceror(0)
 	lenline = nullterm(line,80)
 	if(lenline .EQ. 0)goto 9999
	temp=Cscanf(line,1,lenline,ierr)
	if(ierr .NE. 1)goto 3
 	if (temp.GE.min .AND. temp.LE.max) goto 11
2 	  write(ierrout,1015)
 1015	  format(1x,'Input out of range or something')
 	  go to 10
3 	  write(ierrout,1016)
 1016	  format(1x,'Error on input')
 	  go to 10

11	continue
 	respon= temp

9999	continue
C
C	   THE
	   END




	character function scana(array,loc1,loc2)

c  scana finds the first non-blank character in 'array' and returns
c  it.  If nothing is found, blank is returned.
c
c  input:
c	array - character array; (CHARACTER SCALAR)
c	loc1  - lowest position in array to test;
c	loc2  - hightest position in array to test;
c  output
c	scana - first non-alpha character in array or blank;
c
c$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$ rlr/87

C  "ARRAY" IS NOT REALLY AN ARRAY; IT IS A CHARACTER VARIABLE (SCALAR)
C  WITH LENGTH DETERMINED BY THE CALLING ROUTINE.(Mon Aug 1, 1988 14:05:10)
	character array*(*)
	integer loc1,loc2

	scana = ' '
	
	do 100 i=loc1,loc2
	  if (array(i:i).eq.' ') go to 100
	    scana= array(i:i)
	    go to 120
100	continue

120	continue
C
C	   THE
	   END




	LOGICAL FUNCTION YESNO(QUEST,RESPON)

C  YESNO PROMPTS THE CONSOLE OPERATOR FOR
c  A DOUBLE-PRECISION REAL NUMBER, (that is, for a logical value),
C  AND THEN CHECKS AND RETURNS THE OPERATOR'S RESPONSE.
C
C  ON INPUT --
C	QUEST  - (CHARACTER STRING) PROMPT QUESTION,
C	         Null string => no prompt question.
C	RESPON - (LOGICAL) DEFAULT RESPONSE: .TRUE. or .FALSE..
C
C  ON OUTPUT --
C	QUEST  - UNCHANGED
C	YESNO  - (logical) answer (or default): 'y'->true; 'n'->false
C
C	CALLS -- SCANA (this call removed for the time being
C                       Wed Aug 3, 1988 15:49:26 )
C                                                         RWH 6/82
C$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C  Wed Aug 3, 1988 09:34:27 (sooner, actually)
c  provide for null character.
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
CCC    Mon Aug 8, 1988 10:33:09,  ADD INCLUDE FOR CONSOLE I/O
	INCLUDE 'DEVCOM.INC'
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

	character*1 null
C	parameter (null=char(0))

	character QUEST*(*)
	logical respon, okinput
	character line*80
	character temp*1,reschar*1

	null= char(0)
C  WRITE QUESTION IF ONE WAS PASSED
CC	write(*,*) "[",quest,"]",itrimlen(quest)
 10	continue
	length = len(quest)
c..... Initialize reschar ........................2
	  if(respon) then
	    reschar = 'Y'
	  else
	    reschar = 'N'
	  endif
c.................................................2
c.... Write question if there is one .................1
 	if (quest(1:1) .NE. null) then
 	  write(IOUTDEV,1000) quest,reschar
 1000	  format(1x,A,'(Y or N) <',a,'>: ',$)
	endif
c.....................................................1

C  Obtain, check and return operator response
CC	read(INDEV,1010) line
	read(*,1010) line
 1010	format(a80)
CC 	if (line(1:2).eq.'qq') call ceror(0)
	temp = line(1:1)
	okinput = .FALSE.
	if(temp .EQ. ' ') temp = reschar
	if(temp.EQ.'Y'  .OR.  temp.EQ.'y') then
	  yesno = .TRUE.
	  okinput = .TRUE.
	elseif(temp.EQ.'N'  .OR.  temp.EQ.'n') then
	  yesno = .FALSE.
	  okinput = .TRUE.
	else
	  okinput = .FALSE.
C	  call time(isec)
C	  jsec = isec/8
C	  jsec = jsec*8
C	  jsec = isec - jsec
C 	  if(jsec .NE. 0) write(ierrout,1015)
	  write(ierrout,1015)
 1015	  format(' *** Please enter one of the letters "Y" or "N".')
C 	  if(jsec .EQ. 0) write(ierrout,1016)
C 1016	  format(' *** Now look!  You''ve got to enter ',
C     1		'a "Y" or an "N"!')
 	endif

	if(.NOT.okinput) goto 10
C
C	   THE
	   END



	SUBROUTINE REALIN(QUEST,RESPON,MIN,MAX,FORM)

C  REALIN PROMPTS THE CONSOLE OPERATOR FOR A DOUBLE-PRECISION REAL NUMBER,
C  AND THEN CHECKS AND RETURNS THE OPERATOR'S RESPONSE.
C
C  ON INPUT --
C	QUEST  - (CHARACTER STRING) PROMPT QUESTION,
C	         Null string => no prompt question.
C	RESPON - (DP REAL) DEFAULT RESPONSE
C	MIN    - (DP REAL) LOWEST ACCEPTABLE VALUE
C	MAX    - (DP REAL) HIGHEST ACCEPTABLE VALUE
C	FORM   - (CHARACTER STRING) APPROPRIATE FORMAT FOR RESPON
C
C  ON OUTPUT --
C	QUEST,MIN,MAX - UNCHANGED
C	RESPON - UPDATE W/ OPERATOR RESPONSE IF DEFAULT NOT TAKEN
C
C	CALLS -- SCANF (removed, Dec 89) (put back, 31 May 90)
C                                                                   RWH 6/82
C$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C  Wed Aug 3, 1988 09:34:27
c  provide for null character.
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
c 31 May 90, 1653 - convert back to using improved "cscanf"
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
CCC    Mon Aug 8, 1988 10:35:47,  ADD INCLUDE FOR CONSOLE I/O
	INCLUDE 'DEVCOM.INC'
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

	character*1 null
C	parameter (null= char(0))

	character QUEST*(*),FORM*(*),line*80, rform*80
	double precision RESPON,MIN,MAX,TEMP
	integer length
	double precision CSCANF
	integer itrimlen

	null= char(0)
c  get the format of the expected response (and default)
	rform = '(1x,a,''<'''//form//'''>: '',$)'
C  WRITE QUESTION IF ONE WAS PASSED
 10	length= len(quest)
	if (quest(1:1) .NE. null) then
	  write(IOUTDEV,rform) quest,respon
	endif

C  Obtain, check and return operator response
CC	read(INDEV,1010) line
	read(*,1010) line
 1010	format(a80)
CC 	if (line(1:2).eq.'qq') call ceror(0)
	lenline=nullterm(line,80)
	if(lenline .EQ. 0)goto 9999
 	temp= Cscanf(line,1,lenline,ierr)
 	if (ierr.NE.1) then
 	  write(ierrout,1016)
1016	  format(1x,' Error on input.')
 	  go to 10
 	endif
 	if (temp.LT.min .OR. temp.GT.max) then
 	  write(ierrout,1015)
1015	  format(1x,'Value input is out of bounds.')
 	  go to 10
 	endif
 	respon= temp

9999	return
C
C	   THE
	   END


      integer function itrimlen(array)
C  find the length of an array to the first blank (or end of array)
C
      character array*(*)
      length= len(array)
      itrimlen= length
C      do 10 i=1,length
CC        if (array(i:i).ne.' ') go to 10
C        if (array(i:i).ne.char(0)) go to 10
C          itrimlen= i-1
C          go to 90
      do 10 i=length,1,-1
        if (array(i:i).eq.' ') go to 10
          itrimlen= i
          go to 90
10    continue
90    return
C
C	   THE
	   END

	function nullterm(A,L)
c
c  This function locates trailing blanks in the character
c  variable A and replaces the first trailing blank  with
c  a null. If there are no trailing blanks, then the last
c  character, that is A(L:L), is replaced with a null. If
c  A is all blanks, then A(1:1) is  the  first "trailing"
c  blank  and  is replaced with a null. The length of the
c  result,  not  counting the null, is  returned  as  the
c  value of the function. Imbedded blanks have no effect.
c  If  the value of L (the untrimmed length) is less than
c  or equal to zero, then that  value  is returned as the
c  value  of  the  function,  and  nothing is replaced by
c  a null character. The real purpose of this routine  is
c  to  terminate  a  file  name  with a null instead of a
c  string of blanks.

	character*(*) A

	i1 = L
	if (i1 .LE. 0) goto 3
1	continue
	if(i1 .EQ. 0)goto 2
	if(A(i1:i1) .NE. ' ') goto 2
	i1 = i1-1
	goto 1
2	continue
	if(i1 .EQ. L) i1 = i1-1
	j1 = i1+1
	A(j1:j1) = char(0)
3	continue
	nullterm = i1
c
c	   t
c	   h
	   e n d


	function  nlterm(A,L)
c  This function locates trailing blanks in the character
c  variable A and replaces the first trailing blank  with
c  a new line. If there are no trailing blanks, then the last
c  character, that is A(L:L), is replaced with a new line. If
c  A is all blanks, then A(1:1) is  the  first "trailing"
c  blank  and  is replaced with a new line. The length of the
c  result,  not  counting the new line, is  returned  as  the
c  value of the function. Imbedded blanks have no effect.
c  If  the value of L (the untrimmed length) is less than
c  or equal to zero, then that  value  is returned as the
c  value  of  the  function,  and  nothing is replaced by
c  a new line character. The real purpose of this routine  is
c  to  terminate  a  file  name  with a new line instead of a
c  string of blanks.

	character*(*) A

	i1 = L
	if (i1 .LE. 0) goto 3
1	continue
	if(i1 .EQ. 0)goto 2
	if(A(i1:i1) .NE. ' ') goto 2
	i1 = i1-1
	goto 1
2	continue
	if(i1 .EQ. L) i1 = i1-1
	j1 = i1+1
	A(j1:j1) = char(10)
3	continue
	 nlterm = i1
c
c	   t
c	   h
	   e n d

CC HP and probably other FORTRANs does not recognize '\n' as newline.
C This routine converts those newlines into new output lines explicitly.
C  05/22/00 - rlr
      subroutine questout (dev, quest)
      integer dev
      character *(*) quest
      integer i, itrimlen, j, k
      character *256 lineop
      character null

      null= char(0)
      lineop= ""
      j= 0
      k= 0
      do i= 1,itrimlen (quest)
	if (quest(i:i+1) .NE. "\n") then
C Drop 'n' of '\n' after writing record
	   if (k .EQ. 1) THEN
	      k= 0
           else
C Copy character to output line
	     lineop= lineop(1:j) // quest(i:i)
	     j= j+1
           endif	      
        else
C  Write line to output and re-initialize
CC	    write (dev, "(A)") lineop(1:itrimlen(lineop))
	    write (dev, "(A)") lineop(1:j)
	    lineop= ""
	    j= 0
	    k= 1
        endif
      enddo
C  Write last line...
      if (lineop(1:1) .NE. null) 
     x  write (dev, "(A, $)") lineop(1:j)

      return
      end
