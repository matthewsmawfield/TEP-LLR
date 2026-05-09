      subroutine error(icode,loc,ierr,ios)
      implicit double precision (a-h,o-z)
      INCLUDE 'LICENSE-BSD3.inc'
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c 08 Mar 91, 1643: change to call ceror which calls HP-UX exit.
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
	character*(*) sccsid
	parameter (sccsid = "@(#)error.f	1.2\t03/25/91")

	character*(*) loc

	include 'DEFCMN.INC'

	include 'DEVCOM.INC'

      write(ierrout,1000) loc,ierr,ios
1000	format(1X,A80,/,' Subroutine ERROR: #',I8,',  #',I8)
	call ceror(icode)
c  "stop" should never be executed.
	stop
c
c	   the
	   end

