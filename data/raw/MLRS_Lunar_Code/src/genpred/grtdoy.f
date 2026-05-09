	subroutine grtdoy(iyr,mon,ida,idoy)
	implicit double precision (a-h,o-z)
      INCLUDE 'LICENSE-BSD3.inc'
cccccccccccccccccccccccccccccccccccccccccccccccccccccc
c 13 Jun 90, 1544 - copy day-of-year program from 'c' program.
cccccccccccccccccccccccccccccccccccccccccccccccccccccc
c  subprogram returns day-of-year for given date.

	character*(*) sccsid
	parameter (sccsid = "@(#)grtdoy.f	1.1\t06/13/90")

	if(mon .LE. 2) then
	  jdin=1461.*(iyr-1)/4.
	  jdin=jdin + (153.*(mon+9)+2.)/5.+ida
	else
	  jdin=1461.*iyr/4.
	  jdin=jdin + (153.*(mon-3)+2.)/5.+ida
	endif

	jdjan0=(1461.*(iyr-1)/4.)+306.
	idoy=jdin - jdjan0

c	   the
	   end

