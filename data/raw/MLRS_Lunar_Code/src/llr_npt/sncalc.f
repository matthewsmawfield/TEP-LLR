	subroutine sncalc
      INCLUDE 'LICENSE-BSD3.inc'
c
c  Calculate the signal to noise ratio for this normalpoint.

c  To do this, we accept the o-c residuals read from the data.  A
c  slope and offset of the 'signal' is then computed using lsq.  These 
c  results are used to flatten the slopes of the noise.  To simplify 
c  calculations, the noise from the noise_width nsec immediately below the data
c  is used for the noise part of the computations. Actually, all the noise
c  could be used, but only after compensating for the fact that the top
c  and bottom of the range gate, and hence the noise boundaries, will
c  now be sloped with respect to the 0-line from the signal fit.  
c
c  Another point is that the signal is averaged over +/- 3 sigma from the
c  mean.  It could be argued that the full width at half maximum or some 
c  other metric would give more meaningful results.
c
c  11/26/91
c  history
c  10/21/97 - Don't allow negative S;N ratio. rlr
c
c$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$ rlr
	implicit double precision (a-h,o-z)

	character sccsid*(*)
	parameter (sccsid = "@(#)sncalc.f	1.2\t01/06/92")
	double precision noise_width
	parameter (noise_width = 20.d0)

	INCLUDE 'NOISE.INC'

	double precision newomc
	integer rmode

C  fit the signal using the old o-c residuals
        RMODE= 2
        CALL LSQ(stime,somc,nsignal,RMODE,3.d0,A,B,C,SD,SE,NREJ)
	threesig= 3*sd

C  NORMALPOINT TIME AND RESIDUAL
	nsnoise= 0
	ceiling= -3.d0*sd
	floor  = -3.d0*sd- noise_width
	do 50 i=1,nnoise
	    newomc= nomc(i)- (A+ B*ntime(i)+ C*ntime(i)*ntime(i))
	    if (newomc .lt. ceiling .and. 
     .		newomc .gt. floor) nsnoise= nsnoise+ 1
 50	continue

C  Compute s:n
	if (nsnoise .gt. 0) then
	    snratio= ((nsignal/(6.d0*sd)) / (nsnoise/noise_width)) -1.d0
	    snratio= dmin1(snratio, 99.9d0)
	    if (snratio .lt. 0.d0) snratio= 0.0d0
	else
	    snratio= 99.9d0
	endif

	return
	end
