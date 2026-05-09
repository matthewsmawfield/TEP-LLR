C  INCLUDE FILE; Mon Nov 14, 1988 18:03:06
	COMMON /NPHLD/ JDI,TIM(1500),RANGE(1500),RESID(1500),
     X			NPTS,NREJ,NPSKEL
		DOUBLE PRECISION JDI
		CHARACTER NPSKEL*130
c  max is dimension of above arrays. used for detecting too many normalpoints.
		parameter (max=1500)

