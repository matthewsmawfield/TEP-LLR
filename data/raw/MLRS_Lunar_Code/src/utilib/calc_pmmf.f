CC	double precision x(25), w(25), bsize, pmm
CC	integer n
CC	data x /1,2,3,4,5,6,3,5,4,4,3,2,4,5,3,3,3,4,5,2,1,6,3,4,3/
CC	data w /1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,0,1,1,1,1,1,1,1,1,1/
CC	data bsize /1/
CC	data n /25/  
CC
CC	pmm1= pmm (x,w,n,bsize)
CC	pmm2= pmm (x,w,n,bsize/2)
CC	write (*,*) "bs1 pmm1 bs2 pmm2",bsize,pmm1,bsize/2,pmm2
CC	end
      INCLUDE 'LICENSE-BSD3.inc'

C CALC_PMM computes the peak-mean of a data array.
C  Pmm is a CRD data field for the npt, cal, and session records.
C  Created around 2007. rlr.
	double precision function calc_pmm (x, w, n, bsize, mean)
	double precision x(*), w(*), bsize, mean
        PARAMETER (NH = 5000)
	integer n, h(NH)
	double precision max, min, tmax, tmin
	double precision blo, bhi, maxv
	integer hmax, i, j, nb

	max= -1.e30
	min= 1.e30;
	do i=1,NH
	  h(i)= 0;
        enddo

CC	write(*,*)"x & w"
	tmax= mean+ (NH/2)*bsize
	tmin= mean- (NH/2)*bsize
CC	write(*,*) "tmin, tmax: ", tmin, tmax, mean
	do i=1,n
CC	WRITE(*,*) i,x(i),w(i)
          if (dabs(w(i)) .ge. 1.d-9) then
            if (x(i) .gt. max) max= x(i)
            if (x(i) .lt. min) min= x(i)
          endif
CC	WRITE(*,*) i,x(i),w(i),max,min
        enddo
CC	write(*,*) "min max ",min,max
C	Make sure the range fits into predefined 'h' array.
	if (max .gt. tmax) max= tmax
	if (min .lt. tmin) min= tmin
CC	write(*,*) "min max ",tmin,tmax

	nb= (max- min)/bsize
CC	write(*,*)"max min size nb",max,min,bsize,nb
	if (DABS(DBLE(nb)) .gt. NH) then
CC          write(*,*)"Sorry, too many bins: ",nb
	  calc_pmm= -1.d0
          return
	endif
	blo= min
	bhi= min+ bsize
	hmax= -1
        maxv= 0
	do j=1, nb
	  do i=1, n
            if (dabs(w(i)) .gt. 1.d-9) then
              if (x(i) .ge. blo .and. x(i) .lt. bhi) h(j)= h(j)+ 1
	    endif
	  enddo
CC	write(*,*) "bin: j blo bhi h(j)",j,blo,bhi,h(j)
          if (h(j) .gt. hmax) then
	    hmax= h(j)
	    maxv= (blo+bhi)/2
	  endif
CC	write(*,*) "j hmax, maxv",j,hmax,maxv
          blo= blo+ bsize
          bhi= bhi+ bsize
	enddo
	calc_pmm= maxv
	return
	end
