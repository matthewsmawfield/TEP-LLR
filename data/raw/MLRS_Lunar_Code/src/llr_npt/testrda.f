c	program testrda ()
c  this program tests reading from formatted direct access
c  file.

	double precision a,b,c
5000	format(1x,i8)
5001	format(a80)
5002	format(3(1x,d22.15))
5003	format(2(1x,d22.15))
5004	format(1x,d22.15)

	radeg=3.14159265d0/180.d0
	open(17,iostat=ios,err=999,file='testot',
     1		status='unknown',access='direct',recl=80,
     2		form='formatted')

	read(17,5002,err=998,iostat=ios,rec=2)a,b,c
	type *,' a,b,c=',a,b,c
	read(17,5002,err=998,iostat=ios,rec=6)a,b,c
	type *,' a,b,c=',a,b,c
	read(17,5002,err=998,iostat=ios,rec=18)a,b,c
	type *,' a,b,c=',a,b,c
c------
	read(17,5003,err=998,iostat=ios,rec=16)a,b
	type *,' a,b=',a,b
	read(17,5003,err=998,iostat=ios,rec=8)a,b
	type *,' a,b=',a,b
	read(17,5003,err=998,iostat=ios,rec=11)a,b
	type *,' a,b=',a,b
c------
	read(17,5004,err=998,iostat=ios,rec=13)a
	type *,' a=',a
c------
	read(17,5000,err=998,iostat=ios,rec=1)ix
	type *,' i=',ix
c------

	close(17)
	goto 10000

999	continue
	type *,' Error in opening file "testot".'
	pause
	goto 10000

998	continue
	type *, ' Error on read record ',irecerr
	type *, ' I/O status is ',ios
	pause
	goto 10000

10000	continue
c	   the
	   end

