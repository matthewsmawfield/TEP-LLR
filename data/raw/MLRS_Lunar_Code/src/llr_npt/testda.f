c	program testda ()
	implicit double precision (a-h,o-z)
c  this program tests writing to formatted direct access
c  file.

ccc	double precision a,b,c,d,e,f,g,h,ang1,ang2,x
	character*80 asciirec
	character*8  iu
5000	format(1x,i8,'  Header Record')
5001	format(a80)
5002	format(3(1x,d22.15))
5003	format(2(1x,d22.15))
5004	format(1x,d22.15,' {}')

	pie=3.1415962535897d0
	write(6,5004)pie
	radeg=3.14159265d0/180.d0
	open(17,iostat=ios,err=999,file='testot',
     1		status='unknown',access='direct',recl=80,
     2		form='formatted')

	inum = 0
	iu='asciirec'

	do 1 i1=1,5
	  do 2 i2=1,80
2	  asciirec(i2:i2)=' '
	inum = inum+1
	nrec = 4*(inum-1)+1
	x=i1
	a=0.2345678901234d0+i1
	b=a*1.d-10
	c=a*1.d-100
	ang1=radeg*(i1-1)
	ang2=ang1*10.d0
	d=sin(ang1)
	e=cos(ang1)
	f=sin(ang2)
	g=cos(ang2)
	h=i1
c------
	nrec1=nrec+1
	write(asciirec,5002)a,b,c
cc	asciirec(79:79)=char(13)
	asciirec(80:80)=char(10)
	irecerr=nrec1
	write(17,5001,err=998,iostat=ios,rec=nrec1)asciirec
c------
	nrec2=nrec+2
	write(asciirec,5003) d,e
cc	asciirec(79:79)=char(13)
	asciirec(80:80)=char(10)
	irecerr=nrec2
	write(17,5001,err=998,iostat=ios,rec=nrec2)asciirec
c------
	nrec3=nrec+3
	write(asciirec,5003) f,g
cc	asciirec(79:79)=char(13)
	asciirec(80:80)=char(10)
	irecerr=nrec3
	write(17,5001,err=998,iostat=ios,rec=nrec3)asciirec
c------
	nrec4=nrec+4
	write(asciirec,5004)  h
cc	asciirec(79:79)=char(13)
	asciirec(80:80)=char(10)
	irecerr=nrec4
	write(17,5001,err=998,iostat=ios,rec=nrec4)asciirec
c------
	type *,' records ',nrec,nrec1,nrec2,nrec3,nrec4
	ntrec=nrec+4
1	continue

	write(asciirec,5000) ntrec
cc	asciirec(79:79)=char(13)
	asciirec(80:80)=char(10)
	irecerr=1
	write(17,5001,err=998,iostat=ios,rec=1)asciirec
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

