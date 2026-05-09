cc   Wed Dec  6 14:48:29 CST 1989
c  This include file contains the common block with flags
c  and file names that may be set by command line arguements.
c  (E.g. "npt -r -i mnpifile -l npt-log-13a")
	logical nflg,rflg,iflg,oflg,mflg,lflg,sflg,cflg,tflg,gflg
	character*255 clfilin,clfilout,clfilmini,clfillog,clfilsum,
     1		clfilres
	double precision nptlen,sigmult
        integer ntype
	common /comarg/ 
     1		nflg,rflg,iflg,oflg,mflg,lflg,sflg,cflg,tflg,gflg,
     2		clfilin,clfilout,clfilmini,clfillog,clfilsum,clfilres,
     3		nargs,nptlen,sigmult,ntype

