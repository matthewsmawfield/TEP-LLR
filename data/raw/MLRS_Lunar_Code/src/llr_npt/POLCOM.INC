c 03 Jul 91, 1504: filpol is from user (euler.f) or command line (iprdargs
c  in misclibnp.f). Used by kustner in jeulpkg.f.
	character*255 filpol
	logical ispol
	common /polcom/ filpol, lenpol, ispol

