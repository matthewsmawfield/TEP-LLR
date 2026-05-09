C  INCLUDE FILE; Mon Nov 14, 1988 10:15:30
cc Thu Apr 19 15:48:49 CDT 1990 - tlr
cc Add 7th file type; types 6 & 7 have offset of -1.
cc (Actually, lunio allows only one record to be written
cc  to or read from file types 6 & 7.)
	COMMON/LIOCMN/FNAME(7),RCDLEN(7),OFFSET(7)
	CHARACTER*14 FNAME
	INTEGER RCDLEN,OFFSET
