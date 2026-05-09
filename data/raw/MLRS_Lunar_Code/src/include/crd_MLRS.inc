C 91 - MLRS - pass initial state
	COMMON /d91/ pcayear, pcadoy, pcasec, pcapos, pcavel, progver
	integer pcayear, pcadoy, pcasec
	real*8 pcapos(3), pcavel(3)
	character progver*1

C 92 - MLRS - Point angle corrections
	COMMON /d92/ d92_sec_of_day, azbias, elbias
	double precision d92_sec_of_day, azbias, elbias

C 93 - MLRS - Filtering information
	COMMON /d93/ d93_sec_of_day, d93_sysconfig_id,
     &          full_range_type_ind, k_val, kp_val,
     &		geo_corr, range_refraction, el_refraction, 
     &		range_OmC, range_OmC_post
	real*8 d93_sec_of_day
	character*40 d93_sysconfig_id
        integer full_range_type_ind
	real*8 k_val, kp_val, geo_corr, range_refraction, el_refraction
	real*8 range_OmC, range_OmC_post

C 94 - MLRS - LLR FR separator for normalpoints
	COMMON /d94/ d94_sec_of_day
	double precision d94_sec_of_day

C 95 - MLRS software versions.
	COMMON /d95/ monvers,satvers,decvers,calvers,psnvers,qlkvers,
     &          nptvers,frsvers
        character*40 monvers,satvers,decvers,calvers,psnvers,qlkvers, 
     &          nptvers,frsvers
