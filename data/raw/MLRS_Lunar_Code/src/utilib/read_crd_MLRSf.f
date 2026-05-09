C read_crd_MLRSf.f
C  - reads MLRS-specific CRD records 
C  rlr
C
C initial state vector / tiv for target
C
      INCLUDE 'LICENSE-BSD3.inc'

	SUBROUTINE read_91 (str)
	character*512 str
	INCLUDE '../include/crd_MLRS.inc'

	READ (str(3:120),*,err=100) 
     &		pcayear, pcadoy, pcasec,
     &		pcapos(1), pcapos(2), pcapos(3),
     &		pcavel(1), pcavel(2), pcavel(3), progver

	return

 100	write(*,*) "Error reading CRD record type 91"
	stop 1
	end

	SUBROUTINE read_93 (str)
	character*512 str
	INCLUDE '../include/crd_MLRS.inc'

	READ (str(3:120),*,err=100) d93_sec_of_day, d93_sysconfig_id,
     &          full_range_type_ind, k_val, kp_val, geo_corr, 
     &          range_refraction, el_refraction, 
     &		range_OmC, range_OmC_post

	return

 100	write(*,*) "Error reading CRD record type 93: [",str,"]"
	stop 1
	end

	SUBROUTINE read_94 (str)
	character*512 str
	INCLUDE '../include/crd_MLRS.inc'

	READ (str(3:120),*,err=100) d93_sec_of_day

	return

 100	write(*,*) "Error reading CRD record type 94"
	stop 1
	end

C
C  software versions
C
	SUBROUTINE read_95 (str)
	character*512 str
	INCLUDE '../include/crd_MLRS.inc'

	READ (str(3:120),*,err=100) 
     &          monvers, satvers, decvers, calvers,
     &          psnvers, qlkvers, nptvers, frsvers

	return

 100	write(*,*) "Error reading CRD record type 95"
	stop 1
	end

