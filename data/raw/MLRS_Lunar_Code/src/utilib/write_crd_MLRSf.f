C write_crd_MLRSf.f
C  - writes MLRS-specific CRD records
C rlr
      INCLUDE 'LICENSE-BSD3.inc'

        SUBROUTINE write_93 (str)
	implicit none
        character*512 str
        INCLUDE '../include/crd_MLRS.inc'

        WRITE (str,1000,err=100) d93_sec_of_day, d93_sysconfig_id,
     &          full_range_type_ind, k_val, kp_val, geo_corr, 
     &          range_refraction, el_refraction, 
     &          range_OmC, range_OmC_post
        return

 100    write(*,*) "Error writing CRD record type 93"
        stop 1

 1000	FORMAT ("93 ",f18.12,1x,a4,1x,i1,1x,f7.3,1x,f7.3,1x,f7.3,
     &          1x,f9.5,1x,f9.5,1x,f11.3,1x,f11.3)
        end

C
C  Software versions
C

        SUBROUTINE write_95 (str)
	implicit none
        character*512 str
        INTEGER trimlen
        INTEGER monlen, satlen, declen, callen, psnlen, qlklen
        INTEGER nptlen, frslen
        INCLUDE '../include/crd_MLRS.inc'

        monlen= trimlen(monvers)
        satlen= trimlen(satvers)
        declen= trimlen(decvers)
        callen= trimlen(calvers)
        psnlen= trimlen(psnvers)
        qlklen= trimlen(qlkvers)
        nptlen= trimlen(nptvers)
        frslen= trimlen(frsvers)

        WRITE (str,1000,err=100) monvers(1:monlen), satvers(1:satlen),
     &          decvers(1:declen), calvers(1:callen),
     &          psnvers(1:psnlen), qlkvers(1:qlklen), 
     &          nptvers(1:nptlen), frsvers(1:frslen)
        return

 100    write(*,*) "Error writing CRD record type 95"
        stop 1

 1000	FORMAT ("95 ",a,1x,a,1x,a,1x,a,1x,a,1x,a,1x,a,1x,a)
        end
