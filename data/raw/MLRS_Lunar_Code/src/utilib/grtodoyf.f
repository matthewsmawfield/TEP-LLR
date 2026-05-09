C Convert Gregorian date to day of year (GRTODOY)

      subroutine grtodoy(year,month,day,doy)
      INCLUDE 'LICENSE-BSD3.inc'
      integer	year, month, day, doy

      integer*4	jd_in, jd_jan0

C  calculate # days since noon feburary 29, 1900 (julian date=2415078.0)
 
      if (month .le. 2) then
 
          jd_in = idint(1461.0d0 * (year-1)/4.0d0)
     &      + idint((153.0d0 * (month+9) + 2.0d0)/5.0d0) + day

      else

          jd_in = idint(1461.0d0 * year/4.0d0) 
     &      + idint((153.0d0 * (month-3) + 2.0d0)/5.0d0) + day

      endif

C what is jd of jan 0?
      jd_jan0 = idint(1461.0d0 * (year-1)/4.0d0)+ 306.0d0
    
C  now doy is easy.
      doy= jd_in- jd_jan0

      RETURN
      END
 
