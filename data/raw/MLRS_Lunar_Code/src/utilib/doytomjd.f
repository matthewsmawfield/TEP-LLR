C DOYTOMJD.F
C
C FOR A GIVEN YEAR AND DAY OF YEAR
C ROUTINE CALCULATES THE EQUIVALENT MODIFIED JULIAN-DATE.
C
C INPUT:
C   YEAR   - (INTEGER) GREGORIAN YEAR - 2 DIGITS, since 1900
C   DOY    - (INTEGER) DAY OF YEAR
C
C OUTPUT: 
C   MJD  - DOUBLE PRECISION MODIFIED JULIAN DATE
C
C rlr
      subroutine doytomjd(year,doy,mjd)

      integer*4 year, doy
      double precision mjd

C     calculate # days since noon Jan 0, 1900 (julian date=2415018.0) 
C     Use of the previous year is to get the number of leap days correct. 
      mjd = (1461.0d0 * (year- 1)/4.0d0) + doy + 366.0d0

C     Add mjd for 1/0/1900 */
      mjd = int(mjd) + 15018.0d0

      RETURN
      END
