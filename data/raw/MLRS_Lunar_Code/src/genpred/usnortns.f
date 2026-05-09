C     These routines are from an early version of the United States
C     Naval Observatory (USNO) NAVAL OBSERVATORY VECTOR ASTROMETRY SOFTWARE
C    (NOVAS) software package.
C
      SUBROUTINE ABERAT (POS1,VE,TLIGHT,POS2)
C     THIS SUBROUTINE CORRECTS POSITION VECTOR FOR ABERRATION OF LIGHT.
C     ALGORITHM INCLUDES RELATIVISTIC TERMS.  SEE MURRAY (1981)
C     MON. NOTICES ROYAL AST. SOCIETY 195, 639-648.
C
C          POS1   = POSITION VECTOR, REFERRED TO ORIGIN AT CENTER OF
C                   MASS OF THE EARTH, COMPONENTS IN AU (IN)
C          VE     = VELOCITY VECTOR OF CENTER OF MASS OF THE EARTH,
C                   REFERRED TO ORIGIN AT SOLAR SYSTEM BARYCENTER,
C                   COMPONENTS IN AU/DAY (IN)
C          TLIGHT = LIGHT TIME FROM BODY TO EARTH IN DAYS (IN)
C                   IF TLIGHT = 0.0D0, THIS SUBROUTINE WILL COMPUTE
C          POS2   = POSITION VECTOR, REFERRED TO ORIGIN AT CENTER OF
C                   MASS OF THE EARTH, CORRECTED FOR ABERRATION,
C                   COMPONENTS IN AU (OUT)
C
C NOTE: The term in POS2 involving Q/R is close to that from classic
C physics. In classis physics, GAMMAI/R would = 1. The GAMMAI term leads
C to a POS2 vector very different from the classic formulation. However,
C the resulting HA/Dec are very close to those in the classic case
C (~1 mas). This has been learned twice, now, and I don't want go
C through this again! :-) rlr
C
      DOUBLE PRECISION POS1,VE,TLIGHT,POS2,C,TL,P1MAG,VEMAG,
     .     BETA,DOT,COSD,GAMMAI,P,Q,R,DSQRT
      DIMENSION POS1(3), VE(3), POS2(3)
C
      DATA C / 173.14463348D0 /
C     C = SPEED OF LIGHT IN AU/DAY
C
      TL = TLIGHT
      P1MAG = TL * C
      IF (TL.NE.0.0D0) GO TO 20
      P1MAG = DSQRT(POS1(1)**2 + POS1(2)**2 + POS1(3)**2)
      TL = P1MAG / C
   20 VEMAG = DSQRT(VE(1)**2 + VE(2)**2 + VE(3)**2)
      BETA = VEMAG / C
      DOT = POS1(1)*VE(1) + POS1(2)*VE(2) + POS1(3)*VE(3)
      COSD = DOT / (P1MAG * VEMAG)
      GAMMAI = DSQRT(1.0D0 - BETA**2)
      P = BETA * COSD
      Q = (1.0D0 + P / (1.0D0 + GAMMAI)) * TL
      R = 1.0D0 + P
C
      DO 30 J=1,3
        POS2(J) = (GAMMAI * POS1(J) + Q * VE(J)) / R
CC        write(*,*) "aberat: ",J,POS1(J),VE(J),POS2(J),POS2(J)-POS1(J),
CC     x           (GAMMAI * POS1(J))/R- POS1(J),(Q * VE(J)/R)
   30 CONTINUE
CC      write (*,*) "tl, p1mag ",tl,p1mag
CC      write(*,*)"vemag,beta,dot,cosd,p ",vemag,beta,dot,cosd,p
CC      write(*,*) "aberat gqr:",gammai,q,r
      RETURN
C
      END
 
      SUBROUTINE WOBBLERECT (X,Y,POS1,POS2)
C
C     THIS SUBROUTINE CORRECTS EARTH-FIXED GEOCENTRIC RECTANGULAR
C     COORDINATES FOR POLAR MOTION.  IT TRANSFORMS A VECTOR FROM
C     EARTH-FIXED GEOGRAPHIC SYSTEM TO ROTATING SYSTEM BASED ON
C     ROTATIONAL EQUATOR AND ORTHOGONAL GREENWICH MERIDIAN THROUGH
C     AXIS OF ROTATION.
C
C          X      = CONVENTIONALLY-DEFINED X COORDINATE OF CELESTIAL
C                   EPHEMERIS POLE WITH RESPECT TO IERS REFERENCE
C                   POLE, IN ARCSECONDS (IN)
C          Y      = CONVENTIONALLY-DEFINED Y COORDINATE OF CELESTIAL
C                   EPHEMERIS POLE WITH RESPECT TO IERS REFERENCE
C                   POLE, IN ARCSECONDS (IN)
C          POS1   = VECTOR IN GEOCENTRIC RECTANGULAR
C                   EARTH-FIXED SYSTEM, REFERRED TO GEOGRAPHIC
C                   EQUATOR AND GREENWICH MERIDIAN (IN)
C          POS2   = VECTOR IN GEOCENTRIC RECTANGULAR
C                   ROTATING SYSTEM, REFERRED TO ROTATIONAL EQUATOR
C                   AND ORTHOGONAL GREENWICH MERIDIAN (OUT)
C
C
      DOUBLE PRECISION X,Y,POS1,POS2,SECCON,XPOLE,YPOLE,
     .     XX,YX,ZX,XY,YY,ZY,XZ,YZ,ZZ
      DIMENSION POS1(3), POS2(3)
C
      DATA SECCON / 206264.8062470964D0 /
C
      XPOLE = X / SECCON
      YPOLE = Y / SECCON
C
C     WOBBLE ROTATION MATRIX FOLLOWS
      XX =  1.0D0
      YX =  0.0D0
      ZX = -XPOLE
      XY =  0.0D0
      YY =  1.0D0
      ZY =  YPOLE
      XZ =  XPOLE
      YZ = -YPOLE
      ZZ =  1.0D0
   10 CONTINUE
C
C     PERFORM ROTATION
   20 POS2(1) = XX*POS1(1) + YX*POS1(2) + ZX*POS1(3)
      POS2(2) = XY*POS1(1) + YY*POS1(2) + ZY*POS1(3)
      POS2(3) = XZ*POS1(1) + YZ*POS1(2) + ZZ*POS1(3)
C
   50 RETURN
C
      END
