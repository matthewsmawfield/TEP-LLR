C Note that /ephi/ contains enough space for 8 days' worth of 20 second preds
C if MAXENT= 35000
      PARAMETER (MAXENT=35000)
      DOUBLE PRECISION jdi(MAXENT),jdf(MAXENT),
     x  outv(MAXENT,6),backv(MAXENT,6), c12(MAXENT),c23(MAXENT),
     x  aberoutv(MAXENT,3),aberbackv(MAXENT,3),
     x  oscrelcorr(MAXENT), jds(MAXENT),
     x  offsetjdi(MAXENT),offsetjdf(MAXENT),
     x  offsetoutv(MAXENT,3),offsetbackv(MAXENT,3), offsetjds(MAXENT),
     x  rotjdi(MAXENT),rotjdf(MAXENT),rotangv(MAXENT,3),gast(MAXENT),
     x  rotjds(MAXENT), eopjdi(MAXENT),eopjdf(MAXENT),
     x  eopxy(MAXENT,2),eopdutc(MAXENT), eopjds(MAXENT)
      INTEGER nv, nvoff, nvrot, nveop, leapflag(MAXENT)
      CHARACTER *10 offsetname(MAXENT)
      COMMON /ephi/ jdi,jdf,outv,backv,c12,c23,aberoutv,aberbackv,
     x   oscrelcorr,jds,offsetjdi,offsetjdf,
     x   offsetoutv,offsetbackv, offsetjds, 
     x   rotjdi,rotjdf,rotangv,gast,rotjds,
     x   eopjdi,eopjdf,eopxy, eopdutc, eopjds, nv, nvoff, nvrot, nveop,
     x   leapflag,offsetname 
