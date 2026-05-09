      CHARACTER eph_source*4,eph_notes*10,tar_name*10
      DOUBLE PRECISION prf
      DOUBLE PRECISION xpond_xmit_delay,xpond_utc_off,xpond_osc_drift
      DOUBLE PRECISION cofm_corr
      INTEGER cospar_id,sic,norad,eph_seq
      INTEGER sttpyear,sttpmonth,sttpday,sttphour,sttpmin,sttpsec
      INTEGER endpyear,endpmonth,endpday,endphour,endpmin,endpsec
      INTEGER epyear,epmon,epday,ephour,ephsep,compat,version
      INTEGER tar_type,ref_frame,rot_type,cofm_app,ltcorr
      INTEGER atrk_ro_0,xtrk_ro_0,rtrk_ro_0
      INTEGER atrk_ro_1,xtrk_ro_1,rtrk_ro_1
      INTEGER atrk_ro_2,xtrk_ro_2,rtrk_ro_2
      COMMON /ephihdr/ prf,
     x   xpond_xmit_delay,xpond_utc_off,xpond_osc_drift,cofm_corr,
     x   cospar_id,sic,norad,eph_seq,
     x   sttpyear,sttpmonth,sttpday,sttphour,sttpmin,sttpsec,
     x   endpyear,endpmonth,endpday,endphour,endpmin,endpsec,
     x   epyear,epmon,epday,ephour,ephsep,compat,version,
     x   tar_type,ref_frame,rot_type,cofm_app,ltcorr,
     x   atrk_ro_0,xtrk_ro_0,rtrk_ro_0,atrk_ro_1,xtrk_ro_1,rtrk_ro_1,
     x   atrk_ro_2,xtrk_ro_2,rtrk_ro_2,
     x   eph_source,eph_notes,tar_name
