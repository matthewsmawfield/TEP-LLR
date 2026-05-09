/*
 * MLRS-specific CRD write routines
 */
#include <stdio.h>
#include <string.h>
#include "../include/crd_MLRS.h"

char stro[256];

/* * initial state vector / tiv for target
 */
void
write_91 (FILE * str_out, struct rd91 data_recd)
{
  fprintf (str_out,
	   "91 %2d %3d %5d %11.2f %11.2f %11.2f %10.4f %10.4f %10.4f %1c\n",
	   data_recd.pcayear, data_recd.pcadoy, data_recd.pcasec,
	   data_recd.pcapos[0]/1.e1, data_recd.pcapos[1]/1.e1, 
	   data_recd.pcapos[2]/1.e1, data_recd.pcavel[0]/1.e3,
	   data_recd.pcavel[1]/1.e3, data_recd.pcavel[2]/1.e3, 
	   data_recd.progver);
}

/*
 * pointing biases
 */
void
write_92 (FILE * str_out, struct rd92 data_recd)
{
  fprintf (str_out,
	   "92 %.3Lf %.4f %.4f\n",
	   data_recd.sec_of_day, data_recd.azbias, data_recd.elbias);
}

/*
 * K and K' cals and O-C residual
 */
void
write_93 (FILE * str_out, struct rd93 data_recd)
{
  fprintf (str_out,
	   "93 %.12Lf %s %d %.3f %.3f %.3f %.5f %.5f %.3f %.3f\n",
	   data_recd.sec_of_day, data_recd.sysconfig_id, 
           data_recd.full_range_type_ind, 
	   data_recd.k_val, data_recd.kp_val, 
	   data_recd.geo_corr, data_recd.range_refraction, 
	   data_recd.el_refraction, 
	   data_recd.range_OmC, data_recd.range_OmC_post);
}

/*
 * LLR normalpoint break
 */
void
write_94 (FILE * str_out, struct rd94 data_recd)
{
  fprintf (str_out,
	   "94 %.3Lf\n",
	   data_recd.sec_of_day);
}

/*
 * Software versions
 */
void
write_95 (FILE * str_out, struct rd95 data_recd)
{
  fprintf (str_out,
	   "95 %s %s %s %s %s %s %s %s\n", 
           data_recd.monvers, data_recd.satvers,
           data_recd.decvers, data_recd.calvers,
           data_recd.psnvers, data_recd.qlkvers,
           data_recd.nptvers, data_recd.frsvers);
}
