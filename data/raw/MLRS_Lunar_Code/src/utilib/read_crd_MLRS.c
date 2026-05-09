#include <stdio.h>
#include <string.h>
#include "../include/crd_MLRS.h"

char stro[256];

/*
 * K and K' cals and O-C residual
 */
void
read_93 (char *str, struct rd93 *data_recd)
{
  char recd_id[3];

  sscanf (str,
	  "%s %Lf %s %d %lf %lf %lf %lf %lf %lf %lf",
	  recd_id,
	  &data_recd->sec_of_day, &data_recd->sysconfig_id, 
          &data_recd->full_range_type_ind,
	  &data_recd->k_val, &data_recd->kp_val,
	  &data_recd->geo_corr, &data_recd->range_refraction,
	  &data_recd->el_refraction, &data_recd->range_OmC,
	  &data_recd->range_OmC_post);
}

/*
 * Software versions
 */
void
read_95 (char *str, struct rd95 *data_recd)
{
  char recd_id[3];

  sscanf (str,
          "%s %s %s %s %s %s %s %s %s\n",
	  recd_id,
          &data_recd->monvers, &data_recd->satvers,
          &data_recd->decvers, &data_recd->calvers,
          &data_recd->psnvers, &data_recd->qlkvers,
          &data_recd->nptvers, &data_recd->frsvers);
}
