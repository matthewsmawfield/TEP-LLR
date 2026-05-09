#include <stdio.h>
#include <string.h>
#include <math.h>
#include <stdlib.h>
#include "../include/math_const.h"
#include "../include/crd.h"
#include "../include/crd_MLRS.h"
#define FNLEN   256
#define DATA_DIR_NAME "."
FILE *str_in,
 *str_ls,
 *str_out,
 *str_sum;
int ndata = 0;
int nnoise = 0;
int tnoise = 0;
int ngrp = 0;
double sod[1000];
double resid[1000];
double begint[50];
double endt[50];
int npts[50];
double tresid;
double snr = -10;
double minsnr = 2.0;    /* Signal:noise below which only 1 group is allowed */
int maxgap = 3 * 60;
int mindata = 3;          /* Minimum number of points/normalpoint */
int mindataspan1= 8 * 60; /* Minimum span of data in a normalpoint (>1/run) */
int mindataspan2= 30 * 60;/* Only 1 npt if <= mindata2 over mindataspan2 sec */

void combine ();
void one_grp ();
void read_h4 ();
void read_10 ();
void read_93 ();
void setup_files ();
void write_94 ();

struct rh4 h4;
struct rd10 d10;
struct rd93 d93;
struct rd94 d94;

/*----------------------------------------------------------------------------
**
**   Program Lun_auto_crd      - decide how to break up lunar
**                		 run into normalpoints.
**
******************************************************************************
* 
* Copyright (c) 2017, The University of Texas at Austin
* All rights reserved.
* 
* Redistribution and use in source and binary forms, with or without 
* modification, are permitted provided that the following conditions are met:
* 
* 1. Redistributions of source code must retain the above copyright notice, 
* this list of conditions and the following disclaimer.
*
* 2. Redistributions in binary form must reproduce the above copyright notice, 
* this list of conditions and the following disclaimer in the documentation 
* and/or other materials provided with the distribution.
*
* 3. Neither the name of the copyright holder nor the names of its contributors 
* may be used to endorse or promote products derived from this software without 
* specific prior written permission.
*
* THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" 
* AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE 
* IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE 
* ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE 
* LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR 
* CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF 
* SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS 
* INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN 
* CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) 
* ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF 
* THE POSSIBILITY OF SUCH DAMAGE.
*
******************************************************************************
**
**   NOTE: If any 'printf' statements are left in use, their output will be 
**        interpreted as an indication of a program error by the ldb script. 
**        So, only leave error printf's in operation.
**
** History:
**   03/22/96 - First, experimental version. rlr.
**   08/30/96 - Rounded number of tgrps up (added 0.5) as old version failed
**              to properly break up a group 2.8x the min data span. rlr.
**   07/08/98 - Remove all but the first #1 and #2 headers. This insures that
**              llr_npt will not create normalpoints shorter than dictated
**              by lun_auto. rlr.
**   07/25/07 - CRD format versions. rlr.
**--------------------------------------------------------------------------*/
int
main (argc, argv)
     int argc;
     char *argv[];

{
  int i;

  int good_10= 0,
    nrecd= 0,
    combined_grp = 0,
    status,
    testgrp;

  double dataspan,
    otestsod,
    span,
    t0,
    testsod,
    tgrp;

  char fn_in[FNLEN],
    str[256];

/**==== S E T U P ====**/
/*  get file names and open files  */
  setup_files (argc, argv, fn_in);

/*  report header  */
  fprintf (str_ls, "                       Lun_auto   v 2.1\n\n");
  fprintf (str_ls, "\nAutomated Segmentation Statistics for\n%s:\n\n", fn_in);
  fprintf (str_sum, "                       Lun_auto   v 2.1\n\n");
  fprintf (str_sum, "\nAutomated Segmentation Statistics for\n%s:\n\n", fn_in);

/*  where in the file are we?  */

/**==== R E A D   D A T A ====**/
/*  Read lunar data file, saving time and residual of returns */
  do
    {

      /*   Read every record in file */
      while ((status = fgets (str, 256, str_in)) != NULL)
        {
	  if (strncmp (str, "h4", 2) == 0)
	    {
	      read_h4 (str, &h4);
              fprintf (str_ls, "Date: %02d/%02d/%04d %02d:%02d\n",
		       h4.start_day, h4.start_mon, h4.start_year,
		       h4.start_hour, h4.start_min);
              fprintf (str_sum, "Date: %02d/%02d/%04d %02d:%02d\n",
		       h4.start_day, h4.start_mon, h4.start_year,
		       h4.start_hour, h4.start_min);
	    }
	  else if (strncmp (str, "10", 2) == 0)
	    {
	      read_10(str, &d10);
	      if (d10.filter_flag == 2)
		{
                  if (ndata+ 1 >= 1000)
                    {
                      printf ("Too many returns!!\n");
                      exit (-1);
                    }
                  sod[ndata] = d10.sec_of_day;
                  if (ndata > 0 && sod[ndata] < sod[ndata - 1])
                        sod[ndata] += 86400.;	    /* new day... */
		  good_10= 1;
		}
	      else
		{
		  good_10= -1;
		}
	    }
	  else if (strncmp (str, "93", 2) == 0)
	    {
	      read_93(str, &d93);
	      if (good_10 > 0)
		{
	          resid[ndata]= d93.range_OmC_post;
		  good_10= 0;
                  ndata++;
		}
	      else if (good_10 < 0)
		{
                  if (fabs (d93.range_OmC_post) < 20.)
                    nnoise++;
                  tnoise++;
		}
	  good_10= 0;
	    }
          nrecd++;
        }                       /* end while */
    }
  while (status != NULL);       /*  end do  */
  if (nnoise > 0) snr = (ndata / 2.) / (nnoise / (40. - 2.)) - 1.;
  else snr= 99.;
  /**printf ("nrecd: %6d ndata %6d nnoise %6d tnoise %6d s:r %5.1lf\n",
          nrecd, ndata, nnoise, tnoise, snr);**/
  fprintf (str_ls, "nrecd: %6d ndata %6d nnoise %6d tnoise %6d s:r %5.1lf\n",
           nrecd, ndata, nnoise, tnoise, snr);
  fprintf (str_sum, "nrecd: %6d ndata %6d nnoise %6d tnoise %6d s:r %5.1lf\n",
           nrecd, ndata, nnoise, tnoise, snr);

/**==== F I N D   G A P S   I N   D A T A ====**/
  ngrp = 0;
  begint[ngrp] = sod[0];
  npts[ngrp] = 1;
  for (i = 1; i < ndata; i++)
    {
      if ((sod[i] - sod[i - 1]) > maxgap)
        {
          endt[ngrp] = sod[i - 1];
          ngrp++;
          if (ngrp >= 20)
            {
              printf ("Too many groups!!\n");
              fprintf (str_ls, "Too many groups!!\n");
              fprintf (str_sum, "Too many groups!!\n");
              exit (-1);
            }
          begint[ngrp] = sod[i];
          npts[ngrp] = 1;
        }
      else
        {
          npts[ngrp]++;
        }
    }
  endt[ngrp] = sod[ndata - 1];
  ngrp++;
  /**printf ("Preliminary configuration:\n");**/
  fprintf (str_ls, "Preliminary configuration:\n");
  fprintf (str_sum, "Preliminary configuration:\n");
  for (i = 0; i < ngrp; i++)
    {
      /**printf ("group: %2d begin: %6lf end: %6lf dt: %4lf npts: %6d\n",
              i, begint[i], endt[i], endt[i] - begint[i], npts[i]);**/
      fprintf (str_ls,
               "group: %2d begin: %6lf end: %6lf dt: %4lf npts: %6d\n",
               i, begint[i], endt[i], endt[i] - begint[i], npts[i]);
      fprintf (str_sum,
               "group: %2d begin: %6lf end: %6lf dt: %4lf npts: %6d\n",
               i, begint[i], endt[i], endt[i] - begint[i], npts[i]);
    }
/**==== S E G M E N T   D A T A ====**/
  if (ndata >= mindata)
    {                           /* already required by Poisson??? */
      do
        {
          combined_grp = 0;
          for (i = 0; i < ngrp; i++)
            {
              if (npts[i] < mindata || (endt[i] - begint[i]) < mindataspan1)
                {
                  combine (i);
                  combined_grp = 1;
                  break;
                }
            }
        }
      while (combined_grp);
      /**printf ("Intermediate A configuration:\n");**/
      fprintf (str_ls, "Intermediate A configuration:\n");
      fprintf (str_sum, "Intermediate A configuration:\n");
      for (i = 0; i < ngrp; i++)
        {
          /**printf ("group: %2d begin: %6lf end: %6lf dt: %4lf npts: %6d\n",
                  i, begint[i], endt[i], endt[i] - begint[i], npts[i]);**/
          fprintf (str_ls,
                   "group: %2d begin: %6lf end: %6lf dt: %4lf npts: %6d\n",
                   i, begint[i], endt[i], endt[i] - begint[i], npts[i]);
          fprintf (str_sum,
                   "group: %2d begin: %6lf end: %6lf dt: %4lf npts: %6d\n",
                   i, begint[i], endt[i], endt[i] - begint[i], npts[i]);
        }
      dataspan = sod[ndata - 1] - sod[0];
      /**printf("dataspan %f\n",dataspan);**/
      if (dataspan < mindataspan1)
        {
	  one_grp();
        }
      else if (dataspan < mindataspan2 && ndata <= 6)
        {
	  one_grp();
        }
/**
      if (snr < minsnr)
        {
	  one_grp();
        }
**/
      if (ngrp == 1)
        {
          tgrp = (int) (dataspan / mindataspan1 + 0.5);
          span = dataspan / (tgrp <= 1 ? 1 : tgrp - 1);
          ngrp = 0;
          npts[ngrp] = 0;
          t0 = sod[0];
          for (i = 0; i < ndata; i++)
            {
              if ((sod[i] - t0) > span)
                {
                  endt[ngrp] = sod[i - 1];
                  ngrp++;
                  if (ngrp >= 20)
                    {
                      printf ("Too many groups!!\n");
                      fprintf (str_ls, "Too many groups!!\n");
                      fprintf (str_sum, "Too many groups!!\n");
                      exit (-1);
                    }
                  begint[ngrp] = sod[i];
                  npts[ngrp] = 1;
                  t0 += span;
                }
              else
                {
                  npts[ngrp]++;
                }
            }
          endt[ngrp] = sod[ndata - 1];
          ngrp++;
          /**printf ("Intermediate B configuration:\n");**/
          fprintf (str_ls, "Intermediate B configuration:\n");
          fprintf (str_sum, "Intermediate B configuration:\n");
          for (i = 0; i < ngrp; i++)
            {
              /**printf ("group: %2d begin: %6lf end: %6lf dt: %4lf npts: %6d\n",
                      i, begint[i], endt[i], endt[i] - begint[i], npts[i]);**/
              fprintf (str_ls,
                  "group: %2d begin: %6lf end: %6lf dt: %4lf npts: %6d\n",
                       i, begint[i], endt[i], endt[i] - begint[i], npts[i]);
              fprintf (str_sum,
                  "group: %2d begin: %6lf end: %6lf dt: %4lf npts: %6d\n",
                       i, begint[i], endt[i], endt[i] - begint[i], npts[i]);
            }
        }

      do
        {
          combined_grp = 0;
          for (i = 0; i < ngrp; i++)
            {
              if (npts[i] < mindata || (endt[i] - begint[i]) < mindataspan1)
                {
                  combine (i);
                  combined_grp = 1;
                  break;
                }
            }
        }
      while (combined_grp);
      if (ngrp == 0) ngrp= 1;
      /**printf ("Final configuration:\n");**/
      fprintf (str_ls, "Final configuration:\n");
      fprintf (str_sum, "Final configuration:\n");
      for (i = 0; i < ngrp; i++)
        {
          /**printf ("group: %2d begin: %6lf end: %6lf dt: %4lf npts: %6d\n",
                  i, begint[i], endt[i], endt[i] - begint[i], npts[i]);**/
          fprintf (str_ls,
                   "group: %2d begin: %6lf end: %6lf dt: %4lf npts: %6d\n",
                   i, begint[i], endt[i], endt[i] - begint[i], npts[i]);
          fprintf (str_sum,
                   "group: %2d begin: %6lf end: %6lf dt: %4lf npts: %6d\n",
                   i, begint[i], endt[i], endt[i] - begint[i], npts[i]);
        }
/**==== O U T P U T   D A T A ====**/
      /** write the data back out in separate groups. **/
      /*  now set the file pointer back to the beginning of the file */
      rewind (str_in);
      testgrp = 1;
      /*   Read every record in file */
      while ((status = fgets (str, 256, str_in)) != NULL)
        {
	  if (strncmp (str, "10", 2) == 0)
	    {
	      read_10(str, &d10);
	      if (d10.filter_flag == 2)
                    {
                      testsod = d10.sec_of_day;
                      if (ndata > 0 && testsod < otestsod)
                        testsod += 86400.;	                                        /* new day... */
                      otestsod = testsod;
                      if (testsod >= begint[testgrp] && testgrp < ngrp)
                        {
                          testgrp++;
			  d94.sec_of_day= d10.sec_of_day;
			  write_94 (str_out, d94);
                        }
                    }
              fputs (str, str_out);
	    }
          else
            {
              fputs (str, str_out);
            }
        }                       /* end while */

    }

  fclose (str_in);
  fclose (str_out);
  exit(0);
}

/*-------------------------------------------------------------------------
**
**   one_grp      - combine data group into one, set ngrp = 0;
**
**-----------------------------------------------------------------------*/
void
one_grp ()
{
  int i;

  for (i = 1; i < ngrp; i++)
    {
       npts[0]+= npts[i];
    }
  endt[0] = endt[ngrp-1];

  /* End processing: we don't want to look any further */
  ngrp= 0;
}
/*-------------------------------------------------------------------------
**
**   combine      - combine data group with another on the stack
**
**-----------------------------------------------------------------------*/
void
combine (grp_no)
     int grp_no;
{
  int i,
    indx;

  /* Remove first group in pass */
  if (grp_no == 0)
    {
      endt[0] = endt[1];
      npts[0] += npts[1];
      for (i = 2; i < ngrp; i++)
        {
          begint[i - 1] = begint[i];
          endt[i - 1] = endt[i];
          npts[i - 1] = npts[i];
        }

      /* Remove last group */
    }
  else if (grp_no == ngrp - 1)
    {
      endt[grp_no - 1] = endt[grp_no];
      npts[grp_no - 1] += npts[grp_no];

      /* Remove any other group */
      /** Combine it with smaller of groups on either side */
    }
  else
    {
      if (endt[grp_no - 1] - begint[grp_no - 1] < mindataspan1)
        {
          /* add to previous group */
          endt[grp_no - 1] = endt[grp_no];
          npts[grp_no - 1] += npts[grp_no];
          indx = grp_no;
        }
      else if (endt[grp_no + 1] - begint[grp_no + 1] < mindataspan1)
        {
          /* add to next group */
          endt[grp_no] = endt[grp_no + 1];
          npts[grp_no] += npts[grp_no + 1];
          indx = grp_no + 1;
        }
      else if (npts[grp_no - 1] > npts[grp_no + 1])
        {
          /* add to next group */
          endt[grp_no] = endt[grp_no + 1];
          npts[grp_no] += npts[grp_no + 1];
          indx = grp_no + 1;
        }
      else
        {
          /* add to previous group */
          endt[grp_no - 1] = endt[grp_no];
          npts[grp_no - 1] += npts[grp_no];
          indx = grp_no;
        }
      /* Now, copy groups down by one */
      for (i = indx + 1; i < ngrp; i++)
        {
          begint[i - 1] = begint[i];
          endt[i - 1] = endt[i];
          npts[i - 1] = npts[i];
        }
    }
  ngrp--;
}

/*-------------------------------------------------------------------------
**
**   setup_files      - Open input and output file names
**
**-----------------------------------------------------------------------*/

void
setup_files (argc, argv, fn_in)
     int argc;
     char *argv[];
     char fn_in[FNLEN];
{
  char fn_base[FNLEN],
    fn_ls[FNLEN],
    fn_sum[FNLEN],
    fn_out[FNLEN],
    fn_root[FNLEN];

  int found = 0,
    i;

/*  Get file name  */
  if (argc <= 2)
    {
      if (argc < 2)
        {
          printf ("Enter name of input file: ");
          scanf ("%s", fn_base);
        }
      else
        {
          strcpy (fn_base, argv[1]);
        }
      strcpy (fn_root, DATA_DIR_NAME);
      strcat (fn_root, "/");
      strcat (fn_root, fn_base);
      strcpy (fn_in, fn_root);
      strcat (fn_in, ".lu");
      strcpy (fn_out, fn_root);
      strcat (fn_out, ".lx");
      strcpy (fn_ls, fn_root);
      strcat (fn_ls, ".ll");
      strcpy (fn_sum, fn_root);
      strcat (fn_sum, ".ls");
      found = 15;
    }
  else if (argc == 9)
    {
      for (i = 1; i < argc; i += 2)
        {
          if (argv[i][0] == '-')
            {
              if (argv[i][1] == 'i')
                {
                  strcpy (fn_in, argv[i + 1]);
                  found += 1;
                }
              else if (argv[i][1] == 'o')
                {
                  strcpy (fn_out, argv[i + 1]);
                  found += 2;
                }
              else if (argv[i][1] == 'l')
                {
                  strcpy (fn_ls, argv[i + 1]);
                  found += 4;
                }
              else if (argv[i][1] == 's')
                {
                  strcpy (fn_sum, argv[i + 1]);
                  found += 8;
                }
            }
        }
    }
  if ((argc > 2 && argc != 9) || found != 15)
    {
      printf ("lun_auto [[fnbase] | [-i fnin -o fnout -l fnls -s fnsum]]");
      exit (1);
    }

/*  open input lun file  */
  if ((str_in = fopen (fn_in, "r")) == NULL)
    {
      printf ("Could not open file %s\n", fn_in);
      exit (1);
    }

/*  open output lun file  */
  if ((str_out = fopen (fn_out, "w")) == NULL)
    {
      printf ("Could not open file %s\n", fn_out);
      exit (1);
    }

/*  open output listing file  */
  if ((str_ls = fopen (fn_ls, "a")) == NULL)
    {
      printf ("Could not open file %s\n", fn_ls);
      exit (1);
    }

/*  open output summary file  */
  if ((str_sum = fopen (fn_sum, "a")) == NULL)
    {
      printf ("Could not open file %s\n", fn_sum);
      exit (1);
    }
}
