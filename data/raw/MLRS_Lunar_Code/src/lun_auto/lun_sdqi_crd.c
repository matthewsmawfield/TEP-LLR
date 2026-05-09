#include <stdio.h>
#include <string.h>
#include <math.h>
#include "../include/math_const.h"
#include "../include/crd.h"
#define FNLEN   256
#define DATA_DIR_NAME "."
FILE *str_in, *str_ls, *str_out, *str_sum;

void setup_files (int, char **, char *, char *);
char compute_sdqi (int, int, double, double);
struct rd10 d10;
struct rd50 d50;

/*****************************************************************************
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
*****************************************************************************/

/*----------------------------------------------------------------------------
**
**   Program Lun_sdqi_crd      - computes 'subjective' data quality indicator
**
**   NOTE: If any 'printf' statements are left in use, their output will be 
**        interpreted as an indication of a program error by the ldb script. 
**        So, only leave error printf's in operation.
**
** History:
**   07/06/98 - first, experimental, version. rlr.
**   07/25/07 - CRD format version. rlr.
**--------------------------------------------------------------------------*/

main (argc, argv)
     int argc;
     char *argv[];

{
  static char SccsId[] = "%W%\t%G%";

  char sdqi[20], sdqif = 'e';

  int i = 0, j, tot_obs, exp_obs, good_obs;

  int nrecd = 0, status, got_50 = 0;

  double slope, omc, sd, snr;

  char fn_in[FNLEN], fn_sum[FNLEN], str[256];

/**==== S E T U P ====**/
/*  get file names and open files  */
  setup_files (argc, argv, fn_in, fn_sum);

/**==== R E A D   S U M M A R Y   F I L E ====**/
/*  Read and reformat lunar data file, one shot at a time  */
  do
    {

      /* Read every record in burst */
      while (status = fgets (str, 256, str_sum) != NULL)
	{

	  /*  Get Poisson results
	   */
	  if (strncmp (str, "Poisson: Time bin", 17) == 0)
	    {
	      if (strncmp (&str[25], "Total Points:", 13) == 0)
		{
		  getifield (str, 38, 5, &tot_obs);
		  good_obs = 0;
		}
	      else
		{
		  getifield (str, 39, 5, &tot_obs);
		  getifield (str, 55, 4, &exp_obs);
		  getifield (str, 70, 4, &good_obs);
		}
	      /**printf (" tot_obs %d exp_obs %d good_obs %d\n", 
                 tot_obs, exp_obs, good_obs);**/
	    }
	  if (strncmp (str, "Poisson: Slope", 14) == 0)
	    {
	      getfield (str, 15, 7, &slope);
	      getfield (str, 40, 8, &omc);
	      getfield (str, 59, 6, &sd);
	      getfield (str, 75, 5, &snr);
	      sdqi[i++] = compute_sdqi (good_obs, exp_obs, snr, sd);
	      /**printf ("slope %f omc %f sd %f snr %f sdqi %c\n", 
                 slope, omc, sd, snr, sdqi[i - 1]);**/
	    }
	  if (strncmp (str, "Poisson: Expected points:", 25) == 0)
	    {
	      getifield (str, 25, 4, &exp_obs);
	      sdqi[i++] = compute_sdqi (good_obs, exp_obs, snr, sd);
	      /**printf ("Expected: %d sdqi %c\n", exp_obs, sdqi[i - 1]);**/
	    }
	}
    }
  while (status != NULL);	/*  end do  */

/* For now, choose best sdqi and use for whole group */
/* This may be ultimate solution, since lun_auto can create normalpoints
 * without respect to Poisson boundaries
 */
  for (j = 0; j < i; j++)
    {
      if ((int) sdqif > (int) sdqi[j])
	sdqif = sdqi[j];
      /**printf ("i %d sdqi %c sdqif %c\n", i, sdqi[j], sdqif);**/
    }

/**==== W R I T E   D A T A ====**/
/*  open output summary file for appending (pass 2)  */
  fclose (str_sum);
  if ((str_sum = fopen (fn_sum, "a")) == NULL)
    {
      printf ("Could not open file %s\n", fn_sum);
      exit (1);
    }

/*  report header  */
  fprintf (str_ls, "                       Lun_sdqi   v 2.1  \n\n");
  fprintf (str_ls, "\nAutomated Subjective Data Quality Indicator for\n%s:\n\n", fn_in);
  fprintf (str_sum, "                       Lun_sdqi   v 2.1 \n\n");
  fprintf (str_sum, "\nAutomated Subjective Data Quality Indicator for\n%s:\n\n", fn_in);

/* Record changes to report files */
  fprintf (str_ls, "New quality indicator: %c\n", sdqif);
  fprintf (str_sum, "New quality indicator: %c\n", sdqif);

/*  Copy lunar (.lu) data file, changing the sdqi. */
  do
    {

      /*   Read every record in file */
      while (status = fgets (str, 256, str_in) != NULL)
	{
	  if (strncmp(str, "10", 2) == 0)
	    {
	      read_10 (str, &d10);
            }
	  if (strncmp(str, "50", 2) == 0)
	    {
	      read_50 (str, &d50);
              if (sdqif == ' ') d50.data_qual_ind= 0 ;
              else d50.data_qual_ind= sdqif- 'a' ;
	      write_50 (str_out, d50);
	      got_50= 1;
	    }
	  else
	    {
	      fputs (str, str_out);
	    }
	  nrecd++;
	}			/* end while */
    }
  while (status != NULL);	/*  end do  */
  if (!got_50)
    {
      strncpy(d50.sysconfig_id, d10.sysconfig_id, 5);
      if (sdqif == ' ') d50.data_qual_ind= 0 ;
      else d50.data_qual_ind= sdqif- 'a'+ 1;
      d50.sess_rms= -1;
      d50.sess_skew= -1;
      d50.sess_kurtosis= -1;
      d50.sess_PmM= -1;
      write_50 (str_out, d50);
    }

  fclose (str_in);
  fclose (str_out);
  exit (0);
}

/*-------------------------------------------------------------------------
**
**   setup_files      - Open input and output file names
**
**-----------------------------------------------------------------------*/

void
setup_files (argc, argv, fn_in, fn_sum)
     int argc;
     char *argv[];
     char fn_in[FNLEN];
     char fn_sum[FNLEN];
{
  char fn_base[FNLEN], fn_ls[FNLEN], fn_out[FNLEN], fn_root[FNLEN];

  int found = 0, i;

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
      printf ("lun_adqi [[fnbase] | [-i fnin -o fnout -l fnls -s fnsum]]");
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

/*  open output summary file for reading (pass 1)  */
  if ((str_sum = fopen (fn_sum, "r")) == NULL)
    {
      printf ("Could not open file %s\n", fn_sum);
      exit (1);
    }
}
char 
compute_sdqi (int good_obs, int exp_obs, double snr, double sd)
{
  char sdqi;

/*  Compute sdqi */
  sdqi = 'e';
  if (good_obs - exp_obs >= 1)	/* Has to be something expected! */
    {
      sdqi = 'd';
      if (snr >= 0.4 && snr < 3.0 && sd <= 0.4)
	sdqi = 'c';
      if (snr >= 3.0 && snr < 10.0 && sd <= 0.3)
	sdqi = 'b';
      if (snr >= 10.0 && sd > 0.5)
	sdqi = 'b';
      if (snr >= 10.0 && sd <= 0.5)
	sdqi = 'a';
    }
  return (sdqi);
}
