#include <stdio.h>
#include <string.h>
#include "../include/crd.h"
#include "../include/crd_MLRS.h"
struct rh4 h4;
struct rd10 d10;
struct rd30 d30;
struct rd95 d95;
char *frsvers= "1.03";
void write_vers00 ();

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

/* frd_strip
 *	Remove mlrs-specific records from full rate data, compressing some
 *	other records.
 * rlr - 7/31/07
 *
 * History:
 * 07/01/09 - Set filter flag and time of flight to 0 for combined 1- and 2-way
 *            ranging. v 1.02 rlr.
 * 04/13/10 - Convert '95' (versions) record into 00 record for output fr file.
 *            v 1.03 rlr.
 */
main (argc, argv)
     int argc;
     char *argv[];
{
  char str[256];
  FILE *str_in, *str_out;

  if ((str_in = fopen (argv[1], "r")) == NULL)
    {
      printf ("Could not open file %s\n", argv[1]);
      exit (1);
    }
  if ((str_out = fopen (argv[2], "w")) == NULL)
    {
      printf ("Could not open file %s\n", argv[2]);
      exit (1);
    }

  /* Copy and reformat data */
  while (fgets (str, 256, str_in) != NULL)
    {
      if (strncmp (str, "h4", 2) == 0 || strncmp (str, "H4", 2) == 0)
        {
          read_h4 (str, &h4);	/* for range_type_ind */
        }
      if (strncmp (str, "9", 1) == 0)
	{
          if (strncmp (str, "95", 2) == 0)
            {
              read_95 (str, &d95);	/* for range_type_ind */
              strcpy (d95.frsvers, frsvers);
              write_vers00 (str_out, d95);
            }
          else
            {
	      /* drop it */
            }
	}
      else if (strncmp (str, "10", 2) == 0)
	{
          /* squeeze out some spaces */
          read_10 (str, &d10);
          if (h4.range_type_ind == 4 && d10.filter_flag == 3) 
            {
              d10.filter_flag = 0;
              d10.time_of_flight = 0.e0;
            }
          if (h4.range_type_ind == 4 || 
              (d10.filter_flag >= 0 && d10.filter_flag <= 2))
            write_10 (str_out, d10);
        }
      else if (strncmp (str, "30", 2) == 0)
	{
	  /* drop the extra angle records */
	  read_30 (str, &d30);
	  if (d30.angle_origin_ind >= 0)
	    fputs (str, str_out);
	}
      else
	{
	  fputs (str, str_out);
	}

    }
}

/*
 * Software versions
 */
void
write_vers00 (FILE * str_out, struct rd95 data_recd)
{
  fprintf (str_out,
           "00 Software versions: %s %s %s %s %s %s\n",
           data_recd.monvers, data_recd.satvers,
           data_recd.decvers, data_recd.calvers,
           data_recd.psnvers, data_recd.frsvers);
}

