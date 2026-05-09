/**
 *		BullA
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
*****************************************************************************

 *  convert the USNO mark3.out EOP file from the usno to the MLRS 
 *  processing format.
 *
 * RLR - date unknown
**/
#include <stdio.h>
#include <string.h>

int
main ()
{
  double mjd, jd, x, y, dut;
  double sjd, fjd, djd;
  int try= 1;
  char str[256];

  while (fgets (str, 256, stdin) != NULL)
    {
      if (!isdigit(str[4])) continue;
      sscanf (str, "%*c%*d %*d %*d %lf %lf %*lf %lf %*lf %lf %*lf", 
              &mjd, &x, &y, &dut);
      if (try == 1)
        {
          try= 2;
          sjd= mjd+2400000.5;
        }
      else if (try == 2)
        {
          try= 0;
          djd= mjd+2400000.5- sjd;
        }
      fjd= mjd+2400000.5;
    }
  rewind (stdin);
  fprintf (stdout, "%10.1lf%10.1lf%4.1lf\n", sjd, fjd, djd);
  while (fgets (str, 256, stdin) != NULL)
    {
      if (strlen(str) < 70 || !isdigit(str[4])) continue;
      sscanf (str, "%*c%*d %*d %*d %lf %lf %*lf %lf %*lf %lf %*lf", 
              &mjd, &x, &y, &dut);
      fprintf (stdout, "%10.1lf%8.4lf%8.4lf%9.5lf\n", mjd+2400000.5, x, y, dut); 
    }

}
