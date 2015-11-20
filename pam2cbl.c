/* pam2cbl.c                                http://cbl.sourceforge.net
 * 
 * code to print out a scoring matrix for CBL smith-waterman routines.
 * paste the output of this code into the sslookup definition of your 
 * driver for the cb_sw routine
 *
 * usage: pam2cbl <pam filename>
 *
 * Copyright (C) 2004 University of Alaska Fairbanks
 * Arctic Region Supercomputing Center (ARSC)
 * http://www.arsc.edu
 * 
 *   This library is free software; you can redistribute it and/or modify
 *   it under the terms of the GNU Lesser General Public License as
 *   published by the Free Software Foundation; either version 2.1 of the
 *   License, or (at your option) any later version.
 *  
 *   This library is distributed in the hope that it will be useful, but
 *   WITHOUT ANY WARRANTY; without even the implied warranty of
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 *   Lesser General Public License for more details.
 *  
 *   You should have received a copy of the GNU Lesser General Public
 *   License along with this library; if not, write to the
 *   Free Software Foundation, Inc.
 *   59 Temple Place, Suite 330
 *   Boston, MA  02111-1307 USA
 *
 * $Id: pam2cbl.c,v 1.1 2004/04/27 20:04:42 jlong777 Exp $
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

int main(int argc, char * argv[])
{
  int i, j, index[32], first=1, last_index, max=0, match, mismatch, output[27][32];
  char buffer[128], *p;
  FILE * pam;

  /* initialize output array */
  for(i=0; i<27; i++)
    for(j=0; j<32; j++)
      output[i][j] = -999;

  /* open pam or blosum file */
  if((pam = fopen(argv[1], "r"))==0x0)
  {
    printf("Can not open %s, exiting...\n", argv[1]);
    exit(1);
  }

  while(fgets(buffer, 128, pam))
  {
    if(buffer[0] == '#') continue; /* go down to the data */

    if(first) /* the 1st line is the letter indicies */
    {
      i = 0;
      p = strtok(buffer, " ");
      index[i++] = p[0] - 64;

      while(p = strtok(NULL, " "))
      {
        index[i++] = p[0] - 64;
        /*printf("%d ",index[i-1]);*/
      }
      last_index = i-1;
      i = first = 0;
    }
    else /* the rest of the lines are values, except for the 1st entry */
    {
      j = 0;
      p = strtok(buffer, " ");

      while(p = strtok(NULL, " "))
      {
        if(strlen(p) > max) max = strlen(p);
        if(index[j]>0 && index[j]<27 && i<last_index)
	{
          output[index[i]][index[j++]] = atoi(p);
          /*printf("%3d ",output[index[i]][index[j-1]]);*/
        }
	else
	{
	  if(j==last_index)
	    match = atoi(p);
	  else
	    mismatch = atoi(p);
	  j++;
	}
      }
      i++;
    }
    /*printf("\n");*/
  }
  /*printf("match = %d  mismatch = %d\n", match, mismatch);*/
  /* finalize output array and print */

  for(i=0; i<27; i++)
  {
    for(j=0; j<32; j++)
    {
      if(output[i][j] == -999)
      {
	if(i>0 && j>0 && j<27)
	{
          if(i==j)
	    output[i][j] = match;
	  else
	    output[i][j] = mismatch;
        }
	if(i==0 || j==0 || j>26)
	  output[i][j] = -1;
      }
      if(i==26 && j==31)
      {
        if(max < 4)
          printf("%2d\n",output[i][j]);
	else
	  printf("%3d\n",output[i][j]);
      }
      else
      {
        if(max < 4)
          printf("%2d,",output[i][j]);
	else
	  printf("%3d,",output[i][j]);
      }
    }
    printf("\n");
  }
}

