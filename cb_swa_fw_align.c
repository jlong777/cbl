/* cb_swa_fw_align.c                             http://cbl.sourceforge.net
 *
 * called by cb_swa_fw.
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
 * $Id: cb_swa_fw_align.c,v 1.2 2004/03/10 00:22:26 jlong777 Exp $
 */
 
#include "cb_macro.h"
#include <stdio.h>
#include <stdlib.h>
 
void fill_alg(long *, long, long, long *, long *, long, long **, long **, long **, long);
	      
void cb_swa_fw_align(long *swtab, long *dbl, long dbllen, long *dbs,long dbslen,
                     long smax, long **algl, long **algm, long **algs,long *alglen,
                     long *algstl, long *algsts, long *errno)
{
  long blah, i, index=0, j, savi=0, savj=0;
  
  /* find smax */
  for(i=dbllen;i>0;i--)
  {
    for(j=dbslen;j>0;j--)
    {
#ifdef LONG32
      if((0x0FFFFFFF & swtab[(dbslen+1)*i+j]) == smax)
#endif
#ifdef LONG64
      if((0x0FFFFFFFFFFFFFFF & swtab[(dbslen+1)*i+j]) == smax)
#endif
      {
        savi = i;
        savj = j;
	goto jump;
      }
    }
  }
  
jump:

  /* compute alglen */
  
  i = savi; j = savj; *alglen = 0;
  
#ifdef LONG32
  while(blah=((unsigned long)(0xF0000000 & swtab[(dbslen+1)*i + j]) >> 28))
#endif
#ifdef LONG64
  while(blah=((unsigned long)(0xF000000000000000 & swtab[(dbslen+1)*i + j]) >> 60))
#endif
    switch(blah)
    {
      case 0:  printf("I should never get here!\n");
	       break;
      case 1:  i--;
               (*alglen)++;
	       break;
      case 2:  j--;
               (*alglen)++;
               break;
      case 3:  i--;
               (*alglen)++;
               break;
      default: i--;
	       j--;
               (*alglen)++;
               break;
    }
  
  *algstl = i+1;
  *algsts = j+1;
  
  *algl = (long *) malloc(*alglen+8);
  *algm = (long *) malloc(*alglen+8);
  *algs = (long *) malloc(*alglen+8);
  
  if(*algl==0 || *algm==0 || *algs==0)
  {
    printf("cb_swa_fw_align: malloc error in align, exiting...\n");
    exit(0);
  }
  
  index = *alglen-1;

  /* recursively trace back through scoring matrix for alignment */
  fill_alg(swtab, savi, savj, dbl, dbs, index, algl, algm, algs, dbslen);
}


void fill_alg(long *swtab, long swtab_i, long swtab_j, long *dbl, long *dbs,
              long index, long **algl, long **algm, long **algs, long dbslen)
{
  char *algl_c = (char *)(*algl);
  char *algm_c = (char *)(*algm);
  char *algs_c = (char *)(*algs);
  
  if(swtab_i==0 || swtab_j==0 || index<0) return;
  
#ifdef LONG32
  switch((unsigned long)(0xF0000000&swtab[(dbslen+1)*swtab_i + swtab_j]) >> 28)
#endif
#ifdef LONG64
  switch((unsigned long)(0xF000000000000000&swtab[(dbslen+1)*swtab_i + swtab_j]) >> 60)
#endif
  {
    case 0:  /* end  */
             swtab_i = 0;
	     swtab_j = 0;
	     break;
    case 1:  /* horizontal (fortran), vertical (c) */
             algl_c[index] = ((char *)dbl)[swtab_i-1];
             algm_c[index] = ' ';
             algs_c[index] = '-';
	     swtab_i--;
             break;
    case 2:   /* vertical (fortran), horizontal (c) */
             algl_c[index] = '-';
             algm_c[index] = ' ';
             algs_c[index] = ((char *)dbs)[swtab_j-1];
	     swtab_j--;
	     break;
    case 3:  /* h or v */
             algl_c[index] = ((char *)dbl)[swtab_i-1];
             algm_c[index] = ' ';
             algs_c[index] = '-';
	     swtab_i--;
             break;
    default: /* corner */
             algl_c[index] = ((char *)dbl)[swtab_i-1];
             algs_c[index] = ((char *)dbs)[swtab_j-1];
	     if(algl_c[index] == algs_c[index])
               algm_c[index] = ':';
	     else
	       algm_c[index] = ' ';
	     swtab_i--;
	     swtab_j--;
             break;
  }
  
  index--;

  fill_alg(swtab, swtab_i, swtab_j, dbl, dbs, index, algl, algm, algs, dbslen);
}
