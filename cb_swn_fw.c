/* cb_swn_fw.c                                 http://cbl.sourceforge.net
 *
 * See original man page at the bottom.
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
 * $Id$
 */

#include "cbl.h"
#include "cb_macro.h"
#include <stdio.h>
#include <stdlib.h>

void gotoh_score(long *, long *, long *, long, long, long, long, 
                 long *, long *, long *, long *, long *, long *);

void gotoh_align(long *, long *, long, long *, long, long, long **, long **, 
                 long **, long *, long *, long *, long *, long, long);

void cb_swa_fw(long *dbl, long dbllen, long *dbs, long dbslen, long eg, 
               long og, long *sslookup, long *smax, long **algl, long **algm, 
               long **algs, long *alglen, long *algstl, long *algsts, 
	       long *errno)
{
  /* note - this version has different internal routines than the one 
   *        described in the original Cray version, see the man page
   *        for the portable cb_swa_fw
   */
  
  long *swtab, *gaph, *gapv, eog, sav_i, sav_j, *tb;
  
  tb    = (long *) malloc((dbllen+1)*(dbslen/4+8));
#ifdef LONG32
  gaph  = (long *) calloc(dbslen+1, 4);
  gapv  = (long *) calloc(dbslen+1, 4);
#endif
#ifdef LONG64
  gapv  = (long *) calloc(dbslen+1, 8);
  gaph  = (long *) calloc(dbslen+1, 8);
#endif
  eog = eg + og;

  gotoh_score(tb, gaph, gapv, dbllen, dbslen, eg, eog, smax, sslookup, dbl, dbs,
              &sav_i, &sav_j);
  gotoh_align(tb, dbl, dbllen, dbs, dbslen, *smax, algl, algm, algs, 
              alglen, algstl, algsts, errno, sav_i, sav_j);

  free(gaph);
  free(gapv);
  free(tb);
}

void gotoh_score(long *tb, long *mgaph, long *mgapv, long dbllen, long dbslen, 
		 long eg, long eog, long *smax, long *sslookup, long *dbl, 
		 long *dbs, long *sav_i, long *sav_j)
{
  long i, j;
  register long max, max1, max2, max3, rmax=0, ri=0, rj=0, r;
  register unsigned long t, trace;
  register long di, di1, e, g;

  unsigned long *temp = (unsigned long *) calloc(dbslen+1, sizeof(long));

  e = eg;
  g = eog;
  t = 0x0;
#ifdef LONG32
  r = dbslen/16 +1;
#endif
#ifdef LONG64 
  r = dbslen/32 +1;
#endif

  /* init */
  for(i=0; i<dbslen+1; i++)
  { 
    mgapv[i] = mgaph[i] = -g - (i-1)*e;
    temp[i] = 0x0;
  }
#ifdef LONG32
  for(i=0; i<dbslen/16+1; i++)
#endif
#ifdef LONG64
  for(i=0; i<dbslen/32+1; i++)
#endif
    tb[i] = 0x0;

  for(i=1; i<dbllen+1; i++)
  {
    di = (dbslen+1)*i;
    di1 = (dbslen+1)*(i-1);

    max2 = 0 - g;

    for(j=1; j<dbslen+1; j+=32) /* remember to handle tail cases */
    {
      max =   0x0;
      trace = 0x0;
      max1 = temp[j-1] + 
             sslookup[(4*((unsigned long)dbs[j/32] >> 30))+
	                 ((unsigned long)dbl[i/32] >> 30)];
      temp[j-1] = max2 + g;
      /*max2 = max2 - gaph[j-1];*/
      max3 = temp[j] - g;
      mgaph[j] = (max2 > mgaph[j-1]-e)? max2: mgaph[j-1]-e;
      mgapv[j] = (max3 > mgapv[j]-e)? max3: mgapv[j]-e;

      if(max1>max) max = max1; 
      if(mgaph[j]>max) max = mgaph[j]; 
      if(mgapv[j]>max) max = mgapv[j]; 
      if(max>=rmax)
      {
        rmax = max;
        ri = i;
        rj = j;
      }

      if(max)
      {
        if(max==max2) trace = 0x2; /* vert (fortran), horiz (c) */
        if(max==max3) trace = 0x1; /* horiz (fortran), vert (c) */
        if(max==max1) trace = 0x3; /* came from corner */
#ifdef DEBUG_INFO
#ifdef LONG32
        if(max > 0x0FFFFFFF) 
	  printf("cb_sw_fw_score: cell score exceeds maximum...\n");
#endif
#ifdef LONG64
        if(max > 0x0FFFFFFFFFFFFFFF) 
          printf("cb_sw_fw_score: cell score exceeds maximum...\n");
#endif
#endif
      }
#ifdef LONG32 
      t |= (trace << (30-2*(j%16)));
      if((j%16) == 15)
      {
        tb[r*i + j/16] = t; /* tb is the 2-bit traceback matrix */
	t = 0x0;
      }
#endif
#ifdef LONG64
      t |= (trace << (62-2*(j%32))); /* 0x1F&j is j%32 */
      if((j%32) == 31)
      {
        tb[r*i + j/32] = t; /* (unsigned long)j>>5 is j/32 */
	t = 0x0;
      }
#endif
      max2 = max - g;
    }
#ifdef LONG32 
    tb[r*i + dbslen/16] = t;
#endif
#ifdef LONG64
    tb[r*i + dbslen/32] = t;
#endif
    t = 0x0;
    temp[dbslen] = max;
  }
  *smax  = rmax;
  *sav_i = ri;
  *sav_j = rj;
  free(temp);
}

void backfill(long *tb, long tb_i, long tb_j, long *dbl, long *dbs, long index,
              long **algl, long **algm, long **algs, long dbslen)
{
  char *algl_c = (char *)(*algl);
  char *algm_c = (char *)(*algm);
  char *algs_c = (char *)(*algs);

  if(tb_i==0 || tb_j==0 || index<0) return;

#ifdef LONG32
  switch((unsigned long)(tb[(dbslen/16 +1)*tb_i + tb_j/16] << 2*(tb_j%16)) >> 30)
#endif
#ifdef LONG64
  switch((unsigned long)(tb[(dbslen/32 +1)*tb_i + tb_j/32] << 2*(tb_j%32)) >> 62)
#endif
  {
    case 0x0:  /* end  */
               tb_i = 0;
	       tb_j = 0;
	       break;
    case 0x1:  /* horizontal (fortran), vertical (c) */
               algl_c[index] = ((char *)dbl)[tb_i-1];
               algm_c[index] = ' ';
               algs_c[index] = '-';
	       tb_i--;
               break;
    case 0x2:  /* vertical (fortran), horizontal (c) */
               algl_c[index] = '-';
               algm_c[index] = ' ';
               algs_c[index] = ((char *)dbs)[tb_j-1];
	       tb_j--;
	       break;
    default: /* corner */
             algl_c[index] = ((char *)dbl)[tb_i-1];
             algs_c[index] = ((char *)dbs)[tb_j-1];
	     if(algl_c[index] == algs_c[index])
               algm_c[index] = ':';
	     else
	       algm_c[index] = ' ';
	     tb_i--;
	     tb_j--;
             break;
  }

  index--;

  backfill(tb, tb_i, tb_j, dbl, dbs, index, algl, algm, algs, dbslen);
}

void gotoh_align(long *tb, long *dbl, long dbllen, long *dbs, long dbslen,
                 long smax, long **algl, long **algm, long **algs, long *alglen,
                 long *algstl, long *algsts, long *errno, long sav_i, long sav_j)
{
  long i, index=0, j;
  register unsigned long blah, r;

#ifdef LONG32
  r = dbslen/16 +1;
#endif
#ifdef LONG64 
  r = dbslen/32 +1;
#endif

  i = sav_i;
  j = sav_j;
  *alglen = 0;
#ifdef LONG32
  while(blah=((unsigned long)(tb[r*i + j/16] << 2*(j%16)) >> 30))
#endif
#ifdef LONG64 
  while(blah=((unsigned long)(tb[r*i + j/32] << 2*(j%32)) >> 62))
#endif
    switch(blah)
    {
      case 0x1:  i--;
	         (*alglen)++;
		 break;
      case 0x2:  j--;
                 (*alglen)++;
                 break;
      default:   i--;
                 j--;
                 (*alglen)++;
	         break;
    }
/*printf("alglen = %d  smax = %d\n", *alglen, smax);*/

  *algstl = sav_i - *alglen;
  *algsts = sav_j - *alglen;

  *algl = (long *) malloc(*alglen+8);
  *algm = (long *) malloc(*alglen+8);
  *algs = (long *) malloc(*alglen+8);
  if(*algl==0 || *algm==0 || *algs==0)
  {
    printf("malloc error in align, exiting...\n");
    exit(0);
  }

  index = *alglen-1;

  /* recursively trace back through scoring matrix for alignment */
  backfill(tb, sav_i, sav_j, dbl, dbs, index, algl, algm, algs, dbslen);
}

