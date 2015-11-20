/* cb_sw_fw_score.c                            http://cbl.sourceforge.net
 *
 * called by cb_swa_fw, this is *not* a fully optimized implementation
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
 * $Id: cb_sw_fw_score.c,v 1.2 2004/03/10 00:22:26 jlong777 Exp $
 */
 
#include "cb_macro.h"
#include <stdio.h>
 
void cb_sw_fw_score(long *swtab, long *gaph, long *gapv, long dbllen, 
                    long dbslen, long eg, long eog, long *smax)
{
  int i, j;
  register long max, max1, max2, max3;
  register unsigned long trace;
  register long di, di1, e, g;
  
  e = eg;
  g = eog;
  *smax = 0;
  
  /* init gap penalties */
  for(i=0; i<dbslen+1; i++) 
    gapv[i] = gaph[i] = -g - (i-1)*e;
  
  for(i=1; i<dbllen+1; i++)
  {
    di = (dbslen+1)*i;
    di1 = (dbslen+1)*(i-1);
    max2 = 0 - g;
    for(j=1; j<dbslen+1; j++)
    {
      max =   0x0;
      trace = 0x0;
#ifdef LONG32
      max1 = (0x0FFFFFFF & swtab[di1+j-1]) + swtab[di+j];
      max3 = (0x0FFFFFFF & swtab[di1+j]) - g;
#endif
#ifdef LONG64
      max1 = (0x0FFFFFFFFFFFFFFF & swtab[di1+j-1]) + swtab[di+j];
      max3 = (0x0FFFFFFFFFFFFFFF & swtab[di1+j]) - g;
#endif

      gaph[j] = (max2 > gaph[j-1]-e)? max2: gaph[j-1]-e;
      gapv[j] = (max3 > gapv[j]-e)? max3: gapv[j]-e;

      if(max1>max) max = max1;
      if(gaph[j]>max) max = gaph[j];
      if(gapv[j]>max) max = gapv[j];
      if(max>*smax) *smax = max;

      if(max)
      {
#ifdef LONG32
        if(max==max1)    trace |= 0x40000000; /* came from corner */
        if(max==gaph[j]) trace |= 0x20000000; /* vert (fortran), horiz (c) */
        if(max==gapv[j]) trace |= 0x10000000; /* horiz (fortran), vert (c) */
#ifdef DEBUG_INFO
        if(max > 0x0FFFFFFF) printf("cb_sw_fw_score: cell score exceeds maximum...\n");
#endif
#endif
#ifdef LONG64
        if(max==max1)    trace |= 0x4000000000000000; /* came from corner */
        if(max==gaph[j]) trace |= 0x2000000000000000; /* vert (fortran), horiz (c) */
        if(max==gapv[j]) trace |= 0x1000000000000000; /* horiz (fortran), vert (c) */
#ifdef DEBUG_INFO
        if(max > 0x0FFFFFFFFFFFFFFF) printf("cb_sw_fw_score: cell score exceeds maximum...\n");
#endif
#endif
      }
      swtab[di+j] = max | trace;
      max2 = max - g;
    }
  }
}
