/* cb_searchn.c                                  http://cbl.sourceforge.net
 *
 * search for an approximate substring in a bigger string
 *
 * Copyright (C) 2003 University of Alaska Fairbanks
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
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
 *   Lesser General Public License for more details.
 *  
 *   You should have received a copy of the GNU Lesser General Public
 *   License along with this library; if not, write to the
 *   Free Software Foundation, Inc.
 *   59 Temple Place, Suite 330
 *   Boston, MA  02111-1307 USA
 *
 * $Id: cb_searchn.c,v 1.6 2003/12/03 21:26:02 jlong777 Exp $
 */
 
#include "cb_macro.h"
#include <stdio.h>

void cb_searchn(long *db, long dblen, long *test, long testlen,
                long *found, long *foundlen, long threshold, long *count)
{
  long i, j=0, k=0, foundlen_in, shift, m, m2, mod;
  long j2, k2=0, cutoff;
  unsigned long num[VECTLEN]={0}, xor0[VECTLEN], xor1[VECTLEN], num1, xor2;
  
  struct
  {
    long i, n, s;
  } cand[2*VECTLEN];
  
  foundlen_in = *foundlen;
  *foundlen = 0;
  
  if(dblen==0 || testlen==0 || dblen<=testlen || testlen<=threshold || threshold<0)
  {
    printf("cb_searchn: violation of dblen > testlen > threshold >= 0, returning...\n");
    return;
  }
   
  for(m=0; m<2*VECTLEN; m++)
  {
    cand[m].i = -1;
    cand[m].n = -1;
    cand[m].s = -1;
  }

#ifdef B_ENDIAN

#ifdef LONG32

  mod = testlen%16;
  
  /* cutoff is max index of test for 1st cut, benchmark sensitive.
   * 3*threshold/2 empirically determined with random test string of length
   * 197
   */
  if((((3*threshold/2)/16) > testlen/16) || (testlen/16==0))
    cutoff = testlen/16;
  else
    cutoff =((3*threshold/2)/16)+1;
  
  for(i=0; i<=(dblen-testlen+15)/16 - VECTLEN; i+=VECTLEN)
  {
    /* shift = 0 case (some compilers don't like whole word shifts) */
    if(cutoff==0)
    {
      for(m=0; m<VECTLEN; m++)
      {
        xor0[m] = test[0] ^ db[i+m]; /* shift=0 */
	
	xor1[m] = ((xor0[m] | (xor0[m]>>1)) & 0x55555555) >> (32-2*mod);
        xor0[m] =  (xor1[m] + (xor1[m]>>2)) & 0x33333333;
        xor1[m] =  (xor0[m] + (xor0[m]>>4)) & 0x07070707;
        xor0[m] =  (xor1[m] + (xor1[m]>>8));
        num[m] +=  (xor0[m] + (xor0[m]>>16))& 0x0000001F;
      }
    }
    else
    {
      for(j=0; j<cutoff; j++) /* first cut */
      {
#pragma vdir nodep
        /* compare test[] against db */
	for(m=0; m<VECTLEN; m++)
	{
          xor0[m] = test[j] ^ db[i+m+j]; /* shift=0 */
	  
	  xor1[m] = (xor0[m] | (xor0[m]>>1)) & 0x55555555;
          xor0[m] = (xor1[m] + (xor1[m]>>2)) & 0x33333333;
          xor1[m] = (xor0[m] + (xor0[m]>>4)) & 0x07070707;
          xor0[m] = (xor1[m] + (xor1[m]>>8));
          num[m] += (xor0[m] + (xor0[m]>>16))& 0x0000001F;
	}
      }
    }
    for(m=0; m<VECTLEN; m++) /* save candidates */
    {
      if(num[m] <= threshold)
      {
	cand[k2].i = i+m;
	cand[k2].n = num[m];
	cand[k2++].s = 0;
      }
      num[m] = 0;
    }
      
    if(k2>=VECTLEN) /* there are enough candidates for vector hardware to process */
    {
      k2 = k2 - VECTLEN;
      
      for(j2=j; j2<testlen/16; j2++) /* finish */
      {
	for(m=0; m<VECTLEN; m++)
	{
          if(cand[m].n > threshold) continue;
          if(cand[m].s)
	    xor0[m] = test[j2] ^ ((db[cand[m].i+j2] << cand[m].s) | 
	            (unsigned long)db[cand[m].i+j2+1] >> (32-cand[m].s));
	  else
	    xor0[m] = test[j2] ^ db[cand[m].i+j2];
	  xor1[m] = (xor0[m] | (xor0[m]>>1)) & 0x55555555;
          xor0[m] = (xor1[m] + (xor1[m]>>2)) & 0x33333333;
          xor1[m] = (xor0[m] + (xor0[m]>>4)) & 0x07070707;
          xor0[m] = (xor1[m] + (xor1[m]>>8));
          cand[m].n += (xor0[m] + (xor0[m]>>16))& 0x0000001F;
	}
      }
      if(mod && cutoff!=0) /* now do last partial word of test if there is any */
      {
	for(m=0; m<VECTLEN; m++)
        {
          if(cand[m].n > threshold) continue;
          if(cand[m].s)
	    xor0[m] = test[j2] ^ ((db[cand[m].i+j2+1] << cand[m].s) | 
	            (unsigned long)db[cand[m].i+j2+2] >> (32-cand[m].s));
	  else
	    xor0[m] = test[j2] ^ db[cand[m].i+j2+1];
	  xor1[m] = ((xor0[m] | (xor0[m]>>1)) & 0x55555555) >> (32-2*mod);
          xor0[m] =  (xor1[m] + (xor1[m]>>2)) & 0x33333333;
          xor1[m] =  (xor0[m] + (xor0[m]>>4)) & 0x07070707;
          xor0[m] =  (xor1[m] + (xor1[m]>>8));
          cand[m].n += (xor0[m] + (xor0[m]>>16))& 0x0000001F;
        }
      }
      for(m2=0; m2<VECTLEN; m2++)
      {
        if(cand[m2].n <= threshold)
        {
          if(k < foundlen_in)
          {
            found[k] = 16*cand[m2].i+cand[m2].s/2;
            count[k] = (long)cand[m2].n;
	    *foundlen = ++k;
          }
          else
          {
            *foundlen *= -1;
            return;
          }
        }
      }

      for(m2=0; m2<VECTLEN; m2++)
      {
	if(cand[m2+VECTLEN].i>=0)
	{
	  cand[m2].i = cand[m2+VECTLEN].i;
	  cand[m2+VECTLEN].i = -1;
	  cand[m2].n = cand[m2+VECTLEN].n;
	  cand[m2+VECTLEN].n = -1;
	  cand[m2].s = cand[m2+VECTLEN].s;
	  cand[m2+VECTLEN].s = -1;
	}
	else
	{
	  cand[m2].i = -1;
	  cand[m2].n = -1;
	  cand[m2].s = -1;
	}
      }
    }
  
    for(shift=2; shift<32; shift+=2) /* slide start of test along db */
    {
      if(cutoff==0)
      {
        for(m=0; m<VECTLEN; m++)
	{
          xor0[m] = test[0] ^ ((db[i+m] << shift) | 
	         (unsigned long)db[i+m+1] >> (32-shift));
	  
	  xor1[m] = ((xor0[m] | (xor0[m]>>1)) & 0x55555555) >> (32-2*mod);
          xor0[m] =  (xor1[m] + (xor1[m]>>2)) & 0x33333333;
          xor1[m] =  (xor0[m] + (xor0[m]>>4)) & 0x07070707;
          xor0[m] =  (xor1[m] + (xor1[m]>>8));
          num[m] +=  (xor0[m] + (xor0[m]>>16))& 0x0000001F;
	}
      }
      else
      {
        for(j=0; j<cutoff; j++) /* first cut */
        {
#pragma vdir nodep
          /* compare test[] against db */
	  for(m=0; m<VECTLEN; m++)
	  {
            xor0[m] = test[j] ^ ((db[i+m+j] << shift) | 
	           (unsigned long)db[i+m+j+1] >> (32-shift));
	    
	    xor1[m] = (xor0[m] | (xor0[m]>>1)) & 0x55555555;
            xor0[m] = (xor1[m] + (xor1[m]>>2)) & 0x33333333;
            xor1[m] = (xor0[m] + (xor0[m]>>4)) & 0x07070707;
            xor0[m] = (xor1[m] + (xor1[m]>>8));
            num[m] += (xor0[m] + (xor0[m]>>16))& 0x0000001F;
	  }
        }
      }
      for(m=0; m<VECTLEN; m++) /* save candidates */
      {
        if(num[m] <= threshold)
        {
	  cand[k2].i = i+m;
	  cand[k2].n = num[m];
	  cand[k2++].s = shift;
	}
	num[m] = 0;
      }
      
      if(k2>=VECTLEN) /* there are enough candidates for vector hardware to process */
      {
        k2 = k2 - VECTLEN;
      
	for(j2=j; j2<testlen/16; j2++) /* finish */
        {
	  for(m=0; m<VECTLEN; m++)
	  {
            if(cand[m].n > threshold) continue;
            if(cand[m].s)
	      xor0[m] = test[j2] ^ ((db[cand[m].i+j2] << cand[m].s) | 
	              (unsigned long)db[cand[m].i+j2+1] >> (32-cand[m].s));
	    else
	      xor0[m] = test[j2] ^ db[cand[m].i+j2];
	    xor1[m] = (xor0[m] | (xor0[m]>>1)) & 0x55555555;
            xor0[m] = (xor1[m] + (xor1[m]>>2)) & 0x33333333;
            xor1[m] = (xor0[m] + (xor0[m]>>4)) & 0x07070707;
            xor0[m] = (xor1[m] + (xor1[m]>>8));
            cand[m].n += (xor0[m] + (xor0[m]>>16))& 0x0000001F;
	  }
	}
	if(mod && cutoff!=0) /* now do last partial word of test if there is any */
        {
	  for(m=0; m<VECTLEN; m++)
          {
            if(cand[m].n > threshold) continue;
            if(cand[m].s)
	      xor0[m] = test[j2] ^ ((db[cand[m].i+j2+1] << cand[m].s) | 
	              (unsigned long)db[cand[m].i+j2+2] >> (32-cand[m].s));
	    else
	      xor0[m] = test[j2] ^ db[cand[m].i+j2+1];
	    xor1[m] = ((xor0[m] | (xor0[m]>>1)) & 0x55555555) >> (32-2*mod);
            xor0[m] =  (xor1[m] + (xor1[m]>>2)) & 0x33333333;
            xor1[m] =  (xor0[m] + (xor0[m]>>4)) & 0x07070707;
            xor0[m] =  (xor1[m] + (xor1[m]>>8));
            cand[m].n += (xor0[m] + (xor0[m]>>16))& 0x0000001F;
          }
	}
        for(m2=0; m2<VECTLEN; m2++)
        {
          if(cand[m2].n <= threshold)
          {
            if(k < foundlen_in)
            {
              found[k] = 16*cand[m2].i+cand[m2].s/2;
              count[k] = (long)cand[m2].n;
	      *foundlen = ++k;
            }
            else
            {
              *foundlen *= -1;
              return;
            }
          }
        }

	for(m2=0; m2<VECTLEN; m2++)
	{
	  if(cand[m2+VECTLEN].i>=0)
	  {
	    cand[m2].i = cand[m2+VECTLEN].i;
	    cand[m2+VECTLEN].i = -1;
	    cand[m2].n = cand[m2+VECTLEN].n;
	    cand[m2+VECTLEN].n = -1;
	    cand[m2].s = cand[m2+VECTLEN].s;
	    cand[m2+VECTLEN].s = -1;
	  }
	  else
	  {
	    cand[m2].i = -1;
	    cand[m2].n = -1;
	    cand[m2].s = -1;
	  }
	}
      }
    }
  }
  
  if(cand[0].i>=0 && *foundlen>=0) /* make sure remaining candidates get prosessed */
  {
    for(j2=j; j2<testlen/16; j2++) /* first cut */
    {
      for(m=0; m<VECTLEN; m++)
      {
        if(cand[m].n > threshold) continue;
        if(cand[m].i>=0)
        {
          if(cand[m].s)
	    xor0[m] = test[j2] ^ ((db[cand[m].i+j2] << cand[m].s) | 
	            (unsigned long)db[cand[m].i+j2+1] >> (32-cand[m].s));
	  else
	    xor0[m] = test[j2] ^ db[cand[m].i+j2];
	  xor1[m] = (xor0[m] | (xor0[m]>>1)) & 0x55555555;
          xor0[m] = (xor1[m] + (xor1[m]>>2)) & 0x33333333;
          xor1[m] = (xor0[m] + (xor0[m]>>4)) & 0x07070707;
          xor0[m] = (xor1[m] + (xor1[m]>>8));
          cand[m].n += (xor0[m] + (xor0[m]>>16))& 0x0000001F;
        }
      }
    }
    if(mod && cutoff!=0) /* now do last partial word of test if there is any */
    {
      for(m=0; m<VECTLEN; m++)
      {
        if(cand[m].n > threshold) continue;
        if(cand[m].i>=0)
        {
	  if(cand[m].s)
	    xor0[m] = test[j2] ^ ((db[cand[m].i+j2+1] << cand[m].s) | 
	            (unsigned long)db[cand[m].i+j2+2] >> (32-cand[m].s));
	  else
	    xor0[m] = test[j2] ^ db[cand[m].i+j2+1];
	  xor1[m] = ((xor0[m] | (xor0[m]>>1)) & 0x55555555) >> (32-2*mod);
          xor0[m] =  (xor1[m] + (xor1[m]>>2)) & 0x33333333;
          xor1[m] =  (xor0[m] + (xor0[m]>>4)) & 0x07070707;
          xor0[m] =  (xor1[m] + (xor1[m]>>8));
          cand[m].n += (xor0[m] + (xor0[m]>>16))& 0x0000001F;
	}
      }
    }
    for(m2=0; m2<VECTLEN; m2++)
    {
      if(cand[m2].n<=threshold && cand[m2].i>=0)
      {
        if(k < foundlen_in)
        {
          found[k] = 16*cand[m2].i+cand[m2].s/2;
          count[k] = (long)cand[m2].n;
	  *foundlen = ++k;
        }
        else
        {
          *foundlen *= -1;
          return;
        }
      }
    }
  }

  /* do the final segment of db that was not a multiple of the vector length */
  if(i)
    i = i+1;

  while(i<=(dblen-testlen+15)/16)
  {
    num1=0;
    if((16*i) > (dblen-testlen)) break; /* end of db? */
#pragma vdir nodep
    for(j=0; j<(testlen)/16; j++)
    {
      xor2 = test[j] ^ db[i+j]; /* shift=0 */
      
      xor2 =  (xor2 | (xor2>>1)) & 0x55555555;
      xor2 =  (xor2 + (xor2>>2)) & 0x33333333;
      xor2 =  (xor2 + (xor2>>4)) & 0x07070707;
      xor2 =  (xor2 + (xor2>>8));
      num1 += (xor2 + (xor2>>16))& 0x0000001F;
      
      if(num1 > threshold) break;
    }

    if(mod && (num1 <= threshold))
    {
      xor2 = test[j] ^ db[i+j]; /* shift=0 */
      
      xor2 = ((xor2 | (xor2>>1)) & 0x55555555) >> (32-2*mod); 
      xor2 =  (xor2 + (xor2>>2)) & 0x33333333;
      xor2 =  (xor2 + (xor2>>4)) & 0x07070707;
      xor2 =  (xor2 + (xor2>>8));
      num1 += (xor2 + (xor2>>16))& 0x0000001F;
    }

    if(num1 <= threshold)
    {
      if(k < foundlen_in)
      {
        found[k] = 16*i;
        count[k] = num1;
        *foundlen = ++k;
      }
      else
      {
        *foundlen *= -1;
        return;
      }
    }

    for(shift=2; shift<32; shift+=2) /* slide start of test along db */
    {
      num1=0;
      if((16*i+shift/2) > (dblen-testlen)) {i++; break;} /* end of db? */

      for(j=0; j<testlen/16; j++)
      {
        /* compare test[] against db */
	xor2 = test[j] ^ ((db[i+j] << shift) | 
	    (unsigned long)db[i+j+1] >> (32-shift));
	/* the next line makes number of 1s = number of mismatches */
        xor2 =  (xor2 | (xor2>>1)) & 0x55555555; 
        /* now add up the 1s */
        xor2 =  (xor2 + (xor2>>2)) & 0x33333333;
        xor2 =  (xor2 + (xor2>>4)) & 0x07070707;
        xor2 =  (xor2 + (xor2>>8));
        num1 += (xor2 + (xor2>>16))& 0x0000001F;
          
	if(num1 > threshold) break;
      }

      /* now do last partial word of test if there is any */
      if(mod && (num1 <= threshold))
      {
	xor2 = test[j] ^ ((db[i+j] << shift) | 
	    (unsigned long)db[i+j+1] >> (32-shift));
	xor2 = ((xor2 | (xor2>>1)) & 0x55555555) >> (32-2*mod); 
        xor2 =  (xor2 + (xor2>>2)) & 0x33333333;
        xor2 =  (xor2 + (xor2>>4)) & 0x07070707;
        xor2 =  (xor2 + (xor2>>8));
        num1 += (xor2 + (xor2>>16))& 0x0000001F;
      }

      if(num1 > threshold) continue;
      
      if(k < foundlen_in)
      {
        found[k] = 16*i+shift/2;
        count[k] = (long)num1;
	*foundlen = ++k;
      }
      else
      {
        *foundlen *= -1;
        return;
      }
    }
    i++;
  }

#endif /* 32-bit */

#ifdef LONG64

  mod = testlen%32;
  
  /* cutoff is max index of test for 1st cut, benchmark sensitive.
   * +3*threshold/4 empirically determined with random test string of length
   * 197
   */
  if((((threshold+3*threshold/4)/32) > testlen/32) || (testlen/32==0))
    cutoff = testlen/32;
  else
    cutoff =((threshold+3*threshold/4)/32)+1;
  
  for(i=0; i<=(dblen-testlen+31)/32 - VECTLEN; i+=VECTLEN)
  {
    /* shift = 0 case (some compilers don't like whole word shifts) */
    if(cutoff==0)
    {
      for(m=0; m<VECTLEN; m++)
      {
        xor0[m] = test[0] ^ db[i+m]; /* shift=0 */
	
        xor1[m] = ((xor0[m] | (xor0[m]>>1)) & 0x5555555555555555) >> (64-2*mod); 
	xor0[m] =  (xor1[m] + (xor1[m]>>2)) & 0x3333333333333333;
        xor1[m] =  (xor0[m] + (xor0[m]>>4)) & 0x0707070707070707;
        xor0[m] =  (xor1[m] + (xor1[m]>>8));
        xor1[m] =  (xor0[m] + (xor0[m]>>16));
        num[m] +=  (xor1[m] + (xor1[m]>>32))& 0x000000000000003F;
      }
    }
    else
    {
      for(j=0; j<cutoff; j++) /* first cut */
      {
#pragma vdir nodep
        /* compare test[] against db */
	for(m=0; m<VECTLEN; m++)
	{
          xor0[m] = test[j] ^ db[i+m+j]; /* shift=0 */
	  
          xor1[m] = (xor0[m] | (xor0[m]>>1)) & 0x5555555555555555; 
	  xor0[m] = (xor1[m] + (xor1[m]>>2)) & 0x3333333333333333;
          xor1[m] = (xor0[m] + (xor0[m]>>4)) & 0x0707070707070707;
          xor0[m] = (xor1[m] + (xor1[m]>>8));
          xor1[m] = (xor0[m] + (xor0[m]>>16));
          num[m] += (xor1[m] + (xor1[m]>>32))& 0x000000000000003F;
	}
      }
    }
    for(m=0; m<VECTLEN; m++) /* save candidates */
    {
      if(num[m] <= threshold)
      {
	cand[k2].i = i+m;
	cand[k2].n = num[m];
	cand[k2++].s = 0;
      }
      num[m] = 0;
    }
      
    if(k2>=VECTLEN) /* there are enough candidates for vector hardware to process */
    {
      k2 = k2 - VECTLEN;
      
      for(j2=j; j2<testlen/32; j2++) /* finish */
      {
	for(m=0; m<VECTLEN; m++)
	{
	  if(cand[m].s)
	    xor0[m] = test[j2] ^ ((db[cand[m].i+j2] << cand[m].s) | 
	            (unsigned long)db[cand[m].i+j2+1] >> (64-cand[m].s));
          else
	    xor0[m] = test[j2] ^ db[cand[m].i+j2];
	  xor1[m] = (xor0[m] | (xor0[m]>>1)) & 0x5555555555555555; 
	  xor0[m] = (xor1[m] + (xor1[m]>>2)) & 0x3333333333333333;
          xor1[m] = (xor0[m] + (xor0[m]>>4)) & 0x0707070707070707;
          xor0[m] = (xor1[m] + (xor1[m]>>8));
          xor1[m] = (xor0[m] + (xor0[m]>>16));
          cand[m].n += (xor1[m] + (xor1[m]>>32))& 0x000000000000003F;
	}
      }
      if(mod && cutoff!=0) /* now do last partial word of test if there is any */
      {
	for(m=0; m<VECTLEN; m++)
        {
	  if(cand[m].s)
            xor0[m] = test[j2] ^ ((db[cand[m].i+j2+1] << cand[m].s) | 
	            (unsigned long)db[cand[m].i+j2+2] >> (64-cand[m].s));
          else
	    xor0[m] = test[j2] ^ db[cand[m].i+j2+1];
	  xor1[m] = ((xor0[m] | (xor0[m]>>1)) & 0x5555555555555555) >> (64-2*mod);
	  xor0[m] =  (xor1[m] + (xor1[m]>>2)) & 0x3333333333333333;
          xor1[m] =  (xor0[m] + (xor0[m]>>4)) & 0x0707070707070707;
          xor0[m] =  (xor1[m] + (xor1[m]>>8));
          xor1[m] =  (xor0[m] + (xor0[m]>>16));
          cand[m].n += (xor1[m] + (xor1[m]>>32))& 0x000000000000003F;
        }
      }
      for(m2=0; m2<VECTLEN; m2++)
      {
        if(cand[m2].n <= threshold)
        {
          if(k < foundlen_in)
          {
            found[k] = 32*cand[m2].i+cand[m2].s/2;
            count[k] = (long)cand[m2].n;
	    *foundlen = ++k;
          }
          else
          {
            *foundlen *= -1;
            return;
          }
        }
      }

      for(m2=0; m2<VECTLEN; m2++)
      {
	if(cand[m2+VECTLEN].i>=0)
	{
	  cand[m2].i = cand[m2+VECTLEN].i;
	  cand[m2+VECTLEN].i = -1;
	  cand[m2].n = cand[m2+VECTLEN].n;
	  cand[m2+VECTLEN].n = -1;
	  cand[m2].s = cand[m2+VECTLEN].s;
	  cand[m2+VECTLEN].s = -1;
	}
	else
	{
	  cand[m2].i = -1;
	  cand[m2].n = -1;
	  cand[m2].s = -1;
	}
      }
    }
  
    for(shift=2; shift<64; shift+=2) /* slide start of test along db */
    {
      if(cutoff==0)
      {
        for(m=0; m<VECTLEN; m++)
	{
          xor0[m] = test[0] ^ ((db[i+m] << shift) | 
	         (unsigned long)db[i+m+1] >> (64-shift));
		 
          xor1[m] = ((xor0[m] | (xor0[m]>>1)) & 0x5555555555555555) >> (64-2*mod);
	  xor0[m] =  (xor1[m] + (xor1[m]>>2)) & 0x3333333333333333;
          xor1[m] =  (xor0[m] + (xor0[m]>>4)) & 0x0707070707070707;
          xor0[m] =  (xor1[m] + (xor1[m]>>8));
          xor1[m] =  (xor0[m] + (xor0[m]>>16));
          num[m] +=  (xor1[m] + (xor1[m]>>32))& 0x000000000000003F;
	}
      }
      else
      {
        for(j=0; j<cutoff; j++) /* first cut */
        {
#pragma vdir nodep
          /* compare test[] against db */
	  for(m=0; m<VECTLEN; m++)
	  {
            xor0[m] = test[j] ^ ((db[i+m+j] << shift) | 
	           (unsigned long)db[i+m+j+1] >> (64-shift));
            
	    xor1[m] = (xor0[m] | (xor0[m]>>1)) & 0x5555555555555555; 
	    xor0[m] = (xor1[m] + (xor1[m]>>2)) & 0x3333333333333333;
            xor1[m] = (xor0[m] + (xor0[m]>>4)) & 0x0707070707070707;
            xor0[m] = (xor1[m] + (xor1[m]>>8));
            xor1[m] = (xor0[m] + (xor0[m]>>16));
            num[m] += (xor1[m] + (xor1[m]>>32))& 0x000000000000003F;
	  }
        }
      }
      for(m=0; m<VECTLEN; m++) /* save candidates */
      {
        if(num[m] <= threshold)
        {
	  cand[k2].i = i+m;
	  cand[k2].n = num[m];
	  cand[k2++].s = shift;
	}
	num[m] = 0;
      }
      
      if(k2>=VECTLEN) /* there are enough candidates for vector hardware to process */
      {
        k2 = k2 - VECTLEN;
      
	for(j2=j; j2<testlen/32; j2++) /* finish */
        {
	  for(m=0; m<VECTLEN; m++)
	  {
	    if(cand[m].s)
	      xor0[m] = test[j2] ^ ((db[cand[m].i+j2] << cand[m].s) | 
	              (unsigned long)db[cand[m].i+j2+1] >> (64-cand[m].s));
            else
	      xor0[m] = test[j2] ^ db[cand[m].i+j2];
	    xor1[m] = (xor0[m] | (xor0[m]>>1)) & 0x5555555555555555; 
	    xor0[m] = (xor1[m] + (xor1[m]>>2)) & 0x3333333333333333;
            xor1[m] = (xor0[m] + (xor0[m]>>4)) & 0x0707070707070707;
            xor0[m] = (xor1[m] + (xor1[m]>>8));
            xor1[m] = (xor0[m] + (xor0[m]>>16));
            cand[m].n += (xor1[m] + (xor1[m]>>32))& 0x000000000000003F;
	  }
	}
	if(mod && cutoff!=0) /* now do last partial word of test if there is any */
        {
	  for(m=0; m<VECTLEN; m++)
          {
	    if(cand[m].s)
              xor0[m] = test[j2] ^ ((db[cand[m].i+j2+1] << cand[m].s) | 
	              (unsigned long)db[cand[m].i+j2+2] >> (64-cand[m].s));
            else
	      xor0[m] = test[j2] ^ db[cand[m].i+j2+1];
	    xor1[m] = ((xor0[m] | (xor0[m]>>1)) & 0x5555555555555555) >> (64-2*mod);
	    xor0[m] =  (xor1[m] + (xor1[m]>>2)) & 0x3333333333333333;
            xor1[m] =  (xor0[m] + (xor0[m]>>4)) & 0x0707070707070707;
            xor0[m] =  (xor1[m] + (xor1[m]>>8));
            xor1[m] =  (xor0[m] + (xor0[m]>>16));
            cand[m].n += (xor1[m] + (xor1[m]>>32))& 0x000000000000003F;
          }
	}
        for(m2=0; m2<VECTLEN; m2++)
        {
          if(cand[m2].n <= threshold)
          {
            if(k < foundlen_in)
            {
              found[k] = 32*cand[m2].i+cand[m2].s/2;
              count[k] = (long)cand[m2].n;
	      *foundlen = ++k;
            }
            else
            {
              *foundlen *= -1;
              return;
            }
          }
        }

	for(m2=0; m2<VECTLEN; m2++)
	{
	  if(cand[m2+VECTLEN].i>=0)
	  {
	    cand[m2].i = cand[m2+VECTLEN].i;
	    cand[m2+VECTLEN].i = -1;
	    cand[m2].n = cand[m2+VECTLEN].n;
	    cand[m2+VECTLEN].n = -1;
	    cand[m2].s = cand[m2+VECTLEN].s;
	    cand[m2+VECTLEN].s = -1;
	  }
	  else
	  {
	    cand[m2].i = -1;
	    cand[m2].n = -1;
	    cand[m2].s = -1;
	  }
	}
      }
    }
  }
  
  if(cand[0].i>=0 && *foundlen>=0) /* make sure remaining candidates get prosessed */
  {
    for(j2=j; j2<testlen/32; j2++) /* first cut */
    {
      for(m=0; m<VECTLEN; m++)
      {
        if(cand[m].i>=0)
        {
	  if(cand[m].s)
            xor0[m] = test[j2] ^ ((db[cand[m].i+j2] << cand[m].s) | 
	            (unsigned long)db[cand[m].i+j2+1] >> (64-cand[m].s));
          else
	    xor0[m] = test[j2] ^ db[cand[m].i+j2];
	  xor1[m] =  (xor0[m] | (xor0[m]>>1)) & 0x5555555555555555; 
	  xor0[m] =  (xor1[m] + (xor1[m]>>2)) & 0x3333333333333333;
          xor1[m] =  (xor0[m] + (xor0[m]>>4)) & 0x0707070707070707;
          xor0[m] =  (xor1[m] + (xor1[m]>>8));
          xor1[m] =  (xor0[m] + (xor0[m]>>16));
          cand[m].n += (xor1[m] + (xor1[m]>>32))& 0x000000000000003F;
        }
      }
    }
    if(mod && cutoff!=0) /* now do last partial word of test if there is any */
    {
      for(m=0; m<VECTLEN; m++)
      {
        if(cand[m].i>=0)
        {
	  if(cand[m].s)
            xor0[m] = test[j2] ^ ((db[cand[m].i+j2+1] << cand[m].s) | 
	            (unsigned long)db[cand[m].i+j2+2] >> (64-cand[m].s));
          else
	    xor0[m] = test[j2] ^ db[cand[m].i+j2+1];
	  xor1[m] = ((xor0[m] | (xor0[m]>>1)) & 0x5555555555555555) >> (64-2*mod);
	  xor0[m] =  (xor1[m] + (xor1[m]>>2)) & 0x3333333333333333;
          xor1[m] =  (xor0[m] + (xor0[m]>>4)) & 0x0707070707070707;
          xor0[m] =  (xor1[m] + (xor1[m]>>8));
          xor1[m] =  (xor0[m] + (xor0[m]>>16));
          cand[m].n += (xor1[m] + (xor1[m]>>32))& 0x000000000000003F;
	}
      }
    }
    for(m2=0; m2<VECTLEN; m2++)
    {
      if(cand[m2].n<=threshold && cand[m2].i>=0)
      {
        if(k < foundlen_in)
        {
          found[k] = 32*cand[m2].i+cand[m2].s/2;
          count[k] = (long)cand[m2].n;
	  *foundlen = ++k;
        }
        else
        {
          *foundlen *= -1;
          return;
        }
      }
    }
  }

  /* do the final segment of db that was not a multiple of the vector length */
  if(i)
    i = i+1;
    
  while(i<=(dblen-testlen+31)/32)
  {
    num1=0;
    if((32*i) > (dblen-testlen)) break; /* end of db? */
#pragma vdir nodep
    for(j=0; j<testlen/32; j++)
    {
      xor2 = test[j] ^ db[i+j]; /* shift=0 */
      xor2 =  (xor2 | (xor2>>1)) & 0x5555555555555555; 
      xor2 =  (xor2 + (xor2>>2)) & 0x3333333333333333;
      xor2 =  (xor2 + (xor2>>4)) & 0x0707070707070707;
      xor2 =  (xor2 + (xor2>>8));
      xor2 =  (xor2 + (xor2>>16));
      num1 += (xor2 + (xor2>>32))& 0x000000000000003F;
      
      if(num1 > threshold) break;
    }

    if(mod && (num1 <= threshold))
    {
      xor2 = test[j] ^ db[i+j]; /* shift=0 */
      xor2 = ((xor2 | (xor2>>1)) & 0x5555555555555555) >> (64-2*mod);
      xor2 =  (xor2 + (xor2>>2)) & 0x3333333333333333;
      xor2 =  (xor2 + (xor2>>4)) & 0x0707070707070707;
      xor2 =  (xor2 + (xor2>>8));
      xor2 =  (xor2 + (xor2>>16));
      num1 += (xor2 + (xor2>>32))& 0x000000000000003F;
    }

    if(num1 <= threshold)
    {
      if(k < foundlen_in)
      {
        found[k] = 32*i;
        count[k] = num1;
        *foundlen = ++k;
      }
      else
      {
        *foundlen *= -1;
        return;
      }
    }

    for(shift=2; shift<64; shift+=2) /* slide start of test along db */
    {
      num1=0;
      if((32*i+shift/2) > (dblen-testlen)) {i++; break;} /* end of db? */

      for(j=0; j<testlen/32; j++)
      {
        /* compare test[] against db */
        xor2 = test[j] ^ ((db[i+j] << shift) | 
	    (unsigned long)db[i+j+1] >> (64-shift));
	
	/* the next line makes number of 1s = number of mismatches */
        xor2 = (xor2 | (xor2>>1)) & 0x5555555555555555; 
        /* now add up the 1s */
        xor2 =  (xor2 + (xor2>>2)) & 0x3333333333333333;
        xor2 =  (xor2 + (xor2>>4)) & 0x0707070707070707;
        xor2 =  (xor2 + (xor2>>8));
        xor2 =  (xor2 + (xor2>>16));
        num1 += (xor2 + (xor2>>32))& 0x000000000000003F;
          
	if(num1 > threshold) break;
      }

      /* now do last partial word of test if there is any*/
      if(mod && (num1 <= threshold))
      {
        xor2 = test[j] ^ ((db[i+j] << shift) | 
	    (unsigned long)db[i+j+1] >> (64-shift));
        xor2 = ((xor2 | (xor2>>1)) & 0x5555555555555555) >> (64-2*mod);
        xor2 =  (xor2 + (xor2>>2)) & 0x3333333333333333;
        xor2 =  (xor2 + (xor2>>4)) & 0x0707070707070707;
        xor2 =  (xor2 + (xor2>>8));
        xor2 =  (xor2 + (xor2>>16));
        num1 += (xor2 + (xor2>>32))& 0x000000000000003F;
      }

      if(num1 > threshold) continue;
      
      if(k < foundlen_in)
      {
        found[k] = 32*i+shift/2;
        count[k] = (long)num1;
	*foundlen = ++k;
      }
      else
      {
        *foundlen *= -1;
        return;
      }
    }
    i++;
  }
#endif /* 64-bit */
#endif /* big_endian */

#ifdef L_ENDIAN

#ifdef LONG32
  
  mod = testlen%16;
  
  /* cutoff is max index of test for 1st cut, benchmark sensitive.
   * 3*threshold/2 empirically determined with random test string of length
   * 197
   */
  if((((3*threshold/2)/16) > testlen/16) || (testlen/16==0))
    cutoff = testlen/16;
  else
    cutoff =((3*threshold/2)/16)+1;
    
  for(i=0; i<=(dblen-testlen+15)/16 - VECTLEN; i+=VECTLEN)
  {
    /* shift = 0 case (some compilers don't like whole word shifts) */
    if(cutoff==0)
    {
      for(m=0; m<VECTLEN; m++)
      {
        xor0[m] = test[0] ^ db[i+m]; /* shift=0 */
	
	xor1[m] = ((xor0[m] | (xor0[m]>>1)) & 0x55555555) << (32-2*mod);
        xor0[m] =  (xor1[m] + (xor1[m]>>2)) & 0x33333333;
        xor1[m] =  (xor0[m] + (xor0[m]>>4)) & 0x07070707;
        xor0[m] =  (xor1[m] + (xor1[m]>>8));
        num[m] +=  (xor0[m] + (xor0[m]>>16))& 0x0000001F;
      }
    }
    else
    {
      for(j=0; j<cutoff; j++) /* first cut */
      {
#pragma vdir nodep
        /* compare test[] against db */
	for(m=0; m<VECTLEN; m++)
	{
          xor0[m] = test[j] ^ db[i+m+j]; /* shift=0 */
	  
	  xor1[m] = (xor0[m] | (xor0[m]>>1)) & 0x55555555;
          xor0[m] = (xor1[m] + (xor1[m]>>2)) & 0x33333333;
          xor1[m] = (xor0[m] + (xor0[m]>>4)) & 0x07070707;
          xor0[m] = (xor1[m] + (xor1[m]>>8));
          num[m] += (xor0[m] + (xor0[m]>>16))& 0x0000001F;
	}
      }
    }
    for(m=0; m<VECTLEN; m++) /* save candidates */
    {
      if(num[m] <= threshold)
      {
	cand[k2].i = i+m;
	cand[k2].n = num[m];
	cand[k2++].s = 0;
      }
      num[m] = 0;
    }
      
    if(k2>=VECTLEN) /* there are enough candidates for vector hardware to process */
    {
      k2 = k2 - VECTLEN;
      
      for(j2=j; j2<testlen/16; j2++) /* finish */
      {
	for(m=0; m<VECTLEN; m++)
	{
          if(cand[m].n > threshold) continue;
          if(cand[m].s)
	    xor0[m] = test[j2] ^ 
	              (((unsigned long)db[cand[m].i+j2] >> cand[m].s) | 
	                              (db[cand[m].i+j2+1] << (32-cand[m].s)));
          else
	    xor0[m] = test[j2] ^ db[cand[m].i+j2];
	  xor1[m] = (xor0[m] | (xor0[m]>>1)) & 0x55555555;
          xor0[m] = (xor1[m] + (xor1[m]>>2)) & 0x33333333;
          xor1[m] = (xor0[m] + (xor0[m]>>4)) & 0x07070707;
          xor0[m] = (xor1[m] + (xor1[m]>>8));
          cand[m].n += (xor0[m] + (xor0[m]>>16))& 0x0000001F;
	}
      }
      if(mod && cutoff!=0) /* now do last partial word of test if there is any */
      {
	for(m=0; m<VECTLEN; m++)
        {
          if(cand[m].n > threshold) continue;
          if(cand[m].s)
	    xor0[m] = test[j2] ^ 
	              (((unsigned long)db[cand[m].i+j2+1] >> cand[m].s) | 
	                              (db[cand[m].i+j2+2] << (32-cand[m].s)));
	  else
	    xor0[m] = test[j2] ^ db[cand[m].i+j2+1];
	  xor1[m] = ((xor0[m] | (xor0[m]>>1)) & 0x55555555) << (32-2*mod);
          xor0[m] =  (xor1[m] + (xor1[m]>>2)) & 0x33333333;
          xor1[m] =  (xor0[m] + (xor0[m]>>4)) & 0x07070707;
          xor0[m] =  (xor1[m] + (xor1[m]>>8));
          cand[m].n += (xor0[m] + (xor0[m]>>16))& 0x0000001F;
        }
      }
      for(m2=0; m2<VECTLEN; m2++)
      {
        if(cand[m2].n <= threshold)
        {
          if(k < foundlen_in)
          {
            found[k] = 16*cand[m2].i+cand[m2].s/2;
            count[k] = (long)cand[m2].n;
	    *foundlen = ++k;
          }
          else
          {
            *foundlen *= -1;
            return;
          }
        }
      }

      for(m2=0; m2<VECTLEN; m2++)
      {
	if(cand[m2+VECTLEN].i>=0)
	{
	  cand[m2].i = cand[m2+VECTLEN].i;
	  cand[m2+VECTLEN].i = -1;
	  cand[m2].n = cand[m2+VECTLEN].n;
	  cand[m2+VECTLEN].n = -1;
	  cand[m2].s = cand[m2+VECTLEN].s;
	  cand[m2+VECTLEN].s = -1;
	}
	else
	{
	  cand[m2].i = -1;
	  cand[m2].n = -1;
	  cand[m2].s = -1;
	}
      }
    }
  
    for(shift=2; shift<32; shift+=2) /* slide start of test along db */
    {
      if(cutoff==0)
      {
        for(m=0; m<VECTLEN; m++)
	{
	  xor0[m] = test[0] ^ (((unsigned long)db[i+m] >> shift) | 
	                                      (db[i+m+1] << (32-shift)));
	  
	  xor1[m] = ((xor0[m] | (xor0[m]>>1)) & 0x55555555) << (32-2*mod);
          xor0[m] =  (xor1[m] + (xor1[m]>>2)) & 0x33333333;
          xor1[m] =  (xor0[m] + (xor0[m]>>4)) & 0x07070707;
          xor0[m] =  (xor1[m] + (xor1[m]>>8));
          num[m] +=  (xor0[m] + (xor0[m]>>16))& 0x0000001F;
	}
      }
      else
      {
        for(j=0; j<cutoff; j++) /* first cut */
        {
#pragma vdir nodep
          /* compare test[] against db */
	  for(m=0; m<VECTLEN; m++)
	  {
	    xor0[m] = test[j] ^ (((unsigned long)db[i+m+j] >> shift) | 
	                                        (db[i+m+j+1] << (32-shift)));
						
            xor1[m] = (xor0[m] | (xor0[m]>>1)) & 0x55555555;
            xor0[m] = (xor1[m] + (xor1[m]>>2)) & 0x33333333;
            xor1[m] = (xor0[m] + (xor0[m]>>4)) & 0x07070707;
            xor0[m] = (xor1[m] + (xor1[m]>>8));
            num[m] += (xor0[m] + (xor0[m]>>16))& 0x0000001F;
	  }
        }
      }
      for(m=0; m<VECTLEN; m++) /* save candidates */
      {
        if(num[m] <= threshold)
        {
	  cand[k2].i = i+m;
	  cand[k2].n = num[m];
	  cand[k2++].s = shift;
	}
	num[m] = 0;
      }
      
      if(k2>=VECTLEN) /* there are enough candidates for vector hardware to process */
      {
        k2 = k2 - VECTLEN;

	for(j2=j; j2<testlen/16; j2++) /* finish */
        {
	  for(m=0; m<VECTLEN; m++)
	  {
            if(cand[m].n > threshold) continue;
            if(cand[m].s)
	      xor0[m] = test[j2] ^ 
	               (((unsigned long)db[cand[m].i+j2] >> cand[m].s) | 
	                               (db[cand[m].i+j2+1] << (32-cand[m].s)));
            else
	      xor0[m] = test[j2] ^ db[cand[m].i+j2];
	    xor1[m] = (xor0[m] | (xor0[m]>>1)) & 0x55555555;
            xor0[m] = (xor1[m] + (xor1[m]>>2)) & 0x33333333;
            xor1[m] = (xor0[m] + (xor0[m]>>4)) & 0x07070707;
            xor0[m] = (xor1[m] + (xor1[m]>>8));
            cand[m].n += (xor0[m] + (xor0[m]>>16))& 0x0000001F;
	  }
	}
	if(mod && cutoff!=0) /* now do last partial word of test if there is any */
        {
	  for(m=0; m<VECTLEN; m++)
          {
            if(cand[m].n > threshold) continue;
            if(cand[m].s)
	      xor0[m] = test[j2] ^ 
	                (((unsigned long)db[cand[m].i+j2+1] >> cand[m].s) | 
	                                (db[cand[m].i+j2+2] << (32-cand[m].s)));
            else
	      xor0[m] = test[j2] ^ db[cand[m].i+j2+1];
	    xor1[m] = ((xor0[m] | (xor0[m]>>1)) & 0x55555555) << (32-2*mod);
            xor0[m] =  (xor1[m] + (xor1[m]>>2)) & 0x33333333;
            xor1[m] =  (xor0[m] + (xor0[m]>>4)) & 0x07070707;
            xor0[m] =  (xor1[m] + (xor1[m]>>8));
            cand[m].n += (xor0[m] + (xor0[m]>>16))& 0x0000001F;
          }
	}
        for(m2=0; m2<VECTLEN; m2++)
        {
          if(cand[m2].n <= threshold)
          {
            if(k < foundlen_in)
            {
              found[k] = 16*cand[m2].i+cand[m2].s/2;
              count[k] = (long)cand[m2].n;
	      *foundlen = ++k;
            }
            else
            {
              *foundlen *= -1;
              return;
            }
          }
        }

	for(m2=0; m2<VECTLEN; m2++)
	{
	  if(cand[m2+VECTLEN].i>=0)
	  {
	    cand[m2].i = cand[m2+VECTLEN].i;
	    cand[m2+VECTLEN].i = -1;
	    cand[m2].n = cand[m2+VECTLEN].n;
	    cand[m2+VECTLEN].n = -1;
	    cand[m2].s = cand[m2+VECTLEN].s;
	    cand[m2+VECTLEN].s = -1;
	  }
	  else
	  {
	    cand[m2].i = -1;
	    cand[m2].n = -1;
	    cand[m2].s = -1;
	  }
	}
      }
    }
  }
  
  if(cand[0].i>=0 && *foundlen>=0) /* make sure remaining candidates get prosessed */
  {
    for(j2=j; j2<testlen/16; j2++) /* first cut */
    {
      for(m=0; m<VECTLEN; m++)
      {
	if(cand[m].n > threshold) continue;
	if(cand[m].i>=0)
        {
          if(cand[m].s)
	    xor0[m] = test[j2] ^ 
	           (((unsigned long)db[cand[m].i+j2] >> cand[m].s) | 
	                           (db[cand[m].i+j2+1] << (32-cand[m].s)));
          else
	    xor0[m] = test[j2] ^ db[cand[m].i+j2];
	  xor1[m] =  (xor0[m] | (xor0[m]>>1)) & 0x55555555;
          xor0[m] =  (xor1[m] + (xor1[m]>>2)) & 0x33333333;
          xor1[m] =  (xor0[m] + (xor0[m]>>4)) & 0x07070707;
          xor0[m] =  (xor1[m] + (xor1[m]>>8));
          cand[m].n += (xor0[m] + (xor0[m]>>16))& 0x0000001F;
        }
      }
    }
    if(mod && cutoff!=0) /* now do last partial word of test if there is any */
    {
      for(m=0; m<VECTLEN; m++)
      {
	if(cand[m].n > threshold) continue;
        if(cand[m].i>=0)
        {
	  if(cand[m].s)
	    xor0[m] = test[j2] ^ 
	           (((unsigned long)db[cand[m].i+j2+1] >> cand[m].s) | 
	                           (db[cand[m].i+j2+2] << (32-cand[m].s)));
          else
	    xor0[m] = test[j2] ^ db[cand[m].i+j2+1];
	  xor1[m] = ((xor0[m] | (xor0[m]>>1)) & 0x55555555) << (32-2*mod);
          xor0[m] =  (xor1[m] + (xor1[m]>>2)) & 0x33333333;
          xor1[m] =  (xor0[m] + (xor0[m]>>4)) & 0x07070707;
          xor0[m] =  (xor1[m] + (xor1[m]>>8));
          cand[m].n += (xor0[m] + (xor0[m]>>16))& 0x0000001F;
	}
      }
    }
    for(m2=0; m2<VECTLEN; m2++)
    {
      if(cand[m2].n<=threshold && cand[m2].i>=0)
      {
        if(k < foundlen_in)
        {
          found[k] = 16*cand[m2].i+cand[m2].s/2;
          count[k] = (long)cand[m2].n;
	  *foundlen = ++k;
        }
        else
        {
          *foundlen *= -1;
          return;
        }
      }
    }
  }

  /* do the final segment of db that was not a multiple of the vector length */
  if(i)
    i = i+1;
    
  while(i<=(dblen-testlen+15)/16)
  {
    num1=0;
    if((16*i) > (dblen-testlen)) break; /* end of db? */
#pragma vdir nodep
    for(j=0; j<testlen/16; j++)
    {
      xor2 = test[j] ^ db[i+j]; /* shift=0 */
      
      xor2 =  (xor2 | (xor2>>1)) & 0x55555555;
      xor2 =  (xor2 + (xor2>>2)) & 0x33333333;
      xor2 =  (xor2 + (xor2>>4)) & 0x07070707;
      xor2 =  (xor2 + (xor2>>8));
      num1 += (xor2 + (xor2>>16))& 0x0000001F;
      
      if(num1 > threshold) break;
    }

    if(mod && (num1 <= threshold))
    {
      xor2 = test[j] ^ db[i+j]; /* shift=0 */
      
      xor2 = ((xor2 | (xor2>>1)) & 0x55555555) << (32-2*mod); 
      xor2 =  (xor2 + (xor2>>2)) & 0x33333333;
      xor2 =  (xor2 + (xor2>>4)) & 0x07070707;
      xor2 =  (xor2 + (xor2>>8));
      num1 += (xor2 + (xor2>>16))& 0x0000001F;
    }

    if(num1 <= threshold)
    {
      if(k < foundlen_in)
      {
        found[k] = 16*i;
        count[k] = num1;
        *foundlen = ++k;
      }
      else
      {
        *foundlen *= -1;
        return;
      }
    }

    for(shift=2; shift<32; shift+=2) /* slide start of test along db */
    {
      num1=0;
      if((16*i+shift/2) > (dblen-testlen)) {i++; break;} /* end of db? */

      for(j=0; j<testlen/16; j++)
      {
        /* compare test[] against db */
	xor2 = test[j] ^ (((unsigned long)db[i+j] >> shift) | 
	                                 (db[i+j+1] << (32-shift)));
					 
	/* the next line makes number of 1s = number of mismatches */
        xor2 =  (xor2 | (xor2>>1)) & 0x55555555; 
        /* now add up the 1s */
        xor2 =  (xor2 + (xor2>>2)) & 0x33333333;
        xor2 =  (xor2 + (xor2>>4)) & 0x07070707;
        xor2 =  (xor2 + (xor2>>8));
        num1 += (xor2 + (xor2>>16))& 0x0000001F;
          
	if(num1 > threshold) break;
      }

      /* now do last partial word of test if there is any */
      if(mod && (num1 <= threshold))
      {
        xor2 = test[j] ^ (((unsigned long)db[i+j] >> shift) | 
                                         (db[i+j+1] << (32-shift)));
					 
        xor2 = ((xor2 | (xor2>>1)) & 0x55555555) << (32-2*mod); 
        xor2 =  (xor2 + (xor2>>2)) & 0x33333333;
        xor2 =  (xor2 + (xor2>>4)) & 0x07070707;
        xor2 =  (xor2 + (xor2>>8));
        num1 += (xor2 + (xor2>>16))& 0x0000001F;
      }

      if(num1 > threshold) continue;
      
      if(k < foundlen_in)
      {
        found[k] = 16*i+shift/2;
        count[k] = (long)num1;
	*foundlen = ++k;
      }
      else
      {
        *foundlen *= -1;
        return;
      }
    }
    i++;
  }

#endif /* 32-bit */

#ifdef LONG64

  mod = testlen%32;
  
  /* cutoff is max index of test for 1st cut, benchmark sensitive.
   * +3*threshold/4 empirically determined with random test string of length
   * 197
   */
  if((((threshold+3*threshold/4)/32) > testlen/32) || (testlen/32==0))
    cutoff = testlen/32;
  else
    cutoff =((threshold+3*threshold/4)/32)+1;
  
  for(i=0; i<=(dblen-testlen+31)/32 - VECTLEN; i+=VECTLEN)
  {
    /* shift = 0 case (some compilers don't like whole word shifts) */
    if(cutoff==0)
    {
      for(m=0; m<VECTLEN; m++)
      {
        xor0[m] = test[0] ^ db[i+m]; /* shift=0 */
	
        xor1[m] = ((xor0[m] | (xor0[m]>>1)) & 0x5555555555555555) << (64-2*mod); 
	xor0[m] =  (xor1[m] + (xor1[m]>>2)) & 0x3333333333333333;
        xor1[m] =  (xor0[m] + (xor0[m]>>4)) & 0x0707070707070707;
        xor0[m] =  (xor1[m] + (xor1[m]>>8));
        xor1[m] =  (xor0[m] + (xor0[m]>>16));
        num[m] +=  (xor1[m] + (xor1[m]>>32))& 0x000000000000003F;
      }
    }
    else
    {
      for(j=0; j<cutoff; j++) /* first cut */
      {
#pragma vdir nodep
        /* compare test[] against db */
	for(m=0; m<VECTLEN; m++)
	{
          xor0[m] = test[j] ^ db[i+m+j]; /* shift=0 */
	  
          xor1[m] = (xor0[m] | (xor0[m]>>1)) & 0x5555555555555555; 
	  xor0[m] = (xor1[m] + (xor1[m]>>2)) & 0x3333333333333333;
          xor1[m] = (xor0[m] + (xor0[m]>>4)) & 0x0707070707070707;
          xor0[m] = (xor1[m] + (xor1[m]>>8));
          xor1[m] = (xor0[m] + (xor0[m]>>16));
          num[m] += (xor1[m] + (xor1[m]>>32))& 0x000000000000003F;
	}
      }
    }
    for(m=0; m<VECTLEN; m++) /* save candidates */
    {
      if(num[m] <= threshold)
      {
	cand[k2].i = i+m;
	cand[k2].n = num[m];
	cand[k2++].s = 0;
      }
      num[m] = 0;
    }
      
    if(k2>=VECTLEN) /* there are enough candidates for vector hardware to process */
    {
      k2 = k2 - VECTLEN;
      
      for(j2=j; j2<testlen/32; j2++) /* finish */
      {
	for(m=0; m<VECTLEN; m++)
	{
	  if(cand[m].s)
	    xor0[m] = test[j2] ^ 
	           (((unsigned long)db[cand[m].i+j2] >> cand[m].s) | 
	                           (db[cand[m].i+j2+1] << (64-cand[m].s)));
          else
	    xor0[m] = test[j2] ^ db[cand[m].i+j2];
	  xor1[m] = (xor0[m] | (xor0[m]>>1)) & 0x5555555555555555; 
	  xor0[m] = (xor1[m] + (xor1[m]>>2)) & 0x3333333333333333;
          xor1[m] = (xor0[m] + (xor0[m]>>4)) & 0x0707070707070707;
          xor0[m] = (xor1[m] + (xor1[m]>>8));
          xor1[m] = (xor0[m] + (xor0[m]>>16));
          cand[m].n += (xor1[m] + (xor1[m]>>32))& 0x000000000000003F;
	}
      }
      if(mod && cutoff!=0) /* now do last partial word of test if there is any */
      {
	for(m=0; m<VECTLEN; m++)
        {
	  if(cand[m].s)
            xor0[m] = test[j2] ^ 
	           (((unsigned long)db[cand[m].i+j2+1] >> cand[m].s) | 
	                           (db[cand[m].i+j2+2] << (64-cand[m].s)));
          else
	    xor0[m] = test[j2] ^ db[cand[m].i+j2+1];
	  xor1[m] = ((xor0[m] | (xor0[m]>>1)) & 0x5555555555555555) << (64-2*mod);
	  xor0[m] =  (xor1[m] + (xor1[m]>>2)) & 0x3333333333333333;
          xor1[m] =  (xor0[m] + (xor0[m]>>4)) & 0x0707070707070707;
          xor0[m] =  (xor1[m] + (xor1[m]>>8));
          xor1[m] =  (xor0[m] + (xor0[m]>>16));
          cand[m].n += (xor1[m] + (xor1[m]>>32))& 0x000000000000003F;
        }
      }
      for(m2=0; m2<VECTLEN; m2++)
      {
        if(cand[m2].n <= threshold)
        {
          if(k < foundlen_in)
          {
            found[k] = 32*cand[m2].i+cand[m2].s/2;
            count[k] = (long)cand[m2].n;
	    *foundlen = ++k;
          }
          else
          {
            *foundlen *= -1;
            return;
          }
        }
      }

      for(m2=0; m2<VECTLEN; m2++)
      {
	if(cand[m2+VECTLEN].i>=0)
	{
	  cand[m2].i = cand[m2+VECTLEN].i;
	  cand[m2+VECTLEN].i = -1;
	  cand[m2].n = cand[m2+VECTLEN].n;
	  cand[m2+VECTLEN].n = -1;
	  cand[m2].s = cand[m2+VECTLEN].s;
	  cand[m2+VECTLEN].s = -1;
	}
	else
	{
	  cand[m2].i = -1;
	  cand[m2].n = -1;
	  cand[m2].s = -1;
	}
      }
    }
  
    for(shift=2; shift<64; shift+=2) /* slide start of test along db */
    {
      if(cutoff==0)
      {
        for(m=0; m<VECTLEN; m++)
	{
          xor0[m] = test[0] ^ (((unsigned long)db[i+m] >> shift) | 
	                                      (db[i+m+1] << (64-shift)));
					      
          xor1[m] = ((xor0[m] | (xor0[m]>>1)) & 0x5555555555555555) << (64-2*mod);
	  xor0[m] =  (xor1[m] + (xor1[m]>>2)) & 0x3333333333333333;
          xor1[m] =  (xor0[m] + (xor0[m]>>4)) & 0x0707070707070707;
          xor0[m] =  (xor1[m] + (xor1[m]>>8));
          xor1[m] =  (xor0[m] + (xor0[m]>>16));
          num[m] += (xor1[m] + (xor1[m]>>32))& 0x000000000000003F;
	}
      }
      else
      {
        for(j=0; j<cutoff; j++) /* first cut */
        {
#pragma vdir nodep
          /* compare test[] against db */
	  for(m=0; m<VECTLEN; m++)
	  {
            xor0[m] = test[j] ^ (((unsigned long)db[i+m+j] >> shift) | 
	                                        (db[i+m+j+1] << (64-shift)));
						
            xor1[m] = (xor0[m] | (xor0[m]>>1)) & 0x5555555555555555; 
	    xor0[m] = (xor1[m] + (xor1[m]>>2)) & 0x3333333333333333;
            xor1[m] = (xor0[m] + (xor0[m]>>4)) & 0x0707070707070707;
            xor0[m] = (xor1[m] + (xor1[m]>>8));
            xor1[m] = (xor0[m] + (xor0[m]>>16));
            num[m] += (xor1[m] + (xor1[m]>>32))& 0x000000000000003F;
	  }
        }
      }
      for(m=0; m<VECTLEN; m++) /* save candidates */
      {
        if(num[m] <= threshold)
        {
	  cand[k2].i = i+m;
	  cand[k2].n = num[m];
	  cand[k2++].s = shift;
	}
	num[m] = 0;
      }
      
      if(k2>=VECTLEN) /* there are enough candidates for vector hardware to process */
      {
        k2 = k2 - VECTLEN;
	
	for(j2=j; j2<testlen/32; j2++) /* finish */
        {
	  for(m=0; m<VECTLEN; m++)
	  {
            if(cand[m].s)
	      xor0[m] = test[j2] ^ 
	             (((unsigned long)db[cand[m].i+j2] >> cand[m].s) | 
	                             (db[cand[m].i+j2+1] << (64-cand[m].s)));
            else
	      xor0[m] = test[j2] ^ db[cand[m].i+j2];
	    xor1[m] =  (xor0[m] | (xor0[m]>>1)) & 0x5555555555555555; 
	    xor0[m] =  (xor1[m] + (xor1[m]>>2)) & 0x3333333333333333;
            xor1[m] =  (xor0[m] + (xor0[m]>>4)) & 0x0707070707070707;
            xor0[m] =  (xor1[m] + (xor1[m]>>8));
            xor1[m] =  (xor0[m] + (xor0[m]>>16));
            cand[m].n += (xor1[m] + (xor1[m]>>32))& 0x000000000000003F;
	  }
	}
	if(mod && cutoff!=0) /* now do last partial word of test if there is any */
	{
	  for(m=0; m<VECTLEN; m++)
          {
            if(cand[m].s)
	      xor0[m] = test[j2] ^ 
	             (((unsigned long)db[cand[m].i+j2+1] >> cand[m].s) | 
	                             (db[cand[m].i+j2+2] << (64-cand[m].s)));
            else
	      xor0[m] = test[j2] ^ db[cand[m].i+j2+1];
	    xor1[m] = ((xor0[m] | (xor0[m]>>1)) & 0x5555555555555555) << (64-2*mod);
	    xor0[m] =  (xor1[m] + (xor1[m]>>2)) & 0x3333333333333333;
            xor1[m] =  (xor0[m] + (xor0[m]>>4)) & 0x0707070707070707;
            xor0[m] =  (xor1[m] + (xor1[m]>>8));
            xor1[m] =  (xor0[m] + (xor0[m]>>16));
            cand[m].n += (xor1[m] + (xor1[m]>>32))& 0x000000000000003F;
          }
	}
        for(m2=0; m2<VECTLEN; m2++)
        {
          if(cand[m2].n <= threshold)
          {
            if(k < foundlen_in)
            {
              found[k] = 32*cand[m2].i+cand[m2].s/2;
              count[k] = (long)cand[m2].n;
	      *foundlen = ++k;
            }
	    else
            {
              *foundlen *= -1;
              return;
            }
          }
        }

	for(m2=0; m2<VECTLEN; m2++)
	{
	  if(cand[m2+VECTLEN].i>=0)
	  {
	    cand[m2].i = cand[m2+VECTLEN].i;
	    cand[m2+VECTLEN].i = -1;
	    cand[m2].n = cand[m2+VECTLEN].n;
	    cand[m2+VECTLEN].n = -1;
	    cand[m2].s = cand[m2+VECTLEN].s;
	    cand[m2+VECTLEN].s = -1;
	  }
	  else
	  {
	    cand[m2].i = -1;
	    cand[m2].n = -1;
	    cand[m2].s = -1;
	  }
	}
      }
    }
  }
  
  if(cand[0].i>=0 && *foundlen>=0) /* make sure remaining candidates get prosessed */
  {
    for(j2=j; j2<testlen/32; j2++) /* first cut */
    {
      for(m=0; m<VECTLEN; m++)
      {
        if(cand[m].i>=0)
        {
          if(cand[m].s)
            xor0[m] = test[j2] ^ 
	           (((unsigned long)db[cand[m].i+j2] >> cand[m].s) | 
	                           (db[cand[m].i+j2+1] << (64-cand[m].s)));
          else
	    xor0[m] = test[j2] ^ db[cand[m].i+j2];
	  xor1[m] =  (xor0[m] | (xor0[m]>>1)) & 0x5555555555555555; 
	  xor0[m] =  (xor1[m] + (xor1[m]>>2)) & 0x3333333333333333;
          xor1[m] =  (xor0[m] + (xor0[m]>>4)) & 0x0707070707070707;
          xor0[m] =  (xor1[m] + (xor1[m]>>8));
          xor1[m] =  (xor0[m] + (xor0[m]>>16));
          cand[m].n += (xor1[m] + (xor1[m]>>32))& 0x000000000000003F;
        }
      }
    }
    if(mod && cutoff!=0) /* now do last partial word of test if there is any */
    {
      for(m=0; m<VECTLEN; m++)
      {
        if(cand[m].i>=0)
        {
	  if(cand[m].s)
	    xor0[m] = test[j2] ^ 
	           (((unsigned long)db[cand[m].i+j2+1] >> cand[m].s) | 
	                           (db[cand[m].i+j2+2] << (64-cand[m].s)));
          else
	    xor0[m] = test[j2] ^ db[cand[m].i+j2+1];
	  xor1[m] = ((xor0[m] | (xor0[m]>>1)) & 0x5555555555555555) << (64-2*mod);
	  xor0[m] =  (xor1[m] + (xor1[m]>>2)) & 0x3333333333333333;
          xor1[m] =  (xor0[m] + (xor0[m]>>4)) & 0x0707070707070707;
          xor0[m] =  (xor1[m] + (xor1[m]>>8));
          xor1[m] =  (xor0[m] + (xor0[m]>>16));
          cand[m].n += (xor1[m] + (xor1[m]>>32))& 0x000000000000003F;
	}
      }
    }
    for(m2=0; m2<VECTLEN; m2++)
    {
      if(cand[m2].n<=threshold && cand[m2].i>=0)
      {
        if(k < foundlen_in)
        {
          found[k] = 32*cand[m2].i+cand[m2].s/2;
          count[k] = (long)cand[m2].n;
	  *foundlen = ++k;
        }
        else
        {
          *foundlen *= -1;
          return;
        }
      }
    }
  }

  /* do the final segment of db that was not a multiple of the vector length */
  if(i)
    i = i+1;
    
  while(i<=(dblen-testlen+31)/32)
  {
    num1=0;
    if((32*i) > (dblen-testlen)) break; /* end of db? */
#pragma vdir nodep
    for(j=0; j<testlen/32; j++)
    {
      xor2 = test[j] ^ db[i+j]; /* shift=0 */
      
      xor2 =  (xor2 | (xor2>>1)) & 0x5555555555555555; 
      xor2 =  (xor2 + (xor2>>2)) & 0x3333333333333333;
      xor2 =  (xor2 + (xor2>>4)) & 0x0707070707070707;
      xor2 =  (xor2 + (xor2>>8));
      xor2 =  (xor2 + (xor2>>16));
      num1 += (xor2 + (xor2>>32))& 0x000000000000003F;
      
      if(num1 > threshold) break;
    }

    if(mod && (num1 <= threshold))
    {
      xor2 = test[j] ^ db[i+j]; /* shift=0 */
      
      xor2 = ((xor2 | (xor2>>1)) & 0x5555555555555555) << (64-2*mod);
      xor2 =  (xor2 + (xor2>>2)) & 0x3333333333333333;
      xor2 =  (xor2 + (xor2>>4)) & 0x0707070707070707;
      xor2 =  (xor2 + (xor2>>8));
      xor2 =  (xor2 + (xor2>>16));
      num1 += (xor2 + (xor2>>32))& 0x000000000000003F;
    }

    if(num1 <= threshold)
    {
      if(k < foundlen_in)
      {
        found[k] = 32*i;
        count[k] = num1;
        *foundlen = ++k;
      }
      else
      {
        *foundlen *= -1;
        return;
      }
    }

    for(shift=2; shift<64; shift+=2) /* slide start of test along db */
    {
      num1=0;
      if((32*i+shift/2) > (dblen-testlen)) {i++; break;} /* end of db? */

      for(j=0; j<testlen/32; j++)
      {
        /* compare test[] against db */
	xor2 = test[j] ^ (((unsigned long)db[i+j] >> shift) | 
	                                 (db[i+j+1] << (64-shift)));
			
	/* the next line makes number of 1s = number of mismatches */
	xor2 =  (xor2 | (xor2>>1)) & 0x5555555555555555; 
        /* now add up the 1s */
        xor2 =  (xor2 + (xor2>>2)) & 0x3333333333333333;
        xor2 =  (xor2 + (xor2>>4)) & 0x0707070707070707;
        xor2 =  (xor2 + (xor2>>8));
        xor2 =  (xor2 + (xor2>>16));
        num1 += (xor2 + (xor2>>32))& 0x000000000000003F;
          
	if(num1 > threshold) break;
      }

      /* now do last partial word of test if there is any */
      if(mod && (num1 <= threshold))
      {
        xor2 = test[j] ^ (((unsigned long)db[i+j] >> shift) | 
                                         (db[i+j+1] << (64-shift)));
					 
        xor2 = ((xor2 | (xor2>>1)) & 0x5555555555555555) << (64-2*mod);
        xor2 =  (xor2 + (xor2>>2)) & 0x3333333333333333;
        xor2 =  (xor2 + (xor2>>4)) & 0x0707070707070707;
        xor2 =  (xor2 + (xor2>>8));
        xor2 =  (xor2 + (xor2>>16));
        num1 += (xor2 + (xor2>>32))& 0x000000000000003F;
      }

      if(num1 > threshold) continue;
      
      if(k < foundlen_in)
      {
        found[k] = 32*i+shift/2;
        count[k] = (long)num1;
	*foundlen = ++k;
      }
      else
      {
        *foundlen *= -1;
        return;
      }
    }
    i++;
  }
#endif /* 64-bit */
#endif /* little_endian */
}
