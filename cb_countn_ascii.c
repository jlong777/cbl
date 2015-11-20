/* cb_countn_ascii.c                             http://cbl.sourceforge.net
 *
 * counts the letters A, C, G, N, and T in an ascii string, case insensitive.
 * see original man page at the bottom. This code is the generic version
 * with provision for 32 and 64 bit long ints.
 *
 * thanks to Bill Long of CRAY for hints on this algorithm
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
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 *   Lesser General Public License for more details.
 *  
 *   You should have received a copy of the GNU Lesser General Public
 *   License along with this library; if not, write to the
 *   Free Software Foundation, Inc.
 *   59 Temple Place, Suite 330
 *   Boston, MA  02111-1307 USA
 *
 * $Id: cb_countn_ascii.c,v 1.5 2003/11/20 23:15:43 jlong777 Exp $
 */

#include "cb_macro.h"

void cb_countn_ascii(long *dba, long dblen, long res[])
{
  long i, j, length;
  register unsigned long a=0, c=0, g=0, t=0, n=0;
  register unsigned long r0, r1, r2;
  
  for(i=0; i<5; i++) res[i] = 0;
  
#ifdef LONG32
  length = (dblen+3)/4;
#endif
#ifdef LONG64
  length = (dblen+7)/8;
#endif
  
  j=0;
  while(j<length-25001)
  {
#pragma vdir nodep loopcnt=25000 /* SX6 */
#pragma _CRI ivdep
    for(i=j; i<25000+j; i++)
    {
#ifdef LONG32
      r0 = (unsigned long)dba[i] & 0x5F5F5F5F;

      r1 = ((r0^0x3E3E3E3E) + 0x01010101) & 0x80808080; /* A, a */
      r2 = (r1>>7) + (r1>>23);
      a += (r2 + (r2>>8)) & 0xF;
      
      r1 = ((r0^0x3C3C3C3C) + 0x01010101) & 0x80808080; /* C, c */
      r2 = (r1>>7) + (r1>>23);
      c += (r2 + (r2>>8)) & 0xF;
      
      r1 = ((r0^0x2B2B2B2B) + 0x01010101) & 0x80808080; /* T, t */
      r2 = (r1>>7) + (r1>>23);
      t += (r2 + (r2>>8)) & 0xF;
      
      r1 = ((r0^0x38383838) + 0x01010101) & 0x80808080; /* G, g */
      r2 = (r1>>7) + (r1>>23);
      g += (r2 + (r2>>8)) & 0xF;
      
      r1 = ((r0^0x31313131) + 0x01010101) & 0x80808080; /* N, n */
      r2 = (r1>>7) + (r1>>23);
      n += (r2 + (r2>>8)) & 0xF;
#endif
#ifdef LONG64
      r0 = (unsigned long)dba[i] & 0x5F5F5F5F5F5F5F5F;

      r1 = ((r0^0x3E3E3E3E3E3E3E3E) + 0x0101010101010101) &
                0x8080808080808080;
      r2 = (r1>>7) + (r1>>39);
      r1 = r2 + (r2>>16);
      a += (r1 + (r1>>8)) & 0xF;
      
      r1 = ((r0^0x3C3C3C3C3C3C3C3C) + 0x0101010101010101) &
                0x8080808080808080;
      r2 = (r1>>7) + (r1>>39);
      r1 = r2 + (r2>>16);
      c += (r1 + (r1>>8)) & 0xF;
    
      r1 = ((r0^0x2B2B2B2B2B2B2B2B) + 0x0101010101010101) &
                0x8080808080808080;
      r2 = (r1>>7) + (r1>>39);
      r1 = r2 + (r2>>16);
      t += (r1 + (r1>>8)) & 0xF;
    
      r1 = ((r0^0x3838383838383838) + 0x0101010101010101) & 
                0x8080808080808080;
      r2 = (r1>>7) + (r1>>39);
      r1 = r2 + (r2>>16);
      g += (r1 + (r1>>8)) & 0xF;
    
      r1 = ((r0^0x3131313131313131) + 0x0101010101010101) &
                0x8080808080808080;
      r2 = (r1>>7) + (r1>>39);
      r1 = r2 + (r2>>16);
      n += (r1 + (r1>>8)) & 0xF;
#endif
    }
    j=i;
  }
  
#pragma vdir nodep loopcnt=25000 /* SX6 */
#pragma _CRI ivdep
  for(i=j; i<length-1; i++)
  {
#ifdef LONG32
    r0 = (unsigned long)dba[i] & 0x5F5F5F5F;

    r1 = ((r0^0x3E3E3E3E) + 0x01010101) & 0x80808080; /* A, a */
    r2 = (r1>>7) + (r1>>23);
    a += (r2 + (r2>>8)) & 0xF;
      
    r1 = ((r0^0x3C3C3C3C) + 0x01010101) & 0x80808080; /* C, c */
    r2 = (r1>>7) + (r1>>23);
    c += (r2 + (r2>>8)) & 0xF;
      
    r1 = ((r0^0x2B2B2B2B) + 0x01010101) & 0x80808080; /* T, t */
    r2 = (r1>>7) + (r1>>23);
    t += (r2 + (r2>>8)) & 0xF;
      
    r1 = ((r0^0x38383838) + 0x01010101) & 0x80808080; /* G, g */
    r2 = (r1>>7) + (r1>>23);
    g += (r2 + (r2>>8)) & 0xF;
      
    r1 = ((r0^0x31313131) + 0x01010101) & 0x80808080; /* N, n */
    r2 = (r1>>7) + (r1>>23);
    n += (r2 + (r2>>8)) & 0xF;
#endif
#ifdef LONG64
    r0 = (unsigned long)dba[i] & 0x5F5F5F5F5F5F5F5F;

    r1 = ((r0^0x3E3E3E3E3E3E3E3E) + 0x0101010101010101) &
              0x8080808080808080;
    r2 = (r1>>7) + (r1>>39);
    r1 = r2 + (r2>>16);
    a += (r1 + (r1>>8)) & 0xF;
      
    r1 = ((r0^0x3C3C3C3C3C3C3C3C) + 0x0101010101010101) &
              0x8080808080808080;
    r2 = (r1>>7) + (r1>>39);
    r1 = r2 + (r2>>16);
    c += (r1 + (r1>>8)) & 0xF;
    
    r1 = ((r0^0x2B2B2B2B2B2B2B2B) + 0x0101010101010101) &
              0x8080808080808080;
    r2 = (r1>>7) + (r1>>39);
    r1 = r2 + (r2>>16);
    t += (r1 + (r1>>8)) & 0xF;
    
    r1 = ((r0^0x3838383838383838) + 0x0101010101010101) & 
              0x8080808080808080;
    r2 = (r1>>7) + (r1>>39);
    r1 = r2 + (r2>>16);
    g += (r1 + (r1>>8)) & 0xF;
    
    r1 = ((r0^0x3131313131313131) + 0x0101010101010101) &
              0x8080808080808080;
    r2 = (r1>>7) + (r1>>39);
    r1 = r2 + (r2>>16);
    n += (r1 + (r1>>8)) & 0xF;
#endif
  } 
  
  res[0] = (long)a;
  res[1] = (long)c;
  res[2] = (long)t;
  res[3] = (long)g;
  res[4] = (long)n;
  
/* now do the last word */

#ifdef B_ENDIAN
#ifdef LONG32
  if(4*i <= dblen)
    switch(dba[i] & 0xFF000000)
    {
      case 0x41000000: res[0]++; break;
      case 0x61000000: res[0]++; break;
      case 0x43000000: res[1]++; break;
      case 0x63000000: res[1]++; break;
      case 0x54000000: res[2]++; break;
      case 0x74000000: res[2]++; break;
      case 0x47000000: res[3]++; break;
      case 0x67000000: res[3]++; break;
      case 0x4E000000: res[4]++; break;
      case 0x6E000000: res[4]++; break;
      default: break;
    }
  if(4*i+1 <= dblen)
    switch(dba[i] & 0xFF0000)
    {
      case 0x410000: res[0]++; break;
      case 0x610000: res[0]++; break;
      case 0x430000: res[1]++; break;
      case 0x630000: res[1]++; break;
      case 0x540000: res[2]++; break;
      case 0x740000: res[2]++; break;
      case 0x470000: res[3]++; break;
      case 0x670000: res[3]++; break;
      case 0x4E0000: res[4]++; break;
      case 0x6E0000: res[4]++; break;
      default: break;
    }
  if(4*i+2 <= dblen)
    switch(dba[i] & 0xFF00)
    {
      case 0x4100: res[0]++; break;
      case 0x6100: res[0]++; break;
      case 0x4300: res[1]++; break;
      case 0x6300: res[1]++; break;
      case 0x5400: res[2]++; break;
      case 0x7400: res[2]++; break;
      case 0x4700: res[3]++; break;
      case 0x6700: res[3]++; break;
      case 0x4E00: res[4]++; break;
      case 0x6E00: res[4]++; break;
      default: break;
    }
  if(4*i+3 <= dblen)
    switch(dba[i] & 0xFF)
    {
      case 0x41: res[0]++; break;
      case 0x61: res[0]++; break;
      case 0x43: res[1]++; break;
      case 0x63: res[1]++; break;
      case 0x54: res[2]++; break;
      case 0x74: res[2]++; break;
      case 0x47: res[3]++; break;
      case 0x67: res[3]++; break;
      case 0x4E: res[4]++; break;
      case 0x6E: res[4]++; break;
      default: break;
    }
#endif /* 32-bit */

#ifdef LONG64
if(8*i <= dblen)
    switch(dba[i] & 0xFF00000000000000)
    {
      case 0x4100000000000000: res[0]++; break;
      case 0x6100000000000000: res[0]++; break;
      case 0x4300000000000000: res[1]++; break;
      case 0x6300000000000000: res[1]++; break;
      case 0x5400000000000000: res[2]++; break;
      case 0x7400000000000000: res[2]++; break;
      case 0x4700000000000000: res[3]++; break;
      case 0x6700000000000000: res[3]++; break;
      case 0x4E00000000000000: res[4]++; break;
      case 0x6E00000000000000: res[4]++; break;
      default: break;
    }
  if(8*i+1 <= dblen)
    switch(dba[i] & 0xFF000000000000)
    {
      case 0x41000000000000: res[0]++; break;
      case 0x61000000000000: res[0]++; break;
      case 0x43000000000000: res[1]++; break;
      case 0x63000000000000: res[1]++; break;
      case 0x54000000000000: res[2]++; break;
      case 0x74000000000000: res[2]++; break;
      case 0x47000000000000: res[3]++; break;
      case 0x67000000000000: res[3]++; break;
      case 0x4E000000000000: res[4]++; break;
      case 0x6E000000000000: res[4]++; break;
      default: break;
    }
  if(8*i+2 <= dblen)
    switch(dba[i] & 0xFF0000000000)
    {
      case 0x410000000000: res[0]++; break;
      case 0x610000000000: res[0]++; break;
      case 0x430000000000: res[1]++; break;
      case 0x630000000000: res[1]++; break;
      case 0x540000000000: res[2]++; break;
      case 0x740000000000: res[2]++; break;
      case 0x470000000000: res[3]++; break;
      case 0x670000000000: res[3]++; break;
      case 0x4E0000000000: res[4]++; break;
      case 0x6E0000000000: res[4]++; break;
      default: break;
    }
  if(8*i+3 <= dblen)
    switch(dba[i] & 0xFF00000000)
    {
      case 0x4100000000: res[0]++; break;
      case 0x6100000000: res[0]++; break;
      case 0x4300000000: res[1]++; break;
      case 0x6300000000: res[1]++; break;
      case 0x5400000000: res[2]++; break;
      case 0x7400000000: res[2]++; break;
      case 0x4700000000: res[3]++; break;
      case 0x6700000000: res[3]++; break;
      case 0x4E00000000: res[4]++; break;
      case 0x6E00000000: res[4]++; break;
      default: break;
    }
    if(8*i+4 <= dblen)
    switch(dba[i] & 0xFF000000)
    {
      case 0x41000000: res[0]++; break;
      case 0x61000000: res[0]++; break;
      case 0x43000000: res[1]++; break;
      case 0x63000000: res[1]++; break;
      case 0x54000000: res[2]++; break;
      case 0x74000000: res[2]++; break;
      case 0x47000000: res[3]++; break;
      case 0x67000000: res[3]++; break;
      case 0x4E000000: res[4]++; break;
      case 0x6E000000: res[4]++; break;
      default: break;
    }
  if(8*i+5 <= dblen)
    switch(dba[i] & 0xFF0000)
    {
      case 0x410000: res[0]++; break;
      case 0x610000: res[0]++; break;
      case 0x430000: res[1]++; break;
      case 0x630000: res[1]++; break;
      case 0x540000: res[2]++; break;
      case 0x740000: res[2]++; break;
      case 0x470000: res[3]++; break;
      case 0x670000: res[3]++; break;
      case 0x4E0000: res[4]++; break;
      case 0x6E0000: res[4]++; break;
      default: break;
    }
  if(8*i+6 <= dblen)
    switch(dba[i] & 0xFF00)
    {
      case 0x4100: res[0]++; break;
      case 0x6100: res[0]++; break;
      case 0x4300: res[1]++; break;
      case 0x6300: res[1]++; break;
      case 0x5400: res[2]++; break;
      case 0x7400: res[2]++; break;
      case 0x4700: res[3]++; break;
      case 0x6700: res[3]++; break;
      case 0x4E00: res[4]++; break;
      case 0x6E00: res[4]++; break;
      default: break;
    }
  if(8*i+7 <= dblen)
    switch(dba[i] & 0xFF)
    {
      case 0x41: res[0]++; break;
      case 0x61: res[0]++; break;
      case 0x43: res[1]++; break;
      case 0x63: res[1]++; break;
      case 0x54: res[2]++; break;
      case 0x74: res[2]++; break;
      case 0x47: res[3]++; break;
      case 0x67: res[3]++; break;
      case 0x4E: res[4]++; break;
      case 0x6E: res[4]++; break;
      default: break;
    }
#endif /* 64-bit */
#endif /* big_endian */

#ifdef L_ENDIAN
#ifdef LONG32
  if(4*i <= dblen)
    switch(dba[i] & 0xFF)
    {
      case 0x41: res[0]++; break;
      case 0x61: res[0]++; break;
      case 0x43: res[1]++; break;
      case 0x63: res[1]++; break;
      case 0x54: res[2]++; break;
      case 0x74: res[2]++; break;
      case 0x47: res[3]++; break;
      case 0x67: res[3]++; break;
      case 0x4E: res[4]++; break;
      case 0x6E: res[4]++; break;
      default: break;
    }
  if(4*i+1 <= dblen)
    switch(dba[i] & 0xFF00)
    {
      case 0x4100: res[0]++; break;
      case 0x6100: res[0]++; break;
      case 0x4300: res[1]++; break;
      case 0x6300: res[1]++; break;
      case 0x5400: res[2]++; break;
      case 0x7400: res[2]++; break;
      case 0x4700: res[3]++; break;
      case 0x6700: res[3]++; break;
      case 0x4E00: res[4]++; break;
      case 0x6E00: res[4]++; break;
      default: break;
    }
  if(4*i+2 <= dblen)
    switch(dba[i] & 0xFF0000)
    {
      case 0x410000: res[0]++; break;
      case 0x610000: res[0]++; break;
      case 0x430000: res[1]++; break;
      case 0x630000: res[1]++; break;
      case 0x540000: res[2]++; break;
      case 0x740000: res[2]++; break;
      case 0x470000: res[3]++; break;
      case 0x670000: res[3]++; break;
      case 0x4E0000: res[4]++; break;
      case 0x6E0000: res[4]++; break;
      default: break;
    }
  if(4*i+3 <= dblen)
    switch(dba[i] & 0xFF000000)
    {
      case 0x41000000: res[0]++; break;
      case 0x61000000: res[0]++; break;
      case 0x43000000: res[1]++; break;
      case 0x63000000: res[1]++; break;
      case 0x54000000: res[2]++; break;
      case 0x74000000: res[2]++; break;
      case 0x47000000: res[3]++; break;
      case 0x67000000: res[3]++; break;
      case 0x4E000000: res[4]++; break;
      case 0x6E000000: res[4]++; break;
      default: break;
    }
  
#endif /* 32-bit */

#ifdef LONG64
  if(8*i <= dblen)
    switch(dba[i] & 0xFF)
    {
      case 0x41: res[0]++; break;
      case 0x61: res[0]++; break;
      case 0x43: res[1]++; break;
      case 0x63: res[1]++; break;
      case 0x54: res[2]++; break;
      case 0x74: res[2]++; break;
      case 0x47: res[3]++; break;
      case 0x67: res[3]++; break;
      case 0x4E: res[4]++; break;
      case 0x6E: res[4]++; break;
      default: break;
    }
  if(8*i+1 <= dblen)
    switch(dba[i] & 0xFF00)
    {
      case 0x4100: res[0]++; break;
      case 0x6100: res[0]++; break;
      case 0x4300: res[1]++; break;
      case 0x6300: res[1]++; break;
      case 0x5400: res[2]++; break;
      case 0x7400: res[2]++; break;
      case 0x4700: res[3]++; break;
      case 0x6700: res[3]++; break;
      case 0x4E00: res[4]++; break;
      case 0x6E00: res[4]++; break;
      default: break;
    }
  if(8*i+2 <= dblen)
    switch(dba[i] & 0xFF0000)
    {
      case 0x410000: res[0]++; break;
      case 0x610000: res[0]++; break;
      case 0x430000: res[1]++; break;
      case 0x630000: res[1]++; break;
      case 0x540000: res[2]++; break;
      case 0x740000: res[2]++; break;
      case 0x470000: res[3]++; break;
      case 0x670000: res[3]++; break;
      case 0x4E0000: res[4]++; break;
      case 0x6E0000: res[4]++; break;
      default: break;
    }
  if(8*i+3 <= dblen)
    switch(dba[i] & 0xFF000000)
    {
      case 0x41000000: res[0]++; break;
      case 0x61000000: res[0]++; break;
      case 0x43000000: res[1]++; break;
      case 0x63000000: res[1]++; break;
      case 0x54000000: res[2]++; break;
      case 0x74000000: res[2]++; break;
      case 0x47000000: res[3]++; break;
      case 0x67000000: res[3]++; break;
      case 0x4E000000: res[4]++; break;
      case 0x6E000000: res[4]++; break;
      default: break;
    }
  if(8*i+4 <= dblen)
    switch(dba[i] & 0xFF00000000)
    {
      case 0x4100000000: res[0]++; break;
      case 0x6100000000: res[0]++; break;
      case 0x4300000000: res[1]++; break;
      case 0x6300000000: res[1]++; break;
      case 0x5400000000: res[2]++; break;
      case 0x7400000000: res[2]++; break;
      case 0x4700000000: res[3]++; break;
      case 0x6700000000: res[3]++; break;
      case 0x4E00000000: res[4]++; break;
      case 0x6E00000000: res[4]++; break;
      default: break;
    }
  
  if(8*i+5 <= dblen)
    switch(dba[i] & 0xFF0000000000)
    {
      case 0x410000000000: res[0]++; break;
      case 0x610000000000: res[0]++; break;
      case 0x430000000000: res[1]++; break;
      case 0x630000000000: res[1]++; break;
      case 0x540000000000: res[2]++; break;
      case 0x740000000000: res[2]++; break;
      case 0x470000000000: res[3]++; break;
      case 0x670000000000: res[3]++; break;
      case 0x4E0000000000: res[4]++; break;
      case 0x6E0000000000: res[4]++; break;
      default: break;
    }
  
  if(8*i+6 <= dblen)
    switch(dba[i] & 0xFF000000000000)
    {
      case 0x41000000000000: res[0]++; break;
      case 0x61000000000000: res[0]++; break;
      case 0x43000000000000: res[1]++; break;
      case 0x63000000000000: res[1]++; break;
      case 0x54000000000000: res[2]++; break;
      case 0x74000000000000: res[2]++; break;
      case 0x47000000000000: res[3]++; break;
      case 0x67000000000000: res[3]++; break;
      case 0x4E000000000000: res[4]++; break;
      case 0x6E000000000000: res[4]++; break;
      default: break;
    }
  
  if(8*i+7 <= dblen)
    switch(dba[i] & 0xFF00000000000000)
    {
      case 0x4100000000000000: res[0]++; break;
      case 0x6100000000000000: res[0]++; break;
      case 0x4300000000000000: res[1]++; break;
      case 0x6300000000000000: res[1]++; break;
      case 0x5400000000000000: res[2]++; break;
      case 0x7400000000000000: res[2]++; break;
      case 0x4700000000000000: res[3]++; break;
      case 0x6700000000000000: res[3]++; break;
      case 0x4E00000000000000: res[4]++; break;
      case 0x6E00000000000000: res[4]++; break;
      default: break;
    }
#endif /* 64-bit */
#endif /* little_endian */
}

/*
cb_count_ascii(3B)                                        Last changed: 09-17-02

NAME
        cb_count_ascii - counts A, C, T, G, and N characters in a string

SYNOPSIS

        C/C++:

        #include <cbl.h>
        void cb_count_ascii( long *dba, long dblen, long res[] );

        Fortran:

        call cb_count_ascii( dba, dblen, res )


IMPLEMENTATION

        Cray SV1 series UNICOS systems

DESCRIPTION

        cb_count_ascii counts the number of occurrences of the letters 
        A, C, T, G, and N in a string containing ASCII letters. Lower case 
        letters are treated as equivalent to the corresponding upper case
        letters, and are included in the counts.  The number of letters
        not included in the set of five counted can be obtained by subtracting
        the sum of the five results from dblen.

        dba     (input) sequence of ASCII letters packed 8/word with
                possible null characters in the end of the final word.
                In Fortran, db should be an INTEGER(8) array.

        dblen   (input) number of nucleotides packed into dba.

        res     (output) array of at least 5 elements containing the counts
                of the common nucleotide codes as follows:

                Fortran    C  

                res(1)   res[0] = number of A and a characters in dba.
                res(2)   res[1] = number of C and c characters in dba.
                res(3)   res[2] = number of T and t characters in dba.
                res(4)   res[3] = number of G and g characters in dba.
                res(5)   res[4] = number of N and n characters in dba.

                Memory for res must be allocated before calling
                cb_count_ascii. In Fortran, res should be declared as
                an INTEGER(8) array.

NOTES       

        cb_count_ascii provides a quick way to profile a set of nucleotide
        data to determine how many ambiguous codes are present, and also to
        check for a relative excess or deficiency of C/G letters.

        cb_count_ascii is single-threaded (i.e. not tasked) and may be called 
        from within a parallel region.

        cb_count_ascii replaces the contents of the bmm register.

SEE ALSO

        INTRO_LIBCBL(3B)

*/
