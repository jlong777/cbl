/* cb_copy_bits.c                                http://cbl.sourceforge.net
 *
 * Copy contiguous sequence of memory bits.
 * see original man page at the bottom. This code is the generic version
 * with provision for 32- and 64-bit long ints.
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
 * $Id: cb_copy_bits.c,v 1.4 2003/09/25 00:53:27 jlong777 Exp $
 */

#include "cb_macro.h"

void cb_copy_bits(long *dest, long doffset, long *src, long soffset, long nbits)
{
  long i, length, switcharg, wordsz;
  unsigned long mask_l, mask_r;
  register long r0, r1, r2, r3;
  
#ifdef B_ENDIAN

#ifdef LONG32
  
  wordsz = 32;
  
  r0 = doffset%wordsz;
  
  switcharg = r0;
  switch(switcharg) /* get left mask */
  {
    case  1: mask_l = 0x80000000; break;
    case  2: mask_l = 0xC0000000; break;
    case  3: mask_l = 0xE0000000; break;
    case  4: mask_l = 0xF0000000; break;
    case  5: mask_l = 0xF8000000; break;
    case  6: mask_l = 0xFC000000; break;
    case  7: mask_l = 0xFE000000; break;
    case  8: mask_l = 0xFF000000; break;
    case  9: mask_l = 0xFF800000; break;
    case 10: mask_l = 0xFFC00000; break;
    case 11: mask_l = 0xFFE00000; break;
    case 12: mask_l = 0xFFF00000; break;
    case 13: mask_l = 0xFFF80000; break;
    case 14: mask_l = 0xFFFC0000; break;
    case 15: mask_l = 0xFFFE0000; break;
    case 16: mask_l = 0xFFFF0000; break;
    case 17: mask_l = 0xFFFF8000; break;
    case 18: mask_l = 0xFFFFC000; break;
    case 19: mask_l = 0xFFFFE000; break;
    case 20: mask_l = 0xFFFFF000; break;
    case 21: mask_l = 0xFFFFF800; break;
    case 22: mask_l = 0xFFFFFC00; break;
    case 23: mask_l = 0xFFFFFE00; break;
    case 24: mask_l = 0xFFFFFF00; break;
    case 25: mask_l = 0xFFFFFF80; break;
    case 26: mask_l = 0xFFFFFFC0; break;
    case 27: mask_l = 0xFFFFFFE0; break;
    case 28: mask_l = 0xFFFFFFF0; break;
    case 29: mask_l = 0xFFFFFFF8; break;
    case 30: mask_l = 0xFFFFFFFC; break;
    case 31: mask_l = 0xFFFFFFFE; break;
    default: mask_l = 0x00000000; break;
  }
  
  switcharg = (nbits+r0)%32;
  switch(switcharg) /* get right mask */
  {
    case  1: mask_r = 0x7FFFFFFF; break;
    case  2: mask_r = 0x3FFFFFFF; break;
    case  3: mask_r = 0x1FFFFFFF; break;
    case  4: mask_r = 0x0FFFFFFF; break;
    case  5: mask_r = 0x07FFFFFF; break;
    case  6: mask_r = 0x03FFFFFF; break;
    case  7: mask_r = 0x01FFFFFF; break;
    case  8: mask_r = 0x00FFFFFF; break;
    case  9: mask_r = 0x007FFFFF; break;
    case 10: mask_r = 0x003FFFFF; break;
    case 11: mask_r = 0x001FFFFF; break;
    case 12: mask_r = 0x000FFFFF; break;
    case 13: mask_r = 0x0007FFFF; break;
    case 14: mask_r = 0x0003FFFF; break;
    case 15: mask_r = 0x0001FFFF; break;
    case 16: mask_r = 0x0000FFFF; break;
    case 17: mask_r = 0x00007FFF; break;
    case 18: mask_r = 0x00003FFF; break;
    case 19: mask_r = 0x00001FFF; break;
    case 20: mask_r = 0x00000FFF; break;
    case 21: mask_r = 0x000007FF; break;
    case 22: mask_r = 0x000003FF; break;
    case 23: mask_r = 0x000001FF; break;
    case 24: mask_r = 0x000000FF; break;
    case 25: mask_r = 0x0000007F; break;
    case 26: mask_r = 0x0000003F; break;
    case 27: mask_r = 0x0000001F; break;
    case 28: mask_r = 0x0000000F; break;
    case 29: mask_r = 0x00000007; break;
    case 30: mask_r = 0x00000003; break;
    case 31: mask_r = 0x00000001; break;
    default: mask_r = 0x00000000; break;
  }
  
#endif /* long32 */
  
#ifdef LONG64
  
  wordsz = 64;
  
  r0 = doffset%wordsz;
  
  switcharg = r0;
  switch(switcharg) /* get left mask */
  {
    case  1: mask_l = 0x8000000000000000; break;
    case  2: mask_l = 0xC000000000000000; break;
    case  3: mask_l = 0xE000000000000000; break;
    case  4: mask_l = 0xF000000000000000; break;
    case  5: mask_l = 0xF800000000000000; break;
    case  6: mask_l = 0xFC00000000000000; break;
    case  7: mask_l = 0xFE00000000000000; break;
    case  8: mask_l = 0xFF00000000000000; break;
    case  9: mask_l = 0xFF80000000000000; break;
    case 10: mask_l = 0xFFC0000000000000; break;
    case 11: mask_l = 0xFFE0000000000000; break;
    case 12: mask_l = 0xFFF0000000000000; break;
    case 13: mask_l = 0xFFF8000000000000; break;
    case 14: mask_l = 0xFFFC000000000000; break;
    case 15: mask_l = 0xFFFE000000000000; break;
    case 16: mask_l = 0xFFFF000000000000; break;
    case 17: mask_l = 0xFFFF800000000000; break;
    case 18: mask_l = 0xFFFFC00000000000; break;
    case 19: mask_l = 0xFFFFE00000000000; break;
    case 20: mask_l = 0xFFFFF00000000000; break;
    case 21: mask_l = 0xFFFFF80000000000; break;
    case 22: mask_l = 0xFFFFFC0000000000; break;
    case 23: mask_l = 0xFFFFFE0000000000; break;
    case 24: mask_l = 0xFFFFFF0000000000; break;
    case 25: mask_l = 0xFFFFFF8000000000; break;
    case 26: mask_l = 0xFFFFFFC000000000; break;
    case 27: mask_l = 0xFFFFFFE000000000; break;
    case 28: mask_l = 0xFFFFFFF000000000; break;
    case 29: mask_l = 0xFFFFFFF800000000; break;
    case 30: mask_l = 0xFFFFFFFC00000000; break;
    case 31: mask_l = 0xFFFFFFFE00000000; break;
    case 32: mask_l = 0xFFFFFFFF00000000; break;
    case 33: mask_l = 0xFFFFFFFF80000000; break;
    case 34: mask_l = 0xFFFFFFFFC0000000; break;
    case 35: mask_l = 0xFFFFFFFFE0000000; break;
    case 36: mask_l = 0xFFFFFFFFF0000000; break;
    case 37: mask_l = 0xFFFFFFFFF8000000; break;
    case 38: mask_l = 0xFFFFFFFFFC000000; break;
    case 39: mask_l = 0xFFFFFFFFFE000000; break;
    case 40: mask_l = 0xFFFFFFFFFF000000; break;
    case 41: mask_l = 0xFFFFFFFFFF800000; break;
    case 42: mask_l = 0xFFFFFFFFFFC00000; break;
    case 43: mask_l = 0xFFFFFFFFFFE00000; break;
    case 44: mask_l = 0xFFFFFFFFFFF00000; break;
    case 45: mask_l = 0xFFFFFFFFFFF80000; break;
    case 46: mask_l = 0xFFFFFFFFFFFC0000; break;
    case 47: mask_l = 0xFFFFFFFFFFFE0000; break;
    case 48: mask_l = 0xFFFFFFFFFFFF0000; break;
    case 49: mask_l = 0xFFFFFFFFFFFF8000; break;
    case 50: mask_l = 0xFFFFFFFFFFFFC000; break;
    case 51: mask_l = 0xFFFFFFFFFFFFE000; break;
    case 52: mask_l = 0xFFFFFFFFFFFFF000; break;
    case 53: mask_l = 0xFFFFFFFFFFFFF800; break;
    case 54: mask_l = 0xFFFFFFFFFFFFFC00; break;
    case 55: mask_l = 0xFFFFFFFFFFFFFE00; break;
    case 56: mask_l = 0xFFFFFFFFFFFFFF00; break;
    case 57: mask_l = 0xFFFFFFFFFFFFFF80; break;
    case 58: mask_l = 0xFFFFFFFFFFFFFFC0; break;
    case 59: mask_l = 0xFFFFFFFFFFFFFFE0; break;
    case 60: mask_l = 0xFFFFFFFFFFFFFFF0; break;
    case 61: mask_l = 0xFFFFFFFFFFFFFFF8; break;
    case 62: mask_l = 0xFFFFFFFFFFFFFFFC; break;
    case 63: mask_l = 0xFFFFFFFFFFFFFFFE; break;
    default: mask_l = 0x0000000000000000; break;
  }
  
  switcharg = (nbits+r0)%64;
  switch(switcharg) /* get right mask */
  {
    case  1: mask_r = 0x7FFFFFFFFFFFFFFF; break;
    case  2: mask_r = 0x3FFFFFFFFFFFFFFF; break;
    case  3: mask_r = 0x1FFFFFFFFFFFFFFF; break;
    case  4: mask_r = 0x0FFFFFFFFFFFFFFF; break;
    case  5: mask_r = 0x07FFFFFFFFFFFFFF; break;
    case  6: mask_r = 0x03FFFFFFFFFFFFFF; break;
    case  7: mask_r = 0x01FFFFFFFFFFFFFF; break;
    case  8: mask_r = 0x00FFFFFFFFFFFFFF; break;
    case  9: mask_r = 0x007FFFFFFFFFFFFF; break;
    case 10: mask_r = 0x003FFFFFFFFFFFFF; break;
    case 11: mask_r = 0x001FFFFFFFFFFFFF; break;
    case 12: mask_r = 0x000FFFFFFFFFFFFF; break;
    case 13: mask_r = 0x0007FFFFFFFFFFFF; break;
    case 14: mask_r = 0x0003FFFFFFFFFFFF; break;
    case 15: mask_r = 0x0001FFFFFFFFFFFF; break;
    case 16: mask_r = 0x0000FFFFFFFFFFFF; break;
    case 17: mask_r = 0x00007FFFFFFFFFFF; break;
    case 18: mask_r = 0x00003FFFFFFFFFFF; break;
    case 19: mask_r = 0x00001FFFFFFFFFFF; break;
    case 20: mask_r = 0x00000FFFFFFFFFFF; break;
    case 21: mask_r = 0x000007FFFFFFFFFF; break;
    case 22: mask_r = 0x000003FFFFFFFFFF; break;
    case 23: mask_r = 0x000001FFFFFFFFFF; break;
    case 24: mask_r = 0x000000FFFFFFFFFF; break;
    case 25: mask_r = 0x0000007FFFFFFFFF; break;
    case 26: mask_r = 0x0000003FFFFFFFFF; break;
    case 27: mask_r = 0x0000001FFFFFFFFF; break;
    case 28: mask_r = 0x0000000FFFFFFFFF; break;
    case 29: mask_r = 0x00000007FFFFFFFF; break;
    case 30: mask_r = 0x00000003FFFFFFFF; break;
    case 31: mask_r = 0x00000001FFFFFFFF; break;
    case 32: mask_r = 0x00000000FFFFFFFF; break;
    case 33: mask_r = 0x000000007FFFFFFF; break;
    case 34: mask_r = 0x000000003FFFFFFF; break;
    case 35: mask_r = 0x000000001FFFFFFF; break;
    case 36: mask_r = 0x000000000FFFFFFF; break;
    case 37: mask_r = 0x0000000007FFFFFF; break;
    case 38: mask_r = 0x0000000003FFFFFF; break;
    case 39: mask_r = 0x0000000001FFFFFF; break;
    case 40: mask_r = 0x0000000000FFFFFF; break;
    case 41: mask_r = 0x00000000007FFFFF; break;
    case 42: mask_r = 0x00000000003FFFFF; break;
    case 43: mask_r = 0x00000000001FFFFF; break;
    case 44: mask_r = 0x00000000000FFFFF; break;
    case 45: mask_r = 0x000000000007FFFF; break;
    case 46: mask_r = 0x000000000003FFFF; break;
    case 47: mask_r = 0x000000000001FFFF; break;
    case 48: mask_r = 0x000000000000FFFF; break;
    case 49: mask_r = 0x0000000000007FFF; break;
    case 50: mask_r = 0x0000000000003FFF; break;
    case 51: mask_r = 0x0000000000001FFF; break;
    case 52: mask_r = 0x0000000000000FFF; break;
    case 53: mask_r = 0x00000000000007FF; break;
    case 54: mask_r = 0x00000000000003FF; break;
    case 55: mask_r = 0x00000000000001FF; break;
    case 56: mask_r = 0x00000000000000FF; break;
    case 57: mask_r = 0x000000000000007F; break;
    case 58: mask_r = 0x000000000000003F; break;
    case 59: mask_r = 0x000000000000001F; break;
    case 60: mask_r = 0x000000000000000F; break;
    case 61: mask_r = 0x0000000000000007; break;
    case 62: mask_r = 0x0000000000000003; break;
    case 63: mask_r = 0x0000000000000001; break;    
    default: mask_r = 0x0000000000000000; break;
  }
  
#endif /* long64 */

  r1 = doffset/wordsz;
  r2 = soffset%wordsz;
  r3 = soffset/wordsz;
  
  length = (nbits+r0+(wordsz-1))/wordsz; /* max number of dest words affected */

  if(r0 >= r2)
  { 
    if(nbits > (wordsz-r0)) /* changing at least all of the first word */
    {
      if(r0 == r2)    /* avoid whole word shifts */
      {
        dest[r1] = (dest[r1] & mask_l) | ((unsigned long)(src[r3] << r2) >> r0);
#pragma vdir nodep
#pragma _CRI ivdep
	for(i=1; i<length-1; i++)
          dest[i+r1] = src[i+r3];
      
        dest[i+r1] = (dest[i+r1] & mask_r) | (src[i+r3] & ~mask_r);
      }
      else /* no shifts by whole word */
      {
        dest[r1] = (dest[r1] & mask_l) | ((unsigned long)(src[r3] << r2) >> r0);
#pragma vdir nodep
#pragma _CRI ivdep
        for(i=1; i<length-1; i++)
          dest[i+r1] = (src[i+r3-1] << (wordsz-r0+r2)) | 
	               ((unsigned long)src[i+r3] >> (r0-r2));
    
        dest[i+r1] = (dest[i+r1] & (mask_r))    |
                     (((src[i+r3-1] << (wordsz-r0+r2)) | 
	             ((unsigned long)src[i+r3] >> (r0-r2))) & ~mask_r);
      }
    }
    else /* just changing something in first word */
    {
      if(nbits == (wordsz-r0)) /* avoid shift by whole word */
        dest[r1] = (dest[r1] & mask_l) | ((unsigned long)(src[r3] << r2) >> r0);
      else
        dest[r1] = (dest[r1] & mask_l)                                     |
	           (((unsigned long)((unsigned long)(src[r3] << r2) >> r0) >>
	           (wordsz-r0-nbits)) << (wordsz-r0-nbits))                |
		   ((unsigned long)(dest[r1] << (r0+nbits)) >> (r0+nbits));
    }
  }
  else /* r0 < r2 */
  {
    if(nbits >= (wordsz-r0)) /* changing at least all of the first word */
    {
      dest[r1] = (dest[r1] & mask_l) | ((unsigned long)(src[r3] << r2) >> r0) |
		 (((unsigned long)src[r3+1]) >> (wordsz-r2+r0));

      if(nbits != (wordsz-r0))
      {
#pragma vdir nodep
#pragma _CRI ivdep
        for(i=1; i<length-1; i++)
          dest[i+r1] = (src[i+r3] << (r2-r0)) | 
	               ((unsigned long)src[i+r3+1] >> (wordsz-r2+r0));

        dest[i+r1] = (dest[i+r1] & mask_r)    |
                     (((src[i+r3] << (r2-r0)) |
		     ((unsigned long)src[i+r3+1] >> (wordsz-r2+r0))) & ~mask_r);
      }
    }
    else /* just changing part of the first word */
    {
      dest[r1] = (dest[r1] & mask_l)                                         |
                 (((unsigned long)(((unsigned long)(src[r3] << r2) >> r0)    |
		 ((unsigned long)src[r3+1] >> (wordsz-r2+r0))) >> (wordsz-r0-nbits)) <<
		 (wordsz-r0-nbits))                                          |
		 ((unsigned long)(dest[r1] << (r0+nbits)) >> (r0+nbits));
    }
  }
  
#endif /* big_endian */


#ifdef L_ENDIAN

#ifdef LONG32
  
  wordsz = 32;
  
  r0 = doffset%wordsz;
  
  switcharg = r0;
  switch(switcharg) /* get right mask */
  {
    case  1: mask_r = 0x00000001; break;
    case  2: mask_r = 0x00000003; break;
    case  3: mask_r = 0x00000007; break;
    case  4: mask_r = 0x0000000F; break;
    case  5: mask_r = 0x0000001F; break;
    case  6: mask_r = 0x0000003F; break;
    case  7: mask_r = 0x0000007F; break;
    case  8: mask_r = 0x000000FF; break;
    case  9: mask_r = 0x000001FF; break;
    case 10: mask_r = 0x000003FF; break;
    case 11: mask_r = 0x000007FF; break;
    case 12: mask_r = 0x00000FFF; break;
    case 13: mask_r = 0x00001FFF; break;
    case 14: mask_r = 0x00003FFF; break;
    case 15: mask_r = 0x00007FFF; break;
    case 16: mask_r = 0x0000FFFF; break;
    case 17: mask_r = 0x0001FFFF; break;
    case 18: mask_r = 0x0003FFFF; break;
    case 19: mask_r = 0x0007FFFF; break;
    case 20: mask_r = 0x000FFFFF; break;
    case 21: mask_r = 0x001FFFFF; break;
    case 22: mask_r = 0x003FFFFF; break;
    case 23: mask_r = 0x007FFFFF; break;
    case 24: mask_r = 0x00FFFFFF; break;
    case 25: mask_r = 0x01FFFFFF; break;
    case 26: mask_r = 0x03FFFFFF; break;
    case 27: mask_r = 0x07FFFFFF; break;
    case 28: mask_r = 0x0FFFFFFF; break;
    case 29: mask_r = 0x1FFFFFFF; break;
    case 30: mask_r = 0x3FFFFFFF; break;
    case 31: mask_r = 0x7FFFFFFF; break;
    default: mask_r = 0x00000000; break;
  }
  
  switcharg = (nbits+r0)%32;
  switch(switcharg) /* get left mask */
  {
    case  1: mask_l = 0xFFFFFFFE; break;
    case  2: mask_l = 0xFFFFFFFC; break;
    case  3: mask_l = 0xFFFFFFF8; break;
    case  4: mask_l = 0xFFFFFFF0; break;
    case  5: mask_l = 0xFFFFFFE0; break;
    case  6: mask_l = 0xFFFFFFC0; break;
    case  7: mask_l = 0xFFFFFF80; break;
    case  8: mask_l = 0xFFFFFF00; break;
    case  9: mask_l = 0xFFFFFE00; break;
    case 10: mask_l = 0xFFFFFC00; break;
    case 11: mask_l = 0xFFFFF800; break;
    case 12: mask_l = 0xFFFFF000; break;
    case 13: mask_l = 0xFFFFE000; break;
    case 14: mask_l = 0xFFFFC000; break;
    case 15: mask_l = 0xFFFF8000; break;
    case 16: mask_l = 0xFFFF0000; break;
    case 17: mask_l = 0xFFFE0000; break;
    case 18: mask_l = 0xFFFC0000; break;
    case 19: mask_l = 0xFFF80000; break;
    case 20: mask_l = 0xFFF00000; break;
    case 21: mask_l = 0xFFE00000; break;
    case 22: mask_l = 0xFFC00000; break;
    case 23: mask_l = 0xFF800000; break;
    case 24: mask_l = 0xFF000000; break;
    case 25: mask_l = 0xFE000000; break;
    case 26: mask_l = 0xFC000000; break;
    case 27: mask_l = 0xF8000000; break;
    case 28: mask_l = 0xF0000000; break;
    case 29: mask_l = 0xE0000000; break;
    case 30: mask_l = 0xC0000000; break;
    case 31: mask_l = 0x80000000; break;
    default: mask_l = 0x00000000; break;
  }
  
#endif /* long32 */
  
#ifdef LONG64
  
  wordsz = 64;
  
  r0 = doffset%wordsz;
  
  switcharg = r0;
  switch(switcharg) /* get right mask */
  {
    case 63: mask_r = 0x7FFFFFFFFFFFFFFF; break;
    case 62: mask_r = 0x3FFFFFFFFFFFFFFF; break;
    case 61: mask_r = 0x1FFFFFFFFFFFFFFF; break;
    case 60: mask_r = 0x0FFFFFFFFFFFFFFF; break;
    case 59: mask_r = 0x07FFFFFFFFFFFFFF; break;
    case 58: mask_r = 0x03FFFFFFFFFFFFFF; break;
    case 57: mask_r = 0x01FFFFFFFFFFFFFF; break;
    case 56: mask_r = 0x00FFFFFFFFFFFFFF; break;
    case 55: mask_r = 0x007FFFFFFFFFFFFF; break;
    case 54: mask_r = 0x003FFFFFFFFFFFFF; break;
    case 53: mask_r = 0x001FFFFFFFFFFFFF; break;
    case 52: mask_r = 0x000FFFFFFFFFFFFF; break;
    case 51: mask_r = 0x0007FFFFFFFFFFFF; break;
    case 50: mask_r = 0x0003FFFFFFFFFFFF; break;
    case 49: mask_r = 0x0001FFFFFFFFFFFF; break;
    case 48: mask_r = 0x0000FFFFFFFFFFFF; break;
    case 47: mask_r = 0x00007FFFFFFFFFFF; break;
    case 46: mask_r = 0x00003FFFFFFFFFFF; break;
    case 45: mask_r = 0x00001FFFFFFFFFFF; break;
    case 44: mask_r = 0x00000FFFFFFFFFFF; break;
    case 43: mask_r = 0x000007FFFFFFFFFF; break;
    case 42: mask_r = 0x000003FFFFFFFFFF; break;
    case 41: mask_r = 0x000001FFFFFFFFFF; break;
    case 40: mask_r = 0x000000FFFFFFFFFF; break;
    case 39: mask_r = 0x0000007FFFFFFFFF; break;
    case 38: mask_r = 0x0000003FFFFFFFFF; break;
    case 37: mask_r = 0x0000001FFFFFFFFF; break;
    case 36: mask_r = 0x0000000FFFFFFFFF; break;
    case 35: mask_r = 0x00000007FFFFFFFF; break;
    case 34: mask_r = 0x00000003FFFFFFFF; break;
    case 33: mask_r = 0x00000001FFFFFFFF; break;
    case 32: mask_r = 0x00000000FFFFFFFF; break;
    case 31: mask_r = 0x000000007FFFFFFF; break;
    case 30: mask_r = 0x000000003FFFFFFF; break;
    case 29: mask_r = 0x000000001FFFFFFF; break;
    case 28: mask_r = 0x000000000FFFFFFF; break;
    case 27: mask_r = 0x0000000007FFFFFF; break;
    case 26: mask_r = 0x0000000003FFFFFF; break;
    case 25: mask_r = 0x0000000001FFFFFF; break;
    case 24: mask_r = 0x0000000000FFFFFF; break;
    case 23: mask_r = 0x00000000007FFFFF; break;
    case 22: mask_r = 0x00000000003FFFFF; break;
    case 21: mask_r = 0x00000000001FFFFF; break;
    case 20: mask_r = 0x00000000000FFFFF; break;
    case 19: mask_r = 0x000000000007FFFF; break;
    case 18: mask_r = 0x000000000003FFFF; break;
    case 17: mask_r = 0x000000000001FFFF; break;
    case 16: mask_r = 0x000000000000FFFF; break;
    case 15: mask_r = 0x0000000000007FFF; break;
    case 14: mask_r = 0x0000000000003FFF; break;
    case 13: mask_r = 0x0000000000001FFF; break;
    case 12: mask_r = 0x0000000000000FFF; break;
    case 11: mask_r = 0x00000000000007FF; break;
    case 10: mask_r = 0x00000000000003FF; break;
    case  9: mask_r = 0x00000000000001FF; break;
    case  8: mask_r = 0x00000000000000FF; break;
    case  7: mask_r = 0x000000000000007F; break;
    case  6: mask_r = 0x000000000000003F; break;
    case  5: mask_r = 0x000000000000001F; break;
    case  4: mask_r = 0x000000000000000F; break;
    case  3: mask_r = 0x0000000000000007; break;
    case  2: mask_r = 0x0000000000000003; break;
    case  1: mask_r = 0x0000000000000001; break;    
    default: mask_r = 0x0000000000000000; break;
  }
  
  switcharg = (nbits+r0)%64;
  switch(switcharg) /* get left mask */
  {
    case 63: mask_l = 0x8000000000000000; break;
    case 62: mask_l = 0xC000000000000000; break;
    case 61: mask_l = 0xE000000000000000; break;
    case 60: mask_l = 0xF000000000000000; break;
    case 59: mask_l = 0xF800000000000000; break;
    case 58: mask_l = 0xFC00000000000000; break;
    case 57: mask_l = 0xFE00000000000000; break;
    case 56: mask_l = 0xFF00000000000000; break;
    case 55: mask_l = 0xFF80000000000000; break;
    case 54: mask_l = 0xFFC0000000000000; break;
    case 53: mask_l = 0xFFE0000000000000; break;
    case 52: mask_l = 0xFFF0000000000000; break;
    case 51: mask_l = 0xFFF8000000000000; break;
    case 50: mask_l = 0xFFFC000000000000; break;
    case 49: mask_l = 0xFFFE000000000000; break;
    case 48: mask_l = 0xFFFF000000000000; break;
    case 47: mask_l = 0xFFFF800000000000; break;
    case 46: mask_l = 0xFFFFC00000000000; break;
    case 45: mask_l = 0xFFFFE00000000000; break;
    case 44: mask_l = 0xFFFFF00000000000; break;
    case 43: mask_l = 0xFFFFF80000000000; break;
    case 42: mask_l = 0xFFFFFC0000000000; break;
    case 41: mask_l = 0xFFFFFE0000000000; break;
    case 40: mask_l = 0xFFFFFF0000000000; break;
    case 39: mask_l = 0xFFFFFF8000000000; break;
    case 38: mask_l = 0xFFFFFFC000000000; break;
    case 37: mask_l = 0xFFFFFFE000000000; break;
    case 36: mask_l = 0xFFFFFFF000000000; break;
    case 35: mask_l = 0xFFFFFFF800000000; break;
    case 34: mask_l = 0xFFFFFFFC00000000; break;
    case 33: mask_l = 0xFFFFFFFE00000000; break;
    case 32: mask_l = 0xFFFFFFFF00000000; break;
    case 31: mask_l = 0xFFFFFFFF80000000; break;
    case 30: mask_l = 0xFFFFFFFFC0000000; break;
    case 29: mask_l = 0xFFFFFFFFE0000000; break;
    case 28: mask_l = 0xFFFFFFFFF0000000; break;
    case 27: mask_l = 0xFFFFFFFFF8000000; break;
    case 26: mask_l = 0xFFFFFFFFFC000000; break;
    case 25: mask_l = 0xFFFFFFFFFE000000; break;
    case 24: mask_l = 0xFFFFFFFFFF000000; break;
    case 23: mask_l = 0xFFFFFFFFFF800000; break;
    case 22: mask_l = 0xFFFFFFFFFFC00000; break;
    case 21: mask_l = 0xFFFFFFFFFFE00000; break;
    case 20: mask_l = 0xFFFFFFFFFFF00000; break;
    case 19: mask_l = 0xFFFFFFFFFFF80000; break;
    case 18: mask_l = 0xFFFFFFFFFFFC0000; break;
    case 17: mask_l = 0xFFFFFFFFFFFE0000; break;
    case 16: mask_l = 0xFFFFFFFFFFFF0000; break;
    case 15: mask_l = 0xFFFFFFFFFFFF8000; break;
    case 14: mask_l = 0xFFFFFFFFFFFFC000; break;
    case 13: mask_l = 0xFFFFFFFFFFFFE000; break;
    case 12: mask_l = 0xFFFFFFFFFFFFF000; break;
    case 11: mask_l = 0xFFFFFFFFFFFFF800; break;
    case 10: mask_l = 0xFFFFFFFFFFFFFC00; break;
    case  9: mask_l = 0xFFFFFFFFFFFFFE00; break;
    case  8: mask_l = 0xFFFFFFFFFFFFFF00; break;
    case  7: mask_l = 0xFFFFFFFFFFFFFF80; break;
    case  6: mask_l = 0xFFFFFFFFFFFFFFC0; break;
    case  5: mask_l = 0xFFFFFFFFFFFFFFE0; break;
    case  4: mask_l = 0xFFFFFFFFFFFFFFF0; break;
    case  3: mask_l = 0xFFFFFFFFFFFFFFF8; break;
    case  2: mask_l = 0xFFFFFFFFFFFFFFFC; break;
    case  1: mask_l = 0xFFFFFFFFFFFFFFFE; break;
    default: mask_l = 0x0000000000000000; break;
  }
  
#endif /* long64 */

  r1 = doffset/wordsz;
  r2 = soffset%wordsz;
  r3 = soffset/wordsz;
  
  length = (nbits+r0+(wordsz-1))/wordsz; /* max number of dest words affected */

  if(r0 >= r2)
  { 
    if(nbits > (wordsz-r0)) /* changing at least all of the first word */
    {
      if(r0 == r2)    /* avoid whole word shifts */
      {
        dest[r1] = (dest[r1] & mask_r) | (((unsigned long)src[r3] >> r2) << r0);
	
	for(i=1; i<length-1; i++)
          dest[i+r1] = src[i+r3];
      
        dest[i+r1] = (dest[i+r1] & mask_l) | (src[i+r3] & ~mask_l);
      }
      else /* no shifts by whole word */
      {
        dest[r1] = (dest[r1] & mask_r) | (((unsigned long)src[r3] >> r2) << r0);
	
        for(i=1; i<length-1; i++)
          dest[i+r1] = ((unsigned long)src[i+r3-1] >> (wordsz-r0+r2)) | 
	               (src[i+r3] << (r0-r2));
    
        dest[i+r1] = (dest[i+r1] & (mask_l))    |
                     ((((unsigned long)src[i+r3-1] >> (wordsz-r0+r2)) | 
	             (src[i+r3] << (r0-r2))) & ~mask_l);
      }
    }
    else /* just changing part of the first word */
    {
      if(nbits == (wordsz-r0)) /* avoid shift by whole word */
        dest[r1] = (dest[r1] & mask_r) | (((unsigned long)src[r3] >> r2) << r0);
      else
        dest[r1] = (dest[r1] & mask_r)                                     |
	           ((unsigned long)((((unsigned long)src[r3] >> r2) << r0) <<
	           (wordsz-r0-nbits)) >> (wordsz-r0-nbits))                |
		   (((unsigned long)dest[r1] >> (r0+nbits)) << (r0+nbits));
    }
  }
  else /* r0 < r2 */
  {
    if(nbits >= (wordsz-r0)) /* changing at least all of the first word */
    {
      dest[r1] = (dest[r1] & mask_r) | (((unsigned long)src[r3] >> r2) << r0) |
		 (src[r3+1] << (wordsz-r2+r0));

      if(nbits != (wordsz-r0))
      {
        for(i=1; i<length-1; i++)
          dest[i+r1] = ((unsigned long)src[i+r3] >> (r2-r0)) | 
	               (src[i+r3+1] << (wordsz-r2+r0));

        dest[i+r1] = (dest[i+r1] & mask_l)    |
                     ((((unsigned long)src[i+r3] >> (r2-r0)) |
		     (src[i+r3+1] << (wordsz-r2+r0))) & ~mask_l);
      }
    }
    else /* just changing part of the first word */
    {
      dest[r1] = (dest[r1] & mask_r)                                         |
                 ((unsigned long)(((((unsigned long)src[r3] >> r2) << r0)    |
		 (src[r3+1] << (wordsz-r2+r0))) << (wordsz-r0-nbits)) >>
		 (wordsz-r0-nbits))                                          |
		 (((unsigned long)dest[r1] >> (r0+nbits)) << (r0+nbits));
    }
  }

#endif /* little_endian */
}
 
/*
cb_copy_bits(3B)                                           Last changed: 09-17-02

NAME
        cb_copy_bits - Copy contiguous sequence of memory bits

SYNOPSIS

        C/C++:

        #include <cbl.h>
        void cb_copy_bits( long *dest, long doffset, long *src, long soffset, long nbits );

        Fortran:

        use cb_bits
        call cb_copy_bits( dest, doffset, src, soffset, nbits)


IMPLEMENTATION

        Cray SV1 series UNICOS systems

DESCRIPTION

        cb_copy_bits copies a sequence of bits from one region of memory to 
        another region of memory.  The source and destination regions should 
        not overlap. The starting location for each region may be on an 
        arbitrary bit within a memory word. 

        dest    (output) base address of the destination memory array.
                In Fortran, dest should be an INTEGER(8) array.

        doffset (input) offset in bits from the beginning of dest where
                the destination of the copy actually begins. doffset
                must be zero or positive and may be larger than 64.
                In Fortran, doffset should be an INTEGER(8) variable,
                constant, or expression.

        src     (input) base address of the source memory array.
                In Fortran, src should be an INTEGER(8) array.

        soffset (input) offset in bits from the beginnig of src where
                the source bits of the copy actually begin. soffset 
                must be zero or positive and may be larger than 64.
                In Fortran, soffset should be an INTEGER(8) variable,
                constant, or expression.

        nbits   (input) number of bits to be copied. In Fortran, nbits 
                should be an INTEGER(8) variable, constant, or expression.


NOTES

        cb_copy_bits is optimized for copying large blocks of unaligned memory. 
        If the source and destination fields are each contained within single 
        words of memory, using the Fortran intrinsic MVBITS is more efficient 
        for Fortran programmers.

        The bcopy routine is an alternative to cb_copy_bits for C programmers
        if the source and destination regions are aligned on byte boundaries
        and the number of bits copied is a multiple of 8. cb_move_bits is
        slightly faster than bcopy, but bcopy does not have the restriction
        of non-overlapping regions.

        cb_copy_bits is single-threaded (i.e. not tasked) and may be called 
        from within a parallel region.
               
SEE ALSO

        MVBITS(3I), bstring(3C), INTRO_LIBCBL(3B)


*/
