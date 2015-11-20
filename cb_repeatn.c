/* cb_repeatn.c                                  http://cbl.sourceforge.net
 *
 * find short tandem repeats in a nucleotide string
 * see original man page at the bottom. This code is the generic version
 * with provision for 32 and 64 bit long ints.
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
 * $Id: cb_repeatn.c,v 1.4 2003/09/25 00:53:27 jlong777 Exp $
 */
 
#include "cb_macro.h"
#include <stdio.h>

#ifdef LONG32 /* need a 64-bit array for pattern, use long long */
void cb_repeatn(long *db, long dblen, long repeat_len, long min_repeats,
                long long *pattern, long *location, long *num_found)
#endif
#ifdef LONG64
void cb_repeatn(long *db, long dblen, long repeat_len, long min_repeats,
                long *pattern, long *location, long *num_found)
#endif
{
#ifdef B_ENDIAN
#ifdef LONG32
  long i, k=0, arr_size, shift;
  unsigned long dbshift, nextone, lastone, test;
  register unsigned long  mask, r0, sav;
#endif
#ifdef LONG64
  long i, j, k=0, m, sav[VECTLEN], shift;
  unsigned long  num_rep[VECTLEN], reject[VECTLEN], thisone[VECTLEN],
                 nextone[VECTLEN], lastone[VECTLEN];
  register unsigned long mask, r0;
  unsigned long dbshift[VECTLEN]; /* word containing pattern */
  
  long arr_size, sav2;
  unsigned long num_rep2, reject2, thisone2, nextone2, lastone2;
  unsigned long dbshift2; /* word containing pattern */
#endif
#endif

#ifdef L_ENDIAN /* same as big-endian 32-bit case */
  long i, k=0, arr_size, shift;
  unsigned long dbshift, nextone, lastone, test;
  register unsigned long  mask, r0, sav;
#endif

  if((repeat_len < 2) || (repeat_len > 16))
  {
    fprintf(stderr, "cb_repeatn: repeat_len range error, returning...\n");
    *num_found = -1;
    return;
  }
  
  if(*num_found <= 0)
  {
    fprintf(stderr, "cb_repeatn: num_found range error, returning...\n");
    *num_found = -1;
    return;
  }
  
  if(dblen <= repeat_len*min_repeats)
  {
    fprintf(stderr, "cb_repeatn: dblen <= repeat_len*min_repeats, returning...\n");
    *num_found = -1;
    return;
  }
  
  if(repeat_len*min_repeats > 32) 
    min_repeats = 32/repeat_len;

  arr_size = *num_found;
  *num_found = 0;
  
#ifdef B_ENDIAN
#ifdef LONG32
	
  switch(repeat_len)
  {
    case  2: mask = 0xF0000000; break;
    case  3: mask = 0xFC000000; break;
    case  4: mask = 0xFF000000; break;
    case  5: mask = 0xFFC00000; break;
    case  6: mask = 0xFFF00000; break;
    case  7: mask = 0xFFFC0000; break;
    case  8: mask = 0xFFFF0000; break;
    case  9: mask = 0xFFFFC000; break;
    case 10: mask = 0xFFFFF000; break;
    case 11: mask = 0xFFFFFC00; break;
    case 12: mask = 0xFFFFFF00; break;
    case 13: mask = 0xFFFFFFC0; break;
    case 14: mask = 0xFFFFFFF0; break;
    case 15: mask = 0xFFFFFFFC; break;
    case 16: mask = 0xFFFFFFFF; break;
    default: break;
  }

  r0 = 2*repeat_len;

  for(i=0; i<=(dblen-repeat_len*min_repeats+15)/16; i++) /* db index */
  {
    for(shift=0; shift<32; shift+=2) /* slide start of repeat along db */
    {
      if((16*i+shift/2) > (dblen-repeat_len*min_repeats)) break; /* end of db? */
      
      /* get word from which to construct pattern */
      dbshift = ((db[i]<<shift) | (shift?((unsigned long)db[i+1]>>(32-shift)):0x0));

      /* reject candidate if it is itself a repeat, 0=reject 1=no reject */
      switch(repeat_len)
      {
        case  2: test = (0xC0000000 & dbshift)^((0x30000000 & dbshift)<<2);
		 break;
	case  3: test = ((0xC0000000 & dbshift)^((0x30000000 & dbshift)<<2)) |
	                ((0xC0000000 & dbshift)^((0x0C000000 & dbshift)<<4));
		 break;
	case  4: test = (0xF0000000 & dbshift)^((0x0F000000 & dbshift)<<4);
	         break;
	case  5: test = ((0xC0000000 & dbshift)^((0x30000000 & dbshift)<<2)) |
	                ((0xC0000000 & dbshift)^((0x0C000000 & dbshift)<<4)) |
			((0xC0000000 & dbshift)^((0x03000000 & dbshift)<<6)) |
			((0xC0000000 & dbshift)^((0x00C00000 & dbshift)<<8));
		 break;
	case  6: test = (((0xF0000000 & dbshift)^((0x0F000000 & dbshift)<<4))  |
	                 ((0xF0000000 & dbshift)^((0x00F00000 & dbshift)<<8))) &&
			 ((0xFC000000 & dbshift)^((0x03F00000 & dbshift)<<6));
		 break;
	case  7: test = ((0xC0000000 & dbshift)^((0x30000000 & dbshift)<<2)) |
	                ((0xC0000000 & dbshift)^((0x0C000000 & dbshift)<<4)) |
			((0xC0000000 & dbshift)^((0x03000000 & dbshift)<<6)) |
			((0xC0000000 & dbshift)^((0x00C00000 & dbshift)<<8)) |
			((0xC0000000 & dbshift)^((0x00300000 & dbshift)<<10))|
			((0xC0000000 & dbshift)^((0x000C0000 & dbshift)<<12));
		 break;
	case  8: test = (0xFF000000 & dbshift)^((0x00FF0000 & dbshift)<<8);
		 break;
	case  9: test = ((0xFC000000 & dbshift)^((0x03F00000 & dbshift)<<6)) |
	                ((0xFC000000 & dbshift)^((0x000FC000 & dbshift)<<12));
		 break;
	case 10: test = (((0xF0000000 & dbshift)^((0x0F000000 & dbshift)<<4))  |
	                 ((0xF0000000 & dbshift)^((0x00F00000 & dbshift)<<8))  |
			 ((0xF0000000 & dbshift)^((0x000F0000 & dbshift)<<12)) |
			 ((0xF0000000 & dbshift)^((0x0000F000 & dbshift)<<16)))&&
			 ((0xFFC00000 & dbshift)^((0x003FF000 & dbshift)<<10));
	         break;
	case 11: test = ((0xC0000000 & dbshift)^((0x30000000 & dbshift)<<2)) |
	                ((0xC0000000 & dbshift)^((0x0C000000 & dbshift)<<4)) |
			((0xC0000000 & dbshift)^((0x03000000 & dbshift)<<6)) |
			((0xC0000000 & dbshift)^((0x00C00000 & dbshift)<<8)) |
			((0xC0000000 & dbshift)^((0x00300000 & dbshift)<<10))|
			((0xC0000000 & dbshift)^((0x000C0000 & dbshift)<<12))|
			((0xC0000000 & dbshift)^((0x00030000 & dbshift)<<14))|
			((0xC0000000 & dbshift)^((0x0000C000 & dbshift)<<16))|
			((0xC0000000 & dbshift)^((0x00003000 & dbshift)<<18))|
			((0xC0000000 & dbshift)^((0x00000C00 & dbshift)<<20));
	         break;
	case 12: test = (((0xFF000000 & dbshift)^((0x00FF0000 & dbshift)<<8))  |
	                 ((0xFF000000 & dbshift)^((0x0000FF00 & dbshift)<<16)))&&
			 ((0xFFF00000 & dbshift)^((0x000FFF00 & dbshift)<<12));
	         break;
	case 13: test = ((0xC0000000 & dbshift)^((0x30000000 & dbshift)<<2)) |
	                ((0xC0000000 & dbshift)^((0x0C000000 & dbshift)<<4)) |
			((0xC0000000 & dbshift)^((0x03000000 & dbshift)<<6)) |
			((0xC0000000 & dbshift)^((0x00C00000 & dbshift)<<8)) |
			((0xC0000000 & dbshift)^((0x00300000 & dbshift)<<10))|
			((0xC0000000 & dbshift)^((0x000C0000 & dbshift)<<12))|
			((0xC0000000 & dbshift)^((0x00030000 & dbshift)<<14))|
			((0xC0000000 & dbshift)^((0x0000C000 & dbshift)<<16))|
			((0xC0000000 & dbshift)^((0x00003000 & dbshift)<<18))|
			((0xC0000000 & dbshift)^((0x00000C00 & dbshift)<<20))|
			((0xC0000000 & dbshift)^((0x00000300 & dbshift)<<22))|
			((0xC0000000 & dbshift)^((0x000000C0 & dbshift)<<24));
	         break;
	case 14: test = (((0xF0000000 & dbshift)^((0x0F000000 & dbshift)<<4))  |
	                 ((0xF0000000 & dbshift)^((0x00F00000 & dbshift)<<8))  |
			 ((0xF0000000 & dbshift)^((0x000F0000 & dbshift)<<12)) |
			 ((0xF0000000 & dbshift)^((0x0000F000 & dbshift)<<16)) |
			 ((0xF0000000 & dbshift)^((0x00000F00 & dbshift)<<20)) |
			 ((0xF0000000 & dbshift)^((0x000000F0 & dbshift)<<24)))&&
			 ((0xFFFC0000 & dbshift)^((0x0003FFF0 & dbshift)<<14));
	         break;
	case 15: test = (((0xFC000000 & dbshift)^((0x03F00000 & dbshift)<<6))  |
	                 ((0xFC000000 & dbshift)^((0x000FC000 & dbshift)<<12)) |
			 ((0xFC000000 & dbshift)^((0x00003F00 & dbshift)<<18)) |
			 ((0xFC000000 & dbshift)^((0x000000FC & dbshift)<<24)))&&
			(((0xFFC00000 & dbshift)^((0x003FF000 & dbshift)<<10)) |
			 ((0xFFC00000 & dbshift)^((0x00000FFC & dbshift)<<20)));
		 break;
	case 16: test = (0xFFFF0000 & dbshift)^((0x0000FFFF & dbshift)<<16);
	         break;
      }
      
      if(test)
      {
        /* now save the word just before dbshift */
        if(i==0)
          lastone = shift?((unsigned long)db[i]>>(32-shift)):0x0;
        else
          lastone = (db[i-1]<<shift) | (shift?((unsigned long)db[i]>>(32-shift)):0x0);

        dbshift = dbshift & mask;
      
        /* is this pattern a replicate of the previous rep_len chars? */
        lastone = (lastone << (32-r0)) & mask;
        if((shift>=r0 || i>0) && dbshift == lastone) continue;
      
        test = 1; /* now is number of times pattern is actually repeated */
        sav = shift;
      
        nextone = ((db[i+(sav+r0)/32]<<((sav+r0)%32)) | (((sav+r0)%32)?
                  ((unsigned long)db[i+(sav+r0)/32+1] >> (32-(sav+r0)%32)): 0x0)) & mask;
		 
        while(dbshift == nextone   && 
             ((16*i + sav/2) < (dblen-r0+1))) /* don't overrun db */
        {
	  sav += r0;
	  nextone = ((db[i+(sav+r0)/32]<<((sav+r0)%32)) | (((sav+r0)%32)?
                    ((unsigned long)db[i+(sav+r0)/32+1] >> (32-(sav+r0)%32)): 0x0)) & mask;
          test++;
        }
      
        if(test >= min_repeats)
        {
          if(k < arr_size)
          {
	    location[k] = shift/2 + 16*i + 1;
	    (*num_found)++;
	    pattern[k] = dbshift;
	    pattern[k] = (pattern[k] << 30) | (0x3FFFFFFF & test);
	    k++;
	  }
	  else
	  {
	    *num_found *= -1;
            return;
	  }
        }
      }
    }
  }

#endif /* 32-bit */

#ifdef LONG64
/* vector version implemented only in this section, runs ok on non-vector */

  switch(repeat_len)
  {
    case  2: mask = 0xF000000000000000; break;
    case  3: mask = 0xFC00000000000000; break;
    case  4: mask = 0xFF00000000000000; break;
    case  5: mask = 0xFFC0000000000000; break;
    case  6: mask = 0xFFF0000000000000; break;
    case  7: mask = 0xFFFC000000000000; break;
    case  8: mask = 0xFFFF000000000000; break;
    case  9: mask = 0xFFFFC00000000000; break;
    case 10: mask = 0xFFFFF00000000000; break;
    case 11: mask = 0xFFFFFC0000000000; break;
    case 12: mask = 0xFFFFFF0000000000; break;
    case 13: mask = 0xFFFFFFC000000000; break;
    case 14: mask = 0xFFFFFFF000000000; break;
    case 15: mask = 0xFFFFFFFC00000000; break;
    case 16: mask = 0xFFFFFFFF00000000; break;
    default: break;
  }
	
  r0 = 2*repeat_len;

  for(i=0; i<=(dblen-repeat_len*min_repeats+31)/32 - VECTLEN; i+=VECTLEN) /* db index */
  {
    for(shift=0; shift<64; shift+=2) /* slide start of repeat along db */
    {
      if((32*i+shift/2) > (dblen-repeat_len*min_repeats)) break; /* end of db? */

      /* get word from which to construct pattern */
#pragma _CRI ivdep
      for(m=0; m<VECTLEN; m++)
        dbshift[m] = (db[i+m]<<shift) | (shift?((unsigned long)db[i+m+1]>>(64-shift)):0x0);

      for(m=0; m<VECTLEN; m++)
      {
        /* now save the word just before dbshift */
        if(i==0 && m==0)
          lastone[m] = 0x0 | (shift?((unsigned long)db[i+m]>>(64-shift)):0x0);
        else
          lastone[m] = (db[i+m-1]<<shift) | (shift?((unsigned long)db[i+m]>>(64-shift)):0x0);
      }
      
      switch(repeat_len)
      {
        /* reject candidate if it is itself a repeat, 0=reject 1=no reject */
        case 2:
                for(m=0; m<VECTLEN; m++)
                {
	          reject[m] = (0xC000000000000000 & dbshift[m]) ^
	                     ((0x3000000000000000 & dbshift[m]) << 2);
                } break;
        case 3:
                for(m=0; m<VECTLEN; m++)
                {
	          reject[m] = ((0xC000000000000000 & dbshift[m]) ^
	                      ((0x3000000000000000 & dbshift[m]) << 2)) |
	                      ((0xC000000000000000 & dbshift[m]) ^
                              ((0x0C00000000000000 & dbshift[m]) << 4));
                } break;
        case 4:
                for(m=0; m<VECTLEN; m++)
                {
	          reject[m] = (0xF000000000000000 & dbshift[m]) ^
	                     ((0x0F00000000000000 & dbshift[m]) << 4);
                } break;
        case 5:
                for(m=0; m<VECTLEN; m++)
                {
	          reject[m] = ((0xC000000000000000 & dbshift[m]) ^
	                      ((0x3000000000000000 & dbshift[m]) << 2)) |
	                      ((0xC000000000000000 & dbshift[m]) ^
                              ((0x0C00000000000000 & dbshift[m]) << 4)) |
                              ((0xC000000000000000 & dbshift[m]) ^
                              ((0x0300000000000000 & dbshift[m]) << 6)) |
                              ((0xC000000000000000 & dbshift[m]) ^
                              ((0x00C0000000000000 & dbshift[m]) << 8));
                } break;
        case 6:
                for(m=0; m<VECTLEN; m++)
                {
	          reject[m] = (((0xF000000000000000 & dbshift[m]) ^
                               ((0x0F00000000000000 & dbshift[m]) << 4))  |
                               ((0xF000000000000000 & dbshift[m]) ^
                               ((0x00F0000000000000 & dbshift[m]) << 8))) &&
                               ((0xFC00000000000000 & dbshift[m]) ^
                               ((0x03F0000000000000 & dbshift[m]) << 6));
                } break;
        case 7:
                for(m=0; m<VECTLEN; m++)
                {
	          reject[m] = ((0xC000000000000000 & dbshift[m]) ^
                              ((0x3000000000000000 & dbshift[m]) << 2)) |
                              ((0xC000000000000000 & dbshift[m]) ^
                              ((0x0C00000000000000 & dbshift[m]) << 4)) |
                              ((0xC000000000000000 & dbshift[m]) ^
                              ((0x0300000000000000 & dbshift[m]) << 6)) |
                              ((0xC000000000000000 & dbshift[m]) ^
                              ((0x00C0000000000000 & dbshift[m]) << 8)) |
                              ((0xC000000000000000 & dbshift[m]) ^
                              ((0x0030000000000000 & dbshift[m]) << 10))|
                              ((0xC000000000000000 & dbshift[m]) ^
                              ((0x000C000000000000 & dbshift[m]) << 12));
                } break;
        case 8:
                for(m=0; m<VECTLEN; m++)
                {
	          reject[m] = (0xFF00000000000000 & dbshift[m]) ^
                             ((0x00FF000000000000 & dbshift[m]) << 8);
                } break;
        case 9:
                for(m=0; m<VECTLEN; m++)
                {
	          reject[m] = ((0xFC00000000000000 & dbshift[m]) ^
                              ((0x03F0000000000000 & dbshift[m]) << 6)) |
                              ((0xFC00000000000000 & dbshift[m]) ^
                              ((0x000FC00000000000 & dbshift[m]) << 12));
                } break;
        case 10:
                for(m=0; m<VECTLEN; m++)
                {
	          reject[m] = (((0xF000000000000000 & dbshift[m]) ^
                               ((0x0F00000000000000 & dbshift[m]) << 4))  |
                               ((0xF000000000000000 & dbshift[m]) ^
                               ((0x00F0000000000000 & dbshift[m]) << 8))  |
                               ((0xF000000000000000 & dbshift[m]) ^
                               ((0x000F000000000000 & dbshift[m]) << 12)) |
                               ((0xF000000000000000 & dbshift[m]) ^
                               ((0x0000F00000000000 & dbshift[m]) << 16)))&&
                               ((0xFFC0000000000000 & dbshift[m]) ^
                               ((0x003FF00000000000 & dbshift[m]) << 10));
                } break;
        case 11:
                for(m=0; m<VECTLEN; m++)
                {
	          reject[m] = ((0xC000000000000000 & dbshift[m]) ^
	                      ((0x3000000000000000 & dbshift[m]) << 2)) |
	                      ((0xC000000000000000 & dbshift[m]) ^
		              ((0x0C00000000000000 & dbshift[m]) << 4)) |
		              ((0xC000000000000000 & dbshift[m]) ^
		              ((0x0300000000000000 & dbshift[m]) << 6)) |
		              ((0xC000000000000000 & dbshift[m]) ^
		              ((0x00C0000000000000 & dbshift[m]) << 8)) |
		              ((0xC000000000000000 & dbshift[m]) ^
		              ((0x0030000000000000 & dbshift[m]) << 10))|
		              ((0xC000000000000000 & dbshift[m]) ^
		              ((0x000C000000000000 & dbshift[m]) << 12))|
		              ((0xC000000000000000 & dbshift[m]) ^
		              ((0x0003000000000000 & dbshift[m]) << 14))|
		              ((0xC000000000000000 & dbshift[m]) ^
		              ((0x0000C00000000000 & dbshift[m]) << 16))|
		              ((0xC000000000000000 & dbshift[m]) ^
		              ((0x0000300000000000 & dbshift[m]) << 18))|
		              ((0xC000000000000000 & dbshift[m]) ^
		              ((0x00000C0000000000 & dbshift[m]) << 20));
                } break;
        case 12:
                for(m=0; m<VECTLEN; m++)
                {
	          reject[m] = (((0xFF00000000000000 & dbshift[m]) ^
	                       ((0x00FF000000000000 & dbshift[m]) << 8))  |
	                       ((0xFF00000000000000 & dbshift[m]) ^
		               ((0x0000FF0000000000 & dbshift[m]) << 16)))&&
		               ((0xFFF0000000000000 & dbshift[m]) ^
		               ((0x000FFF0000000000 & dbshift[m]) << 12));
                } break;
        case 13:
                for(m=0; m<VECTLEN; m++)
                {
	          reject[m] = ((0xC000000000000000 & dbshift[m]) ^
	                      ((0x3000000000000000 & dbshift[m]) << 2)) |
	                      ((0xC000000000000000 & dbshift[m]) ^
		              ((0x0C00000000000000 & dbshift[m]) << 4)) |
		              ((0xC000000000000000 & dbshift[m]) ^
		              ((0x0300000000000000 & dbshift[m]) << 6)) |
		              ((0xC000000000000000 & dbshift[m]) ^
		              ((0x00C0000000000000 & dbshift[m]) << 8)) |
		              ((0xC000000000000000 & dbshift[m]) ^
		              ((0x0030000000000000 & dbshift[m]) << 10))|
		              ((0xC000000000000000 & dbshift[m]) ^
		              ((0x000C000000000000 & dbshift[m]) << 12))|
		              ((0xC000000000000000 & dbshift[m]) ^
		              ((0x0003000000000000 & dbshift[m]) << 14))|
		              ((0xC000000000000000 & dbshift[m]) ^
		              ((0x0000C00000000000 & dbshift[m]) << 16))|
		              ((0xC000000000000000 & dbshift[m]) ^
		              ((0x0000300000000000 & dbshift[m]) << 18))|
		              ((0xC000000000000000 & dbshift[m]) ^
		              ((0x00000C0000000000 & dbshift[m]) << 20))|
		              ((0xC000000000000000 & dbshift[m]) ^
		              ((0x0000030000000000 & dbshift[m]) << 22))|
		              ((0xC000000000000000 & dbshift[m]) ^
		              ((0x000000C000000000 & dbshift[m]) << 24));
                } break;
        case 14:
                for(m=0; m<VECTLEN; m++)
                {
	          reject[m] = (((0xF000000000000000 & dbshift[m]) ^
	                       ((0x0F00000000000000 & dbshift[m]) << 4))  |
	                       ((0xF000000000000000 & dbshift[m]) ^
		               ((0x00F0000000000000 & dbshift[m]) << 8))  |
		               ((0xF000000000000000 & dbshift[m]) ^
		               ((0x000F000000000000 & dbshift[m]) << 12)) |
		               ((0xF000000000000000 & dbshift[m]) ^
		               ((0x0000F00000000000 & dbshift[m]) << 16)) |
		               ((0xF000000000000000 & dbshift[m]) ^
		               ((0x00000F0000000000 & dbshift[m]) << 20)) |
		               ((0xF000000000000000 & dbshift[m]) ^
		               ((0x000000F000000000 & dbshift[m]) << 24)))&&
		               ((0xFFFC000000000000 & dbshift[m]) ^
		               ((0x0003FFF000000000 & dbshift[m]) << 14));
                } break;
        case 15:
                for(m=0; m<VECTLEN; m++)
                {
	          reject[m] = (((0xFC00000000000000 & dbshift[m]) ^
	                       ((0x03F0000000000000 & dbshift[m]) << 6))  |
	                       ((0xFC00000000000000 & dbshift[m]) ^
		               ((0x000FC00000000000 & dbshift[m]) << 12)) |
		               ((0xFC00000000000000 & dbshift[m]) ^
		               ((0x00003F0000000000 & dbshift[m]) << 18)) |
		               ((0xFC00000000000000 & dbshift[m]) ^
		               ((0x000000FC00000000 & dbshift[m]) << 24)))&&
		              (((0xFFC0000000000000 & dbshift[m]) ^
		               ((0x003FF00000000000 & dbshift[m]) << 10)) |
		               ((0xFFC0000000000000 & dbshift[m]) ^
		               ((0x00000FFC00000000 & dbshift[m]) << 20)));
                } break;
        case 16:
                for(m=0; m<VECTLEN; m++)
                {
	          reject[m] = (0xFFFF000000000000 & dbshift[m]) ^
	                     ((0x0000FFFF00000000 & dbshift[m]) << 16);
                } break;
        default:  break;
      }
#pragma _CRI ivdep
      for(m=0; m<VECTLEN; m++)
      {	
        num_rep[m] = 1;
      
        if(!reject[m]) num_rep[m] = 0;
        
        thisone[m] = dbshift[m] & mask;
      
        /* is this pattern a replicate of the previous rep_len chars? */
        lastone[m] = (lastone[m] << (64-r0)) & mask;
        if((shift>=r0 || i>0) && thisone[m] == lastone[m]) num_rep[m] = 0;
      
        sav[m] = shift; /* save shift for next loop */
      
        nextone[m] = ((db[i+m+(shift+r0)/64]<<((shift+r0)%64)) | 
	            (((shift+r0)%64)?
                     ((unsigned long)db[i+m+(shift+r0)/64+1] >>
		      (64-(shift+r0)%64)): 0x0)) & mask;
      }
      
      /* is pattern repeated min_repeat times? */
      for(j=1; j<min_repeats && ((32*i + sav[m]/2) < (dblen-r0+1)); j++)
      {
        for(m=0; m<VECTLEN; m++)
        {
	  if(num_rep[m] && (thisone[m] == nextone[m]))
	  {
	    sav[m] += r0;
	    nextone[m] = ((db[i+(sav[m]+r0)/64]<<((sav[m]+r0)%64)) | 
	                (((sav[m]+r0)%64)?
                         ((unsigned long)db[i+(sav[m]+r0)/64+1] >>
		          (64-(sav[m]+r0)%64)): 0x0)) & mask;
            num_rep[m]++;
	  }
	}
      }
      
      /* from qualifiers, count number of times pattern is actually repeated */
      for(m=0; m<VECTLEN; m++)
      {
        if(num_rep[m] >= min_repeats)
	{
          while(thisone[m] == nextone[m]   && 
               ((32*i + sav[m]/2) < (dblen-r0+1))) /* don't overrun db */
          {
	    sav[m] += r0;
	    nextone[m] = ((db[i+(sav[m]+r0)/64]<<((sav[m]+r0)%64)) | 
	                (((sav[m]+r0)%64)?
                         ((unsigned long)db[i+(sav[m]+r0)/64+1] >>
		          (64-(sav[m]+r0)%64)): 0x0)) & mask;
            num_rep[m]++;
          }
	  if(k < arr_size)
          {
	    /* store results */
	    location[k] = shift/2 + 32*i + 1;
	    (*num_found)++;
	    pattern[k] = thisone[m] | (num_rep[m] & 0x000000003FFFFFFF);
	    k++;
	  }
	  else
	  {
	    *num_found *= -1;
            return;
	  }
	}
      }
    }
  }
  
  /* do the final segment of db that was not a multiple of the vector length */
  if(i)
    i = i-VECTLEN+1;
    
  for(; i<=(dblen-repeat_len*min_repeats+31)/32; i++)
  {
    for(shift=0; shift<64; shift+=2) /* slide start of repeat along db */
    {
      if((32*i+shift/2) > (dblen-repeat_len*min_repeats)) break; /* end of db? */

      /* get word from which to construct pattern */
      dbshift2 = (db[i]<<shift) | (shift?((unsigned long)db[i+1]>>(64-shift)):0x0);

      /* now save the word just before dbshift */
      if(i==0)
        lastone2 = 0x0 | (shift?((unsigned long)db[i]>>(64-shift)):0x0);
      else
        lastone2 = (db[i-1]<<shift) | (shift?((unsigned long)db[i]>>(64-shift)):0x0);

      /* reject candidate if it is itself a repeat, 0=reject 1=no reject */
      switch(repeat_len)
      {
        case  2: reject2 = (0xC000000000000000 & dbshift2) ^
	                 ((0x3000000000000000 & dbshift2) << 2);
		 break;
	case  3: reject2 = ((0xC000000000000000 & dbshift2) ^
	                  ((0x3000000000000000 & dbshift2) << 2)) |
	                  ((0xC000000000000000 & dbshift2) ^
			  ((0x0C00000000000000 & dbshift2) << 4));
		 break;
	case  4: reject2 = (0xF000000000000000 & dbshift2) ^
	                 ((0x0F00000000000000 & dbshift2) << 4);
	         break;
	case  5: reject2 = ((0xC000000000000000 & dbshift2) ^
	                  ((0x3000000000000000 & dbshift2) << 2)) |
	                  ((0xC000000000000000 & dbshift2) ^
			  ((0x0C00000000000000 & dbshift2) << 4)) |
			  ((0xC000000000000000 & dbshift2) ^
			  ((0x0300000000000000 & dbshift2) << 6)) |
			  ((0xC000000000000000 & dbshift2) ^
			  ((0x00C0000000000000 & dbshift2) << 8));
		 break;
	case  6: reject2 = (((0xF000000000000000 & dbshift2) ^
	                   ((0x0F00000000000000 & dbshift2) << 4))  |
	                   ((0xF000000000000000 & dbshift2) ^
			   ((0x00F0000000000000 & dbshift2) << 8))) &&
			   ((0xFC00000000000000 & dbshift2) ^
			   ((0x03F0000000000000 & dbshift2) << 6));
		 break;
	case  7: reject2 = ((0xC000000000000000 & dbshift2) ^
	                  ((0x3000000000000000 & dbshift2) << 2)) |
	                  ((0xC000000000000000 & dbshift2) ^
			  ((0x0C00000000000000 & dbshift2) << 4)) |
			  ((0xC000000000000000 & dbshift2) ^
			  ((0x0300000000000000 & dbshift2) << 6)) |
			  ((0xC000000000000000 & dbshift2) ^
			  ((0x00C0000000000000 & dbshift2) << 8)) |
			  ((0xC000000000000000 & dbshift2) ^
			  ((0x0030000000000000 & dbshift2) << 10))|
			  ((0xC000000000000000 & dbshift2) ^
			  ((0x000C000000000000 & dbshift2) << 12));
		 break;
	case  8: reject2 = (0xFF00000000000000 & dbshift2) ^
	                 ((0x00FF000000000000 & dbshift2) << 8);
		 break;
	case  9: reject2 = ((0xFC00000000000000 & dbshift2) ^
	                  ((0x03F0000000000000 & dbshift2) << 6)) |
	                  ((0xFC00000000000000 & dbshift2) ^
			  ((0x000FC00000000000 & dbshift2) << 12));
		 break;
	case 10: reject2 = (((0xF000000000000000 & dbshift2) ^
	                   ((0x0F00000000000000 & dbshift2) << 4))  |
	                   ((0xF000000000000000 & dbshift2) ^
			   ((0x00F0000000000000 & dbshift2) << 8))  |
			   ((0xF000000000000000 & dbshift2) ^
			   ((0x000F000000000000 & dbshift2) << 12)) |
			   ((0xF000000000000000 & dbshift2) ^
			   ((0x0000F00000000000 & dbshift2) << 16)))&&
			   ((0xFFC0000000000000 & dbshift2) ^
			   ((0x003FF00000000000 & dbshift2) << 10));
	         break;
	case 11: reject2 = ((0xC000000000000000 & dbshift2) ^
	                  ((0x3000000000000000 & dbshift2) << 2)) |
	                  ((0xC000000000000000 & dbshift2) ^
			  ((0x0C00000000000000 & dbshift2) << 4)) |
			  ((0xC000000000000000 & dbshift2) ^
			  ((0x0300000000000000 & dbshift2) << 6)) |
			  ((0xC000000000000000 & dbshift2) ^
			  ((0x00C0000000000000 & dbshift2) << 8)) |
			  ((0xC000000000000000 & dbshift2) ^
			  ((0x0030000000000000 & dbshift2) << 10))|
			  ((0xC000000000000000 & dbshift2) ^
			  ((0x000C000000000000 & dbshift2) << 12))|
			  ((0xC000000000000000 & dbshift2) ^
			  ((0x0003000000000000 & dbshift2) << 14))|
			  ((0xC000000000000000 & dbshift2) ^
			  ((0x0000C00000000000 & dbshift2) << 16))|
			  ((0xC000000000000000 & dbshift2) ^
			  ((0x0000300000000000 & dbshift2) << 18))|
			  ((0xC000000000000000 & dbshift2) ^
			  ((0x00000C0000000000 & dbshift2) << 20));
	         break;
	case 12: reject2 = (((0xFF00000000000000 & dbshift2) ^
	                   ((0x00FF000000000000 & dbshift2) << 8))  |
	                   ((0xFF00000000000000 & dbshift2) ^
			   ((0x0000FF0000000000 & dbshift2) << 16)))&&
			   ((0xFFF0000000000000 & dbshift2) ^
			   ((0x000FFF0000000000 & dbshift2) << 12));
	         break;
	case 13: reject2 = ((0xC000000000000000 & dbshift2) ^
	                  ((0x3000000000000000 & dbshift2) << 2)) |
	                  ((0xC000000000000000 & dbshift2) ^
			  ((0x0C00000000000000 & dbshift2) << 4)) |
			  ((0xC000000000000000 & dbshift2) ^
			  ((0x0300000000000000 & dbshift2) << 6)) |
			  ((0xC000000000000000 & dbshift2) ^
			  ((0x00C0000000000000 & dbshift2) << 8)) |
			  ((0xC000000000000000 & dbshift2) ^
			  ((0x0030000000000000 & dbshift2) << 10))|
			  ((0xC000000000000000 & dbshift2) ^
			  ((0x000C000000000000 & dbshift2) << 12))|
			  ((0xC000000000000000 & dbshift2) ^
			  ((0x0003000000000000 & dbshift2) << 14))|
			  ((0xC000000000000000 & dbshift2) ^
			  ((0x0000C00000000000 & dbshift2) << 16))|
			  ((0xC000000000000000 & dbshift2) ^
			  ((0x0000300000000000 & dbshift2) << 18))|
			  ((0xC000000000000000 & dbshift2) ^
			  ((0x00000C0000000000 & dbshift2) << 20))|
			  ((0xC000000000000000 & dbshift2) ^
			  ((0x0000030000000000 & dbshift2) << 22))|
			  ((0xC000000000000000 & dbshift2) ^
			  ((0x000000C000000000 & dbshift2) << 24));
	         break;
	case 14: reject2 = (((0xF000000000000000 & dbshift2) ^
	                   ((0x0F00000000000000 & dbshift2) << 4))  |
	                   ((0xF000000000000000 & dbshift2) ^
			   ((0x00F0000000000000 & dbshift2) << 8))  |
			   ((0xF000000000000000 & dbshift2) ^
			   ((0x000F000000000000 & dbshift2) << 12)) |
			   ((0xF000000000000000 & dbshift2) ^
			   ((0x0000F00000000000 & dbshift2) << 16)) |
			   ((0xF000000000000000 & dbshift2) ^
			   ((0x00000F0000000000 & dbshift2) << 20)) |
			   ((0xF000000000000000 & dbshift2) ^
			   ((0x000000F000000000 & dbshift2) << 24)))&&
			   ((0xFFFC000000000000 & dbshift2) ^
			   ((0x0003FFF000000000 & dbshift2) << 14));
	         break;
	case 15: reject2 = (((0xFC00000000000000 & dbshift2) ^
	                   ((0x03F0000000000000 & dbshift2) << 6))  |
	                   ((0xFC00000000000000 & dbshift2) ^
			   ((0x000FC00000000000 & dbshift2) << 12)) |
			   ((0xFC00000000000000 & dbshift2) ^
			   ((0x00003F0000000000 & dbshift2) << 18)) |
			   ((0xFC00000000000000 & dbshift2) ^
			   ((0x000000FC00000000 & dbshift2) << 24)))&&
			  (((0xFFC0000000000000 & dbshift2) ^
			   ((0x003FF00000000000 & dbshift2) << 10)) |
			   ((0xFFC0000000000000 & dbshift2) ^
			   ((0x00000FFC00000000 & dbshift2) << 20)));
		 break;
	case 16: reject2 = (0xFFFF000000000000 & dbshift2) ^
	                  ((0x0000FFFF00000000 & dbshift2) << 16);
	         break;
      }
      
      if(!reject2) continue;

      thisone2 = dbshift2 & mask;
      
      /* is this pattern a replicate of the previous rep_len chars? */
      lastone2 = (lastone2 << (64-r0)) & mask;
      if((shift>=r0 || i>0) && thisone2 == lastone2) continue;
      
      /* count number of times pattern is actually repeated */
      
      num_rep2 = 1;
      sav2 = shift;
      
      nextone2 = ((db[i+(sav2+r0)/64]<<((sav2+r0)%64)) | 
	       (((sav2+r0)%64)?
                ((unsigned long)db[i+(sav2+r0)/64+1] >>
		 (64-(sav2+r0)%64)): 0x0)) & mask;
		 
      while(thisone2 == nextone2   && 
           ((32*i + sav2/2) < (dblen-r0+1))) /* don't overrun db */
      {
	sav2 += r0;
	nextone2 = ((db[i+(sav2+r0)/64]<<((sav2+r0)%64)) | 
	         (((sav2+r0)%64)?
                  ((unsigned long)db[i+(sav2+r0)/64+1] >>
		   (64-(sav2+r0)%64)): 0x0)) & mask;
        num_rep2++;
      }
      
      if(num_rep2 >= min_repeats)
      {
	if(k < arr_size)
        {
	  location[k] = shift/2 + 32*i + 1;
	  (*num_found)++;
	  pattern[k] = thisone2 | (num_rep2 & 0x000000003FFFFFFF);
	  k++;
	}
	else
	{
	  *num_found *= -1;
          return;
	}
      }
    }
  }

#endif /* 64-bit */
#endif /* big_endian */

#ifdef L_ENDIAN
/*
  pattern	output array of packed words containing 
  		the repeated pattern and the number of
		times it was repeated, in this format:
                            
		------------------------------------------
		| count (30 bits) | pattern (32 bits) |00|
		------------------------------------------
*/
#ifdef LONG32
	
  switch(repeat_len)
  {
    case  2: mask = 0x0000000F; break;
    case  3: mask = 0x0000003F; break;
    case  4: mask = 0x000000FF; break;
    case  5: mask = 0x000003FF; break;
    case  6: mask = 0x00000FFF; break;
    case  7: mask = 0x00003FFF; break;
    case  8: mask = 0x0000FFFF; break;
    case  9: mask = 0x0003FFFF; break;
    case 10: mask = 0x000FFFFF; break;
    case 11: mask = 0x003FFFFF; break;
    case 12: mask = 0x00FFFFFF; break;
    case 13: mask = 0x03FFFFFF; break;
    case 14: mask = 0x0FFFFFFF; break;
    case 15: mask = 0x3FFFFFFF; break;
    case 16: mask = 0xFFFFFFFF; break;
    default: break;
  }
  
  r0 = 2*repeat_len;

  for(i=0; i<=(dblen-repeat_len*min_repeats+15)/16; i++) /* db index */
  {
    for(shift=0; shift<32; shift+=2) /* slide start of repeat along db */
    {
      if((16*i+shift/2) > (dblen-repeat_len*min_repeats)) break; /* end of db? */
      
      /* get word from which to construct pattern */
      dbshift = (((unsigned long)db[i]>>shift) | (shift?(db[i+1]<<(32-shift)):0x0));

      /* reject candidate if it is itself a repeat, 0=reject 1=no reject */
      switch(repeat_len)
      {
        case  2: test = (0x00000003 & dbshift)^((0x0000000C & dbshift)>>2);
		 break;
	case  3: test = ((0x00000003 & dbshift)^((0x0000000C & dbshift)>>2)) |
	                ((0x00000003 & dbshift)^((0x00000030 & dbshift)>>4));
		 break;
	case  4: test = (0x0000000F & dbshift)^((0x000000F0 & dbshift)>>4);
	         break;
	case  5: test = ((0x00000003 & dbshift)^((0x0000000C & dbshift)>>2)) |
	                ((0x00000003 & dbshift)^((0x00000030 & dbshift)>>4)) |
			((0x00000003 & dbshift)^((0x000000C0 & dbshift)>>6)) |
			((0x00000003 & dbshift)^((0x00000300 & dbshift)>>8));
		 break;
	case  6: test = (((0x0000000F & dbshift)^((0x000000F0 & dbshift)>>4))  |
	                 ((0x0000000F & dbshift)^((0x00000F00 & dbshift)>>8))) &&
			 ((0x0000003F & dbshift)^((0x00000FC0 & dbshift)>>6));
		 break;
	case  7: test = ((0x00000003 & dbshift)^((0x0000000C & dbshift)>>2)) |
	                ((0x00000003 & dbshift)^((0x00000030 & dbshift)>>4)) |
			((0x00000003 & dbshift)^((0x000000C0 & dbshift)>>6)) |
			((0x00000003 & dbshift)^((0x00000300 & dbshift)>>8)) |
			((0x00000003 & dbshift)^((0x00000C00 & dbshift)>>10))|
			((0x00000003 & dbshift)^((0x00003000 & dbshift)>>12));
		 break;
	case  8: test = (0x000000FF & dbshift)^((0x0000FF00 & dbshift)>>8);
		 break;
	case  9: test = ((0x0000003F & dbshift)^((0x00000FC0 & dbshift)>>6)) |
	                ((0x0000003F & dbshift)^((0x0003F000 & dbshift)>>12));
		 break;
	case 10: test = (((0x0000000F & dbshift)^((0x000000F0 & dbshift)>>4))  |
	                 ((0x0000000F & dbshift)^((0x00000F00 & dbshift)>>8))  |
			 ((0x0000000F & dbshift)^((0x0000F000 & dbshift)>>12)) |
			 ((0x0000000F & dbshift)^((0x000F0000 & dbshift)>>16)))&&
			 ((0x000003FF & dbshift)^((0x000FFC00 & dbshift)>>10));
	         break;
	case 11: test = ((0x00000003 & dbshift)^((0x0000000C & dbshift)>>2)) |
	                ((0x00000003 & dbshift)^((0x00000030 & dbshift)>>4)) |
			((0x00000003 & dbshift)^((0x000000C0 & dbshift)>>6)) |
			((0x00000003 & dbshift)^((0x00000300 & dbshift)>>8)) |
			((0x00000003 & dbshift)^((0x00000C00 & dbshift)>>10))|
			((0x00000003 & dbshift)^((0x00003000 & dbshift)>>12))|
			((0x00000003 & dbshift)^((0x0000C000 & dbshift)>>14))|
			((0x00000003 & dbshift)^((0x00030000 & dbshift)>>16))|
			((0x00000003 & dbshift)^((0x000C0000 & dbshift)>>18))|
			((0x00000003 & dbshift)^((0x00300000 & dbshift)>>20));
	         break;
	case 12: test = (((0x000000FF & dbshift)^((0x0000FF00 & dbshift)>>8))  |
	                 ((0x000000FF & dbshift)^((0x00FF0000 & dbshift)>>16)))&&
			 ((0x00000FFF & dbshift)^((0x00FFF000 & dbshift)>>12));
	         break;
	case 13: test = ((0x00000003 & dbshift)^((0x0000000C & dbshift)>>2)) |
	                ((0x00000003 & dbshift)^((0x00000030 & dbshift)>>4)) |
			((0x00000003 & dbshift)^((0x000000C0 & dbshift)>>6)) |
			((0x00000003 & dbshift)^((0x00000300 & dbshift)>>8)) |
			((0x00000003 & dbshift)^((0x00000C00 & dbshift)>>10))|
			((0x00000003 & dbshift)^((0x00003000 & dbshift)>>12))|
			((0x00000003 & dbshift)^((0x0000C000 & dbshift)>>14))|
			((0x00000003 & dbshift)^((0x00030000 & dbshift)>>16))|
			((0x00000003 & dbshift)^((0x000C0000 & dbshift)>>18))|
			((0x00000003 & dbshift)^((0x00300000 & dbshift)>>20))|
			((0x00000003 & dbshift)^((0x00C00000 & dbshift)>>22))|
			((0x00000003 & dbshift)^((0x03000000 & dbshift)>>24));
	         break;
	case 14: test = (((0x0000000F & dbshift)^((0x000000F0 & dbshift)>>4))  |
	                 ((0x0000000F & dbshift)^((0x00000F00 & dbshift)>>8))  |
			 ((0x0000000F & dbshift)^((0x0000F000 & dbshift)>>12)) |
			 ((0x0000000F & dbshift)^((0x000F0000 & dbshift)>>16)) |
			 ((0x0000000F & dbshift)^((0x00F00000 & dbshift)>>20)) |
			 ((0x0000000F & dbshift)^((0x0F000000 & dbshift)>>24)))&&
			 ((0x00003FFF & dbshift)^((0x0FFFC000 & dbshift)>>14));
	         break;
	case 15: test = (((0x0000003F & dbshift)^((0x00000FC0 & dbshift)>>6))  |
	                 ((0x0000003F & dbshift)^((0x0003F000 & dbshift)>>12)) |
			 ((0x0000003F & dbshift)^((0x00FC0000 & dbshift)>>18)) |
			 ((0x0000003F & dbshift)^((0x3F000000 & dbshift)>>24)))&&
			(((0x000003FF & dbshift)^((0x000FFC00 & dbshift)>>10)) |
			 ((0x000003FF & dbshift)^((0x3FF00000 & dbshift)>>20)));
		 break;
	case 16: test = (0x0000FFFF & dbshift)^((0xFFFF0000 & dbshift)>>16);
	         break;
      }
      
      if(test)
      {
        /* now save the word just before dbshift */
        if(i==0)
          lastone = shift?(db[i]<<(32-shift)):0x0;
        else
          lastone = ((unsigned long)db[i-1]>>shift) | (shift?(db[i]<<(32-shift)):0x0);

        dbshift = dbshift & mask;
      
        /* is this pattern a replicate of the previous rep_len chars? */
        lastone = (lastone >> (32-r0)) & mask;
        if((shift>=r0 || i>0) && (dbshift == lastone)) continue;
      
        test = 1; /* now is number of times pattern is actually repeated */
        sav = shift;
      
        nextone = (((unsigned long)db[i+(sav+r0)/32]>>((sav+r0)%32)) | (((sav+r0)%32)?
                  (db[i+(sav+r0)/32+1] << (32-(sav+r0)%32)): 0x0)) & mask;

        while(dbshift == nextone   && 
             ((16*i + sav/2) < (dblen-r0+1))) /* don't overrun db */
        {
	  sav += r0;
	  nextone = (((unsigned long)db[i+(sav+r0)/32]>>((sav+r0)%32)) | (((sav+r0)%32)?
                    (db[i+(sav+r0)/32+1] << (32-(sav+r0)%32)): 0x0)) & mask;
          test++;
        }
      
        if(test >= min_repeats)
        {
	  if(k < arr_size)
          {
	    location[k] = shift/2 + 16*i + 1;
	    (*num_found)++;
	    pattern[k] = test;
	    pattern[k] = ((pattern[k] << 32) | dbshift) << 2;
	    k++;
	  }
	  else
	  {
	    *num_found *= -1;
            return;
	  }
        }
      }
    }
  }

#endif /* 32-bit */

#ifdef LONG64

  switch(repeat_len)
  {
    case  2: mask = 0x000000000000000F; break;
    case  3: mask = 0x000000000000003F; break;
    case  4: mask = 0x00000000000000FF; break;
    case  5: mask = 0x00000000000003FF; break;
    case  6: mask = 0x0000000000000FFF; break;
    case  7: mask = 0x0000000000003FFF; break;
    case  8: mask = 0x000000000000FFFF; break;
    case  9: mask = 0x000000000003FFFF; break;
    case 10: mask = 0x00000000000FFFFF; break;
    case 11: mask = 0x00000000003FFFFF; break;
    case 12: mask = 0x0000000000FFFFFF; break;
    case 13: mask = 0x0000000003FFFFFF; break;
    case 14: mask = 0x000000000FFFFFFF; break;
    case 15: mask = 0x000000003FFFFFFF; break;
    case 16: mask = 0x00000000FFFFFFFF; break;
    default: break;
  }
  
  r0 = 2*repeat_len;

  for(i=0; i<=(dblen-repeat_len*min_repeats+31)/32; i++) /* db index */
  {
    for(shift=0; shift<64; shift+=2) /* slide start of repeat along db */
    {
      if((32*i+shift/2) > (dblen-repeat_len*min_repeats)) break; /* end of db? */
      
      /* get word from which to construct pattern */
      dbshift = (((unsigned long)db[i]>>shift) | (shift?(db[i+1]<<(64-shift)):0x0));

      /* reject candidate if it is itself a repeat, 0=reject 1=no reject */
      switch(repeat_len)
      {
        case  2: test = (0x0000000000000003 & dbshift) ^
	               ((0x000000000000000C & dbshift) >> 2);
		 break;
	case  3: test = ((0x0000000000000003 & dbshift) ^
	                ((0x000000000000000C & dbshift) >> 2)) |
	                ((0x0000000000000003 & dbshift) ^
			((0x0000000000000030 & dbshift) >> 4));
		 break;
	case  4: test = (0x000000000000000F & dbshift) ^
	               ((0x00000000000000F0 & dbshift) >> 4);
	         break;
	case  5: test = ((0x0000000000000003 & dbshift) ^
	                ((0x000000000000000C & dbshift) >> 2)) |
	                ((0x0000000000000003 & dbshift) ^
			((0x0000000000000030 & dbshift) >> 4)) |
			((0x0000000000000003 & dbshift) ^
			((0x00000000000000C0 & dbshift) >> 6)) |
			((0x0000000000000003 & dbshift) ^
			((0x0000000000000300 & dbshift) >> 8));
		 break;
	case  6: test = (((0x000000000000000F & dbshift) ^
	                 ((0x00000000000000F0 & dbshift) >> 4))  |
	                 ((0x000000000000000F & dbshift) ^
			 ((0x0000000000000F00 & dbshift) >> 8))) &&
			 ((0x000000000000003F & dbshift) ^
			 ((0x0000000000000FC0 & dbshift) >> 6));
		 break;
	case  7: test = ((0x0000000000000003 & dbshift) ^
	                ((0x000000000000000C & dbshift) >> 2)) |
	                ((0x0000000000000003 & dbshift) ^
			((0x0000000000000030 & dbshift) >> 4)) |
			((0x0000000000000003 & dbshift) ^
			((0x00000000000000C0 & dbshift) >> 6)) |
			((0x0000000000000003 & dbshift) ^
			((0x0000000000000300 & dbshift) >> 8)) |
			((0x0000000000000003 & dbshift) ^
			((0x0000000000000C00 & dbshift) >> 10))|
			((0x0000000000000003 & dbshift) ^
			((0x0000000000003000 & dbshift) >> 12));
		 break;
	case  8: test = (0x00000000000000FF & dbshift) ^
	               ((0x000000000000FF00 & dbshift) >> 8);
		 break;
	case  9: test = ((0x000000000000003F & dbshift) ^
	                ((0x0000000000000FC0 & dbshift) >> 6)) |
	                ((0x000000000000003F & dbshift) ^
			((0x000000000003F000 & dbshift) >> 12));
		 break;
	case 10: test = (((0x000000000000000F & dbshift) ^
	                 ((0x00000000000000F0 & dbshift) >> 4))  |
	                 ((0x000000000000000F & dbshift) ^
			 ((0x0000000000000F00 & dbshift) >> 8))  |
			 ((0x000000000000000F & dbshift) ^
			 ((0x000000000000F000 & dbshift) >> 12)) |
			 ((0x000000000000000F & dbshift) ^
			 ((0x00000000000F0000 & dbshift) >> 16)))&&
			 ((0x00000000000003FF & dbshift) ^
			 ((0x00000000000FFC00 & dbshift) >> 10));
	         break;
	case 11: test = ((0x0000000000000003 & dbshift) ^
	                ((0x000000000000000C & dbshift) >> 2)) |
	                ((0x0000000000000003 & dbshift) ^
			((0x0000000000000030 & dbshift) >> 4)) |
			((0x0000000000000003 & dbshift) ^
			((0x00000000000000C0 & dbshift) >> 6)) |
			((0x0000000000000003 & dbshift) ^
			((0x0000000000000300 & dbshift) >> 8)) |
			((0x0000000000000003 & dbshift) ^
			((0x0000000000000C00 & dbshift) >> 10))|
			((0x0000000000000003 & dbshift) ^
			((0x0000000000003000 & dbshift) >> 12))|
			((0x0000000000000003 & dbshift) ^
			((0x000000000000C000 & dbshift) >> 14))|
			((0x0000000000000003 & dbshift) ^
			((0x0000000000030000 & dbshift) >> 16))|
			((0x0000000000000003 & dbshift) ^
			((0x00000000000C0000 & dbshift) >> 18))|
			((0x0000000000000003 & dbshift) ^
			((0x0000000000300000 & dbshift) >> 20));
	         break;
	case 12: test = (((0x00000000000000FF & dbshift) ^
	                 ((0x000000000000FF00 & dbshift) >> 8))  |
	                 ((0x00000000000000FF & dbshift) ^
			 ((0x0000000000FF0000 & dbshift) >> 16)))&&
			 ((0x0000000000000FFF & dbshift) ^
			 ((0x0000000000FFF000 & dbshift) >> 12));
	         break;
	case 13: test = ((0x0000000000000003 & dbshift) ^
	                ((0x000000000000000C & dbshift) >> 2)) |
	                ((0x0000000000000003 & dbshift) ^
			((0x0000000000000030 & dbshift) >> 4)) |
			((0x0000000000000003 & dbshift) ^
			((0x00000000000000C0 & dbshift) >> 6)) |
			((0x0000000000000003 & dbshift) ^
			((0x0000000000000300 & dbshift) >> 8)) |
			((0x0000000000000003 & dbshift) ^
			((0x0000000000000C00 & dbshift) >> 10))|
			((0x0000000000000003 & dbshift) ^
			((0x0000000000003000 & dbshift) >> 12))|
			((0x0000000000000003 & dbshift) ^
			((0x000000000000C000 & dbshift) >> 14))|
			((0x0000000000000003 & dbshift) ^
			((0x0000000000030000 & dbshift) >> 16))|
			((0x0000000000000003 & dbshift) ^
			((0x00000000000C0000 & dbshift) >> 18))|
			((0x0000000000000003 & dbshift) ^
			((0x0000000000300000 & dbshift) >> 20))|
			((0x0000000000000003 & dbshift) ^
			((0x0000000000C00000 & dbshift) >> 22))|
			((0x0000000000000003 & dbshift) ^
			((0x0000000003000000 & dbshift) >> 24));
	         break;
	case 14: test = (((0x000000000000000F & dbshift) ^
	                 ((0x00000000000000F0 & dbshift) >> 4))  |
	                 ((0x000000000000000F & dbshift) ^
			 ((0x0000000000000F00 & dbshift) >> 8))  |
			 ((0x000000000000000F & dbshift) ^
			 ((0x000000000000F000 & dbshift) >> 12)) |
			 ((0x000000000000000F & dbshift) ^
			 ((0x00000000000F0000 & dbshift) >> 16)) |
			 ((0x000000000000000F & dbshift) ^
			 ((0x0000000000F00000 & dbshift) >> 20)) |
			 ((0x000000000000000F & dbshift) ^
			 ((0x000000000F000000 & dbshift) >> 24)))&&
			 ((0x0000000000003FFF & dbshift) ^
			 ((0x000000000FFFC000 & dbshift) >> 14));
	         break;
	case 15: test = (((0x000000000000003F & dbshift) ^
	                 ((0x0000000000000FC0 & dbshift) >> 6))  |
	                 ((0x000000000000003F & dbshift) ^
			 ((0x000000000003F000 & dbshift) >> 12)) |
			 ((0x000000000000003F & dbshift) ^
			 ((0x0000000000FC0000 & dbshift) >> 18)) |
			 ((0x000000000000003F & dbshift) ^
			 ((0x000000003F000000 & dbshift) >> 24)))&&
			(((0x00000000000003FF & dbshift) ^
			 ((0x00000000000FFC00 & dbshift) >> 10)) |
			 ((0x00000000000003FF & dbshift) ^
			 ((0x000000003FF00000 & dbshift) >> 20)));
		 break;
	case 16: test = (0x000000000000FFFF & dbshift) ^
	               ((0x00000000FFFF0000 & dbshift) >> 16);
	         break;
      }
      
      if(test)
      {
        /* now save the word just before dbshift */
        if(i==0)
          lastone = shift?(db[i]<<(64-shift)):0x0;
        else
          lastone = ((unsigned long)db[i-1]>>shift) | (shift?(db[i]<<(64-shift)):0x0);

        dbshift = dbshift & mask;
      
        /* is this pattern a replicate of the previous rep_len chars? */
        lastone = (lastone >> (64-r0)) & mask;
        if((shift>=r0 || i>0) && (dbshift == lastone)) continue;
      
        test = 1; /* now is number of times pattern is actually repeated */
        sav = shift;
      
        nextone = (((unsigned long)db[i+(sav+r0)/64]>>((sav+r0)%64)) | (((sav+r0)%64)?
                  (db[i+(sav+r0)/64+1] << (64-(sav+r0)%64)): 0x0)) & mask;

        while(dbshift == nextone   && 
             ((32*i + sav/2) < (dblen-r0+1))) /* don't overrun db */
        {
	  sav += r0;
	  nextone = (((unsigned long)db[i+(sav+r0)/64]>>((sav+r0)%64)) | (((sav+r0)%64)?
                    (db[i+(sav+r0)/64+1] << (64-(sav+r0)%64)): 0x0)) & mask;
          test++;
        }
      
        if(test >= min_repeats)
        {
	  if(k < arr_size)
          {
	    location[k] = shift/2 + 32*i + 1;
	    (*num_found)++;
	    pattern[k] = test;
	    pattern[k] = ((pattern[k] << 32) | dbshift) << 2;
	    k++;
	  }
	  else
	  {
	    *num_found *= -1;
            return;
	  }
        }
      }
    }
  }

#endif /* 64-bit */
#endif /* little_endian */
}

/*
cb_repeatn(3B)                                           Last changed: 09-17-02

NAME
        cb_repeatn - find short tandem repeats in a nucleotide string

SYNOPSIS

        C/C++:

        #include <cbl.h>
        void cb_repeatn(long *db, long dblen, long repeat_len, long min_repeats,
                long *pattern, long *location, long *num_found);

        Fortran:

        call cb_repeatn( db, dblen, repeat_len, min_repeats, pattern, &
                  location, num_found)
                

IMPLEMENTATION

        Cray SV1 series UNICOS systems

DESCRIPTION

        cb_repeatn finds regions of a nucleotide string that contain
        a short pattern of nucleotides repeated consecutively in the
        region.  For each nucleotide starting location, N, in db
        between the beginning of db and repeat_len*min_repeats locations 
        from the end of db, repeat_len nucleotides in db are replicated 
        min_repeats times and compared to the the original nucleotide 
        sequence. If there is an exact match, the current pattern becomes
        a candidate. The candidate is checked to determine if it itself is 
        a repeat sequence. For example, if repeat_len is 6 and the candidate 
        is ATATAT, the candidate is itself a repeat of a length 2 pattern.
        If the candidate is an overlay of a smaller size repeat, the candidate 
        is rejected. The list of recent patterns is checked to be sure that 
        the candidate is not the tail end of an already recorded pattern. 
        For example, if at position N=100 the sequence is 
        ACTACTACTACTACTACTACTACT, and repeat_len = 3, a pattern with 8 repeats 
        of ACT is recorded. The candidate of 7 repeats starting at location 
        N=103 is rejected.  If the candidate survives these filters, then the
        actual number of times it is repeated is computed. This value, count, 
        is placed in the lower 30 bits of the next free location of the pattern 
        array. The actual repeated nucleotide pattern is placed in bits 30-61 of 
        the same word, and the pattern is left justified in this field.  
        The starting location, N, is placed in the corresponding element of the 
        location array.
                            
        db      input nucleotide data for database string, packed using 
                2-bit compression (see cb_compress, mode=2). In Fortran
                db should be an INTEGER(8) array.
                          
        dblen   input number of nucleotides packed into db. In Fortran
                dblen should be an INTEGER(8) variable, constant, or
                expression.
                          
        repeat_len   input length, in number of nucleotides, of repeat 
                pattern. In Fortran repeat_len should be an INTEGER(8)
                variable, constant, or expression.
              
        min_repeats  input minimum number of times the pattern must be
                repeated to qualify as a find. In Fortran min_repeats
                should be an INTEGER(8) variable, constant, or 
                expression.
                            
        pattern output array of packed words containing the repeated pattern 
                and the number of times it was repeated, in this format:
                            
                       ------------------------------------------
                       |00| pattern (32 bits) | count (30 bits) |               
                       ------------------------------------------

                In Fortran pattern should be an INTEGER(8) array. The
                memory for pattern must be allocated before calling
                cb_repeatn.
                            
        location   output array of starting locations in db for repeated 
                patterns in the corresponding locations of the pattern array. 
                The first location of db is number 1. In Fortran location
                should be an INTEGER(8) array. The memory for location
                must be allocated before calling cb_repeatn.
                            
        num_found  on input - size of pattern and location arrays.
                   on output - number of valid entries in the pattern array. 
                If the number of patterns found is larger than the size of 
                the pattern array, num_found is returned as -(number found 
                before overflow), and the remainder of db is not searched.
                 Three error conditions cause the special value of -1 to be 
                returned. See the NOTES below.
                            
NOTES   

        cb_repeatn assumes the following conditions:

                repeat_len is in the range 2..16

                num_found > 0 on entry

                dblen > repeat_len*min_repeats

        If any of the three conditions is false on entry, num_found
        is set to -1 and no search is done.

        If repeat_len*min_repeats > 32 the internal value for min_repeats
        is reset to 32/repeat_len.
    
        cb_repeatn is single-threaded (i.e. not tasked) and may be called 
        from within a parallel region.

        cb_repeatn replaces the contents of the bmm register.

SEE ALSO

        cb_compress(3B), INTRO_LIBCBL(3B)

*/
