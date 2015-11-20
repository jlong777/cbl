/* cb_revcompl.c                                 http://cbl.sourceforge.net
 *
 * See original man page at the bottom. This code is the generic version
 * with provision for 32 and 64 bit long ints.
 * cb_revcompl works on both big & little endian boxes
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
 * $Id: cb_revcompl.c,v 1.5 2003/09/25 00:53:27 jlong777 Exp $
 */

#include "cb_macro.h"
#include <stdio.h>

void cb_revcompl( long *db, long *dbrc, long dblen, long mode)
{
  int i, length;
  register int lshf, rshf;
  register unsigned long thisone;

#ifdef NOT_VECTOR_BOX
  register unsigned long lastone;
#endif
#ifdef VECTOR_BOX
  unsigned long lastone[VECTLEN];
  register unsigned long mask;
#endif
  
#ifdef LONG32
  
  switch(mode)
  {
    case 2:
      length = (dblen+15)/16; /* length of db in words */
    
      /* algorithm is as follows: find out number of used 2-bit slots in
       * the last word of db so we know how much to shift each previous
       * word of db. reverse the 2-bit ordering, complement, save part that 
       * will get shifted off, shift and | with saved part from last iteration
       */
      
      switch(dblen%16) /* find amount to shift all words */
      {
        case  1: lshf = 30; rshf =  2; break;
	case  2: lshf = 28; rshf =  4; break;
	case  3: lshf = 26; rshf =  6; break;
	case  4: lshf = 24; rshf =  8; break;
	case  5: lshf = 22; rshf = 10; break;
	case  6: lshf = 20; rshf = 12; break;
	case  7: lshf = 18; rshf = 14; break;
	case  8: lshf = 16; rshf = 16; break;
	case  9: lshf = 14; rshf = 18; break;
	case 10: lshf = 12; rshf = 20; break;
	case 11: lshf = 10; rshf = 22; break;
	case 12: lshf =  8; rshf = 24; break;
	case 13: lshf =  6; rshf = 26; break;
	case 14: lshf =  4; rshf = 28; break;
	case 15: lshf =  2; rshf = 30; break;
	default: lshf =  0; rshf = 32; 
      }
      if(dblen%16)
      {
        /* got this reversal idea at http://remus.rutgers.edu/~rhoads/Code/rev.long.c */
        lastone = (((unsigned long)(0xCCCCCCCC & db[length-1])) >> 2) |
	          (                (0x33333333 & db[length-1])  << 2);
        lastone = ((0xF0F0F0F0 & lastone) >> 4) | ((0x0F0F0F0F & lastone) << 4);
        lastone = ((0xFF00FF00 & lastone) >> 8) | ((0x00FF00FF & lastone) << 8);
        lastone = ((lastone >> 16) | (lastone << 16)) ^ 0xAAAAAAAA; /* complement */ 
      
        for(i=1; i<length; i++)
        {
           /* reverse */
	  thisone = (((unsigned long)(0xCCCCCCCC & db[length-1-i])) >> 2) |
	            (                (0x33333333 & db[length-1-i])  << 2);
          thisone = ((0xF0F0F0F0 & thisone) >> 4) | ((0x0F0F0F0F & thisone) << 4);
          thisone = ((0xFF00FF00 & thisone) >> 8) | ((0x00FF00FF & thisone) << 8);
          thisone = ((thisone >> 16) | (thisone << 16)) ^ 0xAAAAAAAA; /* complement */   
#ifdef B_ENDIAN
	  dbrc[i-1] = (lastone << lshf) | (thisone >> rshf);
#endif
#ifdef L_ENDIAN
          dbrc[i-1] = (lastone >> lshf) | (thisone << rshf);
#endif
          lastone = thisone;
        }
#ifdef B_ENDIAN
        dbrc[length-1] = lastone << lshf;
#endif
#ifdef L_ENDIAN
        dbrc[length-1] = lastone >> lshf;
#endif
      }
      else
      {
        for(i=0; i<length; i++)
        {
           /* reverse */
	  thisone = (((unsigned long)(0xCCCCCCCC & db[length-1-i])) >> 2) |
	            (                (0x33333333 & db[length-1-i])  << 2);
          thisone = ( (unsigned long)(0xF0F0F0F0 & thisone) >> 4) | ((0x0F0F0F0F & thisone) << 4);
          thisone = ( (unsigned long)(0xFF00FF00 & thisone) >> 8) | ((0x00FF00FF & thisone) << 8);
          dbrc[i] = (((unsigned long)(thisone) >> 16) | (thisone << 16)) ^ 0xAAAAAAAA; /* complement */ 
        }
      }
      break;
    
    case 4:
      length = (dblen+7)/8; /* length of db in words */
    
      switch(dblen%8) /* amount to shift all words */
      {
        case  1: lshf = 28; rshf =  4; break;
	case  2: lshf = 24; rshf =  8; break;
	case  3: lshf = 20; rshf = 12; break;
	case  4: lshf = 16; rshf = 16; break;
	case  5: lshf = 12; rshf = 20; break;
	case  6: lshf =  8; rshf = 24; break;
	case  7: lshf =  4; rshf = 28; break;
	default: lshf =  0; rshf = 32; 
      }
      if(dblen%8)
      {
        lastone = (((unsigned long)(0xAAAAAAAA & db[length-1])) >> 1) |
	          (                (0x55555555 & db[length-1])  << 1);
        lastone = ((0xCCCCCCCC & lastone) >> 2) | ((0x33333333 & lastone) << 2);
	lastone = ((0xF0F0F0F0 & lastone) >> 4) | ((0x0F0F0F0F & lastone) << 4);
        lastone = ((0xFF00FF00 & lastone) >> 8) | ((0x00FF00FF & lastone) << 8);
        lastone = (lastone >> 16) | (lastone << 16);
      
        for(i=1; i<length; i++)
        {
           /* reverse */
	  thisone  = (((unsigned long)(0xAAAAAAAA & db[length-1-i])) >> 1) |
	             (                (0x55555555 & db[length-1-i])  << 1);
          thisone = ((0xCCCCCCCC & thisone) >> 2) | ((0x33333333 & thisone) << 2);
          thisone = ((0xF0F0F0F0 & thisone) >> 4) | ((0x0F0F0F0F & thisone) << 4);
          thisone = ((0xFF00FF00 & thisone) >> 8) | ((0x00FF00FF & thisone) << 8);
          thisone = (thisone >> 16) | (thisone << 16);
#ifdef B_ENDIAN
	  dbrc[i-1] = (lastone << lshf) | (thisone >> rshf);
#endif
#ifdef L_ENDIAN
	  dbrc[i-1] = (lastone >> lshf) | (thisone << rshf);
#endif
          lastone = thisone;
        }
#ifdef B_ENDIAN
        dbrc[length-1] = lastone << lshf;
#endif
#ifdef L_ENDIAN
        dbrc[length-1] = lastone >> lshf;
#endif
      }
      else
      {
        for(i=0; i<length; i++)
        {
           /* reverse */
	  thisone = (((unsigned long)(0xAAAAAAAA & db[length-1-i])) >> 1) |
	            (                (0x55555555 & db[length-1-i])  << 1);
	  thisone = ((unsigned long)(0xCCCCCCCC & thisone) >>  2) | ((0x33333333 & thisone) << 2);
	  thisone = ((unsigned long)(0xF0F0F0F0 & thisone) >>  4) | ((0x0F0F0F0F & thisone) << 4);
	  thisone = ((unsigned long)(0xFF00FF00 & thisone) >>  8) | ((0x00FF00FF & thisone) << 8);
	  dbrc[i] = ((unsigned long)(thisone) >> 16) | (thisone << 16);
        }
      }
      break;
    
    default:
      fprintf(stderr, "cb_revcompl: Invalid mode parameter %d, returning...\n", mode);
      return;
  }
#endif /* 32-bit */

#ifdef LONG64

  switch(mode)
  {
    case 2:
      length = (dblen+31)/32; /* length of db in words */
      
      switch(dblen%32) /* find amount to shift all words */
      {
        case  1: lshf = 62; rshf =  2; break;
	case  2: lshf = 60; rshf =  4; break;
	case  3: lshf = 58; rshf =  6; break;
	case  4: lshf = 56; rshf =  8; break;
	case  5: lshf = 54; rshf = 10; break;
	case  6: lshf = 52; rshf = 12; break;
	case  7: lshf = 50; rshf = 14; break;
	case  8: lshf = 48; rshf = 16; break;
	case  9: lshf = 46; rshf = 18; break;
	case 10: lshf = 44; rshf = 20; break;
	case 11: lshf = 42; rshf = 22; break;
	case 12: lshf = 40; rshf = 24; break;
	case 13: lshf = 38; rshf = 26; break;
	case 14: lshf = 36; rshf = 28; break;
	case 15: lshf = 34; rshf = 30; break;
	case 16: lshf = 32; rshf = 32; break;
	case 17: lshf = 30; rshf = 34; break;
	case 18: lshf = 28; rshf = 36; break;
	case 19: lshf = 26; rshf = 38; break;
	case 20: lshf = 24; rshf = 40; break;
	case 21: lshf = 22; rshf = 42; break;
	case 22: lshf = 20; rshf = 44; break;
	case 23: lshf = 18; rshf = 46; break;
	case 24: lshf = 16; rshf = 48; break;
	case 25: lshf = 14; rshf = 50; break;
	case 26: lshf = 12; rshf = 52; break;
	case 27: lshf = 10; rshf = 54; break;
	case 28: lshf =  8; rshf = 56; break;
	case 29: lshf =  6; rshf = 58; break;
	case 30: lshf =  4; rshf = 60; break;
	case 31: lshf =  2; rshf = 62; break;
	default: lshf =  0; rshf = 64; 
      }
      if(dblen%32)
      {
#ifdef NOT_VECTOR_BOX
        lastone  = (((unsigned long)(0xCCCCCCCCCCCCCCCC & db[length-1])) >> 2) |
	           (                (0x3333333333333333 & db[length-1])  << 2);
        lastone = ((0xF0F0F0F0F0F0F0F0 & lastone) >>  4) | ((0x0F0F0F0F0F0F0F0F & lastone)  <<  4);
        lastone = ((0xFF00FF00FF00FF00 & lastone) >>  8) | ((0x00FF00FF00FF00FF & lastone) <<  8);
        lastone = ((0xFFFF0000FFFF0000 & lastone) >> 16) | ((0x0000FFFF0000FFFF & lastone)  << 16);
	lastone = ((lastone >> 32) | (lastone << 32)) ^ 0xAAAAAAAAAAAAAAAA;
	
	for(i=1; i<length; i++)
        {
	  thisone   = (((unsigned long)(0xCCCCCCCCCCCCCCCC & db[length-1-i])) >> 2) |
	              (                (0x3333333333333333 & db[length-1-i])  << 2);
          thisone = ((0xF0F0F0F0F0F0F0F0 & thisone) >>  4) | ((0x0F0F0F0F0F0F0F0F & thisone) <<  4);
          thisone = ((0xFF00FF00FF00FF00 & thisone) >>  8) | ((0x00FF00FF00FF00FF & thisone) <<  8);
          thisone = ((0xFFFF0000FFFF0000 & thisone) >> 16) | ((0x0000FFFF0000FFFF & thisone) << 16);
	  thisone = ((thisone >> 32) | (thisone << 32)) ^ 0xAAAAAAAAAAAAAAAA; /* complement */   

#ifdef B_ENDIAN
	  dbrc[i-1] = (lastone << lshf) | (thisone >> rshf);
#endif
#ifdef L_ENDIAN
          dbrc[i-1] = (lastone >> lshf) | (thisone << rshf);
#endif
          lastone = thisone;
        }
	
#ifdef B_ENDIAN
        dbrc[length-1] = lastone << lshf;
#endif
#ifdef L_ENDIAN
        dbrc[length-1] = lastone >> lshf;
#endif

#endif /* not a vector machine */


#ifdef VECTOR_BOX /* big endian only */
#pragma vdir nodep
#pragma _CRI ivdep
	for(i=0; i<VECTLEN; i++)
        {
	  lastone[i] = (((unsigned long)(0xCCCCCCCCCCCCCCCC & db[length-1-i])) >> 2) |
	               (                (0x3333333333333333 & db[length-1-i])  << 2);
	  lastone[i] = ((0xF0F0F0F0F0F0F0F0 & lastone[i]) >>  4) | ((0x0F0F0F0F0F0F0F0F & lastone[i]) <<  4);
          lastone[i] = ((0xFF00FF00FF00FF00 & lastone[i]) >>  8) | ((0x00FF00FF00FF00FF & lastone[i]) <<  8);
          lastone[i] = ((0xFFFF0000FFFF0000 & lastone[i]) >> 16) | ((0x0000FFFF0000FFFF & lastone[i]) << 16);
	  lastone[i] = ((lastone[i] >> 32) | (lastone[i] << 32)) ^ 0xAAAAAAAAAAAAAAAA;
	}
	
	switch(VECTLEN) /* instead of i%VECTLEN, i&mask is faster */
	{
	  case 2:    mask = 0x0000000000000001; break;
	  case 4:    mask = 0x0000000000000003; break;
	  case 8:    mask = 0x0000000000000007; break;
	  case 16:   mask = 0x000000000000000F; break;
	  case 32:   mask = 0x000000000000001F; break;
	  case 64:   mask = 0x000000000000003F; break;
	  case 128:  mask = 0x000000000000007F; break;
	  case 256:  mask = 0x00000000000000FF; break;
	  case 512:  mask = 0x00000000000001FF; break;
	  case 1024: mask = 0x00000000000003FF; break;
	  default:   mask = 0xFFFFFFFFFFFFFFFF; break;
	}
	
#pragma vdir nodep
#pragma _CRI ivdep
	for(i=1; i<length; i++)
        {
	  thisone   = (((unsigned long)(0xCCCCCCCCCCCCCCCC & db[length-1-i])) >> 2) |
	              (                (0x3333333333333333 & db[length-1-i])  << 2);
          thisone = ((0xF0F0F0F0F0F0F0F0 & thisone) >>  4) | ((0x0F0F0F0F0F0F0F0F & thisone) <<  4);
          thisone = ((0xFF00FF00FF00FF00 & thisone) >>  8) | ((0x00FF00FF00FF00FF & thisone) <<  8);
          thisone = ((0xFFFF0000FFFF0000 & thisone) >> 16) | ((0x0000FFFF0000FFFF & thisone) << 16);
	  thisone = ((thisone >> 32) | (thisone << 32)) ^ 0xAAAAAAAAAAAAAAAA; /* complement */   

	  dbrc[i-1] = (lastone[(i-1)&mask] << lshf) | (thisone >> rshf);
	  
          lastone[i&mask] = thisone;
        }
	
        dbrc[length-1] = lastone[(i-1)&mask] << lshf;
	  
#endif /* a vector machine */
      }
      else
      {
#pragma vdir nodep
#pragma _CRI ivdep
        for(i=0; i<length; i++)
        {
           /* reverse */
	  thisone = (((unsigned long)(0xCCCCCCCCCCCCCCCC & db[length-1-i])) >> 2) |
	            (                (0x3333333333333333 & db[length-1-i])  << 2);
          thisone = ( (unsigned long)(0xF0F0F0F0F0F0F0F0 & thisone) >>  4) | ((0x0F0F0F0F0F0F0F0F & thisone) <<  4);
          thisone = ( (unsigned long)(0xFF00FF00FF00FF00 & thisone) >>  8) | ((0x00FF00FF00FF00FF & thisone) <<  8);
	  thisone = ( (unsigned long)(0xFFFF0000FFFF0000 & thisone) >> 16) | ((0x0000FFFF0000FFFF & thisone) << 16);
          dbrc[i] = (((unsigned long)(thisone) >> 32) | (thisone << 32)) ^ 0xAAAAAAAAAAAAAAAA; /* complement */ 
        }
      }
      break;
    
    case 4:
      length = (dblen+15)/16; /* length of db in words */
    
      switch(dblen%16) /* amount to shift all words */
      {
        case  1: lshf = 60; rshf =  4; break;
	case  2: lshf = 56; rshf =  8; break;
	case  3: lshf = 52; rshf = 12; break;
	case  4: lshf = 48; rshf = 16; break;
	case  5: lshf = 44; rshf = 20; break;
	case  6: lshf = 40; rshf = 24; break;
	case  7: lshf = 36; rshf = 28; break;
	case  8: lshf = 32; rshf = 32; break;
	case  9: lshf = 28; rshf = 36; break;
	case 10: lshf = 24; rshf = 40; break;
	case 11: lshf = 20; rshf = 44; break;
	case 12: lshf = 16; rshf = 48; break;
	case 13: lshf = 12; rshf = 52; break;
	case 14: lshf =  8; rshf = 56; break;
	case 15: lshf =  4; rshf = 60; break;
	default: lshf =  0; rshf = 64; 
      }
      if(dblen%16)
      {
#ifdef NOT_VECTOR_BOX
	lastone = (((unsigned long)(0xAAAAAAAAAAAAAAAA & db[length-1])) >> 1) |
	          (                (0x5555555555555555 & db[length-1])  << 1);
        lastone = ((0xCCCCCCCCCCCCCCCC & lastone) >>  2) | ((0x3333333333333333 & lastone) <<  2);
	lastone = ((0xF0F0F0F0F0F0F0F0 & lastone) >>  4) | ((0x0F0F0F0F0F0F0F0F & lastone) <<  4);
        lastone = ((0xFF00FF00FF00FF00 & lastone) >>  8) | ((0x00FF00FF00FF00FF & lastone) <<  8);
        lastone = ((0xFFFF0000FFFF0000 & lastone) >> 16) | ((0x0000FFFF0000FFFF & lastone) << 16);
	lastone = (lastone >> 32) | (lastone << 32);

        for(i=1; i<length; i++)
        {
	  thisone = (((unsigned long)(0xAAAAAAAAAAAAAAAA & db[length-1-i])) >> 1) |
	            (                (0x5555555555555555 & db[length-1-i])  << 1);
          thisone = ((0xCCCCCCCCCCCCCCCC & thisone) >>  2) | ((0x3333333333333333 & thisone) <<  2);
          thisone = ((0xF0F0F0F0F0F0F0F0 & thisone) >>  4) | ((0x0F0F0F0F0F0F0F0F & thisone) <<  4);
          thisone = ((0xFF00FF00FF00FF00 & thisone) >>  8) | ((0x00FF00FF00FF00FF & thisone) <<  8);
          thisone = ((0xFFFF0000FFFF0000 & thisone) >> 16) | ((0x0000FFFF0000FFFF & thisone) << 16);
	  thisone = (thisone >> 32) | (thisone << 32);
#ifdef B_ENDIAN
	  dbrc[i-1] = (lastone << lshf) | (thisone >> rshf);
#endif
#ifdef L_ENDIAN
	  dbrc[i-1] = (lastone >> lshf) | (thisone << rshf);
#endif
          lastone = thisone;
        }
#ifdef B_ENDIAN
        dbrc[length-1] = lastone << lshf;
#endif
#ifdef L_ENDIAN
        dbrc[length-1] = lastone >> lshf;
#endif

#endif /* not a vector machine */


#ifdef VECTOR_BOX /* big endian only */
#pragma vdir nodep
#pragma _CRI ivdep
        for(i=0; i<VECTLEN; i++)
        {
	  lastone[i] = (((unsigned long)(0xAAAAAAAAAAAAAAAA & db[length-1-i])) >> 1) |
	               (                (0x5555555555555555 & db[length-1-i])  << 1);
          lastone[i] = ((0xCCCCCCCCCCCCCCCC & lastone[i]) >>  2) | ((0x3333333333333333 & lastone[i]) <<  2);
	  lastone[i] = ((0xF0F0F0F0F0F0F0F0 & lastone[i]) >>  4) | ((0x0F0F0F0F0F0F0F0F & lastone[i]) <<  4);
          lastone[i] = ((0xFF00FF00FF00FF00 & lastone[i]) >>  8) | ((0x00FF00FF00FF00FF & lastone[i]) <<  8);
          lastone[i] = ((0xFFFF0000FFFF0000 & lastone[i]) >> 16) | ((0x0000FFFF0000FFFF & lastone[i]) << 16);
	  lastone[i] = (lastone[i] >> 32) | (lastone[i] << 32);
	}
	
	switch(VECTLEN) /* instead of i%VECTLEN, i&mask is faster */
	{
	  case 2:    mask = 0x0000000000000001; break;
	  case 4:    mask = 0x0000000000000003; break;
	  case 8:    mask = 0x0000000000000007; break;
	  case 16:   mask = 0x000000000000000F; break;
	  case 32:   mask = 0x000000000000001F; break;
	  case 64:   mask = 0x000000000000003F; break;
	  case 128:  mask = 0x000000000000007F; break;
	  case 256:  mask = 0x00000000000000FF; break;
	  case 512:  mask = 0x00000000000001FF; break;
	  case 1024: mask = 0x00000000000003FF; break;
	  default:   mask = 0xFFFFFFFFFFFFFFFF; break;
	}
	
#pragma vdir nodep
#pragma _CRI ivdep
	for(i=1; i<length; i++)
        {
	  thisone = (((unsigned long)(0xAAAAAAAAAAAAAAAA & db[length-1-i])) >> 1) |
	            (                (0x5555555555555555 & db[length-1-i])  << 1);
          thisone = ((0xCCCCCCCCCCCCCCCC & thisone) >>  2) | ((0x3333333333333333 & thisone) <<  2);
          thisone = ((0xF0F0F0F0F0F0F0F0 & thisone) >>  4) | ((0x0F0F0F0F0F0F0F0F & thisone) <<  4);
          thisone = ((0xFF00FF00FF00FF00 & thisone) >>  8) | ((0x00FF00FF00FF00FF & thisone) <<  8);
          thisone = ((0xFFFF0000FFFF0000 & thisone) >> 16) | ((0x0000FFFF0000FFFF & thisone) << 16);
	  thisone = (thisone >> 32) | (thisone << 32);
	  
	  dbrc[i-1] = (lastone[(i-1)&mask] << lshf) | (thisone >> rshf);
	  
          lastone[i&mask] = thisone;
        }
	
        dbrc[length-1] = lastone[(i-1)&mask] << lshf;
#endif
      }
      else
      {
#pragma vdir nodep
#pragma _CRI ivdep
        for(i=0; i<length; i++)
        {
           /* reverse */
	  thisone = (((unsigned long)(0xAAAAAAAAAAAAAAAA & db[length-1-i])) >> 1) |
	            (                (0x5555555555555555 & db[length-1-i])  << 1);
	  thisone = ( (unsigned long)(0xCCCCCCCCCCCCCCCC & thisone) >>  2) | ((0x3333333333333333 & thisone) <<  2);
	  thisone = ( (unsigned long)(0xF0F0F0F0F0F0F0F0 & thisone) >>  4) | ((0x0F0F0F0F0F0F0F0F & thisone) <<  4);
	  thisone = ( (unsigned long)(0xFF00FF00FF00FF00 & thisone) >>  8) | ((0x00FF00FF00FF00FF & thisone) <<  8);
	  thisone = ( (unsigned long)(0xFFFF0000FFFF0000 & thisone) >> 16) | ((0x0000FFFF0000FFFF & thisone) << 16);
	  dbrc[i] = ( (unsigned long)(thisone) >> 32) | (thisone << 32);
        }
      }
      break;
    
    default:
      fprintf(stderr, "cb_revcompl: Invalid mode parameter %d, returning...\n", mode);
      return;
  }
#endif /* 64-bit */
}


/*
cb_revcompl(3B)                                           Last changed: 09-17-02

NAME
        cb_revcompl - Reverse complements compressed nucleotide data.

SYNOPSIS

        C/C++:

        #include <cbl.h>
        void cb_revcompl( long *db, long *dbrc, long dblen, long mode);

        Fortran:

        use cb_revc
        call cb_revcompl( db, dbrc, dblen, mode)


IMPLEMENTATION

        Cray SV1 series UNICOS systems

DESCRIPTION

        cb_revcompl reverses the order of the nucleotides in db and
        complements each nucleotide, and stores the result in dbrc.
        The input data is assumed to be stored in compressed form. See
        the man page for cb_compress(3B) for details on nucleotide
        compression. The complement operation changes each nucleotide 
        as follows:
                    
                A -> T
                C -> G
                T -> A
                G -> C               
                    
        For 4-bit encoded data, additional letters are 
        transformed according to these same rules. For 
        example,
                    
                R (A or G) -> V (T or C)
                    

        db      (input) nucleotide sequence in either 2-bit or 4-bit 
                compressed format. See cb_compress(3B) for details on
                compression. Length of db, in 64-bit words is
                (dblen+31)/32 for mode=2, and (dblen+15)/16 for mode=4.
                In Fortran, db should be an INTEGER(8) array.   
                    
        dbrc    (output) nucleotide sequence in same compression
                format as db.  The result in dbrc is the data from
                db in reverse order, and also each nucleotide is
                complemented.  The data in dbrc is shifted so that
                the initial character is left justified in the first
                word of dbrc.  Length of dbrc is the same as length
                of db, and dbrc must be allocated by the caller.
                In Fortran, dbrc should be an INTEGER(8) array.
     
        dblen   (input) number of nucleotide codes packed into db.
                In Fortran, db should be an INTEGER(8) variable,
                constant, or expression.
                    
        mode    (input) allowed values are 2 and 4, and indicates the
                number of bits used to represent each nucleotide
                in the db and dbrc arrays. See the cb_compress(3B) 
                man page for a detailed explanation of the compression
                mappings for mode values 2 and 4. In Fortran, mode should
                be an INTEGER(8) variable, constant, or expression.

NOTES       
        cb_revcompl is single-threaded (i.e. not tasked) and may be called 
        from within a parallel region.

        cb_revcompl replaces the contents of the bmm register.

SEE ALSO

        cb_compress(3B), INTRO_LIBCBL(3B)

*/
