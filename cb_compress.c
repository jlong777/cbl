/* cb_compress.c                                 http://cbl.sourceforge.net
 *
 * See original man page at the bottom, also see two articles I wrote for the
 * ARSC newsletter, issues 255 & 256 at http://www.arsc.edu/pubs/HPCnews.shtml.
 * This code is the generic version with provision for 32 and 64 bit long ints.
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
 * $Id: cb_compress.c,v 1.4 2003/09/25 00:53:27 jlong777 Exp $
 */

#include "cb_macro.h"
#include <stdio.h>
#include <stdlib.h>

void cb_compress(long *db, long *dbc, long dblen, long mode)
{
  long i, length;
  
  static unsigned long bit4[27]; /* 4-bit compression lookup table */
  
  bit4[0]  = 0x00; /* null */
  bit4[1]  = 0x08; /* A, a */
  bit4[2]  = 0x07; /* B, b */
  bit4[3]  = 0x04; /* C, c */
  bit4[4]  = 0x0B; /* D, d */
  bit4[5]  = 0x0F; /* E, e */
  bit4[6]  = 0x0F; /* F, f */
  bit4[7]  = 0x02; /* G, g */
  bit4[8]  = 0x0D; /* H, h */
  bit4[9]  = 0x0F; /* I, i */
  bit4[10] = 0x0F; /* J, j */
  bit4[11] = 0x03; /* K, k */
  bit4[12] = 0x0F; /* L, l */
  bit4[13] = 0x0C; /* M, m */
  bit4[14] = 0x0F; /* N, n */
  bit4[15] = 0x0F; /* O, o */
  bit4[16] = 0x0F; /* P, p */
  bit4[17] = 0x0F; /* Q, q */
  bit4[18] = 0x0A; /* R, r */
  bit4[19] = 0x06; /* S, s */
  bit4[20] = 0x01; /* T, t */
  bit4[21] = 0x01; /* U, u */
  bit4[22] = 0x0E; /* V, v */
  bit4[23] = 0x09; /* W, w */
  bit4[24] = 0x0F; /* X, x */
  bit4[25] = 0x05; /* Y, y */
  bit4[26] = 0x0F; /* Z, z */
	  

#ifdef LONG32
  length = (dblen+3)/4;
  switch (mode)
  {
#ifdef B_ENDIAN
    case 2:
      for(i=0; i<length-4; i+=4)
      {
        dbc[i/4] = ((unsigned long)(0x06000000 & db[i])   <<  5) |
	           ((unsigned long)(0x00060000 & db[i])   << 11) |
	           ((unsigned long)(0x00000600 & db[i])   << 17) |
	           ((unsigned long)(0x00000006 & db[i])   << 23) |
	           ((unsigned long)(0x06000000 & db[i+1]) >>  3) |
	           ((unsigned long)(0x00060000 & db[i+1]) <<  3) |
	           ((unsigned long)(0x00000600 & db[i+1]) <<  9) |
	           ((unsigned long)(0x00000006 & db[i+1]) << 15) |
                   ((unsigned long)(0x06000000 & db[i+2]) >> 11) |
                   ((unsigned long)(0x00060000 & db[i+2]) >>  5) |
		   ((unsigned long)(0x00000600 & db[i+2]) <<  1) |
		   ((unsigned long)(0x00000006 & db[i+2]) <<  7) |
		   ((unsigned long)(0x06000000 & db[i+3]) >> 19) |
		   ((unsigned long)(0x00060000 & db[i+3]) >> 13) |
		   ((unsigned long)(0x00000600 & db[i+3]) >>  7) |
		   ((unsigned long)(0x00000006 & db[i+3]) >>  1);
      }
      dbc[i/4] = ((4*i   >= dblen)? 0x00 :         /* out of chars? */
                 (  (unsigned long)(0x06000000 & db[i])   <<  5))  |
		 ((4*i+1 >= dblen)? 0x00 :
	         (  (unsigned long)(0x00060000 & db[i])   << 11))  |
		 ((4*i+2 >= dblen)? 0x00 :
	         (  (unsigned long)(0x00000600 & db[i])   << 17))  |
		 ((4*i+3 >= dblen)? 0x00 :
	         (  (unsigned long)(0x00000006 & db[i])   << 23))  |
	         ((4*i+4 >= dblen)? 0x00 :
		 (  (unsigned long)(0x06000000 & db[i+1]) >>  3))  |
		 ((4*i+5 >= dblen)? 0x00 :
	         (  (unsigned long)(0x00060000 & db[i+1]) <<  3))  |
		 ((4*i+6 >= dblen)? 0x00 :
	         (  (unsigned long)(0x00000600 & db[i+1]) <<  9))  |
		 ((4*i+7 >= dblen)? 0x00 :
	         (  (unsigned long)(0x00000006 & db[i+1]) << 15))  |
                 ((4*i+8 >= dblen)? 0x00 :
		 (  (unsigned long)(0x06000000 & db[i+2]) >> 11))  |
		 ((4*i+9 >= dblen)? 0x00 :
                 (  (unsigned long)(0x00060000 & db[i+2]) >>  5))  |
		 ((4*i+10>= dblen)? 0x00 :
		 (  (unsigned long)(0x00000600 & db[i+2]) <<  1))  |
		 ((4*i+11>= dblen)? 0x00 :
		 (  (unsigned long)(0x00000006 & db[i+2]) <<  7))  |
		 ((4*i+12>= dblen)? 0x00 :
		 (  (unsigned long)(0x06000000 & db[i+3]) >> 19))  |
		 ((4*i+13>= dblen)? 0x00 :
		 (  (unsigned long)(0x00060000 & db[i+3]) >> 13))  |
		 ((4*i+14>= dblen)? 0x00 :
		 (  (unsigned long)(0x00000600 & db[i+3]) >>  7))  |
		 ((4*i+15>= dblen)? 0x00 :
		 (  (unsigned long)(0x00000006 & db[i+3]) >>  1));
      break;
	  
    case 4:
      for(i=0; i<length-2; i+=2)
      {
        dbc[i/2] = (bit4[(unsigned long)(0x1F000000 & db[i]) >> 24] << 28) |
	           (bit4[(unsigned long)(0x001F0000 & db[i]) >> 16] << 24) |
	           (bit4[(unsigned long)(0x00001F00 & db[i]) >>  8] << 20) |
	           (bit4[                0x0000001F & db[i]]        << 16) |
	           (bit4[(unsigned long)(0x1F000000 & db[i+1])>>24] << 12) |
	           (bit4[(unsigned long)(0x001F0000 & db[i+1])>>16] <<  8) |
	           (bit4[(unsigned long)(0x00001F00 & db[i+1])>> 8] <<  4) |
	           (bit4[                0x0000001F & db[i+1]]);
      }
      dbc[i/2] = ((4*i  >= dblen)?     0x00 :      /* out of chars? */
                 (bit4[(unsigned long)(0x1F000000 & db[i]) >> 24] << 28))  |
		 ((4*i+1>= dblen)?     0x00 :
	         (bit4[(unsigned long)(0x001F0000 & db[i]) >> 16] << 24))  |
		 ((4*i+2>= dblen)?     0x00 :
	         (bit4[(unsigned long)(0x00001F00 & db[i]) >>  8] << 20))  |
		 ((4*i+3>= dblen)?     0x00 :
	         (bit4[                0x0000001F & db[i]]        << 16))  |
		 ((4*i+4>= dblen)?     0x00 :
	         (bit4[(unsigned long)(0x1F000000 & db[i+1])>>24] << 12))  |
		 ((4*i+5>= dblen)?     0x00 :
	         (bit4[(unsigned long)(0x001F0000 & db[i+1])>>16] <<  8))  |
		 ((4*i+6>= dblen)?     0x00 :
	         (bit4[(unsigned long)(0x00001F00 & db[i+1])>> 8] <<  4))  |
		 ((4*i+7>= dblen)?     0x00 :
	          bit4[                0x0000001F & db[i+1]]);
      break;
	  
    case 5:
      for(i=0; i<length-3; i+=3)
      {
        dbc[(2*i)/3]   = ((unsigned long)(0x1F000000 & db[i])   <<  1) |
	                 ((unsigned long)(0x001F0000 & db[i])   <<  4) |
	                 ((unsigned long)(0x00001F00 & db[i])   <<  7) |
	                 ((unsigned long)(0x0000001F & db[i])   << 10) |
	                 ((unsigned long)(0x1F000000 & db[i+1]) >> 19) |
	                 ((unsigned long)(0x001F0000 & db[i+1]) >> 16);
		     
        dbc[(2*i)/3+1] = ((unsigned long)(0x00001F00 & db[i+1]) << 17) |
	                 ((unsigned long)(0x0000001F & db[i+1]) << 20) |
		         ((unsigned long)(0x1F000000 & db[i+2]) >>  9) |
	                 ((unsigned long)(0x001F0000 & db[i+2]) >>  6) |
	                 ((unsigned long)(0x00001F00 & db[i+2]) >>  3) |
	                                 (0x0000001F & db[i+2]);
      }
      dbc[(2*i)/3]   = ((4*i   >= dblen)? 0x00 :      /* out of chars? */
                       (  (unsigned long)(0x1F000000 & db[i])   <<  1))  |
		       ((4*i+1 >= dblen)? 0x00 :
	               (  (unsigned long)(0x001F0000 & db[i])   <<  4))  |
		       ((4*i+2 >= dblen)? 0x00 :
	               (  (unsigned long)(0x00001F00 & db[i])   <<  7))  |
		       ((4*i+3 >= dblen)? 0x00 :
	               (  (unsigned long)(0x0000001F & db[i])   << 10))  |
		       ((4*i+4 >= dblen)? 0x00 :
	               (  (unsigned long)(0x1F000000 & db[i+1]) >> 19))  |
		       ((4*i+5 >= dblen)? 0x00 :
	               (  (unsigned long)(0x001F0000 & db[i+1]) >> 16));
		     
      dbc[(2*i)/3+1] = ((4*i+6 >= dblen)? 0x00 :
                       (  (unsigned long)(0x00001F00 & db[i+1]) << 17))  |
		       ((4*i+7 >= dblen)? 0x00 :
	               (  (unsigned long)(0x0000001F & db[i+1]) << 20))  |
		       ((4*i+8 >= dblen)? 0x00 :
		       (  (unsigned long)(0x1F000000 & db[i+2]) >>  9))  |
		       ((4*i+9 >= dblen)? 0x00 :
	               (  (unsigned long)(0x001F0000 & db[i+2]) >>  6))  |
		       ((4*i+10>= dblen)? 0x00 :
	               (  (unsigned long)(0x00001F00 & db[i+2]) >>  3))  |
		       ((4*i+11>= dblen)? 0x00 :
	                                 (0x0000001F & db[i+2]));

      break;

    case 6:
      for(i=0; i<length-2; i+=2)
      {
        dbc[i/2] = ((0x80000000 &  db[i])? 0x00 : /*  is it an A=1xxx ? */
	           (((unsigned long)(0x40000000 & ~db[i])  << 1)   |
		    ((unsigned long)(0x10000000 & ~db[i])  << 2))) |
	           ((0x08000000 &  db[i])? 0x00 :
		   (((unsigned long)(0x04000000 & ~db[i])  << 3)   |
		    ((unsigned long)(0x01000000 & ~db[i])  << 4))) |
	           ((0x00800000 &  db[i])? 0x00 :
		   (((unsigned long)(0x00400000 & ~db[i])  << 5)   |
		    ((unsigned long)(0x00100000 & ~db[i])  << 6))) |
	           ((0x00080000 &  db[i])? 0x00 :
		   (((unsigned long)(0x00040000 & ~db[i])  << 7)   |
		    ((unsigned long)(0x00010000 & ~db[i])  << 8))) |
	           ((0x00008000 &  db[i])? 0x00 :
		   (((unsigned long)(0x00004000 & ~db[i])  << 9)   |
		    ((unsigned long)(0x00001000 & ~db[i])  <<10))) |
	           ((0x00000800 &  db[i])? 0x00 :
		   (((unsigned long)(0x00000400 & ~db[i])  <<11)   |
		    ((unsigned long)(0x00000100 & ~db[i])  <<12))) |
	           ((0x00000080 &  db[i])? 0x00 :
		   (((unsigned long)(0x00000040 & ~db[i])  <<13)   |
		    ((unsigned long)(0x00000010 & ~db[i])  <<14))) |
	           ((0x00000008 &  db[i])? 0x00 :
		   (((unsigned long)(0x00000004 & ~db[i])  <<15)   |
		    ((unsigned long)(0x00000001 & ~db[i])  <<16))) |
		   ((0x80000000 &  db[i+1])? 0x00 :
	           (((unsigned long)(0x40000000 & ~db[i+1])>>15)   |
		    ((unsigned long)(0x10000000 & ~db[i+1])>>14))) |
	           ((0x08000000 &  db[i+1])? 0x00 :
		   (((unsigned long)(0x04000000 & ~db[i+1])>>13)   |
		    ((unsigned long)(0x01000000 & ~db[i+1])>>12))) |
	           ((0x00800000 &  db[i+1])? 0x00 :
		   (((unsigned long)(0x00400000 & ~db[i+1])>>11)   |
		    ((unsigned long)(0x00100000 & ~db[i+1])>>10))) |
	           ((0x00080000 &  db[i+1])? 0x00 :
		   (((unsigned long)(0x00040000 & ~db[i+1])>> 9)   |
		    ((unsigned long)(0x00010000 & ~db[i+1])>> 8))) |
	           ((0x00008000 &  db[i+1])? 0x00 :
		   (((unsigned long)(0x00004000 & ~db[i+1])>> 7)   |
		    ((unsigned long)(0x00001000 & ~db[i+1])>> 6))) |
	           ((0x00000800 &  db[i+1])? 0x00 :
		   (((unsigned long)(0x00000400 & ~db[i+1])>> 5)   |
		    ((unsigned long)(0x00000100 & ~db[i+1])>> 4))) |
	           ((0x00000080 &  db[i+1])? 0x00 :
		   (((unsigned long)(0x00000040 & ~db[i+1])>> 3)   |
		    ((unsigned long)(0x00000010 & ~db[i+1])>> 2))) |
	           ((0x00000008 &  db[i+1])? 0x00 :
		   (((unsigned long)(0x00000004 & ~db[i+1])>> 1)   |
		    ((unsigned long)(0x00000001 & ~db[i+1]))));
      }
      dbc[i/2] = ((4*i   >= dblen)? 0x00 :     /* out of chars? */
                 (  (unsigned long)(0x80000000 &  db[i])?   0x00 :/* is it an A=1xxx ? */
	         (( (unsigned long)(0x40000000 & ~db[i])   <<  1)      |
		   ((unsigned long)(0x10000000 & ~db[i])   <<  2))))   |
		 ((4*i+1 >= dblen)? 0x00 :
	         (  (unsigned long)(0x08000000 &  db[i])?   0x00 :
		 (( (unsigned long)(0x04000000 & ~db[i])   <<  3)      |
		   ((unsigned long)(0x01000000 & ~db[i])   <<  4))))   |
		 ((4*i+2 >= dblen)? 0x00 :
	         (  (unsigned long)(0x00800000 & db[i])?    0x00 :
		 (( (unsigned long)(0x00400000 & ~db[i])   <<  5)      |
		   ((unsigned long)(0x00100000 & ~db[i])   <<  6))))   |
		 ((4*i+3 >= dblen)? 0x00 :
	         (  (unsigned long)(0x00080000 &  db[i])?   0x00 :
		 (( (unsigned long)(0x00040000 & ~db[i])   <<  7)      |
		   ((unsigned long)(0x00010000 & ~db[i])   <<  8))))   |
		 ((4*i+4 >= dblen)? 0x00 :
	         (  (unsigned long)(0x00008000 &  db[i])?   0x00 :
		 (( (unsigned long)(0x00004000 & ~db[i])   <<  9)      |
		   ((unsigned long)(0x00001000 & ~db[i])   << 10))))   |
		 ((4*i+5 >= dblen)? 0x00 :
	         (  (unsigned long)(0x00000800 &  db[i])?   0x00 :
		 (( (unsigned long)(0x00000400 & ~db[i])   << 11)      |
		   ((unsigned long)(0x00000100 & ~db[i])   << 12))))   |
		 ((4*i+6 >= dblen)? 0x00 :
	         (  (unsigned long)(0x00000080 &  db[i])?   0x00 :
		 (( (unsigned long)(0x00000040 & ~db[i])   << 13)      |
		   ((unsigned long)(0x00000010 & ~db[i])   << 14))))   |
		 ((4*i+7 >= dblen)? 0x00 :
	         (  (unsigned long)(0x00000008 &  db[i])?   0x00 :
		 (( (unsigned long)(0x00000004 & ~db[i])   << 15)      |
		   ((unsigned long)(0x00000001 & ~db[i])   << 16))))   | 
		 ((4*i+8 >= dblen)? 0x00 :
		 (  (unsigned long)(0x80000000 &  db[i+1])? 0x00 :
	         (( (unsigned long)(0x40000000 & ~db[i+1]) >> 15)      |
		   ((unsigned long)(0x10000000 & ~db[i+1]) >> 14))))   |
		 ((4*i+9 >= dblen)? 0x00 :
	         (  (unsigned long)(0x08000000 &  db[i+1])? 0x00 :
		 (( (unsigned long)(0x04000000 & ~db[i+1]) >> 13)      |
		   ((unsigned long)(0x01000000 & ~db[i+1]) >> 12))))   |
		 ((4*i+10>= dblen)? 0x00 :
	         (  (unsigned long)(0x00800000 &  db[i+1])? 0x00 :
		 (( (unsigned long)(0x00400000 & ~db[i+1]) >> 11)      |
		   ((unsigned long)(0x00100000 & ~db[i+1]) >> 10))))   |
		 ((4*i+11>= dblen)? 0x00 :
	         (  (unsigned long)(0x00080000 &  db[i+1])? 0x00 :
		 (( (unsigned long)(0x00040000 & ~db[i+1]) >> 9)       |
		   ((unsigned long)(0x00010000 & ~db[i+1]) >> 8))))    |
		 ((4*i+12>= dblen)? 0x00 :
	         (  (unsigned long)(0x00008000 &  db[i+1])? 0x00 :
		 (( (unsigned long)(0x00004000 & ~db[i+1]) >> 7)       |
		   ((unsigned long)(0x00001000 & ~db[i+1]) >> 6))))    |
		 ((4*i+13>= dblen)? 0x00 :
	         (  (unsigned long)(0x00000800 &  db[i+1])? 0x00 :
		 (( (unsigned long)(0x00000400 & ~db[i+1]) >> 5)       |
		   ((unsigned long)(0x00000100 & ~db[i+1]) >> 4))))    |
		 ((4*i+14>= dblen)? 0x00 :
	         (  (unsigned long)(0x00000080 &  db[i+1])? 0x00 :
		 (( (unsigned long)(0x00000040 & ~db[i+1]) >> 3)       |
		   ((unsigned long)(0x00000010 & ~db[i+1]) >> 2))))    |
		 ((4*i+15>= dblen)? 0x00 :
	         (  (unsigned long)(0x00000008 &  db[i+1])? 0x00 :
		 (( (unsigned long)(0x00000004 & ~db[i+1]) >> 1)       |
		   ((unsigned long)(0x00000001 & ~db[i+1])))));
      break;
#endif /* big endian */
#ifdef L_ENDIAN
    case 2:
      for(i=0; i<length-4; i+=4)
      {
        dbc[i/4] = ((unsigned long)(0x00000006 & db[i])   >>  1) |
	           ((unsigned long)(0x00000600 & db[i])   >>  7) |
	           ((unsigned long)(0x00060000 & db[i])   >> 13) |
	           ((unsigned long)(0x06000000 & db[i])   >> 19) |
	           ((unsigned long)(0x00000006 & db[i+1]) <<  7) |
	           ((unsigned long)(0x00000600 & db[i+1]) <<  1) |
	           ((unsigned long)(0x00060000 & db[i+1]) >>  5) |
	           ((unsigned long)(0x06000000 & db[i+1]) >> 11) |
                   ((unsigned long)(0x00000006 & db[i+2]) << 15) |
                   ((unsigned long)(0x00000600 & db[i+2]) <<  9) |
		   ((unsigned long)(0x00060000 & db[i+2]) <<  3) |
		   ((unsigned long)(0x06000000 & db[i+2]) >>  3) |
		   ((unsigned long)(0x00000006 & db[i+3]) << 23) |
		   ((unsigned long)(0x00000600 & db[i+3]) << 17) |
		   ((unsigned long)(0x00060000 & db[i+3]) << 11) |
		   ((unsigned long)(0x06000000 & db[i+3]) <<  5);
      }
      dbc[i/4] = ((4*i   >= dblen)? 0x00 :         /* out of chars? */
                 (  (unsigned long)(0x00000006 & db[i])   >>  1)) |
		 ((4*i+1 >= dblen)? 0x00 :
		 (  (unsigned long)(0x00000600 & db[i])   >>  7)) |
		 ((4*i+2 >= dblen)? 0x00 :
	         (  (unsigned long)(0x00060000 & db[i])   >> 13)) |
		 ((4*i+3 >= dblen)? 0x00 :
	         (  (unsigned long)(0x06000000 & db[i])   >> 19)) |
	         ((4*i+4 >= dblen)? 0x00 :
	         (  (unsigned long)(0x00000006 & db[i+1]) <<  7)) |
		 ((4*i+5 >= dblen)? 0x00 :
	         (  (unsigned long)(0x00000600 & db[i+1]) <<  1)) |
		 ((4*i+6 >= dblen)? 0x00 :
	         (  (unsigned long)(0x00060000 & db[i+1]) >>  5)) |
		 ((4*i+7 >= dblen)? 0x00 :
	         (  (unsigned long)(0x06000000 & db[i+1]) >> 11)) |
                 ((4*i+8 >= dblen)? 0x00 :
                 (  (unsigned long)(0x00000006 & db[i+2]) << 15)) |
		 ((4*i+9 >= dblen)? 0x00 :
                 (  (unsigned long)(0x00000600 & db[i+2]) <<  9)) |
		 ((4*i+10>= dblen)? 0x00 :
		 (  (unsigned long)(0x00060000 & db[i+2]) <<  3)) |
		 ((4*i+11>= dblen)? 0x00 :
		 (  (unsigned long)(0x06000000 & db[i+2]) >>  3)) |
		 ((4*i+12>= dblen)? 0x00 :
		 (  (unsigned long)(0x00000006 & db[i+3]) << 23)) |
		 ((4*i+13>= dblen)? 0x00 :
		 (  (unsigned long)(0x00000600 & db[i+3]) << 17)) |
		 ((4*i+14>= dblen)? 0x00 :
		 (  (unsigned long)(0x00060000 & db[i+3]) << 11)) |
		 ((4*i+15>= dblen)? 0x00 :
		 (  (unsigned long)(0x06000000 & db[i+3]) <<  5));
      break;
		 	  
    case 4:
      for(i=0; i<length-2; i+=2)
      {
        dbc[i/2] = (bit4[                0x0000001F & db[i]])              |
	           (bit4[(unsigned long)(0x00001F00 & db[i]) >>  8] <<  4) |
	           (bit4[(unsigned long)(0x001F0000 & db[i]) >> 16] <<  8) |
	           (bit4[(unsigned long)(0x1F000000 & db[i]) >> 24] << 12) |
	           (bit4[                0x0000001F & db[i+1]]      << 16) |
	           (bit4[(unsigned long)(0x00001F00 & db[i+1])>> 8] << 20) |
	           (bit4[(unsigned long)(0x001F0000 & db[i+1])>>16] << 24) |
	           (bit4[(unsigned long)(0x1F000000 & db[i+1])>>24] << 28);
      }
      dbc[i/2] = ((4*i   >= dblen)?    0x00 :       /* out of chars? */
                 (bit4[                0x0000001F & db[i]]))              |
		 ((4*i+1 >= dblen)?    0x00 :
	         (bit4[(unsigned long)(0x00001F00 & db[i]) >>  8] <<  4)) |
		 ((4*i+2 >= dblen)?    0x00 :
	         (bit4[(unsigned long)(0x001F0000 & db[i]) >> 16] <<  8)) |
		 ((4*i+3 >= dblen)?    0x00 :
	         (bit4[(unsigned long)(0x1F000000 & db[i]) >> 24] << 12)) |
		 ((4*i+4 >= dblen)?    0x00 :
	         (bit4[                0x0000001F & db[i+1]]      << 16)) |
		 ((4*i+5 >= dblen)?    0x00 :
	         (bit4[(unsigned long)(0x00001F00 & db[i+1])>> 8] << 20)) |
		 ((4*i+6 >= dblen)?    0x00 :
	         (bit4[(unsigned long)(0x001F0000 & db[i+1])>>16] << 24)) |
		 ((4*i+7 >= dblen)?    0x00 :
	         (bit4[(unsigned long)(0x1F000000 & db[i+1])>>24] << 28));
      break;
	  
    case 5:
      for(i=0; i<length-3; i+=3)
      {
        dbc[(2*i)/3]   = ((unsigned long)(0x0000001F & db[i])   <<  2) |
	                 ((unsigned long)(0x00001F00 & db[i])   >>  1) |
	                 ((unsigned long)(0x001F0000 & db[i])   >>  4) |
	                 ((unsigned long)(0x1F000000 & db[i])   >>  7) |
	                 ((unsigned long)(0x0000001F & db[i+1]) << 22) |
	                 ((unsigned long)(0x00001F00 & db[i+1]) << 19);
		     
        dbc[(2*i)/3+1] = ((unsigned long)(0x001F0000 & db[i+1]) >> 14) |
		         ((unsigned long)(0x1F000000 & db[i+1]) >> 17) |
	                 ((unsigned long)(0x0000001F & db[i+2]) << 12) |
	                 ((unsigned long)(0x00001F00 & db[i+2]) <<  9) |
	                 ((unsigned long)(0x001F0000 & db[i+2]) <<  6) |
		         ((unsigned long)(0x1F000000 & db[i+2]) <<  3);
      }
      dbc[(2*i)/3]   = ((4*i   >= dblen)? 0x00 :      /* out of chars? */
                       (  (unsigned long)(0x0000001F & db[i])   <<  2)) |
		       ((4*i+1 >= dblen)? 0x00 :
	               (  (unsigned long)(0x00001F00 & db[i])   >>  1)) |
		       ((4*i+2 >= dblen)? 0x00 :
	               (  (unsigned long)(0x001F0000 & db[i])   >>  4)) |
		       ((4*i+3 >= dblen)? 0x00 :
	               (  (unsigned long)(0x1F000000 & db[i])   >>  7)) |
		       ((4*i+4 >= dblen)? 0x00 :
	               (  (unsigned long)(0x0000001F & db[i+1]) << 22)) |
		       ((4*i+5 >= dblen)? 0x00 :
	               (  (unsigned long)(0x00001F00 & db[i+1]) << 19));
		     
      dbc[(2*i)/3+1] = ((4*i+6 >= dblen)? 0x00 :
                       (  (unsigned long)(0x001F0000 & db[i+1]) >> 14)) |
		       ((4*i+7 >= dblen)? 0x00 :
		       (  (unsigned long)(0x1F000000 & db[i+1]) >> 17)) |
		       ((4*i+8 >= dblen)? 0x00 :
	               (  (unsigned long)(0x0000001F & db[i+2]) << 12)) |
		       ((4*i+9 >= dblen)? 0x00 :
	               (  (unsigned long)(0x00001F00 & db[i+2]) <<  9)) |
		       ((4*i+10>= dblen)? 0x00 :
	               (  (unsigned long)(0x001F0000 & db[i+2]) <<  6)) |
		       ((4*i+11>= dblen)? 0x00 :
		       (  (unsigned long)(0x1F000000 & db[i+2]) <<  3));

      break;
      
    case 6:
      for(i=0; i<length-2; i+=2)
      {
        dbc[i/2] = ((0x00000008 & db[i])?   0x00 : /*  is it an A=1xxx ? */
		   (((unsigned long)(0x00000004 & ~db[i])   >>  1)   | 
		    ((unsigned long)(0x00000001 & ~db[i]))))         |
		   ((0x00000080 & db[i])?   0x00 :
		   (((unsigned long)(0x00000040 & ~db[i])   >>  3)   |
		    ((unsigned long)(0x00000010 & ~db[i])   >>  2))) |
		   ((0x00000800 & db[i])?   0x00 :
		   (((unsigned long)(0x00000400 & ~db[i])   >>  5)   |
		    ((unsigned long)(0x00000100 & ~db[i])   >>  4))) |
		   ((0x00008000 & db[i])?   0x00 :
		   (((unsigned long)(0x00004000 & ~db[i])   >>  7)   |
		    ((unsigned long)(0x00001000 & ~db[i])   >>  6))) |
		   ((0x00080000 & db[i])?   0x00 :
		   (((unsigned long)(0x00040000 & ~db[i])   >>  9)   |
		    ((unsigned long)(0x00010000 & ~db[i])   >>  8))) |
		   ((0x00800000 & db[i])?   0x00 :
		   (((unsigned long)(0x00400000 & ~db[i])   >> 11)   |
		    ((unsigned long)(0x00100000 & ~db[i])   >> 10))) |
		   ((0x08000000 & db[i])?   0x00 :
		   (((unsigned long)(0x04000000 & ~db[i])   >> 13)   |
		    ((unsigned long)(0x01000000 & ~db[i])   >> 12))) |
		   ((0x80000000 & db[i])?   0x00 :
	           (((unsigned long)(0x40000000 & ~db[i])   >> 15)   |
		    ((unsigned long)(0x10000000 & ~db[i])   >> 14))) |
		   ((0x00000008 & db[i+1])? 0x00 :
		   (((unsigned long)(0x00000004 & ~db[i+1]) << 15)   |
		    ((unsigned long)(0x00000001 & ~db[i+1]) << 16))) |
	           ((0x00000080 & db[i+1])? 0x00 :
		   (((unsigned long)(0x00000040 & ~db[i+1]) << 13)   |
		    ((unsigned long)(0x00000010 & ~db[i+1]) << 14))) |
		   ((0x00000800 & db[i+1])? 0x00 :
		   (((unsigned long)(0x00000400 & ~db[i+1]) << 11)   |
		    ((unsigned long)(0x00000100 & ~db[i+1]) << 12))) |
		   ((0x00008000 & db[i+1])? 0x00 :
		   (((unsigned long)(0x00004000 & ~db[i+1]) <<  9)   |
		    ((unsigned long)(0x00001000 & ~db[i+1]) << 10))) |
		   ((0x00080000 & db[i+1])? 0x00 :
		   (((unsigned long)(0x00040000 & ~db[i+1]) <<  7)   |
		    ((unsigned long)(0x00010000 & ~db[i+1]) <<  8))) |
		   ((0x00800000 & db[i+1])? 0x00 :
		   (((unsigned long)(0x00400000 & ~db[i+1]) <<  5)   |
		    ((unsigned long)(0x00100000 & ~db[i+1]) <<  6))) |
		   ((0x08000000 & db[i+1])? 0x00 :
		   (((unsigned long)(0x04000000 & ~db[i+1]) <<  3)   |
		    ((unsigned long)(0x01000000 & ~db[i+1]) <<  4))) |
		   ((0x80000000 & db[i+1])? 0x00 :
	           (((unsigned long)(0x40000000 & ~db[i+1]) <<  1)   |
		    ((unsigned long)(0x10000000 & ~db[i+1]) <<  2)));
      }
      dbc[i/2] = ((4*i   >= dblen)? 0x00 :     /* out of chars? */
	         (  (unsigned long)(0x00000008 &  db[i])?   0x00 : /*  is it an A=1xxx ? */
		 (( (unsigned long)(0x00000004 & ~db[i])   >>  1)      |
		   ((unsigned long)(0x00000001 & ~db[i])))))           |
		 ((4*i+1 >= dblen)? 0x00 :
	         (  (unsigned long)(0x00000080 &  db[i])?   0x00 :
		 (( (unsigned long)(0x00000040 & ~db[i])   >>  3)      |
		   ((unsigned long)(0x00000010 & ~db[i])   >>  2))))   |
		 ((4*i+2 >= dblen)? 0x00 :
	         (  (unsigned long)(0x00000800 &  db[i])?   0x00 :
		 (( (unsigned long)(0x00000400 & ~db[i])   >>  5)      |
		   ((unsigned long)(0x00000100 & ~db[i])   >>  4))))   |
		 ((4*i+3 >= dblen)? 0x00 :
	         (  (unsigned long)(0x00008000 &  db[i])?   0x00 :
		 (( (unsigned long)(0x00004000 & ~db[i])   >>  7)      |
		   ((unsigned long)(0x00001000 & ~db[i])   >>  6))))   |
		 ((4*i+4 >= dblen)? 0x00 :
	         (  (unsigned long)(0x00080000 &  db[i])?   0x00 :
		 (( (unsigned long)(0x00040000 & ~db[i])   >>  9)      |
		   ((unsigned long)(0x00010000 & ~db[i])   >>  8))))   |
		 ((4*i+5 >= dblen)? 0x00 :
	         (  (unsigned long)(0x00800000 &  db[i])?   0x00 :
		 (( (unsigned long)(0x00400000 & ~db[i])   >> 11)      |
		   ((unsigned long)(0x00100000 & ~db[i])   >> 10))))   |
		 ((4*i+6 >= dblen)? 0x00 :
	         (  (unsigned long)(0x08000000 &  db[i])?   0x00 :
		 (( (unsigned long)(0x04000000 & ~db[i])   >> 13)      |
		   ((unsigned long)(0x01000000 & ~db[i])   >> 12))))   |
		 ((4*i+7 >= dblen)? 0x00 :
		 (  (unsigned long)(0x80000000 &  db[i])?   0x00 :
	         (( (unsigned long)(0x40000000 & ~db[i])   >> 15)      |
		   ((unsigned long)(0x10000000 & ~db[i])   >> 14))))   |
		 ((4*i+8 >= dblen)? 0x00 :
	         (  (unsigned long)(0x00000008 &  db[i+1])? 0x00 :
		 (( (unsigned long)(0x00000004 & ~db[i+1]) << 15)      |
		   ((unsigned long)(0x00000001 & ~db[i+1]) << 16))))   | 
		 ((4*i+9 >= dblen)? 0x00 :
	         (  (unsigned long)(0x00000080 &  db[i+1])? 0x00 :
		 (( (unsigned long)(0x00000040 & ~db[i+1]) << 13)      |
		   ((unsigned long)(0x00000010 & ~db[i+1]) << 14))))   |
		 ((4*i+10>= dblen)? 0x00 :
	         (  (unsigned long)(0x00000800 &  db[i+1])? 0x00 :
		 (( (unsigned long)(0x00000400 & ~db[i+1]) << 11)      |
		   ((unsigned long)(0x00000100 & ~db[i+1]) << 12))))   |
		 ((4*i+11>= dblen)? 0x00 :
	         (  (unsigned long)(0x00008000 &  db[i+1])? 0x00 :
		 (( (unsigned long)(0x00004000 & ~db[i+1]) <<  9)      |
		   ((unsigned long)(0x00001000 & ~db[i+1]) << 10))))   |
		 ((4*i+12>= dblen)? 0x00 :
	         (  (unsigned long)(0x00080000 &  db[i+1])? 0x00 :
		 (( (unsigned long)(0x00040000 & ~db[i+1]) <<  7)      |
		   ((unsigned long)(0x00010000 & ~db[i+1]) <<  8))))   |
		 ((4*i+13>= dblen)? 0x00 :
	         (  (unsigned long)(0x00800000 &  db[i+1])? 0x00 :
		 (( (unsigned long)(0x00400000 & ~db[i+1]) <<  5)      |
		   ((unsigned long)(0x00100000 & ~db[i+1]) <<  6))))   |
                 ((4*i+14>= dblen)? 0x00 :
	         (  (unsigned long)(0x08000000 &  db[i+1])? 0x00 :
		 (( (unsigned long)(0x04000000 & ~db[i+1]) <<  3)      |
		   ((unsigned long)(0x01000000 & ~db[i+1]) <<  4))))   |
                 ((4*i+15>= dblen)? 0x00 :
                 (  (unsigned long)(0x80000000 &  db[i+1])? 0x00 :
	         (( (unsigned long)(0x40000000 & ~db[i+1]) <<  1)      |
		   ((unsigned long)(0x10000000 & ~db[i+1]) <<  2))));
      break;
#endif /* little-endian */

    default: 
      fprintf(stderr, "cb_compress: Invalid mode parameter %d, returning...\n", mode);
      return;
  }
#endif /* 32 bit */

#ifdef LONG64
  length = (dblen+7)/8;
  switch (mode)
  {
#ifdef B_ENDIAN
    case 2:
#pragma vdir nodep altcode /* SX6 */
#pragma _CRI ivdep /* SV1 */
      for(i=0; i<length-4; i+=4)
      {
        dbc[i/4] = ((unsigned long)(0x0600000000000000 & db[i])   <<  5) |
	           ((unsigned long)(0x0006000000000000 & db[i])   << 11) |
	           ((unsigned long)(0x0000060000000000 & db[i])   << 17) |
	           ((unsigned long)(0x0000000600000000 & db[i])   << 23) |
		   ((unsigned long)(0x0000000006000000 & db[i])   << 29) |
	           ((unsigned long)(0x0000000000060000 & db[i])   << 35) |
	           ((unsigned long)(0x0000000000000600 & db[i])   << 41) |
	           ((unsigned long)(0x0000000000000006 & db[i])   << 47) |
	           ((unsigned long)(0x0600000000000000 & db[i+1]) >> 11) |
	           ((unsigned long)(0x0006000000000000 & db[i+1]) >>  5) |
	           ((unsigned long)(0x0000060000000000 & db[i+1]) <<  1) |
	           ((unsigned long)(0x0000000600000000 & db[i+1]) <<  7) |
		   ((unsigned long)(0x0000000006000000 & db[i+1]) << 13) |
	           ((unsigned long)(0x0000000000060000 & db[i+1]) << 19) |
	           ((unsigned long)(0x0000000000000600 & db[i+1]) << 25) |
	           ((unsigned long)(0x0000000000000006 & db[i+1]) << 31) |
                   ((unsigned long)(0x0600000000000000 & db[i+2]) >> 27) |
                   ((unsigned long)(0x0006000000000000 & db[i+2]) >> 21) |
		   ((unsigned long)(0x0000060000000000 & db[i+2]) >> 15) |
		   ((unsigned long)(0x0000000600000000 & db[i+2]) >>  9) |
		   ((unsigned long)(0x0000000006000000 & db[i+2]) >>  3) |
	           ((unsigned long)(0x0000000000060000 & db[i+2]) <<  3) |
	           ((unsigned long)(0x0000000000000600 & db[i+2]) <<  9) |
	           ((unsigned long)(0x0000000000000006 & db[i+2]) << 15) |
		   ((unsigned long)(0x0600000000000000 & db[i+3]) >> 43) |
		   ((unsigned long)(0x0006000000000000 & db[i+3]) >> 37) |
		   ((unsigned long)(0x0000060000000000 & db[i+3]) >> 31) |
		   ((unsigned long)(0x0000000600000000 & db[i+3]) >> 25) |
		   ((unsigned long)(0x0000000006000000 & db[i+3]) >> 19) |
	           ((unsigned long)(0x0000000000060000 & db[i+3]) >> 13) |
	           ((unsigned long)(0x0000000000000600 & db[i+3]) >>  7) |
	           ((unsigned long)(0x0000000000000006 & db[i+3]) >>  1);
      }
      dbc[i/4] = ((8*i   >= dblen)? 0x00 :   /* out of chars? */
                 (  (unsigned long)(0x0600000000000000 & db[i])   <<  5))  |
		 ((8*i+1 >= dblen)? 0x00 :
	         (  (unsigned long)(0x0006000000000000 & db[i])   << 11))  |
		 ((8*i+2 >= dblen)? 0x00 :
	         (  (unsigned long)(0x0000060000000000 & db[i])   << 17))  |
		 ((8*i+3 >= dblen)? 0x00 :
	         (  (unsigned long)(0x0000000600000000 & db[i])   << 23))  |
		 ((8*i+4 >= dblen)? 0x00 :
		 (  (unsigned long)(0x0000000006000000 & db[i])   << 29))  |
		 ((8*i+5 >= dblen)? 0x00 :
	         (  (unsigned long)(0x0000000000060000 & db[i])   << 35))  |
		 ((8*i+6 >= dblen)? 0x00 :
	         (  (unsigned long)(0x0000000000000600 & db[i])   << 41))  |
		 ((8*i+7 >= dblen)? 0x00 :
	         (  (unsigned long)(0x0000000000000006 & db[i])   << 47))  |
		 ((8*i+8 >= dblen)? 0x00 :
	         (  (unsigned long)(0x0600000000000000 & db[i+1]) >> 11))  |
		 ((8*i+9 >= dblen)? 0x00 :
	         (  (unsigned long)(0x0006000000000000 & db[i+1]) >>  5))  |
		 ((8*i+10>= dblen)? 0x00 :
	         (  (unsigned long)(0x0000060000000000 & db[i+1]) <<  1))  |
		 ((8*i+11>= dblen)? 0x00 :
	         (  (unsigned long)(0x0000000600000000 & db[i+1]) <<  7))  |
		 ((8*i+12>= dblen)? 0x00 :
		 (  (unsigned long)(0x0000000006000000 & db[i+1]) << 13))  |
		 ((8*i+13>= dblen)? 0x00 :
	         (  (unsigned long)(0x0000000000060000 & db[i+1]) << 19))  |
		 ((8*i+14>= dblen)? 0x00 :
	         (  (unsigned long)(0x0000000000000600 & db[i+1]) << 25))  |
		 ((8*i+15>= dblen)? 0x00 :
	         (  (unsigned long)(0x0000000000000006 & db[i+1]) << 31))  |
		 ((8*i+16>= dblen)? 0x00 :
                 (  (unsigned long)(0x0600000000000000 & db[i+2]) >> 27))  |
		 ((8*i+17>= dblen)? 0x00 :
                 (  (unsigned long)(0x0006000000000000 & db[i+2]) >> 21))  |
		 ((8*i+18>= dblen)? 0x00 :
		 (  (unsigned long)(0x0000060000000000 & db[i+2]) >> 15))  |
		 ((8*i+19>= dblen)? 0x00 :
		 (  (unsigned long)(0x0000000600000000 & db[i+2]) >>  9))  |
		 ((8*i+20>= dblen)? 0x00 :
		 (  (unsigned long)(0x0000000006000000 & db[i+2]) >>  3))  |
		 ((8*i+21>= dblen)? 0x00 :
	         (  (unsigned long)(0x0000000000060000 & db[i+2]) <<  3))  |
		 ((8*i+22>= dblen)? 0x00 :
	         (  (unsigned long)(0x0000000000000600 & db[i+2]) <<  9))  |
		 ((8*i+23>= dblen)? 0x00 :
	         (  (unsigned long)(0x0000000000000006 & db[i+2]) << 15))  |
		 ((8*i+24>= dblen)? 0x00 :
		 (  (unsigned long)(0x0600000000000000 & db[i+3]) >> 43))  |
		 ((8*i+25>= dblen)? 0x00 :
		 (  (unsigned long)(0x0006000000000000 & db[i+3]) >> 37))  |
		 ((8*i+26>= dblen)? 0x00 :
		 (  (unsigned long)(0x0000060000000000 & db[i+3]) >> 31))  |
		 ((8*i+27>= dblen)? 0x00 :
		 (  (unsigned long)(0x0000000600000000 & db[i+3]) >> 25))  |
		 ((8*i+28>= dblen)? 0x00 :
		 (  (unsigned long)(0x0000000006000000 & db[i+3]) >> 19))  |
		 ((8*i+29>= dblen)? 0x00 :
	         (  (unsigned long)(0x0000000000060000 & db[i+3]) >> 13))  |
		 ((8*i+30>= dblen)? 0x00 :
	         (  (unsigned long)(0x0000000000000600 & db[i+3]) >>  7))  |
		 ((8*i+31>= dblen)? 0x00 :
	         (  (unsigned long)(0x0000000000000006 & db[i+3]) >>  1));
      break;
	  
    case 4:
#pragma vdir nodep altcode
#pragma _CRI ivdep /* SV1 */
      for(i=0; i<length-2; i+=2)
      {
        dbc[i/2] = (bit4[(unsigned long)(0x1F00000000000000 & db[i]) >> 56] << 60) |
	           (bit4[(unsigned long)(0x001F000000000000 & db[i]) >> 48] << 56) |
	           (bit4[(unsigned long)(0x00001F0000000000 & db[i]) >> 40] << 52) |
	           (bit4[(unsigned long)(0x0000001F00000000 & db[i]) >> 32] << 48) |
		   (bit4[(unsigned long)(0x000000001F000000 & db[i]) >> 24] << 44) |
	           (bit4[(unsigned long)(0x00000000001F0000 & db[i]) >> 16] << 40) |
	           (bit4[(unsigned long)(0x0000000000001F00 & db[i]) >>  8] << 36) |
	           (bit4[                0x000000000000001F & db[i]]        << 32) |
	           (bit4[(unsigned long)(0x1F00000000000000 & db[i+1])>>56] << 28) |
	           (bit4[(unsigned long)(0x001F000000000000 & db[i+1])>>48] << 24) |
	           (bit4[(unsigned long)(0x00001F0000000000 & db[i+1])>>40] << 20) |
	           (bit4[(unsigned long)(0x0000001F00000000 & db[i+1])>>32] << 16) |
                   (bit4[(unsigned long)(0x000000001F000000 & db[i+1])>>24] << 12) |
	           (bit4[(unsigned long)(0x00000000001F0000 & db[i+1])>>16] <<  8) |
                   (bit4[(unsigned long)(0x0000000000001F00 & db[i+1])>> 8] <<  4) |
                   (bit4[                0x000000000000001F & db[i+1]]);
      }
      dbc[i/2] = ((8*i   >= dblen)?    0x00 :    /* out of chars? */
                 (bit4[(unsigned long)(0x1F00000000000000 & db[i]) >> 56] << 60))  |
                 ((8*i+1 >= dblen)?    0x00 :
	         (bit4[(unsigned long)(0x001F000000000000 & db[i]) >> 48] << 56))  |
		 ((8*i+2 >= dblen)?    0x00 :
	         (bit4[(unsigned long)(0x00001F0000000000 & db[i]) >> 40] << 52))  |
		 ((8*i+3 >= dblen)?    0x00 :
	         (bit4[(unsigned long)(0x0000001F00000000 & db[i]) >> 32] << 48))  |
		 ((8*i+4 >= dblen)?    0x00 :
		 (bit4[(unsigned long)(0x000000001F000000 & db[i]) >> 24] << 44))  |
		 ((8*i+5 >= dblen)?    0x00 :
	         (bit4[(unsigned long)(0x00000000001F0000 & db[i]) >> 16] << 40))  |
		 ((8*i+6 >= dblen)?    0x00 :
	         (bit4[(unsigned long)(0x0000000000001F00 & db[i]) >>  8] << 36))  |
		 ((8*i+7 >= dblen)?    0x00 :
	         (bit4[                0x000000000000001F & db[i]]        << 32))  |
		 ((8*i+8 >= dblen)?    0x00 :
	         (bit4[(unsigned long)(0x1F00000000000000 & db[i+1])>>56] << 28))  |
		 ((8*i+9 >= dblen)?    0x00 :
	         (bit4[(unsigned long)(0x001F000000000000 & db[i+1])>>48] << 24))  |
		 ((8*i+10>= dblen)?    0x00 :
	         (bit4[(unsigned long)(0x00001F0000000000 & db[i+1])>>40] << 20))  |
		 ((8*i+11>= dblen)?    0x00 :
	         (bit4[(unsigned long)(0x0000001F00000000 & db[i+1])>>32] << 16))  |
		 ((8*i+12>= dblen)?    0x00 :
                 (bit4[(unsigned long)(0x000000001F000000 & db[i+1])>>24] << 12))  |
		 ((8*i+13>= dblen)?    0x00 :
	         (bit4[(unsigned long)(0x00000000001F0000 & db[i+1])>>16] <<  8))  |
		 ((8*i+14>= dblen)?    0x00 :
                 (bit4[(unsigned long)(0x0000000000001F00 & db[i+1])>> 8] <<  4))  |
		 ((8*i+15>= dblen)?    0x00 :
                  bit4[                0x000000000000001F & db[i+1]]);
      break;
	  
    case 5:
#pragma vdir nodep altcode
#pragma _CRI ivdep
      for(i=0; i<length-3; i+=3)
      {
        dbc[(2*i)/3]   = ((unsigned long)(0x1F00000000000000 & db[i])   >>  1) |
	                 ((unsigned long)(0x001F000000000000 & db[i])   <<  2) |
	                 ((unsigned long)(0x00001F0000000000 & db[i])   <<  5) |
	                 ((unsigned long)(0x0000001F00000000 & db[i])   <<  8) |
	                 ((unsigned long)(0x000000001F000000 & db[i])   << 11) |
	                 ((unsigned long)(0x00000000001F0000 & db[i])   << 14) |
	                 ((unsigned long)(0x0000000000001F00 & db[i])   << 17) |
	                 ((unsigned long)(0x000000000000001F & db[i])   << 20) |
	                 ((unsigned long)(0x1F00000000000000 & db[i+1]) >> 41) |
	                 ((unsigned long)(0x001F000000000000 & db[i+1]) >> 38) |
			 ((unsigned long)(0x00001F0000000000 & db[i+1]) >> 35) |
	                 ((unsigned long)(0x0000001F00000000 & db[i+1]) >> 32);
		     
        dbc[(2*i)/3+1] = ((unsigned long)(0x000000001F000000 & db[i+1]) << 31) |
	                 ((unsigned long)(0x00000000001F0000 & db[i+1]) << 34) |
	                 ((unsigned long)(0x0000000000001F00 & db[i+1]) << 37) |
	                 ((unsigned long)(0x000000000000001F & db[i+1]) << 40) |
	                 ((unsigned long)(0x1F00000000000000 & db[i+2]) >> 21) |
	                 ((unsigned long)(0x001F000000000000 & db[i+2]) >> 18) |
	                 ((unsigned long)(0x00001F0000000000 & db[i+2]) >> 15) |
	                 ((unsigned long)(0x0000001F00000000 & db[i+2]) >> 12) |
		         ((unsigned long)(0x000000001F000000 & db[i+2]) >>  9) |
	                 ((unsigned long)(0x00000000001F0000 & db[i+2]) >>  6) |
	                 ((unsigned long)(0x0000000000001F00 & db[i+2]) >>  3) |
	                                 (0x000000000000001F & db[i+2]);
      }
      dbc[(2*i)/3]   = ((8*i   >= dblen)? 0x00 :/* out of chars? */
                       (  (unsigned long)(0x1F00000000000000 & db[i])   >>  1))  |
		       ((8*i+1 >= dblen)? 0x00 :
	               (  (unsigned long)(0x001F000000000000 & db[i])   <<  2))  |
		       ((8*i+2 >= dblen)? 0x00 :
	               (  (unsigned long)(0x00001F0000000000 & db[i])   <<  5))  |
		       ((8*i+3 >= dblen)? 0x00 :
	               (  (unsigned long)(0x0000001F00000000 & db[i])   <<  8))  |
		       ((8*i+4 >= dblen)? 0x00 :
	               (  (unsigned long)(0x000000001F000000 & db[i])   << 11))  |
		       ((8*i+5 >= dblen)? 0x00 :
	               (  (unsigned long)(0x00000000001F0000 & db[i])   << 14))  |
		       ((8*i+6 >= dblen)? 0x00 :
	               (  (unsigned long)(0x0000000000001F00 & db[i])   << 17))  |
		       ((8*i+7 >= dblen)? 0x00 :
	               (  (unsigned long)(0x000000000000001F & db[i])   << 20))  |
		       ((8*i+8 >= dblen)? 0x00 :
	               (  (unsigned long)(0x1F00000000000000 & db[i+1]) >> 41))  |
		       ((8*i+9 >= dblen)? 0x00 :
	               (  (unsigned long)(0x001F000000000000 & db[i+1]) >> 38))  |
		       ((8*i+10>= dblen)? 0x00 :
	               (  (unsigned long)(0x00001F0000000000 & db[i+1]) >> 35))  |
		       ((8*i+11>= dblen)? 0x00 :
	               (  (unsigned long)(0x0000001F00000000 & db[i+1]) >> 32));
		     
      dbc[(2*i)/3+1] = ((8*i+12>= dblen)? 0x00 :
                       (  (unsigned long)(0x000000001F000000 & db[i+1]) << 31))  |
		       ((8*i+13>= dblen)? 0x00 :
	               (  (unsigned long)(0x00000000001F0000 & db[i+1]) << 34))  |
		       ((8*i+14>= dblen)? 0x00 :
		       (  (unsigned long)(0x0000000000001F00 & db[i+1]) << 37))  |
		       ((8*i+15>= dblen)? 0x00 :
	               (  (unsigned long)(0x000000000000001F & db[i+1]) << 40))  |
		       ((8*i+16>= dblen)? 0x00 :
                       (  (unsigned long)(0x1F00000000000000 & db[i+2]) >> 21))  |
		       ((8*i+17>= dblen)? 0x00 :
	               (  (unsigned long)(0x001F000000000000 & db[i+2]) >> 18))  |
		       ((8*i+18>= dblen)? 0x00 :
	               (  (unsigned long)(0x00001F0000000000 & db[i+2]) >> 15))  |
		       ((8*i+19>= dblen)? 0x00 :
	               (  (unsigned long)(0x0000001F00000000 & db[i+2]) >> 12))  |
		       ((8*i+20>= dblen)? 0x00 :
	               (  (unsigned long)(0x000000001F000000 & db[i+2]) >>  9))  |
		       ((8*i+21>= dblen)? 0x00 :
	               (  (unsigned long)(0x00000000001F0000 & db[i+2]) >>  6))  |
		       ((8*i+22>= dblen)? 0x00 :
	               (  (unsigned long)(0x0000000000001F00 & db[i+2]) >>  3))  |
		       ((8*i+23>= dblen)? 0x00 :
	                                 (0x000000000000001F & db[i+2]));

      break;
      
    case 6:
#pragma vdir nodep altcode
#pragma _CRI ivdep
      for(i=0; i<length-2; i+=2)
      {
        dbc[i/2] = ((0x8000000000000000 & db[i])? 0x00 : /*  is it an A=1xxx ? */
	           (((unsigned long)(0x4000000000000000 & ~db[i])  << 1)   |
		    ((unsigned long)(0x1000000000000000 & ~db[i])  << 2))) |
	           ((0x0800000000000000 & db[i])? 0x00 :
		   (((unsigned long)(0x0400000000000000 & ~db[i])  << 3)   |
		    ((unsigned long)(0x0100000000000000 & ~db[i])  << 4))) |
	           ((0x0080000000000000 & db[i])? 0x00 :
		   (((unsigned long)(0x0040000000000000 & ~db[i])  << 5)   |
		    ((unsigned long)(0x0010000000000000 & ~db[i])  << 6))) |
	           ((0x0008000000000000 & db[i])? 0x00 :
		   (((unsigned long)(0x0004000000000000 & ~db[i])  << 7)   |
		    ((unsigned long)(0x0001000000000000 & ~db[i])  << 8))) |
	           ((0x0000800000000000 & db[i])? 0x00 :
		   (((unsigned long)(0x0000400000000000 & ~db[i])  << 9)   |
		    ((unsigned long)(0x0000100000000000 & ~db[i])  <<10))) |
	           ((0x0000080000000000 & db[i])? 0x00 :
		   (((unsigned long)(0x0000040000000000 & ~db[i])  <<11)   |
		    ((unsigned long)(0x0000010000000000 & ~db[i])  <<12))) |
	           ((0x0000008000000000 & db[i])? 0x00 :
		   (((unsigned long)(0x0000004000000000 & ~db[i])  <<13)   |
		    ((unsigned long)(0x0000001000000000 & ~db[i])  <<14))) |
	           ((0x0000000800000000 & db[i])? 0x00 :
		   (((unsigned long)(0x0000000400000000 & ~db[i])  <<15)   |
		    ((unsigned long)(0x0000000100000000 & ~db[i])  <<16))) |
		   ((0x0000000080000000 & db[i])? 0x00 :
	           (((unsigned long)(0x0000000040000000 & ~db[i])  <<17)   |
		    ((unsigned long)(0x0000000010000000 & ~db[i])  <<18))) |
	           ((0x0000000008000000 & db[i])? 0x00 :
		   (((unsigned long)(0x0000000004000000 & ~db[i])  <<19)   |
		    ((unsigned long)(0x0000000001000000 & ~db[i])  <<20))) |
	           ((0x0000000000800000 & db[i])? 0x00 :
		   (((unsigned long)(0x0000000000400000 & ~db[i])  <<21)   |
		    ((unsigned long)(0x0000000000100000 & ~db[i])  <<22))) |
	           ((0x0000000000080000 & db[i])? 0x00 :
		   (((unsigned long)(0x0000000000040000 & ~db[i])  <<23)   |
		    ((unsigned long)(0x0000000000010000 & ~db[i])  <<24))) |
	           ((0x0000000000008000 & db[i])? 0x00 :
		   (((unsigned long)(0x0000000000004000 & ~db[i])  <<25)   |
		    ((unsigned long)(0x0000000000001000 & ~db[i])  <<26))) |
	           ((0x0000000000000800 & db[i])? 0x00 :
		   (((unsigned long)(0x0000000000000400 & ~db[i])  <<27)   |
		    ((unsigned long)(0x0000000000000100 & ~db[i])  <<28))) |
	           ((0x0000000000000080 & db[i])? 0x00 :
		   (((unsigned long)(0x0000000000000040 & ~db[i])  <<29)   |
		    ((unsigned long)(0x0000000000000010 & ~db[i])  <<30))) |
	           ((0x0000000000000008 & db[i])? 0x00 :
		   (((unsigned long)(0x0000000000000004 & ~db[i])  <<31)   |
		    ((unsigned long)(0x0000000000000001 & ~db[i])  <<32))) |
		   ((0x8000000000000000 & db[i+1])? 0x00 :
	           (((unsigned long)(0x4000000000000000 & ~db[i+1])>>31)   |
		    ((unsigned long)(0x1000000000000000 & ~db[i+1])>>30))) |
	           ((0x0800000000000000 & db[i+1])? 0x00 :
		   (((unsigned long)(0x0400000000000000 & ~db[i+1])>>29)   |
		    ((unsigned long)(0x0100000000000000 & ~db[i+1])>>28))) |
	           ((0x0080000000000000 & db[i+1])? 0x00 :
		   (((unsigned long)(0x0040000000000000 & ~db[i+1])>>27)   |
		    ((unsigned long)(0x0010000000000000 & ~db[i+1])>>26))) |
	           ((0x0008000000000000 & db[i+1])? 0x00 :
		   (((unsigned long)(0x0004000000000000 & ~db[i+1])>>25)   |
		    ((unsigned long)(0x0001000000000000 & ~db[i+1])>>24))) |
	           ((0x0000800000000000 & db[i+1])? 0x00 :
		   (((unsigned long)(0x0000400000000000 & ~db[i+1])>>23)   |
		    ((unsigned long)(0x0000100000000000 & ~db[i+1])>>22))) |
	           ((0x0000080000000000 & db[i+1])? 0x00 :
		   (((unsigned long)(0x0000040000000000 & ~db[i+1])>>21)   |
		    ((unsigned long)(0x0000010000000000 & ~db[i+1])>>20))) |
	           ((0x0000008000000000 & db[i+1])? 0x00 :
		   (((unsigned long)(0x0000004000000000 & ~db[i+1])>>19)   |
		    ((unsigned long)(0x0000001000000000 & ~db[i+1])>>18))) |
	           ((0x0000000800000000 & db[i+1])? 0x00 :
		   (((unsigned long)(0x0000000400000000 & ~db[i+1])>>17)   |
		    ((unsigned long)(0x0000000100000000 & ~db[i+1])>>16))) |
		   ((0x0000000080000000 & db[i+1])? 0x00 :
	           (((unsigned long)(0x0000000040000000 & ~db[i+1])>>15)   |
		    ((unsigned long)(0x0000000010000000 & ~db[i+1])>>14))) |
	           ((0x0000000008000000 & db[i+1])? 0x00 :
		   (((unsigned long)(0x0000000004000000 & ~db[i+1])>>13)   |
		    ((unsigned long)(0x0000000001000000 & ~db[i+1])>>12))) |
	           ((0x0000000000800000 & db[i+1])? 0x00 :
		   (((unsigned long)(0x0000000000400000 & ~db[i+1])>>11)   |
		    ((unsigned long)(0x0000000000100000 & ~db[i+1])>>10))) |
	           ((0x0000000000080000 & db[i+1])? 0x00 :
		   (((unsigned long)(0x0000000000040000 & ~db[i+1])>> 9)   |
		    ((unsigned long)(0x0000000000010000 & ~db[i+1])>> 8))) |
	           ((0x0000000000008000 & db[i+1])? 0x00 :
		   (((unsigned long)(0x0000000000004000 & ~db[i+1])>> 7)   |
		    ((unsigned long)(0x0000000000001000 & ~db[i+1])>> 6))) |
	           ((0x0000000000000800 & db[i+1])? 0x00 :
		   (((unsigned long)(0x0000000000000400 & ~db[i+1])>> 5)   |
		    ((unsigned long)(0x0000000000000100 & ~db[i+1])>> 4))) |
	           ((0x0000000000000080 & db[i+1])? 0x00 :
		   (((unsigned long)(0x0000000000000040 & ~db[i+1])>> 3)   |
		    ((unsigned long)(0x0000000000000010 & ~db[i+1])>> 2))) |
	           ((0x0000000000000008 & db[i+1])? 0x00 :
		   (((unsigned long)(0x0000000000000004 & ~db[i+1])>> 1)   |
		    ((unsigned long)(0x0000000000000001 & ~db[i+1]))));
      }
      dbc[i/2] = ((8*i   >= dblen)?0x00 :              /* out of chars? */
                 ((0x8000000000000000 & db[i])? 0x00 : /*  is it an A=1xxx ? */
	         (((unsigned long)(0x4000000000000000 & ~db[i])  << 1)    |
		  ((unsigned long)(0x1000000000000000 & ~db[i])  << 2)))) |
		 ((8*i+1 >= dblen)?0x00 :
	         ((0x0800000000000000 & db[i])? 0x00 :
		 (((unsigned long)(0x0400000000000000 & ~db[i])  << 3)    |
		  ((unsigned long)(0x0100000000000000 & ~db[i])  << 4)))) |
		 ((8*i+2 >= dblen)?0x00 :
	         ((0x0080000000000000 & db[i])? 0x00 :
		 (((unsigned long)(0x0040000000000000 & ~db[i])  << 5)    |
		  ((unsigned long)(0x0010000000000000 & ~db[i])  << 6)))) |
		 ((8*i+3 >= dblen)?0x00 :
	         ((0x0008000000000000 & db[i])? 0x00 :
		 (((unsigned long)(0x0004000000000000 & ~db[i])  << 7)    |
		  ((unsigned long)(0x0001000000000000 & ~db[i])  << 8)))) |
		 ((8*i+4 >= dblen)?0x00 :
	         ((0x0000800000000000 & db[i])? 0x00 :
		 (((unsigned long)(0x0000400000000000 & ~db[i])  << 9)    |
		  ((unsigned long)(0x0000100000000000 & ~db[i])  <<10)))) |
		 ((8*i+5 >= dblen)?0x00 :
	         ((0x0000080000000000 & db[i])? 0x00 :
		 (((unsigned long)(0x0000040000000000 & ~db[i])  <<11)    |
		  ((unsigned long)(0x0000010000000000 & ~db[i])  <<12)))) |
		 ((8*i+6 >= dblen)?0x00 :
	         ((0x0000008000000000 & db[i])? 0x00 :
		 (((unsigned long)(0x0000004000000000 & ~db[i])  <<13)    |
		  ((unsigned long)(0x0000001000000000 & ~db[i])  <<14)))) |
		 ((8*i+7 >= dblen)?0x00 :
	         ((0x0000000800000000 & db[i])? 0x00 :
		 (((unsigned long)(0x0000000400000000 & ~db[i])  <<15)    |
		  ((unsigned long)(0x0000000100000000 & ~db[i])  <<16)))) |
		 ((8*i+8 >= dblen)?0x00 :
		 ((0x0000000080000000 & db[i])? 0x00 :
	         (((unsigned long)(0x0000000040000000 & ~db[i])  <<17)    |
		  ((unsigned long)(0x0000000010000000 & ~db[i])  <<18)))) |
		 ((8*i+9 >= dblen)?0x00 :
	         ((0x0000000008000000 & db[i])? 0x00 :
		 (((unsigned long)(0x0000000004000000 & ~db[i])  <<19)    |
		  ((unsigned long)(0x0000000001000000 & ~db[i])  <<20)))) |
		 ((8*i+10>= dblen)?0x00 :
	         ((0x0000000000800000 & db[i])? 0x00 :
		 (((unsigned long)(0x0000000000400000 & ~db[i])  <<21)    |
		  ((unsigned long)(0x0000000000100000 & ~db[i])  <<22)))) |
		 ((8*i+11>= dblen)?0x00 :
	         ((0x0000000000080000 & db[i])? 0x00 :
		 (((unsigned long)(0x0000000000040000 & ~db[i])  <<23)    |
		  ((unsigned long)(0x0000000000010000 & ~db[i])  <<24)))) |
		 ((8*i+12>= dblen)?0x00 :
	         ((0x0000000000008000 & db[i])? 0x00 :
		 (((unsigned long)(0x0000000000004000 & ~db[i])  <<25)    |
		  ((unsigned long)(0x0000000000001000 & ~db[i])  <<26)))) |
		 ((8*i+13>= dblen)?0x00 :
	         ((0x0000000000000800 & db[i])? 0x00 :
		 (((unsigned long)(0x0000000000000400 & ~db[i])  <<27)    |
		  ((unsigned long)(0x0000000000000100 & ~db[i])  <<28)))) |
		 ((8*i+14>= dblen)?0x00 :
	         ((0x0000000000000080 & db[i])? 0x00 :
		 (((unsigned long)(0x0000000000000040 & ~db[i])  <<29)    |
		  ((unsigned long)(0x0000000000000010 & ~db[i])  <<30)))) |
		 ((8*i+15>= dblen)?0x00 :
	         ((0x0000000000000008 & db[i])? 0x00 :
		 (((unsigned long)(0x0000000000000004 & ~db[i])  <<31)    |
		  ((unsigned long)(0x0000000000000001 & ~db[i])  <<32)))) |
		 ((8*i+16>= dblen)?0x00 :
		 ((0x8000000000000000 & db[i+1])? 0x00 :
	         (((unsigned long)(0x4000000000000000 & ~db[i+1])>>31)    |
		  ((unsigned long)(0x1000000000000000 & ~db[i+1])>>30)))) |
		 ((8*i+17>= dblen)?0x00 :
	         ((0x0800000000000000 & db[i+1])? 0x00 :
		 (((unsigned long)(0x0400000000000000 & ~db[i+1])>>29)    |
		  ((unsigned long)(0x0100000000000000 & ~db[i+1])>>28)))) |
		 ((8*i+18>= dblen)?0x00 :
	         ((0x0080000000000000 & db[i+1])? 0x00 :
		 (((unsigned long)(0x0040000000000000 & ~db[i+1])>>27)    |
		  ((unsigned long)(0x0010000000000000 & ~db[i+1])>>26)))) |
		 ((8*i+19>= dblen)?0x00 :
	         ((0x0008000000000000 & db[i+1])? 0x00 :
		 (((unsigned long)(0x0004000000000000 & ~db[i+1])>>25)    |
		  ((unsigned long)(0x0001000000000000 & ~db[i+1])>>24)))) |
		 ((8*i+20>= dblen)?0x00 :
	         ((0x0000800000000000 & db[i+1])? 0x00 :
		 (((unsigned long)(0x0000400000000000 & ~db[i+1])>>23)    |
		  ((unsigned long)(0x0000100000000000 & ~db[i+1])>>22)))) |
		 ((8*i+21>= dblen)?0x00 :
	         ((0x0000080000000000 & db[i+1])? 0x00 :
		 (((unsigned long)(0x0000040000000000 & ~db[i+1])>>21)    |
		  ((unsigned long)(0x0000010000000000 & ~db[i+1])>>20)))) |
		 ((8*i+22>= dblen)?0x00 :
	         ((0x0000008000000000 & db[i+1])? 0x00 :
		 (((unsigned long)(0x0000004000000000 & ~db[i+1])>>19)    |
		  ((unsigned long)(0x0000001000000000 & ~db[i+1])>>18)))) |
		 ((8*i+23>= dblen)?0x00 :
	         ((0x0000000800000000 & db[i+1])? 0x00 :
		 (((unsigned long)(0x0000000400000000 & ~db[i+1])>>17)    |
		  ((unsigned long)(0x0000000100000000 & ~db[i+1])>>16)))) |
		 ((8*i+24>= dblen)?0x00 :
		 ((0x0000000080000000 & db[i+1])? 0x00 :
	         (((unsigned long)(0x0000000040000000 & ~db[i+1])>>15)    |
		  ((unsigned long)(0x0000000010000000 & ~db[i+1])>>14)))) |
		 ((8*i+25>= dblen)?0x00 :
	         ((0x0000000008000000 & db[i+1])? 0x00 :
		 (((unsigned long)(0x0000000004000000 & ~db[i+1])>>13)    |
		  ((unsigned long)(0x0000000001000000 & ~db[i+1])>>12)))) |
		 ((8*i+26>= dblen)?0x00 :
	         ((0x0000000000800000 & db[i+1])? 0x00 :
		 (((unsigned long)(0x0000000000400000 & ~db[i+1])>>11)    |
		  ((unsigned long)(0x0000000000100000 & ~db[i+1])>>10)))) |
		 ((8*i+27>= dblen)?0x00 :
	         ((0x0000000000080000 & db[i+1])? 0x00 :
		 (((unsigned long)(0x0000000000040000 & ~db[i+1])>> 9)    |
		  ((unsigned long)(0x0000000000010000 & ~db[i+1])>> 8)))) |
		 ((8*i+28>= dblen)?0x00 :
	         ((0x0000000000008000 & db[i+1])? 0x00 :
		 (((unsigned long)(0x0000000000004000 & ~db[i+1])>> 7)    |
		  ((unsigned long)(0x0000000000001000 & ~db[i+1])>> 6)))) |
		 ((8*i+29>= dblen)?0x00 :
	         ((0x0000000000000800 & db[i+1])? 0x00 :
		 (((unsigned long)(0x0000000000000400 & ~db[i+1])>> 5)    |
		  ((unsigned long)(0x0000000000000100 & ~db[i+1])>> 4)))) |
		 ((8*i+30>= dblen)?0x00 :
	         ((0x0000000000000080 & db[i+1])? 0x00 :
		 (((unsigned long)(0x0000000000000040 & ~db[i+1])>> 3)    |
		  ((unsigned long)(0x0000000000000010 & ~db[i+1])>> 2)))) |
		 ((8*i+31>= dblen)?0x00 :
	         ((0x0000000000000008 & db[i+1])? 0x00 :
		 (((unsigned long)(0x0000000000000004 & ~db[i+1])>> 1)    |
		  (               (0x0000000000000001 & ~db[i+1])))));
      break;
#endif /* big-endian */
#ifdef L_ENDIAN
    case 2:
      for(i=0; i<length-4; i+=4)
      {
        dbc[i/4] = ((unsigned long)(0x0000000000000006 & db[i])   >>  1) |
	           ((unsigned long)(0x0000000000000600 & db[i])   >>  7) |
	           ((unsigned long)(0x0000000000060000 & db[i])   >> 13) |
		   ((unsigned long)(0x0000000006000000 & db[i])   >> 19) |
		   ((unsigned long)(0x0000000600000000 & db[i])   >> 25) |
		   ((unsigned long)(0x0000060000000000 & db[i])   >> 31) |
		   ((unsigned long)(0x0006000000000000 & db[i])   >> 37) |
		   ((unsigned long)(0x0600000000000000 & db[i])   >> 43) |
	           ((unsigned long)(0x0000000000000006 & db[i+1]) << 15) |
	           ((unsigned long)(0x0000000000000600 & db[i+1]) <<  9) |
	           ((unsigned long)(0x0000000000060000 & db[i+1]) <<  3) |
		   ((unsigned long)(0x0000000006000000 & db[i+1]) >>  3) |
		   ((unsigned long)(0x0000000600000000 & db[i+1]) >>  9) |
		   ((unsigned long)(0x0000060000000000 & db[i+1]) >> 15) |
                   ((unsigned long)(0x0006000000000000 & db[i+1]) >> 21) |
                   ((unsigned long)(0x0600000000000000 & db[i+1]) >> 27) |
	           ((unsigned long)(0x0000000000000006 & db[i+2]) << 31) |
	           ((unsigned long)(0x0000000000000600 & db[i+2]) << 25) |
	           ((unsigned long)(0x0000000000060000 & db[i+2]) << 19) |
		   ((unsigned long)(0x0000000006000000 & db[i+2]) << 13) |
	           ((unsigned long)(0x0000000600000000 & db[i+2]) <<  7) |
	           ((unsigned long)(0x0000060000000000 & db[i+2]) <<  1) |
	           ((unsigned long)(0x0006000000000000 & db[i+2]) >>  5) |
	           ((unsigned long)(0x0600000000000000 & db[i+2]) >> 11) |
	           ((unsigned long)(0x0000000000000006 & db[i+3]) << 47) |
	           ((unsigned long)(0x0000000000000600 & db[i+3]) << 41) |
	           ((unsigned long)(0x0000000000060000 & db[i+3]) << 35) |
		   ((unsigned long)(0x0000000006000000 & db[i+3]) << 29) |
	           ((unsigned long)(0x0000000600000000 & db[i+3]) << 23) |
	           ((unsigned long)(0x0000060000000000 & db[i+3]) << 17) |
	           ((unsigned long)(0x0006000000000000 & db[i+3]) << 11) |
		   ((unsigned long)(0x0600000000000000 & db[i+3]) <<  5);
      }
      dbc[i/4] = ((8*i   >= dblen)? 0x00 :  /* out of chars? */
	         (  (unsigned long)(0x0000000000000006 & db[i])   >>  1)) |
		 ((8*i+1 >= dblen)? 0x00 :
	         (  (unsigned long)(0x0000000000000600 & db[i])   >>  7)) |
		 ((8*i+2 >= dblen)? 0x00 :
	         (  (unsigned long)(0x0000000000060000 & db[i])   >> 13)) |
		 ((8*i+3 >= dblen)? 0x00 :
		 (  (unsigned long)(0x0000000006000000 & db[i])   >> 19)) |
		 ((8*i+4 >= dblen)? 0x00 :
		 (  (unsigned long)(0x0000000600000000 & db[i])   >> 25)) |
		 ((8*i+5 >= dblen)? 0x00 :
		 (  (unsigned long)(0x0000060000000000 & db[i])   >> 31)) |
		 ((8*i+6 >= dblen)? 0x00 :
		 (  (unsigned long)(0x0006000000000000 & db[i])   >> 37)) |
		 ((8*i+7 >= dblen)? 0x00 :
		 (  (unsigned long)(0x0600000000000000 & db[i])   >> 43)) |
		 ((8*i+8 >= dblen)? 0x00 :
	         (  (unsigned long)(0x0000000000000006 & db[i+1]) << 15)) |
		 ((8*i+9 >= dblen)? 0x00 :
	         (  (unsigned long)(0x0000000000000600 & db[i+1]) <<  9)) |
		 ((8*i+10>= dblen)? 0x00 :
	         (  (unsigned long)(0x0000000000060000 & db[i+1]) <<  3)) |
		 ((8*i+11>= dblen)? 0x00 :
		 (  (unsigned long)(0x0000000006000000 & db[i+1]) >>  3)) |
		 ((8*i+12>= dblen)? 0x00 :
		 (  (unsigned long)(0x0000000600000000 & db[i+1]) >>  9)) |
		 ((8*i+13>= dblen)? 0x00 :
		 (  (unsigned long)(0x0000060000000000 & db[i+1]) >> 15)) |
		 ((8*i+14>= dblen)? 0x00 :
                 (  (unsigned long)(0x0006000000000000 & db[i+1]) >> 21)) |
		 ((8*i+15>= dblen)? 0x00 :
                 (  (unsigned long)(0x0600000000000000 & db[i+1]) >> 27)) |
		 ((8*i+16>= dblen)? 0x00 :
	         (  (unsigned long)(0x0000000000000006 & db[i+2]) << 31)) |
		 ((8*i+17>= dblen)? 0x00 :
	         (  (unsigned long)(0x0000000000000600 & db[i+2]) << 25)) |
		 ((8*i+18>= dblen)? 0x00 :
	         (  (unsigned long)(0x0000000000060000 & db[i+2]) << 19)) |
		 ((8*i+19>= dblen)? 0x00 :
		 (  (unsigned long)(0x0000000006000000 & db[i+2]) << 13)) |
		 ((8*i+20>= dblen)? 0x00 :
	         (  (unsigned long)(0x0000000600000000 & db[i+2]) <<  7)) |
		 ((8*i+21>= dblen)? 0x00 :
	         (  (unsigned long)(0x0000060000000000 & db[i+2]) <<  1)) |
		 ((8*i+22>= dblen)? 0x00 :
	         (  (unsigned long)(0x0006000000000000 & db[i+2]) >>  5)) |
		 ((8*i+23>= dblen)? 0x00 :
	         (  (unsigned long)(0x0600000000000000 & db[i+2]) >> 11)) |
		 ((8*i+24>= dblen)? 0x00 :
	         (  (unsigned long)(0x0000000000000006 & db[i+3]) << 47)) |
		 ((8*i+25>= dblen)? 0x00 :
	         (  (unsigned long)(0x0000000000000600 & db[i+3]) << 41)) |
		 ((8*i+26>= dblen)? 0x00 :
	         (  (unsigned long)(0x0000000000060000 & db[i+3]) << 35)) |
		 ((8*i+27>= dblen)? 0x00 :
		 (  (unsigned long)(0x0000000006000000 & db[i+3]) << 29)) |
		 ((8*i+28>= dblen)? 0x00 :
	         (  (unsigned long)(0x0000000600000000 & db[i+3]) << 23)) |
		 ((8*i+29>= dblen)? 0x00 :
	         (  (unsigned long)(0x0000060000000000 & db[i+3]) << 17)) |
		 ((8*i+30>= dblen)? 0x00 :
	         (  (unsigned long)(0x0006000000000000 & db[i+3]) << 11)) |
		 ((8*i+31>= dblen)? 0x00 : 
                 (  (unsigned long)(0x0600000000000000 & db[i+3]) <<  5));
      break;
	  
    case 4:
      for(i=0; i<length-2; i+=2)
      {
        dbc[i/2] = (bit4[                0x000000000000001F & db[i]])                |
                   (bit4[(unsigned long)(0x0000000000001F00 & db[i])   >>  8] <<  4) |
	           (bit4[(unsigned long)(0x00000000001F0000 & db[i])   >> 16] <<  8) |
                   (bit4[(unsigned long)(0x000000001F000000 & db[i])   >> 24] << 12) |
	           (bit4[(unsigned long)(0x0000001F00000000 & db[i])   >> 32] << 16) |
	           (bit4[(unsigned long)(0x00001F0000000000 & db[i])   >> 40] << 20) |
	           (bit4[(unsigned long)(0x001F000000000000 & db[i])   >> 48] << 24) |
	           (bit4[(unsigned long)(0x1F00000000000000 & db[i])   >> 56] << 28) |
	           (bit4[                0x000000000000001F & db[i+1]]        << 32) |
	           (bit4[(unsigned long)(0x0000000000001F00 & db[i+1]) >>  8] << 36) |
	           (bit4[(unsigned long)(0x00000000001F0000 & db[i+1]) >> 16] << 40) |
		   (bit4[(unsigned long)(0x000000001F000000 & db[i+1]) >> 24] << 44) |
	           (bit4[(unsigned long)(0x0000001F00000000 & db[i+1]) >> 32] << 48) |
	           (bit4[(unsigned long)(0x00001F0000000000 & db[i+1]) >> 40] << 52) |
	           (bit4[(unsigned long)(0x001F000000000000 & db[i+1]) >> 48] << 56) |
		   (bit4[(unsigned long)(0x1F00000000000000 & db[i+1]) >> 56] << 60);
      }
      dbc[i/2] = ((8*i   >= dblen)?    0x00 :   /* out of chars? */
                  bit4[                0x000000000000001F & db[i]])                |
		 ((8*i+1 >= dblen)?    0x00 :
                 (bit4[(unsigned long)(0x0000000000001F00 & db[i])  >>  8] <<  4)) |
		 ((8*i+2 >= dblen)?    0x00 :
	         (bit4[(unsigned long)(0x00000000001F0000 & db[i])  >> 16] <<  8)) |
		 ((8*i+3 >= dblen)?    0x00 :
                 (bit4[(unsigned long)(0x000000001F000000 & db[i])  >> 24] << 12)) |
		 ((8*i+4 >= dblen)?    0x00 :
	         (bit4[(unsigned long)(0x0000001F00000000 & db[i])  >> 32] << 16)) |
		 ((8*i+5 >= dblen)?    0x00 :
	         (bit4[(unsigned long)(0x00001F0000000000 & db[i])  >> 40] << 20)) |
		 ((8*i+6 >= dblen)?    0x00 :
	         (bit4[(unsigned long)(0x001F000000000000 & db[i])  >> 48] << 24)) |
		 ((8*i+7 >= dblen)?    0x00 :
	         (bit4[(unsigned long)(0x1F00000000000000 & db[i])  >> 56] << 28)) |
		 ((8*i+8 >= dblen)?    0x00 :
	         (bit4[                0x000000000000001F & db[i+1]]       << 32)) |
		 ((8*i+9 >= dblen)?    0x00 :
	         (bit4[(unsigned long)(0x0000000000001F00 & db[i+1])>>  8] << 36)) |
		 ((8*i+10>= dblen)?    0x00 :
	         (bit4[(unsigned long)(0x00000000001F0000 & db[i+1])>> 16] << 40)) |
		 ((8*i+11>= dblen)?    0x00 :
		 (bit4[(unsigned long)(0x000000001F000000 & db[i+1])>> 24] << 44)) |
		 ((8*i+12>= dblen)?    0x00 :
	         (bit4[(unsigned long)(0x0000001F00000000 & db[i+1])>> 32] << 48)) |
		 ((8*i+13>= dblen)?    0x00 :
	         (bit4[(unsigned long)(0x00001F0000000000 & db[i+1])>> 40] << 52)) |
                 ((8*i+14>= dblen)?    0x00 :
	         (bit4[(unsigned long)(0x001F000000000000 & db[i+1])>> 48] << 56)) |
		 ((8*i+15>= dblen)?    0x00 : 
                 (bit4[(unsigned long)(0x1F00000000000000 & db[i+1])>> 56] << 60));
      break;
	  
    case 5:
      for(i=0; i<length-3; i+=3)
      {
        dbc[(2*i)/3]   = ((unsigned long)(0x000000000000001F & db[i])   <<  4) |
	                 ((unsigned long)(0x0000000000001F00 & db[i])   <<  1) |
	                 ((unsigned long)(0x00000000001F0000 & db[i])   >>  2) |
		         ((unsigned long)(0x000000001F000000 & db[i])   >>  5) |
	                 ((unsigned long)(0x0000001F00000000 & db[i])   >>  8) |
	                 ((unsigned long)(0x00001F0000000000 & db[i])   >> 11) |
	                 ((unsigned long)(0x001F000000000000 & db[i])   >> 14) |
	                 ((unsigned long)(0x1F00000000000000 & db[i])   >> 17) |
	                 ((unsigned long)(0x000000000000001F & db[i+1]) << 44) |
	                 ((unsigned long)(0x0000000000001F00 & db[i+1]) << 41) |
	                 ((unsigned long)(0x00000000001F0000 & db[i+1]) << 38) |
	                 ((unsigned long)(0x000000001F000000 & db[i+1]) << 35);
		     
        dbc[(2*i)/3+1] = ((unsigned long)(0x0000001F00000000 & db[i+1]) >> 28) |
			 ((unsigned long)(0x00001F0000000000 & db[i+1]) >> 31) |
	                 ((unsigned long)(0x001F000000000000 & db[i+1]) >> 34) |
	                 ((unsigned long)(0x1F00000000000000 & db[i+1]) >> 37) |
	                 ((unsigned long)(0x000000000000001F & db[i+2]) << 24) |
	                 ((unsigned long)(0x0000000000001F00 & db[i+2]) << 21) |
	                 ((unsigned long)(0x00000000001F0000 & db[i+2]) << 18) |
	                 ((unsigned long)(0x000000001F000000 & db[i+2]) << 15) |
	                 ((unsigned long)(0x0000001F00000000 & db[i+2]) << 12) |
	                 ((unsigned long)(0x00001F0000000000 & db[i+2]) <<  9) |
	                 ((unsigned long)(0x001F000000000000 & db[i+2]) <<  6) |
	                 ((unsigned long)(0x1F00000000000000 & db[i+2]) <<  3);
      }
      dbc[(2*i)/3]   = ((8*i   >= dblen)? 0x00 :/* out of chars? */
	               (  (unsigned long)(0x000000000000001F & db[i])   <<  4)) |
		       ((8*i+1 >= dblen)? 0x00 :
	               (  (unsigned long)(0x0000000000001F00 & db[i])   <<  1)) |
		       ((8*i+2 >= dblen)? 0x00 :
	               (  (unsigned long)(0x00000000001F0000 & db[i])   >>  2)) |
		       ((8*i+3 >= dblen)? 0x00 :
	               (  (unsigned long)(0x000000001F000000 & db[i])   >>  5)) |
		       ((8*i+4 >= dblen)? 0x00 :
	               (  (unsigned long)(0x0000001F00000000 & db[i])   >>  8)) |
		       ((8*i+5 >= dblen)? 0x00 :
	               (  (unsigned long)(0x00001F0000000000 & db[i])   >> 11)) |
		       ((8*i+6 >= dblen)? 0x00 :
	               (  (unsigned long)(0x001F000000000000 & db[i])   >> 14)) |
		       ((8*i+7 >= dblen)? 0x00 :
                       (  (unsigned long)(0x1F00000000000000 & db[i])   >> 17)) |
		       ((8*i+8 >= dblen)? 0x00 :
	               (  (unsigned long)(0x000000000000001F & db[i+1]) << 44)) |
		       ((8*i+9 >= dblen)? 0x00 :
		       (  (unsigned long)(0x0000000000001F00 & db[i+1]) << 41)) |
		       ((8*i+10>= dblen)? 0x00 :
	               (  (unsigned long)(0x00000000001F0000 & db[i+1]) << 38)) |
		       ((8*i+11>= dblen)? 0x00 :
                       (  (unsigned long)(0x000000001F000000 & db[i+1]) << 35));
		     
      dbc[(2*i)/3+1] = ((8*i+12>= dblen)? 0x00 :
	               (  (unsigned long)(0x0000001F00000000 & db[i+1]) >> 28)) |
		       ((8*i+13>= dblen)? 0x00 :
	               (  (unsigned long)(0x00001F0000000000 & db[i+1]) >> 31)) |
		       ((8*i+14>= dblen)? 0x00 :
	               (  (unsigned long)(0x001F000000000000 & db[i+1]) >> 34)) |
		       ((8*i+15>= dblen)? 0x00 :
	               (  (unsigned long)(0x1F00000000000000 & db[i+1]) >> 37)) |
		       ((8*i+16>= dblen)? 0x00 :
	               (  (unsigned long)(0x000000000000001F & db[i+2]) << 24)) |
		       ((8*i+17>= dblen)? 0x00 :
	               (  (unsigned long)(0x0000000000001F00 & db[i+2]) << 21)) |
		       ((8*i+18>= dblen)? 0x00 :
	               (  (unsigned long)(0x00000000001F0000 & db[i+2]) << 18)) |
		       ((8*i+19>= dblen)? 0x00 :
	               (  (unsigned long)(0x000000001F000000 & db[i+2]) << 15)) |
		       ((8*i+20>= dblen)? 0x00 :
	               (  (unsigned long)(0x0000001F00000000 & db[i+2]) << 12)) |
		       ((8*i+21>= dblen)? 0x00 :
	               (  (unsigned long)(0x00001F0000000000 & db[i+2]) <<  9)) |
		       ((8*i+22>= dblen)? 0x00 :
	               (  (unsigned long)(0x001F000000000000 & db[i+2]) <<  6)) |
		       ((8*i+23>= dblen)? 0x00 :
                       (  (unsigned long)(0x1F00000000000000 & db[i+2]) <<  3));
      break;
      
    case 6:
      for(i=0; i<length-2; i+=2)
      {
        dbc[i/2] = ((0x0000000000000008 & db[i])?   0x00 : /*  is it an A=1xxx ? */
		   (((unsigned long)(0x0000000000000004 & ~db[i])   >> 1)    |
		    ((unsigned long)(0x0000000000000001 & ~db[i]))))         |
	           ((0x0000000000000080 & db[i])?   0x00 :
		   (((unsigned long)(0x0000000000000040 & ~db[i])   >>  3)   |
		    ((unsigned long)(0x0000000000000010 & ~db[i])   >>  2))) |
	           ((0x0000000000000800 & db[i])?   0x00 :
		   (((unsigned long)(0x0000000000000400 & ~db[i])   >>  5)   |
		    ((unsigned long)(0x0000000000000100 & ~db[i])   >>  4))) |
	           ((0x0000000000008000 & db[i])?   0x00 :
		   (((unsigned long)(0x0000000000004000 & ~db[i])   >>  7)   |
		    ((unsigned long)(0x0000000000001000 & ~db[i])   >>  6))) |
	           ((0x0000000000080000 & db[i])?   0x00 :
		   (((unsigned long)(0x0000000000040000 & ~db[i])   >>  9)   |
		    ((unsigned long)(0x0000000000010000 & ~db[i])   >>  8))) |
	           ((0x0000000000800000 & db[i])?   0x00 :
		   (((unsigned long)(0x0000000000400000 & ~db[i])   >> 11)   |
		    ((unsigned long)(0x0000000000100000 & ~db[i])   >> 10))) |
	           ((0x0000000008000000 & db[i])?   0x00 :
		   (((unsigned long)(0x0000000004000000 & ~db[i])   >> 13)   |
		    ((unsigned long)(0x0000000001000000 & ~db[i])   >> 12))) |
		   ((0x0000000080000000 & db[i])?   0x00 :
	           (((unsigned long)(0x0000000040000000 & ~db[i])   >> 15)   |
		    ((unsigned long)(0x0000000010000000 & ~db[i])   >> 14))) |
	           ((0x0000000800000000 & db[i])?   0x00 :
		   (((unsigned long)(0x0000000400000000 & ~db[i])   >> 17)   |
		    ((unsigned long)(0x0000000100000000 & ~db[i])   >> 16))) |
	           ((0x0000008000000000 & db[i])?   0x00 :
		   (((unsigned long)(0x0000004000000000 & ~db[i])   >> 19)   |
		    ((unsigned long)(0x0000001000000000 & ~db[i])   >> 18))) |
	           ((0x0000080000000000 & db[i])?   0x00 :
		   (((unsigned long)(0x0000040000000000 & ~db[i])   >> 21)   |
		    ((unsigned long)(0x0000010000000000 & ~db[i])   >> 20))) |
	           ((0x0000800000000000 & db[i])?   0x00 :
		   (((unsigned long)(0x0000400000000000 & ~db[i])   >> 23)   |
		    ((unsigned long)(0x0000100000000000 & ~db[i])   >> 22))) |
	           ((0x0008000000000000 & db[i])?   0x00 :
		   (((unsigned long)(0x0004000000000000 & ~db[i])   >> 25)   |
		    ((unsigned long)(0x0001000000000000 & ~db[i])   >> 24))) |
	           ((0x0080000000000000 & db[i])?   0x00 :
		   (((unsigned long)(0x0040000000000000 & ~db[i])   >> 27)   |
		    ((unsigned long)(0x0010000000000000 & ~db[i])   >> 26))) |
	           ((0x0800000000000000 & db[i])?   0x00 :
		   (((unsigned long)(0x0400000000000000 & ~db[i])   >> 29)   |
		    ((unsigned long)(0x0100000000000000 & ~db[i])   >> 28))) |
		   ((0x8000000000000000 & db[i])?   0x00 :
	           (((unsigned long)(0x4000000000000000 & ~db[i])   >> 31)   |
		    ((unsigned long)(0x1000000000000000 & ~db[i])   >> 30))) |
	           ((0x0000000000000008 & db[i+1])? 0x00 :
		   (((unsigned long)(0x0000000000000004 & ~db[i+1]) << 31)   |
		    ((unsigned long)(0x0000000000000001 & ~db[i+1]) << 32))) |
	           ((0x0000000000000080 & db[i+1])? 0x00 :
		   (((unsigned long)(0x0000000000000040 & ~db[i+1]) << 29)   |
		    ((unsigned long)(0x0000000000000010 & ~db[i+1]) << 30))) |
	           ((0x0000000000000800 & db[i+1])? 0x00 :
		   (((unsigned long)(0x0000000000000400 & ~db[i+1]) << 27)   |
		    ((unsigned long)(0x0000000000000100 & ~db[i+1]) << 28))) |
	           ((0x0000000000008000 & db[i+1])? 0x00 :
		   (((unsigned long)(0x0000000000004000 & ~db[i+1]) << 25)   |
		    ((unsigned long)(0x0000000000001000 & ~db[i+1]) << 26))) |
	           ((0x0000000000080000 & db[i+1])? 0x00 :
		   (((unsigned long)(0x0000000000040000 & ~db[i+1]) << 23)   |
		    ((unsigned long)(0x0000000000010000 & ~db[i+1]) << 24))) |
	           ((0x0000000000800000 & db[i+1])? 0x00 :
		   (((unsigned long)(0x0000000000400000 & ~db[i+1]) << 21)   |
		    ((unsigned long)(0x0000000000100000 & ~db[i+1]) << 22))) |
	           ((0x0000000008000000 & db[i+1])? 0x00 :
		   (((unsigned long)(0x0000000004000000 & ~db[i+1]) << 19)   |
		    ((unsigned long)(0x0000000001000000 & ~db[i+1]) << 20))) |
		   ((0x0000000080000000 & db[i+1])? 0x00 :
	           (((unsigned long)(0x0000000040000000 & ~db[i+1]) << 17)   |
		    ((unsigned long)(0x0000000010000000 & ~db[i+1]) << 18))) |
	           ((0x0000000800000000 & db[i+1])? 0x00 :
		   (((unsigned long)(0x0000000400000000 & ~db[i+1]) << 15)   |
		    ((unsigned long)(0x0000000100000000 & ~db[i+1]) << 16))) |
	           ((0x0000008000000000 & db[i+1])? 0x00 :
		   (((unsigned long)(0x0000004000000000 & ~db[i+1]) << 13)   |
		    ((unsigned long)(0x0000001000000000 & ~db[i+1]) << 14))) |
	           ((0x0000080000000000 & db[i+1])? 0x00 :
		   (((unsigned long)(0x0000040000000000 & ~db[i+1]) << 11)   |
		    ((unsigned long)(0x0000010000000000 & ~db[i+1]) << 12))) |
	           ((0x0000800000000000 & db[i+1])? 0x00 :
		   (((unsigned long)(0x0000400000000000 & ~db[i+1]) <<  9)   |
		    ((unsigned long)(0x0000100000000000 & ~db[i+1]) << 10))) |
	           ((0x0008000000000000 & db[i+1])? 0x00 :
		   (((unsigned long)(0x0004000000000000 & ~db[i+1]) <<  7)   |
		    ((unsigned long)(0x0001000000000000 & ~db[i+1]) <<  8))) |
	           ((0x0080000000000000 & db[i+1])? 0x00 :
		   (((unsigned long)(0x0040000000000000 & ~db[i+1]) <<  5)   |
		    ((unsigned long)(0x0010000000000000 & ~db[i+1]) <<  6))) |
	           ((0x0800000000000000 & db[i+1])? 0x00 :
		   (((unsigned long)(0x0400000000000000 & ~db[i+1]) <<  3)   |
		    ((unsigned long)(0x0100000000000000 & ~db[i+1]) <<  4))) |
	           ((0x8000000000000000 & db[i+1])? 0x00 :
	           (((unsigned long)(0x4000000000000000 & ~db[i+1]) <<  1)   |
		    ((unsigned long)(0x1000000000000000 & ~db[i+1]) <<  2)));
      }
      dbc[i/2] = ((8*i   >= dblen)?0x00 :               /* out of chars? */
	         ((0x0000000000000008 & db[i])?  0x00 : /*  is it an A=1xxx ? */
		 (((unsigned long)(0x0000000000000004 & ~db[i])   >>  1)    |
		  (               (0x0000000000000001 & ~db[i])))))         |
		 ((8*i+1 >= dblen)?0x00 :
	         ((0x0000000000000080 & db[i])?  0x00 :
		 (((unsigned long)(0x0000000000000040 & ~db[i])   >>  3)    |
		  ((unsigned long)(0x0000000000000010 & ~db[i])   >>  2)))) |
		 ((8*i+2 >= dblen)?0x00 :
	         ((0x0000000000000800 & db[i])?  0x00 :
		 (((unsigned long)(0x0000000000000400 & ~db[i])   >>  5)    |
		  ((unsigned long)(0x0000000000000100 & ~db[i])   >>  4)))) |
		 ((8*i+3 >= dblen)?0x00 :
	         ((0x0000000000008000 & db[i])?  0x00 :
		 (((unsigned long)(0x0000000000004000 & ~db[i])   >>  7)    |
		  ((unsigned long)(0x0000000000001000 & ~db[i])   >>  6)))) |
		 ((8*i+4 >= dblen)?0x00 :
	         ((0x0000000000080000 & db[i])?  0x00 :
		 (((unsigned long)(0x0000000000040000 & ~db[i])   >>  9)    |
		  ((unsigned long)(0x0000000000010000 & ~db[i])   >>  8)))) |
		 ((8*i+5 >= dblen)?0x00 :
	         ((0x0000000000800000 & db[i])?  0x00 :
		 (((unsigned long)(0x0000000000400000 & ~db[i])   >> 11)    |
		  ((unsigned long)(0x0000000000100000 & ~db[i])   >> 10)))) |
		 ((8*i+6 >= dblen)?0x00 :
	         ((0x0000000008000000 & db[i])?  0x00 :
		 (((unsigned long)(0x0000000004000000 & ~db[i])   >> 13)    |
		  ((unsigned long)(0x0000000001000000 & ~db[i])   >> 12)))) |
		 ((8*i+7 >= dblen)?0x00 :
		 ((0x0000000080000000 & db[i])?  0x00 :
	         (((unsigned long)(0x0000000040000000 & ~db[i])   >> 15)    |
		  ((unsigned long)(0x0000000010000000 & ~db[i])   >> 14)))) |
		 ((8*i+8 >= dblen)?0x00 :
	         ((0x0000000800000000 & db[i])?  0x00 :
		 (((unsigned long)(0x0000000400000000 & ~db[i])   >> 17)    |
		  ((unsigned long)(0x0000000100000000 & ~db[i])   >> 16)))) |
		 ((8*i+9 >= dblen)?0x00 :
	         ((0x0000008000000000 & db[i])?  0x00 :
		 (((unsigned long)(0x0000004000000000 & ~db[i])   >> 19)    |
		  ((unsigned long)(0x0000001000000000 & ~db[i])   >> 18)))) |
		 ((8*i+10>= dblen)?0x00 :
	         ((0x0000080000000000 & db[i])?  0x00 :
		 (((unsigned long)(0x0000040000000000 & ~db[i])   >> 21)    |
		  ((unsigned long)(0x0000010000000000 & ~db[i])   >> 20)))) |
		 ((8*i+11>= dblen)?0x00 :
	         ((0x0000800000000000 & db[i])?  0x00 :
		 (((unsigned long)(0x0000400000000000 & ~db[i])   >> 23)    |
		  ((unsigned long)(0x0000100000000000 & ~db[i])   >> 22)))) |
		 ((8*i+12>= dblen)?0x00 :
	         ((0x0008000000000000 & db[i])?  0x00 :
		 (((unsigned long)(0x0004000000000000 & ~db[i])   >> 25)    |
		  ((unsigned long)(0x0001000000000000 & ~db[i])   >> 24)))) |
		 ((8*i+13>= dblen)?0x00 :
	         ((0x0080000000000000 & db[i])?  0x00 :
		 (((unsigned long)(0x0040000000000000 & ~db[i])   >> 27)    |
		  ((unsigned long)(0x0010000000000000 & ~db[i])   >> 26)))) |
		 ((8*i+14>= dblen)?0x00 :
	         ((0x0800000000000000 & db[i])?  0x00 :
		 (((unsigned long)(0x0400000000000000 & ~db[i])   >> 29)    |
		  ((unsigned long)(0x0100000000000000 & ~db[i])   >> 28)))) |
		 ((8*i+15>= dblen)?0x00 :
		 ((0x8000000000000000 & db[i])?  0x00 :
	         (((unsigned long)(0x4000000000000000 & ~db[i])   >> 31)    |
		  ((unsigned long)(0x1000000000000000 & ~db[i])   >> 30)))) |
		 ((8*i+16>= dblen)?0x00 :
	         ((0x0000000000000008 & db[i+1])?0x00 :
		 (((unsigned long)(0x0000000000000004 & ~db[i+1]) << 31)    |
		  ((unsigned long)(0x0000000000000001 & ~db[i+1]) << 32)))) |
		 ((8*i+17>= dblen)?0x00 :
	         ((0x0000000000000080 & db[i+1])?0x00 :
		 (((unsigned long)(0x0000000000000040 & ~db[i+1]) << 29)    |
		  ((unsigned long)(0x0000000000000010 & ~db[i+1]) << 30)))) |
		 ((8*i+18>= dblen)?0x00 :
	         ((0x0000000000000800 & db[i+1])?0x00 :
		 (((unsigned long)(0x0000000000000400 & ~db[i+1]) << 27)    |
		  ((unsigned long)(0x0000000000000100 & ~db[i+1]) << 28)))) |
		 ((8*i+19>= dblen)?0x00 :
	         ((0x0000000000008000 & db[i+1])?0x00 :
		 (((unsigned long)(0x0000000000004000 & ~db[i+1]) << 25)    |
		  ((unsigned long)(0x0000000000001000 & ~db[i+1]) << 26)))) |
		 ((8*i+20>= dblen)?0x00 :
	         ((0x0000000000080000 & db[i+1])?0x00 :
		 (((unsigned long)(0x0000000000040000 & ~db[i+1]) << 23)    |
		  ((unsigned long)(0x0000000000010000 & ~db[i+1]) << 24)))) |
		 ((8*i+21>= dblen)?0x00 :
	         ((0x0000000000800000 & db[i+1])?0x00 :
		 (((unsigned long)(0x0000000000400000 & ~db[i+1]) << 21)    |
		  ((unsigned long)(0x0000000000100000 & ~db[i+1]) << 22)))) |
		 ((8*i+22>= dblen)?0x00 :
	         ((0x0000000008000000 & db[i+1])?0x00 :
		 (((unsigned long)(0x0000000004000000 & ~db[i+1]) << 19)    |
		  ((unsigned long)(0x0000000001000000 & ~db[i+1]) << 20)))) |
		 ((8*i+23>= dblen)?0x00 :
		 ((0x0000000080000000 & db[i+1])?0x00 :
	         (((unsigned long)(0x0000000040000000 & ~db[i+1]) << 17)    |
		  ((unsigned long)(0x0000000010000000 & ~db[i+1]) << 18)))) |
		 ((8*i+24>= dblen)?0x00 :
	         ((0x0000000800000000 & db[i+1])?0x00 :
		 (((unsigned long)(0x0000000400000000 & ~db[i+1]) << 15)    |
		  ((unsigned long)(0x0000000100000000 & ~db[i+1]) << 16)))) |
		 ((8*i+25>= dblen)?0x00 :
	         ((0x0000008000000000 & db[i+1])?0x00 :
		 (((unsigned long)(0x0000004000000000 & ~db[i+1]) << 13)    |
		  ((unsigned long)(0x0000001000000000 & ~db[i+1]) << 14)))) |
		 ((8*i+26>= dblen)?0x00 :
	         ((0x0000080000000000 & db[i+1])?0x00 :
		 (((unsigned long)(0x0000040000000000 & ~db[i+1]) << 11)    |
		  ((unsigned long)(0x0000010000000000 & ~db[i+1]) << 12)))) |
		 ((8*i+27>= dblen)?0x00 :
	         ((0x0000800000000000 & db[i+1])?0x00 :
		 (((unsigned long)(0x0000400000000000 & ~db[i+1]) <<  9)    |
		  ((unsigned long)(0x0000100000000000 & ~db[i+1]) << 10)))) |
		 ((8*i+28>= dblen)?0x00 :
	         ((0x0008000000000000 & db[i+1])?0x00 :
		 (((unsigned long)(0x0004000000000000 & ~db[i+1]) <<  7)    |
		  ((unsigned long)(0x0001000000000000 & ~db[i+1]) <<  8)))) |
		 ((8*i+29>= dblen)?0x00 :
	         ((0x0080000000000000 & db[i+1])?0x00 :
		 (((unsigned long)(0x0040000000000000 & ~db[i+1]) <<  5)    |
		  ((unsigned long)(0x0010000000000000 & ~db[i+1]) <<  6)))) |
		 ((8*i+30>= dblen)?0x00 :
	         ((0x0800000000000000 & db[i+1])?0x00 :
		 (((unsigned long)(0x0400000000000000 & ~db[i+1]) <<  3)    |
		  ((unsigned long)(0x0100000000000000 & ~db[i+1]) <<  4)))) |
                 ((8*i+31>= dblen)?0x00 :
                 ((0x8000000000000000 & db[i+1])?0x00 :
	         (((unsigned long)(0x4000000000000000 & ~db[i+1]) <<  1)    |
		  ((unsigned long)(0x1000000000000000 & ~db[i+1]) <<  2))));
      break;
#endif /* little-endian */

    default: 
      fprintf(stderr, "cb_compress: Invalid mode parameter %d, returning...\n", mode);
      return;
  }
#endif /* 64-bit */
}


/*
cb_compress(3B)                                           Last changed: 09-17-02

NAME
        cb_compress - Compresses nucleotide or amino acid ASCII data

SYNOPSIS

        C/C++:

        #include <cbl.h>
        void cb_compress( long *db, long *dbc, long dblen, long mode);

        Fortran:

        use cb_comp
        call cb_compress( db, dbc, dblen, mode)


IMPLEMENTATION

        Cray SV1 series UNICOS systems

DESCRIPTION

        cb_compress compresses nucleotide or amino acid ASCII data. The input 
        array is assumed to be a sequence of upper case or lower case ASCII 
        letters (A..Z or a..z).  Each letter is translated into a code with
        the number of bits determined by the mode parameter. Lower case
        letters are translated into the same code as the corresponding upper
        case letters.  The codes are packed into the output array as specified
        below. A special mode, 6, is supplied to convert data packed in 4-bit
        fields into data packed in 2-bit fields.
        
        db      (input) If mode = 2, 4, or 5, db contains ASCII text, packed 
                8 letters/word. Length of db is (dblen+7)/8 (64-bit words). 
                if mode = 6, db contains letters packed in the 4-bit encoding
                scheme described below, packed 16 letters/word. Length of db 
                is (dblen+15)/16 (64-bit words). The final word may be filled 
                with trailing nulls. In Fortran, db should be an INTEGER(8)
                array.
                    
        dbc     (output) packed codes. The number of bits/code and the amount
                of memory required depends on the mode value. The dbc array 
                must be allocated before calling cb_compress. In Fortran, 
                dbc should be an INTEGER(8) array.
                    
        dblen   (input) Number of letters represented in db. In Fortran, 
                dblen should be an INTEGER(8) variable, constant, or
                expression.
            
        mode    (input) Allowed values are 2, 4, 5, and 6. In Fortran, mode
                should be an INTEGER(8) variable, constant, or sxpression.
                    
                mode = 2:  Each ASCII letter is converted to a 2-bit code as 
                        follows:
               
                        00 <- A  (also H, I, P, Q, X, Y)
               
                        01 <- C  (also B, J, K, R, S, Z)
               
                        10 <- T  (also D, E, L, M, U)
               
                        11 <- G  (also F, N, O, V, W)
               
                        Length of dbc must be >= (dblen+31)/32 (64-bit words).
               
               
                mode = 4: Each ASCII letter is converted to a 4-bit code as 
                        follows:
               
                        1000 <- A
                        0100 <- C
                        0010 <- G
                        0001 <- T, U
                        0111 <- B = C or G or T
                        1011 <- D = A or G or T
                        1101 <- H = A or C or T
                        0011 <- K = G or T
                        1100 <- M = A or C
                        1111 <- N = A or C or G or T
                        1010 <- R = A or G
                        0110 <- S = C or G
                        1110 <- V = A or C or G
                        1001 <- W = A or T
                        0101 <- Y = C or T
                        1111 <- X (same as N)
                        1111 <- All other letters
               
                        Length of dbc must be >= (dblen+15)/16 (64-bit words).
               
               
                mode = 5:  Each ASCII letter is converted to a 5-bit code 
                        containing the integer sequence number of  the letter 
                        in the alphabet:
               
                        {A,B,C,...,Z} -> {1,2,3,...,26}
               
                        {a,b,c,...,z} -> {1,2,3,...,26}
               
                        Packing is 12 letters in each 64-bit word, with the 
                        left 4 bits set to zero followed by 5*12 = 60 bits of 
                        data.
               
                        Length of dbc must be >= (dblen+11)/12 (64-bit words).

                mode = 6: Each input letter is assumed to be encoded using the
                        4-bit code described under mode=4 above. The result
                        is the 2-bit code described under mode=2 above.

                        Length of dbc must be >= (dblen+31)/32 (64-bit words).

NOTES       
        cb_compress is single-threaded (i.e. not tasked) and may be called 
        from within a parallel region.  Before cb_compress is called in a       
        parallel region with mode=4, it must have been previously called
        outside a parallel region with mode=4 to initialize internal tables.
        This initial call can be made with dblen=0.
    
        cb_compress replaces the contents of the bmm register.
               
SEE ALSO

        cb_uncompress(3B), INTRO_LIBCBL(3B)
*/
