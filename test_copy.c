/* test_copy.c                                   http://cbl.sourceforge.net
 *
 * test and benchmark cb_copy_bits
 *
 * test:
 * generate 2 random equal length strings up to 10 words long, and copy bits
 * from second string to first. copy_bits will be invoked for each bit position
 * of first string and number of bits to be copied will be determined by a loop 
 * over second string offset that copies the lesser of (length of src - offset)
 * or (length of dest - offset) bits.
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
 * $Id: test_copy.c,v 1.5 2004/02/04 20:00:12 jlong777 Exp $
 */

#include "cbl.h"
#include "cb_macro.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/times.h>
#include <time.h>

#define REP 200  /* repetitions for benchmark */

int compare_bits(long *, long, long *, long, long);

int main()
{
  int i, j, error=0, len, wsz;
  long k, nbits;
  char ch[52], temp[137], version[6];
  long *first, *second; /* destination and source strings */
  
  /* timing stuff */
  struct tms mytime;               /* struct filled by calls to times() */
  float users, userf, syss, sysf, utot=0.0, stot=0.0; /* start and finish */
  clock_t stime, ftime;
  long totsec = 0;
  
  ch[0]  = 'A';  ch[13] = 'N';  ch[26] = 'a';  ch[39] = 'n';
  ch[1]  = 'B';  ch[14] = 'O';  ch[27] = 'b';  ch[40] = 'o';
  ch[2]  = 'C';  ch[15] = 'P';  ch[28] = 'c';  ch[41] = 'p';
  ch[3]  = 'D';  ch[16] = 'Q';  ch[29] = 'd';  ch[42] = 'q';
  ch[4]  = 'E';  ch[17] = 'R';  ch[30] = 'e';  ch[43] = 'r';
  ch[5]  = 'F';  ch[18] = 'S';  ch[31] = 'f';  ch[44] = 's';
  ch[6]  = 'G';  ch[19] = 'T';  ch[32] = 'g';  ch[45] = 't';
  ch[7]  = 'H';  ch[20] = 'U';  ch[33] = 'h';  ch[46] = 'u';
  ch[8]  = 'I';  ch[21] = 'V';  ch[34] = 'i';  ch[47] = 'v';
  ch[9]  = 'J';  ch[22] = 'W';  ch[35] = 'j';  ch[48] = 'w';
  ch[10] = 'K';  ch[23] = 'X';  ch[36] = 'k';  ch[49] = 'x';
  ch[11] = 'L';  ch[24] = 'Y';  ch[37] = 'l';  ch[50] = 'y';
  ch[12] = 'M';  ch[25] = 'Z';  ch[38] = 'm';  ch[51] = 'z';
  
  wsz = sizeof(long);

  printf("running test of cb_copy_bits...\n");
  
  for(i=0; i<17*wsz; i++) /* random strings up to 17 words long */
  {
    first  = (long *) calloc(i+4, 1);
    second = (long *) calloc(i+4, 1);
    if(first==NULL || second==NULL)
    {
      if(DEBUG_INFO)
        printf("test_copy: calloc error, exiting...\n");
      error = 1;
      exit(0);
    }
     
    /* create random text strings */
    for(j=0; j<i; j++)
    {
      ((char *)first)[j] = ch[rand()%52];
      ((char *)second)[j]= ch[rand()%52];
    }
    ((char *)first)[j] = ((char *)second)[j] = 0x0;

    /* i=number of chars, so j will =number of bits */
    for(j=0; j<8*i; j++) /* loop over every bit position in first */
    {
      strcpy(temp, (char *)first); len = strlen((char *)first);
      
      for(k=0; k<i; k++) /* loop over every bit position in second */
      {
	nbits = ((8*strlen((char *)second)-k) < (8*strlen((char *)first)-j))? 
	         (8*strlen((char *)second)-k):(8*strlen((char *)first)-j);
      
	cb_copy_bits(first, j, second, k, nbits);
        if(compare_bits(first, j, second, k, nbits))
	{
          if(DEBUG_INFO)
          {
	    printf("dest offset = %d of \"%s\", copy %d bits from\n", j, temp, nbits);
	    printf(" src offset = %d of \"%s\" (bit length = %d)\n\n", k, (char *)second, 8*len);
          }
          error = 1;
          strcpy((char *)first, temp);
	  break;
	}
	strcpy((char *)first, temp);
      }
    }
    free(first);
    free(second);
  }

  if(error) printf("**failure**\n");
  else      printf("**success**\n");

  if(BMARK) /* benchmark */
  {
    printf("times() in clock ticks to do %d cb_copy_bits:\n", REP);
    
    for(i=0; i<18; i++)
    {
    k = 0x00000100 << i;
    
    first  = (long *) calloc(k+1, 1);
    second = (long *) calloc(k+1, 1);
    if(first==NULL || second==NULL)
    {
      printf("test_copy: calloc error, exiting...\n");
      exit(0);
    }
    
    /* create random text strings */
    for(j=0; j<k; j++)
    {
      ((char *)second)[j]= ch[rand()%52];
    }
    ((char *)first)[j] = ((char *)second)[j] = 0x0;
    
    printf("n = %d\n", k);
    
    nbits = 8*k;
    
    stime = clock()/CLOCKS_PER_SEC;
    times(&mytime);
    users = (float) mytime.tms_utime;
    syss  = (float) mytime.tms_stime;
    
    for(j=0; j<REP; j++)
    {
      cb_copy_bits(first, 20, second, 34, nbits);
    }
    
    times(&mytime);
    userf = (float) mytime.tms_utime;
    sysf  = (float) mytime.tms_stime;
    ftime = clock()/CLOCKS_PER_SEC;
    
    printf("user = %7.1f  sys = %7.1f  secs (approx) = %d\n", 
	    userf-users, sysf-syss, ftime-stime);
    utot += userf-users;
    stot += sysf-syss;
    totsec += ftime-stime;
    
    free(first);
    free(second);
  }
  printf("total clock ticks:user = %9.1f sys = %9.1f\n", utot, stot);
  printf("total approx secs = %d for bit copy\n\n", totsec);
  }

  cb_version(version);
  printf("CBL version: %s\n\n", version);
  
  return 0;
}

/* function to compare memory bit-by-bit between two buffers where copied bit
 * string is not necessarily word or byte aligned */

int compare_bits(long *dest, long doffset, long *src, long soffset, long nbits)
{
  long bit=0, wordsz, mask;
  
#ifdef LONG32

  wordsz = 32;
  
#ifdef B_ENDIAN
    mask = 0x80000000;
#endif
#ifdef L_ENDIAN
    mask = 0x00000001;
#endif
  
#endif /* 32-bit */

#ifdef LONG64

  wordsz = 64;
  
#ifdef B_ENDIAN
    mask = 0x8000000000000000;
#endif
#ifdef L_ENDIAN
    mask = 0x0000000000000001;
#endif

#endif /* 64-bit */

#ifdef B_ENDIAN
  while(bit < nbits) /* compare leftmost bit each time */
  {
    if((mask & (dest[(doffset+bit)/wordsz] << (doffset+bit)%wordsz)) !=
       (mask & ( src[(soffset+bit)/wordsz] << (soffset+bit)%wordsz)))
    {
      printf("copy_bits error at bit %d of %d, exiting further comparison.\n", bit+1, nbits);
      return 1;
    }
    bit++;
  }
#endif
#ifdef L_ENDIAN
  while(bit < nbits) /* compare rightmost bit each time */
  {
    if((mask & (dest[(doffset+bit)/wordsz] >> (doffset+bit)%wordsz)) !=
       (mask & ( src[(soffset+bit)/wordsz] >> (soffset+bit)%wordsz)))
    {
      printf("copy_bits error at bit %d of %d, exiting further comparison.\n", bit+1, nbits);
      return 1;
    }
    bit++;
  }
#endif
  return 0;
}

