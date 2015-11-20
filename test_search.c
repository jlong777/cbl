/* test_search.c                                 http://cbl.sourceforge.net
 * 
 * test cb_searchn for random strings that incrementally cross word 
 * boundries
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
 * $Id: test_search.c,v 1.4 2003/09/25 00:53:27 jlong777 Exp $
 */
 
#include "cbl.h"
#include "cb_macro.h"
#include <ctype.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/times.h>
#include <time.h>

#define REP   40  /* max threshold test size, also number of benchmark reps */
#define SIZE 240  /* change size of test cases */

void slow_search(char *db, int dblen, char *test, int testlen, long *sfound,
                 long *sfoundlen, int threshold, long *scount);

int main()
{
  int i, j, k, l, error=0, mode, n, t, lastthreshold;
  long in[SIZE], test[SIZE], *in2;
  char lastin[SIZE], lasttest[SIZE];
  char ch2[4], version[6];
  long count[SIZE], scount[SIZE], found[SIZE], sfound[SIZE], 
        pack[SIZE],  packt[SIZE], *count2, *found2, *pack2;
  long lastpack[SIZE],  lastpackt[SIZE];
  long foundlen=SIZE, sfoundlen=SIZE, tempc, tempf;
  
  /* timing stuff */
  struct tms mytime;               /* struct filled by calls to times() */
  float users, userf, syss, sysf, utot, stot; /* start and finish times */
  clock_t stime, ftime;
  long totsec = 0;
  
  ch2[0] = 'A';
  ch2[1] = 'C';
  ch2[2] = 'G';
  ch2[3] = 'T';

  printf("running test of cb_searchn...\n");
  
  mode = 2;
  /* test legal thresholds from 0 to REP */
  for(t=0; t<REP+1; t++)
  {
    /* test random strings up to size SIZE-1 */
    for(i=1; i<SIZE-1; i++) /* dblen must be > 0 */
    {
      if(i <= t) continue; /* strlen(in) must be > threshold */
      
strcpy(lastin, (char *)in); /* save last in */
      
      for(j=0; j<i; j++)
        ((char *)in)[j] = ch2[rand()%4];
      ((char *)in)[j] = 0x0;
    
      for(j=1; j<i; j++) /* j will be length of the test string > 0 & < i */
      {
        /* manufacture test string from random substring of "in"
         * this insures we will have at least one hit each time
         */
strcpy(lasttest, (char *)test);/* save last test */
        strncpy((char *)test, ((char *)in)+(rand()%(i-j)), j);
       ((char *)test)[j] = 0x0;

        for(k=0; k<SIZE; k++) /* initialize arrays */
        {
          found[k] = sfound[k] = count[k] = scount[k] = 0;
lastpack[k] = pack[k]; lastpackt[k] = packt[k]; 
        }
      
	if(j <= t) continue; /* testlen must be > t (& > 0) */
lastthreshold = t;

        foundlen = sfoundlen = SIZE;
	
        cb_compress(in, pack, i, mode);

        cb_compress(test, packt, j, mode);

        cb_searchn(pack, i, packt, j, found, &foundlen, t, count);

        slow_search((char *)in, i, (char *)test, j, sfound, &sfoundlen, t, scount);

        /* bubble sort outputs if needed */
	for(k=0; k<foundlen-1; k++)
        {
	  for(l=k+1; l<foundlen; l++)
	    if(found[l]<found[k])
	    {
	      tempf = found[l];
	      tempc = count[l];
	      found[l] = found[k];
	      count[l] = count[k];
	      found[k] = tempf;
	      count[k] = tempc;
	    }
	}

        if(foundlen != sfoundlen)
        {
          if(DEBUG_INFO)
          {
            printf("arg!!************************\n");
            printf("foundlen = %d  sfoundlen = %d\n", foundlen, sfoundlen);
            printf("in = %s  test=%s\n\n", (char *)in, (char *)test);
          }
          error = 1;
        }
      
        for(k=0; k<SIZE; k++)
        {
	  if(count[k] != scount[k])
	  {
            if(DEBUG_INFO)                                                                  {
	      printf("count[%d] = %d  scount[%d] = %d\n",k,count[k],k,scount[k]);
	      printf("threshold = %d\n", t);
	      printf("in = %s  test=%s\n\n", (char *)in, (char *)test);
printf("lastin = %s  lasttest=%s\n\n", lastin, lasttest);
printf("threshold = %d lastthreshold = %d\n", t, lastthreshold);
printf("in = ");
for(k=0; k<SIZE/8; k++) printf("x%X",in[k]);
printf("\n");
printf("test = ");
for(k=0; k<SIZE/8; k++) printf("x%X",test[k]);
printf("\n");
printf("lastpack = ");
for(k=0; k<SIZE; k++) printf("x%X",lastpack[k]);
printf("\n");
printf("lastpackt = ");
for(k=0; k<SIZE; k++) printf("x%X",lastpackt[k]);
printf("\n");
printf("found = ");
for(k=0; k<SIZE; k++) printf("x%X",found[k]);
printf("\n");
printf("count = ");
for(k=0; k<SIZE; k++) printf("x%X",count[k]);
printf("\n");
exit(0);
            }
            error = 1;
  	  }
	  if(found[k]!=sfound[k])
          {
            if(DEBUG_INFO)
            {
	      printf("found[%d] = %d  sfound[%d] = %d\n",k,found[k],k,sfound[k]);
	      printf("threshold = %d\n", t);
	      printf("in = %s  test=%s\n\n\n", (char *)in, (char *)test);
printf("lastin = %s  lasttest=%s\n\n", lastin, lasttest);
printf("threshold = %d lastthreshold = %d\n", t, lastthreshold);
printf("in = ");
for(k=0; k<SIZE/8; k++) printf("x%X",in[k]);
printf("\n");
printf("test = ");
for(k=0; k<SIZE/8; k++) printf("x%X",test[k]);
printf("\n");
printf("lastpack = ");
for(k=0; k<SIZE; k++) printf("x%X",lastpack[k]);
printf("\n");
printf("lastpackt = ");
for(k=0; k<SIZE; k++) printf("x%X",lastpackt[k]);
printf("\n");
printf("found = ");
for(k=0; k<SIZE; k++) printf("x%X",found[k]);
printf("\n");
printf("count = ");
for(k=0; k<SIZE; k++) printf("x%X",count[k]);
printf("\n");
exit(0);
            }
            error = 1;
	  }
        }
      }
    }
  }

  if(error) printf("**failure**\n");
  else      printf("**success**\n");
  
  /* =================================================================== */
  /* benchmark */
  if(BMARK)
  {
    printf("times() in clock ticks to do %d cb_searchns:\n", REP);
  
    utot = stot = 0.0;
    
    printf("length of string to find = %d\n", strlen((char*)test));
  
    for(i=0; i<18; i++)
    {
      n = 0x00000100 << i;
      
      in2 = (long *) calloc(n/(sizeof(long))+9, sizeof(long));
      if(in2==NULL)
      {
        printf("test_search: calloc error, exiting...\n");
        exit(0);
      }
      
      count2 = (long *) calloc(n/(sizeof(long))+9, sizeof(long));
      if(count2==NULL)
      {
        printf("test_search: calloc error, exiting...\n");
        exit(0);
      }
  
      found2 = (long *) calloc(n/(sizeof(long))+9, sizeof(long));
      if(found2==NULL)
      {
        printf("test_search: calloc error, exiting...\n");
        exit(0);
      }
  
      pack2 = (long *) calloc(n/(sizeof(long))+9, sizeof(long));
      if(pack2==NULL)
      {
        printf("test_search: calloc error, exiting...\n");
        exit(0);
      }
      
      /* create random string */
      for(j=0; j<n; j++)
        ((char *)in2)[j] = ch2[rand()%4];
      ((char *)in2)[j] = 0x0;

      printf("n = %d\n", n);
      
      cb_compress(in2, pack2, n, mode);
      
      j = strlen((char*)test);

      stime = clock()/CLOCKS_PER_SEC;
      times(&mytime);
      users = (float) mytime.tms_utime;
      syss  = (float) mytime.tms_stime;
      
      /* search REP times */
      for(k=0; k<REP; k++)
      {
        foundlen = n;
        cb_searchn(pack2, n, packt, j, found2, &foundlen, k, count2);
      }
	
      times(&mytime);
      userf = (float) mytime.tms_utime;
      sysf  = (float) mytime.tms_stime;
      ftime = clock()/CLOCKS_PER_SEC;
      
      printf("mode = %d user = %7.1f  sys = %7.1f  secs (approx) = %d\n", 
	      mode, userf-users, sysf-syss, ftime-stime);
	
      utot += userf-users;
      stot += sysf-syss;
      totsec += ftime-stime;
      
      free(in2);
      free(count2);
      free(found2);
      free(pack2);
    }
    printf("total clock ticks:user = %9.1f sys = %9.1f\n", utot, stot);
    printf("total approx secs = %d for searchn\n\n", totsec);
  
  }
  
  cb_version(version);
  printf("CBL version: %s\n\n", version);
  
  return 0;
}

void slow_search(char *db, int dblen, char *test, int testlen, long *sfound,
                 long *sfoundlen, int threshold, long *scount)
{
  int i, j, k=0, num;
  long foundlen_in;
  
  foundlen_in = *sfoundlen;
  *sfoundlen = 0;
  
  for(i=0; i<dblen-testlen+1; i++)
  {
    num = 0;
    for(j=0; j<testlen; j++)
    {
      if(test[j] != db[i+j]) num++;
      if(num > threshold) break;
    }
    if(num > threshold) continue;
    
    if(k < foundlen_in)
    {
      sfound[k] = i;
      scount[k++] = num;
      *sfoundlen = k;
    }
    else
    {
      *sfoundlen *= -1;
      break;
    }
  }
}
