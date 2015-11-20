/* test_repeat.c                                 http://cbl.sourceforge.net
 * 
 * test cb_repeatn for random strings 
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
 * $Id: test_repeat.c,v 1.4 2003/09/25 00:53:27 jlong777 Exp $
 */
 
#include "cbl.h"
#include "cb_macro.h"
#include <ctype.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/times.h>
#include <time.h>

#define SIZE 200  /* change max size of test cases */

void slow_repeat(char *db, int dblen, int rep_len, int min_rep, char (*srpattern)[17],
                 long *srfound, long *srfoundlen, long *srcount);
void slow_search(char *db, int dblen, char *test, int testlen, long *sfound,
                 long *sfoundlen, int threshold, long *scount);
void write_word(unsigned long word, int);

int main()
{
  int i, j, k, l, error=0, mode, n, t, insert, rnd;
  long in[SIZE], cat[SIZE/4], *in2;
  char test[SIZE];
  char ch2[4], srpattern[SIZE][17], version[6], tempp[17];
  long srcount[SIZE], found[SIZE], srfound[SIZE], pack[SIZE], *found2, *pack2;
  long foundlen=SIZE, srfoundlen=SIZE, tempc, tempf;
#ifdef LONG32
  long long pattern[SIZE], *pattern2, templ;
#endif
#ifdef LONG64
  long pattern[SIZE], *pattern2, templ;
#endif
  /* timing stuff */
  struct tms mytime;               /* struct filled by calls to times() */
  float users, userf, syss, sysf, utot, stot; /* start and finish times */
  clock_t stime, ftime;
  long totsec = 0;
  
  ch2[0] = 'A';
  ch2[1] = 'C';
  ch2[2] = 'G';
  ch2[3] = 'T';

  printf("running test of cb_repeatn...\n");
  
  mode = 2;
  for(j=2; j<17; j++) /* j will be length of the pattern (repeat_len) */
  {
    /* test min_repeats from 2 to 16 */
    for(t=2; t<17; t++)
    {
      /* test random strings up to size SIZE-1 */
      for(i=1; i<SIZE-1; i++)
      {
        if(i <= j*t) continue; /* dblen must be > repeat_len*min_repeats */
	if(j*t > 32) continue; /* min_rep is internally modified so that j*t<=32 */
    
        for(k=0; k<i; k++)
          ((char *)in)[k] = ch2[rand()%4];
        ((char *)in)[k] = 0x0;
	((char *)cat)[0]= 0x0;
      
        /* manufacture pattern string from random substring of "in",
	 * replicate t times, and randomly re-insert into "in" and "pack".
         * this insures we will have at least one hit each time
         */
        strncpy(test, ((char *)in)+(rand()%(i-j)), j);
        test[j] = 0x0; /* pattern */
	
	/* replicate t times */
        for(k=0; k<t; k++)
	  strcat((char *)cat, test);
	  
	for(k=0; k<SIZE; k++) /* initialize arrays */
        {
          found[k] = srfound[k] = srcount[k] = 0;
        }
	
	foundlen = srfoundlen = SIZE;
	
	rnd = ((i-j*t)>1)?(rand()%((i-j*t))):0;
	
	insert = 8*(rnd/8);
	/* reinsert into in */
	cb_copy_bits(in, insert, cat, 0, 8*j*t);
	
        cb_compress(in, pack, strlen((char *)in), mode);

        cb_repeatn(pack, i, j, t, pattern, found, &foundlen);
        slow_repeat((char *)in, strlen((char *)in), strlen(test), t, srpattern,
	             srfound, &srfoundlen, srcount);

        /* bubble sort the outputs */
	for(k=0; k<foundlen-1; k++)
        {
	  for(l=k+1; l<foundlen; l++)
	    if(found[l]<found[k])
	    {
	      tempf = found[l];
	      templ = pattern[l];
	      found[l] = found[k];
	      found[k] = tempf;
	      pattern[l] = pattern[k];
	      pattern[k] = templ;
	    }
	}
	
	for(k=0; k<srfoundlen-1; k++)
        {
	  for(l=k+1; l<srfoundlen; l++)
	    if(srfound[l]<srfound[k])
	    {
	      tempf = srfound[l];
	      tempc = srcount[l];
	      strcpy(tempp, srpattern[l]);
	      srfound[l] = srfound[k];
	      srfound[k] = tempf;
	      srcount[l] = srcount[k];
	      srcount[k] = tempc;
	      strcpy(srpattern[l], srpattern[k]);
	      strcpy(srpattern[k], tempp);
	    }
	}

        if(foundlen != srfoundlen && foundlen>=0)
        {
          if(DEBUG_INFO)
          {
            printf("arg!!************************\n");
            printf("foundlen = %d  srfoundlen = %d\n", foundlen, srfoundlen);
            printf("in = %s  test=%s, cat=%s inserted at = %d\n\n", (char *)in, test, (char *)cat, rnd/8);
          }
          error = 1;
        }
	
        for(k=0; k<SIZE; k++)
        {
	  if(found[k]==0 && srfound[k]==0) break;
	  if(found[k]!=srfound[k])
          {
            if(DEBUG_INFO)
            {
	      printf("found[%d] = %d srfound[%d] = %d srcount[%d] = %d\n",
	              k,found[k],k,srfound[k],k,srcount[k]);
	      printf("min_repeats = %d\n", t);
	      printf("in = %s  srpattern[%d]=%s\n", (char *)in, k, srpattern[k]);
	      printf("cbpattern[%d] = ",k);
#ifdef LONG32
	      write_word(pattern[k]>>30, j);
#endif
#ifdef LONG64
	      write_word(pattern[k]<<2, j);
#endif
	      printf("\n");
            }
            error = 1;
	  }
#ifdef B_ENDIAN
	  if((found[k]==srfound[k]) &&((0x3FFFFFFF & pattern[k]) != srcount[k]))
	  {
            if(DEBUG_INFO)
            {
  	      printf("count mismatch: srcount[%d] = %d\n",k,srcount[k]);
  	      printf("                cbcount[%d] = %d\n",k,0x3FFFFFFF & pattern[k]);
  	      printf("found[%d] = %d srfound[%d] = %d\n",k,found[k],k,srfound[k]);
  	      printf("min_repeats = %d\n", t);
	      printf("in = %s  srpattern[%d]=%s\n", (char *)in, k, srpattern[k]);
	      printf("cbpattern[%d] = ",k);
#ifdef LONG32
	      write_word(pattern[k]>>30, j);
#endif
#ifdef LONG64
	      write_word(pattern[k]<<2, j);
#endif
	      printf("\n");
            }
            error = 1;
	  }
#endif
#ifdef L_ENDIAN
	  if((found[k]==srfound[k]) &&((0x3FFFFFFF & (pattern[k]>>34)) != srcount[k]))
	  {
            if(DEBUG_INFO)
            {
	      printf("count mismatch: srcount[%d] = %d\n",k,srcount[k]);
	      printf("                cbcount[%d] = %d\n",k,0x3FFFFFFF & (pattern[k]>>34));
	      printf("found[%d] = %d srfound[%d] = %d\n",k,found[k],k,srfound[k]);
	      printf("min_repeats = %d\n", t);
	      printf("in = %s  srpattern[%d]=%s\n", (char *)in, k, srpattern[k]);
	      printf("cbpattern[%d] = ",k);
#ifdef LONG32
	      write_word(pattern[k]>>30, j);
#endif
#ifdef LONG64
	      write_word(pattern[k]<<2, j);
#endif
	      printf("\n");
            }
            error = 1;
	  }
#endif
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
    printf("times() in clock ticks to do 15 cb_repeatns:\n");
  
    utot = stot = 0.0;
  
    for(i=0; i<18; i++)
    {
      n = 0x00000100 << i;
      
      in2 = (long *) calloc(n/(sizeof(long))+9, sizeof(long));
      if(in2==NULL)
      {
        printf("test_search: calloc error, exiting...\n");
        exit(0);
      }
#ifdef LONG32
      pattern2 = (long long *) calloc(n/(sizeof(long long))+9, sizeof(long long));
#endif
#ifdef LONG64
      pattern2 = (long *) calloc(n/(sizeof(long))+9, sizeof(long));
#endif
      if(pattern2==NULL)
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

      stime = clock()/CLOCKS_PER_SEC;
      times(&mytime);
      users = (float) mytime.tms_utime;
      syss  = (float) mytime.tms_stime;
      
      /* search 15 times */
      for(k=2; k<17; k++)
      {
        foundlen = n;
	cb_repeatn(pack2, n, k, 2, pattern2, found2, &foundlen);
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
      free(pattern2);
      free(found2);
      free(pack2);
    }
    printf("total clock ticks:user = %9.1f sys = %9.1f\n", utot, stot);
    printf("total approx secs = %d for repeatn\n\n", totsec);
    
  }

  cb_version(version);
  printf("CBL version: %s\n\n", version);
  
  return 0;
}

void slow_repeat(char *db, int dblen, int rep_len, int min_rep, char (*srpattern)[17],
                 long *srfound, long *srfoundlen, long *srcount)
{
  int i, j, k, m, reject;
  char pattern[17]={0x0}, tandem[33];
  /*long foundlen_in;*/
  long sfound[1], sfoundlen, scount[1];
  
  if((rep_len < 2) || (rep_len > 16))
  {
    fprintf(stderr, "slow_repeat: rep_len range error, returning...\n");
    *srfoundlen = -1;
    return;
  }
  
  if(*srfoundlen <= 0)
  {
    fprintf(stderr, "slow_repeatn: num_found range error, returning...\n");
    *srfoundlen = -1;
    return;
  }
  
  if(dblen <= rep_len*min_rep)
  {
    fprintf(stderr, "slow_repeatn: dblen <= rep_len*min_rep, returning...\n");
    *srfoundlen = -1;
    return;
  }
  
  *srfoundlen = 0;
  
  if(rep_len*min_rep > 32) 
    min_rep = 32/rep_len;
    
  for(i=0; i<=dblen-rep_len*min_rep; i++)
  {
    /* get candidate pattern & replicate it */
    strncpy(pattern, db+i, rep_len);

    tandem[0] = 0x0;
    for(j=0; j<min_rep; j++)
      strcat(tandem, pattern);

    /* look for match and filter */
    k = i;
    while(k<=dblen-rep_len*min_rep)
    {
      sfoundlen = 1; scount[0] = sfound[0] = 0; 
      slow_search(db+k, dblen-k, tandem, rep_len*min_rep, sfound, &sfoundlen, 0, scount);

      if(sfoundlen!=0) /* found a candidate */
      {
        /* is the pattern itself a repeat? */
	reject = 0;
	switch(rep_len)
	{
	  case 2: if(pattern[0] == pattern[1]) reject = 1; break; /* 1s */
	  case 3: if((pattern[0]==pattern[1]) && (pattern[0]==pattern[2]))
	            reject = 1; /* 1s */
		  break;
	  case 4: if((pattern[0]==pattern[1]) && (pattern[0]==pattern[2]) &&
	             (pattern[0]==pattern[3]))
	            reject = 1; /* 1s */
		  if((pattern[0]==pattern[2]) && (pattern[1]==pattern[3]))
		    reject = 1; /* 2s */
		  break;
	  case 5: if((pattern[0]==pattern[1]) && (pattern[0]==pattern[2]) &&
	             (pattern[0]==pattern[3]) && (pattern[0]==pattern[4]))
	            reject = 1; /* 1s */
		  break;
	  case 6: if((pattern[0]==pattern[2]) && (pattern[0]==pattern[4]) &&
		     (pattern[1]==pattern[3]) && (pattern[1]==pattern[5]))
		    reject = 1; /* 2s */
		  if((pattern[0]==pattern[3]) && (pattern[1]==pattern[4]) &&
		     (pattern[2]==pattern[5]))
		    reject = 1; /* 3s */
		  break;
	  case 7: if((pattern[0]==pattern[1]) && (pattern[0]==pattern[2]) &&
	             (pattern[0]==pattern[3]) && (pattern[0]==pattern[4]) &&
		     (pattern[0]==pattern[5]) && (pattern[0]==pattern[6]))
	            reject = 1; /* 1s */
		  break;
	  case 8: if((pattern[0]==pattern[4]) && (pattern[1]==pattern[5]) &&
		     (pattern[2]==pattern[6]) && (pattern[3]==pattern[7]))
		    reject = 1; /* 4s */
		  break;
	  case 9: if((pattern[0]==pattern[3]) && (pattern[0]==pattern[6]) &&
		     (pattern[1]==pattern[4]) && (pattern[1]==pattern[7]) &&
		     (pattern[2]==pattern[5]) && (pattern[2]==pattern[8]))
		    reject = 1; /* 3s */
		  break;
	  case 10:if((pattern[0]==pattern[2]) && (pattern[0]==pattern[4]) &&
		     (pattern[0]==pattern[6]) && (pattern[0]==pattern[8]) &&
		     (pattern[1]==pattern[3]) && (pattern[1]==pattern[5]) &&
		     (pattern[1]==pattern[7]) && (pattern[1]==pattern[9]))
		    reject = 1; /* 2s */
	          if((pattern[0]==pattern[5]) && (pattern[1]==pattern[6]) &&
		     (pattern[2]==pattern[7]) && (pattern[3]==pattern[8]) &&
		     (pattern[4]==pattern[9]))
		    reject = 1; /* 5s */
		  break;
	  case 11:if((pattern[0]==pattern[1]) && (pattern[0]==pattern[2]) &&
	             (pattern[0]==pattern[3]) && (pattern[0]==pattern[4]) &&
		     (pattern[0]==pattern[5]) && (pattern[0]==pattern[6]) &&
		     (pattern[0]==pattern[7]) && (pattern[0]==pattern[8]) &&
		     (pattern[0]==pattern[9]) && (pattern[0]==pattern[10]))
	            reject = 1; /* 1s */
		  break;
	  case 12:if((pattern[0]==pattern[4]) && (pattern[0]==pattern[8]) &&
		     (pattern[1]==pattern[5]) && (pattern[1]==pattern[9]) &&
		     (pattern[2]==pattern[6]) && (pattern[2]==pattern[10])&&
		     (pattern[3]==pattern[7]) && (pattern[3]==pattern[11]))
		    reject = 1; /* 4s */
		  if((pattern[0]==pattern[6]) && (pattern[1]==pattern[7]) &&
		     (pattern[2]==pattern[8]) && (pattern[3]==pattern[9]) &&
		     (pattern[4]==pattern[10])&& (pattern[5]==pattern[11]))
		    reject = 1; /* 6s */
		  break;
	  case 13:if((pattern[0]==pattern[1]) && (pattern[0]==pattern[2]) &&
	             (pattern[0]==pattern[3]) && (pattern[0]==pattern[4]) &&
		     (pattern[0]==pattern[5]) && (pattern[0]==pattern[6]) &&
		     (pattern[0]==pattern[7]) && (pattern[0]==pattern[8]) &&
		     (pattern[0]==pattern[9]) && (pattern[0]==pattern[10])&&
		     (pattern[0]==pattern[11])&& (pattern[0]==pattern[12]))
	            reject = 1; /* 1s */
		  break;
	  case 14:if((pattern[0]==pattern[2]) && (pattern[0]==pattern[4]) &&
		     (pattern[0]==pattern[6]) && (pattern[0]==pattern[8]) &&
		     (pattern[0]==pattern[10])&& (pattern[0]==pattern[12])&&
		     (pattern[1]==pattern[3]) && (pattern[1]==pattern[5]) &&
		     (pattern[1]==pattern[7]) && (pattern[1]==pattern[9]) &&
		     (pattern[1]==pattern[11])&& (pattern[1]==pattern[13]))
		    reject = 1; /* 2s */
		  if((pattern[0]==pattern[7]) && (pattern[1]==pattern[8]) &&
		     (pattern[2]==pattern[9]) && (pattern[3]==pattern[10])&&
		     (pattern[4]==pattern[11])&& (pattern[5]==pattern[12])&&
		     (pattern[6]==pattern[13]))
		    reject = 1; /* 7s */
		  break;
	  case 15:if((pattern[0]==pattern[3]) && (pattern[0]==pattern[6]) &&
		     (pattern[0]==pattern[9]) && (pattern[0]==pattern[12])&&
		     (pattern[1]==pattern[4]) && (pattern[1]==pattern[7]) &&
		     (pattern[1]==pattern[10])&& (pattern[1]==pattern[13])&&
		     (pattern[2]==pattern[5]) && (pattern[2]==pattern[8]) &&
		     (pattern[2]==pattern[11])&& (pattern[2]==pattern[14]))
		    reject = 1; /* 3s */
	          if((pattern[0]==pattern[5]) && (pattern[0]==pattern[10])&&
		     (pattern[1]==pattern[6]) && (pattern[1]==pattern[11])&&
		     (pattern[2]==pattern[7]) && (pattern[2]==pattern[12])&&
		     (pattern[3]==pattern[8]) && (pattern[3]==pattern[13])&&
		     (pattern[4]==pattern[9]) && (pattern[4]==pattern[14]))
		    reject = 1; /* 5s */
	  case 16:if((pattern[0]==pattern[8]) && (pattern[1]==pattern[9]) &&
		     (pattern[2]==pattern[10])&& (pattern[3]==pattern[11])&&
		     (pattern[4]==pattern[12])&& (pattern[5]==pattern[13])&&
		     (pattern[6]==pattern[14])&& (pattern[7]==pattern[15]))
		    reject = 1; /* 8s */
		  break;
	  default:break;
	}
	if(reject) {k++; continue;}
	else
          k += sfound[0]; /* k is now the location of the candidate in db */
      }
      else
      {
        k++;
	continue;
      }
      
      /* make sure candidate is not already found 
       * or is not tail end of match already found
       */
      for(m=0; m<(*srfoundlen); m++)
      {
	if(!strcmp(pattern, srpattern[m]) && !((k - srfound[m]+1)%rep_len) && 
	   (k <= srfound[m]+rep_len*srcount[m]-1))
	{
	  reject = 1;
	}
      }
      
      if(reject) {k++; continue;}
      
      /* candidate survived filters, now count actual reps & store results */
      m = k + rep_len;
      srcount[*srfoundlen] = 1;
      
      while(!strncmp(pattern, db+m, rep_len) && (m<=dblen-rep_len)) 
      {
        m+=rep_len;
	srcount[*srfoundlen]++;
      }
      
      if(srcount[*srfoundlen] >= min_rep)
      {
        strcpy(srpattern[*srfoundlen], pattern);
        srfound[*srfoundlen] = k+1; /* first location in db is 1 for output*/
        (*srfoundlen)++;
      }
      k++;
    }
  }
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

void write_word(unsigned long word, int j)
{
  int i;
  char string[17], lookup[4];
  
  string[16] = 0x0;
  lookup[0] = 'A';
  lookup[1] = 'C';
  lookup[2] = 'T';
  lookup[3] = 'G';

  for(i=0; i<j; i++)
  {
#ifdef B_ENDIAN
#ifdef LONG32
    string[i] = lookup[(word << 2*i) >> 30];
#endif
#ifdef LONG64
    string[i] = lookup[(word << 2*i) >> 62];
#endif
#endif

#ifdef L_ENDIAN
#ifdef LONG32
    string[i] = lookup[(word << (30-2*i)) >> 30];
#endif
#ifdef LONG64
    string[i] = lookup[(word << (62-2*i)) >> 62];
#endif
#endif
  }
  string[j] = 0x0;
  printf("%s\n", string);
}
