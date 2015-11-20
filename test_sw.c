/* test_sw.c                                    http://cbl.sourceforge.net
 * 
 * test cb_sw for random strings that incrementally cross word 
 * boundries
 *
 * Copyright (C) 2004 University of Alaska Fairbanks
 * Institute of Arctic Biology
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
 * $Id: test_sw.c,v 1.2 2005/01/13 22:51:01 jlong777 Exp $
 */
 
#include "cbl.h"
#include "cb_macro.h"
#include <ctype.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/times.h>
#include <time.h>

#define REP 30000  /* number of benchmark reps */
#define SIZE   64  /* change size of test cases */

void slow_sw(long *dbl, long dbllen, long *dbs, long dbslen, long eg, 
             long og, long *sslookup, long *smax, long **algl, long **algm, 
             long **algs, long *alglen, long *algstl, long *algsts, long *errno);

int main()
{
  int i, j, k, m, n, failure=0;
  long in[SIZE], test[3*SIZE], *in2, og, eg=2, smax, smax2;
  long *algl, *algl2, *algm, *algm2, *algs, *algs2;
  long  alglen, alglen2, algstl, algstl2, algsts, algsts2, error=0, error2=0;
  char lastin[SIZE], lasttest[SIZE];
  char ch2[4], version[6];
  
  /* timing stuff */
  struct tms mytime;               /* struct filled by calls to times() */
  float users, userf, syss, sysf, utot, stot; /* start and finish times */
  clock_t stime, ftime;
  long totsec = 0;
  
  long  sslookup[27][32]= /* PAM 1 */
         /* A  B  C  D  E  F  G  H  I  J  K  L  M  N  O  P  Q  R  S  T  U  V  W  X  Y  Z */
       {-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,
  /*A*/ -1, 2, 0,-6, 0, 0, 0,-6, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,-6, 0, 0, 0, 0, 0, 0,-1,-1,-1,-1,-1,
  /*B*/ -1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,-1,-1,-1,-1,-1,
  /*C*/ -1,-6, 0, 2, 0, 0, 0,-6, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,-6, 0, 0, 0, 0, 0, 0,-1,-1,-1,-1,-1,
  /*D*/ -1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,-1,-1,-1,-1,-1,
  /*E*/ -1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,-1,-1,-1,-1,-1,
  /*F*/ -1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,-1,-1,-1,-1,-1,
  /*G*/ -1,-6, 0,-6, 0, 0, 0, 2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,-6, 0, 0, 0, 0, 0, 0,-1,-1,-1,-1,-1,
  /*H*/ -1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,-1,-1,-1,-1,-1,
  /*I*/ -1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,-1,-1,-1,-1,-1,
  /*J*/ -1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,-1,-1,-1,-1,-1,
  /*K*/ -1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,-1,-1,-1,-1,-1,
  /*L*/ -1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,-1,-1,-1,-1,-1,
  /*M*/ -1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,-1,-1,-1,-1,-1,
  /*N*/ -1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,-1,-1,-1,-1,-1,
  /*O*/ -1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,-1,-1,-1,-1,-1,
  /*P*/ -1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,-1,-1,-1,-1,-1,
  /*Q*/ -1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,-1,-1,-1,-1,-1,
  /*R*/ -1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,-1,-1,-1,-1,-1,
  /*S*/ -1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,-1,-1,-1,-1,-1,
  /*T*/ -1,-6, 0,-6, 0, 0, 0,-6, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 2, 0, 0, 0, 0, 0, 0,-1,-1,-1,-1,-1,
  /*U*/ -1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,-1,-1,-1,-1,-1,
  /*V*/ -1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,-1,-1,-1,-1,-1,
  /*W*/ -1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,-1,-1,-1,-1,-1,
  /*X*/ -1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,-1,-1,-1,-1,-1,
  /*Y*/ -1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,-1,-1,-1,-1,-1,
  /*Z*/ -1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,-1,-1,-1,-1,-1};
  
  
  ch2[0] = 'A';
  ch2[1] = 'C';
  ch2[2] = 'G';
  ch2[3] = 'T';

  printf("running test of cb_swa...\n");
  
  for(og=0; og<30; og++)
  {
    /* test random strings up to size SIZE-1 */
    for(i=1; i<SIZE-1; i++) /* dblen must be > 0 */
    {
      for(k=0; k<i; k++)
        ((char *)in)[k] = ch2[rand()%4];
      ((char *)in)[k] = 0x0;
    
      for(j=1; j<3*SIZE-1; j++) /* j will be length of the test string > 0 & < i */
      {
        for(k=0; k<j; k++)
          ((char *)test)[k] = ch2[rand()%4];
        ((char *)test)[k] = 0x0;

        cb_swa_fw(in, i, test, j, eg, og, &sslookup[0][0], &smax, &algl, &algm, 
	          &algs, &alglen, &algstl, &algsts, &error);

        slow_sw(in, i, test, j, eg, og, &sslookup[0][0], &smax2, &algl2, &algm2, 
	        &algs2, &alglen2, &algstl2, &algsts2, &error2);

        if(smax!=smax2 || alglen!=alglen2 || algstl!=algstl2 || (algsts)!=algsts2 ||
	   error!=error2)
        {
          if(DEBUG_INFO)
          {
	    printf("og = %d\n", og);
            printf("smax = %d  smax2 = %d\n", smax, smax2);
	    printf("alglen = %d  alglen2 = %d\n",alglen, alglen2);
	    printf("algstl = %d  algstl2 = %d\n",algstl, algstl2);
	    printf("algsts = %d  algsts2 = %d\n",algsts, algsts2);
	    printf("error = %d  error2 = %d\n",error, error2);
            printf("in = %s\ntest=%s\n\n", (char *)in, (char *)test);
	    
	    printf("%ld\n", algstl);
            for(m=0; m<alglen; m++)
              printf("%c", *(((char *) algl)+m));
            printf("\n");
            for(m=0; m<alglen; m++)
              printf("%c", *(((char *) algm)+m));
            printf("\n");
            for(m=0; m<alglen; m++)
              printf("%c", *(((char *) algs)+m));
            printf("\n");
            printf("%ld\n\n\n", algsts);
	    
	    printf("%ld\n", algstl2);
            for(m=0; m<alglen2; m++)
              printf("%c", *(((char *) algl2)+m));
            printf("\n");
            for(m=0; m<alglen2; m++)
              printf("%c", *(((char *) algm2)+m));
            printf("\n");
            for(m=0; m<alglen2; m++)
              printf("%c", *(((char *) algs2)+m));
            printf("\n");
            printf("%ld\n\n\n", algsts2);
          }
          failure = 1;
        }
      
      }
    }
  }

  if(failure) printf("**failure**\n");
  else        printf("**success**\n");
  
  /* =================================================================== */
  /* benchmark */
  if(BMARK)
  {
    printf("times() in clock ticks to do %d cb_swa:\n", REP);
  
    utot = stot = 0.0;
    
    printf("length of string to align = %d\n", strlen((char*)test));
  
    for(i=0; i<18; i++)
    {
      n = 0x00000100 << i;
      
      in2 = (long *) calloc(n/(sizeof(long))+9, sizeof(long));
      if(in2==NULL)
      {
        printf("test_sw: calloc error, exiting...\n");
        exit(0);
      }
      
      /* create random string */
      for(j=0; j<n; j++)
        ((char *)in2)[j] = ch2[rand()%4];
      ((char *)in2)[j] = 0x0;

      printf("n = %d\n", n);
      
      j = strlen((char*)test);

      stime = clock()/CLOCKS_PER_SEC;
      times(&mytime);
      users = (float) mytime.tms_utime;
      syss  = (float) mytime.tms_stime;

      /* search REP times */
      for(k=0; k<REP; k++)
      {
        cb_swa_fw(in2, i, test, j, eg, 10, &sslookup[0][0], &smax, &algl, &algm, 
	          &algs, &alglen, &algstl, &algsts, &error);
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
      
      free(in2);
    }
    printf("total clock ticks:user = %9.1f sys = %9.1f\n", utot, stot);
    printf("total approx secs = %d for swa\n\n", totsec);
  
  }
  
  cb_version(version);
  printf("CBL version: %s\n\n", version);
  
  return 0;
}

void slow_init(long *dbl, long dbllen, long *dbs, long dbslen, 
               long *sslookup, long *swtab, long *gaph, long *gapv)
{
  register int i, i1, j, dbsi;
  
  for(j=0; j<dbslen+1; j++)
    swtab[j] = 0;               /* row 0 is 0 */
  
  for(i=1; i<dbllen+1; i++)
  {
    dbsi = (dbslen+1)*i;
    i1 = 0x1F & ((char *)dbl)[i-1];
    swtab[(dbslen+1)*i] = 0;    /* col 0 is 0 */
    for(j=1; j<dbslen+1; j++)
      swtab[dbsi+j] = sslookup[(32*(0x1F & ((char *)dbs)[j-1]))+ i1];
  }
}

void slow_score(long *swtab, long *gaph, long *gapv, long dbllen, 
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

void slow_fill(long *swtab, long swtab_i, long swtab_j, long *dbl, long *dbs,
               long index, long **algl, long **algm, long **algs, long dbslen)
{
  char *algl_c = (char *)(*algl);
  char *algm_c = (char *)(*algm);
  char *algs_c = (char *)(*algs);
  
  if(swtab_i==0 || swtab_j==0 || index<0) return;
  
#ifdef LONG32
  switch((unsigned long)(0xF0000000&swtab[(dbslen+1)*swtab_i + swtab_j]) >> 28)
#endif
#ifdef LONG64
  switch((unsigned long)(0xF000000000000000&swtab[(dbslen+1)*swtab_i + swtab_j]) >> 60)
#endif
  {
    case 0:  /* end  */
             swtab_i = 0;
	     swtab_j = 0;
	     break;
    case 1:  /* horizontal (fortran), vertical (c) */
             algl_c[index] = ((char *)dbl)[swtab_i-1];
             algm_c[index] = ' ';
             algs_c[index] = '-';
	     swtab_i--;
             break;
    case 2:   /* vertical (fortran), horizontal (c) */
             algl_c[index] = '-';
             algm_c[index] = ' ';
             algs_c[index] = ((char *)dbs)[swtab_j-1];
	     swtab_j--;
	     break;
    case 3:  /* h or v */
             algl_c[index] = ((char *)dbl)[swtab_i-1];
             algm_c[index] = ' ';
             algs_c[index] = '-';
	     swtab_i--;
             break;
    default: /* corner */
             algl_c[index] = ((char *)dbl)[swtab_i-1];
             algs_c[index] = ((char *)dbs)[swtab_j-1];
	     if(algl_c[index] == algs_c[index])
               algm_c[index] = ':';
	     else
	       algm_c[index] = ' ';
	     swtab_i--;
	     swtab_j--;
             break;
  }
  
  index--;

  slow_fill(swtab, swtab_i, swtab_j, dbl, dbs, index, algl, algm, algs, dbslen);
}

void slow_align(long *swtab, long *dbl, long dbllen, long *dbs,long dbslen,
                     long smax, long **algl, long **algm, long **algs,long *alglen,
                     long *algstl, long *algsts, long *errno)
{
  long blah, i, index=0, j, savi=0, savj=0;
  
  /* find smax */
  for(i=dbllen;i>0;i--)
  {
    for(j=dbslen;j>0;j--)
    {
#ifdef LONG32
      if((0x0FFFFFFF & swtab[(dbslen+1)*i+j]) == smax)
#endif
#ifdef LONG64
      if((0x0FFFFFFFFFFFFFFF & swtab[(dbslen+1)*i+j]) == smax)
#endif
      {
        savi = i;
        savj = j;
	goto jump;
      }
    }
  }
  
jump:

  /* compute alglen */
  
  i = savi; j = savj; *alglen = 0;
  
#ifdef LONG32
  while(blah=((unsigned long)(0xF0000000 & swtab[(dbslen+1)*i + j]) >> 28))
#endif
#ifdef LONG64
  while(blah=((unsigned long)(0xF000000000000000 & swtab[(dbslen+1)*i + j]) >> 60))
#endif
    switch(blah)
    {
      case 0:  printf("I should never get here!\n");
	       break;
      case 1:  i--;
               (*alglen)++;
	       break;
      case 2:  j--;
               (*alglen)++;
               break;
      case 3:  i--;
               (*alglen)++;
               break;
      default: i--;
	       j--;
               (*alglen)++;
               break;
    }
  
  *algstl = i+1;
  *algsts = j+1;
  
  *algl = (long *) malloc(*alglen+8);
  *algm = (long *) malloc(*alglen+8);
  *algs = (long *) malloc(*alglen+8);
  
  if(*algl==0 || *algm==0 || *algs==0)
  {
    printf("cb_swa_fw_align: malloc error in align, exiting...\n");
    exit(0);
  }
  
  index = *alglen-1;

  /* recursively trace back through scoring matrix for alignment */
  slow_fill(swtab, savi, savj, dbl, dbs, index, algl, algm, algs, dbslen);
}


void slow_sw(long *dbl, long dbllen, long *dbs, long dbslen, long eg, 
             long og, long *sslookup, long *smax, long **algl, long **algm, 
             long **algs, long *alglen, long *algstl, long *algsts, long *errno)
{
  long *swtab, *gaph, *gapv, eog;
#ifdef LONG32
  swtab = (long *) calloc((dbllen+1)*(dbslen+1), 4);
  if(swtab == 0x0)
  {
    printf("test_sw: swtab allocation error, returning...\n");
    *smax = 0;
    return;
  }
  gaph  = (long *) calloc(dbslen+1, 4);
  if(gaph == 0x0)
  {
    printf("test_sw: gapv allocation error, returning...\n");
    *smax = 0;
    free(swtab);
    return;
  }
  gapv  = (long *) calloc(dbslen+1, 4);
  if(gapv == 0x0)
  {
    printf("test_sw: gapv allocation error, returning...\n");
    *smax = 0;
    free(swtab);
    free(gaph);
    return;
  }
#endif
#ifdef LONG64
  swtab = (long *) calloc((dbllen+1)*(dbslen+1), 8);
  if(swtab == 0x0)
  {
    printf("test_sw: swtab allocation error, returning...\n");
    *smax = 0;
    return;
  }
  gaph  = (long *) calloc(dbslen+1, 8);
  if(gaph == 0x0)
  {
    printf("test_sw: gapv allocation error, returning...\n");
    *smax = 0;
    free(swtab);
    return;
  }
  gapv  = (long *) calloc(dbslen+1, 8);
  if(gapv == 0x0)
  {
    printf("test_sw: gaph allocation error, returning...\n");
    *smax = 0;
    free(swtab);
    free(gaph);
    return;
  }
#endif
  eog = eg + og;
  
  slow_init(dbl, dbllen, dbs, dbslen, sslookup, swtab, gaph, gapv);
  slow_score(swtab, gaph, gapv, dbllen, dbslen, eg, eog, smax);
  slow_align(swtab, dbl, dbllen, dbs, dbslen, *smax, algl, algm, algs, 
             alglen, algstl, algsts, errno);
  free(swtab);
  free(gaph);
  free(gapv);
}
