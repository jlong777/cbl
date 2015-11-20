/* test_suite.c                                  http://cbl.sourceforge.net
 * 
 * test cb_compress, cb_uncompress, cb_revcompl, cb_amino_translate_ascii, 
 * and cb_count_ascii for random strings that incrementally cross word 
 * boundries for all modes.
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
 * $Id: test_suite.c,v 1.6 2004/02/04 18:05:37 jlong777 Exp $
 */
 
#include "cbl.h"
#include "cb_macro.h"
#include <ctype.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/times.h>
#include <time.h>

#define REP   40  /* number of times to repeat benchmark case */
#define TSZ  1000  /* set the number of words for testing */

char codon_trans(char *, int);

int main()
{
  char ch2[4], ch4[15], ch5[26], rc2[256]={0x0};
  char test[5], version[6];

  int i, j, k, mode, n, failure=0, size, counter=0;
  char *ptc;
  long *in, *out, *pack, *four, *rc, *a0, *a1, *a2;
  long res[5], a_cnt, c_cnt, g_cnt, n_cnt, t_cnt; /* cb_count_ascii stuff */
  
  /* timing stuff */
  struct tms mytime;               /* struct filled by calls to times() */
  float users, userf, syss, sysf, utot, stot; /* start and finish times */
  clock_t stime, ftime;
  long totsec = 0;
  
  size = sizeof(long);
  
  ch2[0] = 'A';
  ch2[1] = 'C';
  ch2[2] = 'G';
  ch2[3] = 'T';
  
  ch4[0] = 'A';  ch4[4] = 'G';  ch4[8]  = 'N';  ch4[12] = 'V';
  ch4[1] = 'B';  ch4[5] = 'H';  ch4[9]  = 'R';  ch4[13] = 'W';
  ch4[2] = 'C';  ch4[6] = 'K';  ch4[10] = 'S';  ch4[14] = 'Y';
  ch4[3] = 'D';  ch4[7] = 'M';  ch4[11] = 'T'; 
  
  ch5[0] = 'A';  ch5[7]  = 'H';  ch5[14] = 'O';  ch5[21] = 'V';
  ch5[1] = 'B';  ch5[8]  = 'I';  ch5[15] = 'P';  ch5[22] = 'W';
  ch5[2] = 'C';  ch5[9]  = 'J';  ch5[16] = 'Q';  ch5[23] = 'X';
  ch5[3] = 'D';  ch5[10] = 'K';  ch5[17] = 'R';  ch5[24] = 'Y';
  ch5[4] = 'E';  ch5[11] = 'L';  ch5[18] = 'S';  ch5[25] = 'Z';
  ch5[5] = 'F';  ch5[12] = 'M';  ch5[19] = 'T';
  ch5[6] = 'G';  ch5[13] = 'N';  ch5[20] = 'U';
  
  rc2['A'] = 'T'; rc2['C'] = 'G'; rc2['T'] = 'A'; rc2['G'] = 'C';

  printf("\n\n\n");
  printf("=============\n");
  printf("RUNNING TESTS\n");
  printf("=============\n\n");

  printf("running tests of cb_compress, cb_uncompress,\n");
  printf("                 cb_count_ascii, cb_revcompl,\n");
  printf("                 and cb_amino_translate_ascii...\n");  

  in = (long *) calloc((TSZ+2)*size,1);
  if(in==NULL)
  {
    printf("test_compress: calloc error, exiting...\n");
    exit(0);
  }
  ((char *)in)[TSZ*size+1] = 0x0;
  
  out = (long *) calloc((TSZ+2)*size,1);
  if(out==NULL)
  {
    printf("test_compress: calloc error, exiting...\n");
    exit(0);
  }
  
  pack = (long *) calloc((TSZ+1)*size,1);
  if(pack==NULL)
  {
    printf("test_compress: calloc error, exiting...\n");
    exit(0);
  }
  
  rc = (long *) calloc((TSZ+1)*size,1);
  if(rc==NULL)
  {
    printf("test_compress: calloc error, exiting...\n");
    exit(0);
  }
  
  four = (long *) calloc((TSZ+1)*size,1);
  if(four==NULL)
  {
    printf("test_compress: calloc error, exiting...\n");
    exit(0);
  }
  
  a0 = (long *) calloc((TSZ+1)*size,1);
  if(a0==NULL)
  {
    printf("test_compress: calloc error, exiting...\n");
    exit(0);
  }
  
  a1 = (long *) calloc((TSZ+1)*size,1);
  if(a1==NULL)
  {
    printf("test_compress: calloc error, exiting...\n");
    exit(0);
  }
  
  a2 = (long *) calloc((TSZ+1)*size,1);
  if(a2==NULL)
  {
    printf("test_compress: calloc error, exiting...\n");
    exit(0);
  }
  
  mode = 2;
  while(mode)
  {
    /* test random strings up to TSZ compressed words*/
    for(i=0; i<TSZ*size; i++)
    {
      for(j=0; j<i; j++)
      {
        if(mode==2)
          ((char *)in)[j] = ch2[rand()%4];
	else if(mode==4)
          ((char *)in)[j] = ch4[rand()%15];
	else
	  ((char *)in)[j] = ch5[rand()%26];
      }
      ((char *)in)[j] = 0x0; 
      cb_compress(in, pack, strlen((char *)in), mode);
      cb_uncompress(pack, out, strlen((char *)in), mode);
      
      if(strncmp((char *)in, (char *)out, strlen((char *)in)))
      {
        if(DEBUG_INFO)                                                                  {
          printf("test_compress: comparison test failure\n");
	  printf("mode = %d strlen(in) = %d strlen(out) = %d\n", mode, strlen((char *)in), strlen((char *)out));
	  printf("original:     %s\n", (char *)in);
	  printf("uncompr pack: %s\n\n",(char *)out);
        }
        failure = 1;
      }
      
      if(mode==2)
      {
        /* test cb_revcompl for mode 2 */
        cb_revcompl(pack, rc, strlen((char *)in), mode);
	cb_uncompress(rc, out, strlen((char *)in), mode);
	ptc = (char *)out + strlen((char *)in);
	/* explicitly lookup rev compl for each input */
	for(k=0; ((char *)in)[k]!=0x0; k++)
	{
	  if(((char *)in)[k] != rc2[*(ptc-1-k)])
	  {
            if(DEBUG_INFO)                                                                  {
	      printf("test_revcompl: mode 2 test failure\n");
	      printf("in[%d] = %c revcomp[strlen(in)-k-1] = %c\n", k, ((char *)in)[k], rc2[*(ptc-1-k)]);
            }
            failure = 1;
	    break;
	  }
	}
	
        /* test cb_amino_translate_ascii */
	cb_amino_translate_ascii(pack, strlen((char *)in), a0, a1, a2);
	
        /*
        printf("original:     %s\n", in);
        printf("amino a0:  %s\n\n",(char *)a0);
        printf("amino a1:  %s\n\n",(char *)a1);
        printf("amino a2:  %s\n\n",(char *)a2);
        */
	
	if(strlen((char *)in) > 5) /* wait till there is enough input to do all 3 */
	{
	  for(k=0; k<strlen((char *)in)-5; k+=3)
          {
            test[k%3]   = ((char *)in)[k];
            test[k%3+1] = ((char *)in)[k+1];
            test[k%3+2] = ((char *)in)[k+2];
            test[k%3+3] = ((char *)in)[k+3];
            test[k%3+4] = ((char *)in)[k+4];
    
            /* explicitly compute traslation at a higher level for comparison */
            if(codon_trans(test, k%3) != ((char *)a0)[k/3])
            {
              if(DEBUG_INFO)
	      {
                printf("failure:\n");
                printf("codon_trans = %c a0[%d] = %c\n",codon_trans(test, k%3),k/3,((char *)a0)[k/3]);
              printf("original:     %s\n", in);
        printf("amino a0:  %s\n\n",(char *)a0);
        printf("amino a1:  %s\n\n",(char *)a1);
        printf("amino a2:  %s\n\n",(char *)a2);
	      }
              failure = 1;
            }
            if(codon_trans(test, k%3+1) != ((char *)a1)[k/3])
            {
              if(DEBUG_INFO)
	      {
                printf("failure:\n");
                printf("codon_trans = %c a1[%d] = %c\n",codon_trans(test, k%3+1),k/3,((char *)a1)[k/3]);
              }
              failure = 1;
            }
            if(codon_trans(test, k%3+2) != ((char *)a2)[k/3])
            {
              if(DEBUG_INFO)
	      {
                printf("failure:\n");
                printf("codon_trans = %c a2[%d] = %c\n",codon_trans(test, k%3+2),k/3,((char *)a2)[k/3]);
              }
              failure = 1;
            }
          }
	}
	
	/* test cb_count_ascii */
        a_cnt = c_cnt = g_cnt = n_cnt = t_cnt = 0;
	for(k=0; k<5; k++)
          res[k] = 0;
	
	while(j%size) ((char *)in)[j++] = 0x0; /* pad final word with nulls */
	
	cb_countn_ascii(in, strlen((char *)in), res);
	
	/* count the hard way and compare against res */
	for(k=0; k<strlen((char *)in); k++)
        {
	  switch(((char *)in)[k])
	  {
	    case 'A': a_cnt++; break;
	    case 'a': a_cnt++; break;
	    case 'C': c_cnt++; break;
	    case 'c': c_cnt++; break;
	    case 'G': g_cnt++; break;
	    case 'g': g_cnt++; break;
	    case 'N': n_cnt++; break;
	    case 'n': n_cnt++; break;
	    case 'T': t_cnt++; break;
	    case 't': t_cnt++; break;
	    default:  break;
	  }
	}
	if(a_cnt!=res[0] || c_cnt!=res[1] || g_cnt!=res[3] || t_cnt!=res[2] || n_cnt!=res[4])
	{
          if(DEBUG_INFO)                                                                  {
	    printf("failure:\n");
	    printf("input string: %s\n", (char *)in);
            printf("countn_ascii: a_cnt = %d c_cnt = %d g_cnt = %d n_cnt = %d t_cnt = %d\n",a_cnt,c_cnt,g_cnt,n_cnt,t_cnt);
	    printf("              res_A = %d res_C = %d res_G = %d res_N = %d res_T = %d\n",res[0],res[1],res[3],res[4],res[2]);
          }
          failure = 1;
	}
      }
    }
    switch(mode)
    {
      case 2: mode=4; break;
      case 4: mode=5; break;
      case 5: mode=0; break;
      default:mode=0;
    }
  }
  
  free(a0);
  free(a1);
  free(a2);
  
  /* test mode 6 */
  for(i=0; i<TSZ*size; i++)
  {
    for(j=0; j<i; j++)
      ((char *)in)[j] = ch2[rand()%4];
   ((char *)in)[j] = 0x0;
    
    cb_compress(in, pack, strlen((char *)in), 4);
    cb_compress(pack, four, strlen((char *)in), 6);
    cb_uncompress(four, out, strlen((char *)in), 2);

    if(strncmp((char *)in, (char *)out, strlen((char *)in)))
    {
      if(DEBUG_INFO)                                                                  {
        printf("test_compress: mode 6 test failure\n");
        printf("mode = %d strlen(in) = %d\n", mode, strlen((char *)in));
        printf("original:     %s\n", (char *)in);
        printf("uncompr pack: %s\n\n",(char *)out);
      }
      failure = 1;
    }
    
    /* test cb_revcompl for mode 4 */
    cb_revcompl(pack, rc, strlen((char *)in), 4);
    cb_uncompress(rc, out, strlen((char *)in), 4);
    ptc = (char *)out + strlen((char *)in);
    for(k=0; ((char *)in)[k]!=0x0; k++)
    {
      if(((char *)in)[k] != rc2[*(ptc-1-k)])
      {
        if(DEBUG_INFO)                                                                  {
          printf("original:     %s\n", (char *)in);
  	  printf("rev comp:     %s\n\n",(char *)out);
          printf("test_revcompl: mode 4 test failure\n");
          printf("in[%d] = %c revcomp[strlen(in)-k-1] = %c\n", k, ((char *)in)[k], rc2[*(ptc-1-k)]);
        }
        failure = 1;
	break;
      }
    }
  }
  
  free(in);
  free(out);
  free(pack);
  free(four);
  free(rc);
    
  if(failure)
    printf("**failure**\n");
  else
    printf("**success**\n");

  /* =================================================================== */
  /* benchmark */
  if(BMARK)
  {
    printf("times() in clock ticks to do %d compress/uncompresses for modes 2, 4, & 5:\n", REP);
  
    utot = stot = 0.0;
  
    for(i=0; i<18; i++)
    {
      n = 0x00000100 << i;
      in = (long *) calloc(n/(sizeof(long))+1, sizeof(long));
      if(in==NULL)
      {
        printf("test_compress: calloc error, exiting...\n");
        exit(0);
      }
      ((char *)in)[n] = 0x0;
  
      out = (long *) calloc(n/(sizeof(long))+9, sizeof(long));
      if(out==NULL)
      {
        printf("test_compress: calloc error, exiting...\n");
        exit(0);
      }
  
      pack = (long *) calloc(n/(sizeof(long))+1, sizeof(long));
      if(pack==NULL)
      {
        printf("test_compress: calloc error, exiting...\n");
        exit(0);
      }
      
      /* create random string */
      for(j=0; j<n; j++)
        ((char *)in)[j] = ch2[rand()%4];
      ((char *)in)[j] = 0x0;

      printf("n = %d\n", n);
      
      mode = 2;
      while(mode)
      {
        /* compress/uncompress */
	stime = clock()/CLOCKS_PER_SEC;
        times(&mytime);
        users = (float) mytime.tms_utime;
        syss  = (float) mytime.tms_stime;
      
        /* compress/uncompress REP times */
	for(k=0; k<REP; k++)
	{
          cb_compress(in, pack, n, mode);
          cb_uncompress(pack, out, n, mode);
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
	
	switch(mode)
        {
          case 2: mode=4; break;
          case 4: mode=5; break;
          case 5: mode=0; break;
          default:mode=0;
        }
      }
      
      free(in);
      free(out);
      free(pack);
    }
    printf("total clock ticks:user = %9.1f sys = %9.1f\n", utot, stot);
    printf("total approx secs = %d for compression/uncompression\n\n", totsec);
  
    /* =================================================================== */
    /* reverse complement */
    
    printf("times() in clock ticks to do %d reverse complements for modes 2 & 4:\n", 10*REP);
  
    utot = stot = 0.0;
    totsec = 0;
  
    for(i=0; i<18; i++)
    {
      n = (0x00000100 << i)+1; /* string inputs are doubled in size each time */
      in = (long *) calloc(n/(sizeof(long))+1, sizeof(long));
      if(in==NULL)
      {
        printf("test_revcompl: calloc error, exiting...\n");
        exit(0);
      }
      ((char *)in)[n] = 0x0;
  
      pack = (long *) calloc(n/(sizeof(long))+1, sizeof(long));
      if(pack==NULL)
      {
        printf("test_revcompl: calloc error, exiting...\n");
        exit(0);
      }
  
      rc = (long *) calloc(n/(sizeof(long))+1, sizeof(long));
      if(rc==NULL)
      {
        printf("test_revcompl: calloc error, exiting...\n");
        exit(0);
      }
      
      /* create random string */
      for(j=0; j<n; j++)
        ((char *)in)[j] = ch2[rand()%4];
      ((char *)in)[j] = 0x0;

      printf("n = %d\n", n);
      
      mode = 2;
      while(mode)
      {
        cb_compress(in, pack, n, mode);
	stime = clock()/CLOCKS_PER_SEC;
        times(&mytime);
        users = (float) mytime.tms_utime;
        syss  = (float) mytime.tms_stime;
      
        /* revcompl 10*REP times */
	for(k=0; k<10*REP; k++)
	  cb_revcompl(pack, rc, n, mode);
	
	times(&mytime);
        userf = (float) mytime.tms_utime;
        sysf  = (float) mytime.tms_stime;
	ftime = clock()/CLOCKS_PER_SEC;
	
	printf("mode = %d user = %7.1f  sys = %7.1f  secs (approx) = %d\n", 
	        mode, userf-users, sysf-syss, ftime-stime);
	
	
	utot += userf-users;
	stot += sysf-syss;
	totsec += ftime-stime;
	
	switch(mode)
        {
          case 2: mode=4; break;
          case 4: mode=0; break;
          default:mode=0;
        }
      }
      
      free(in);
      free(pack);
      free(rc);
    }
    printf("total clock ticks:user = %9.1f sys = %9.1f\n", utot, stot);
    printf("total approx secs = %d for reverse complement\n\n", totsec);
    
    /* =================================================================== */
    /* amino translate ascii */
    
    printf("times() in clock ticks to do %d amino_translate_asciis:\n", 5*REP);
  
    utot = stot = 0.0;
    totsec = 0;
  
    for(i=0; i<18; i++)
    {
      n = 0x00000100 << i;
      in = (long *) calloc(n/(sizeof(long))+1, sizeof(long));
      if(in==NULL)
      {
        printf("test_amino_translate_ascii: calloc error, exiting...\n");
        exit(0);
      }
      ((char *)in)[n] = 0x0;
  
      pack = (long *) calloc(n/(sizeof(long))+1, sizeof(long));
      if(pack==NULL)
      {
        printf("test_amino_translate_ascii: calloc error, exiting...\n");
        exit(0);
      }
  
      a0 = (long *) calloc(n/(sizeof(long))+1, sizeof(long));
      if(a0==NULL)
      {
        printf("test_amino_translate_ascii: calloc error, exiting...\n");
        exit(0);
      }
      
      a1 = (long *) calloc(n/(sizeof(long))+1, sizeof(long));
      if(a1==NULL)
      {
        printf("test_amino_translate_ascii: calloc error, exiting...\n");
        exit(0);
      }
      
      a2 = (long *) calloc(n/(sizeof(long))+1, sizeof(long));
      if(a2==NULL)
      {
        printf("test_amino_translate_ascii: calloc error, exiting...\n");
        exit(0);
      }
      
      /* create random string */
      for(j=0; j<n; j++)
        ((char *)in)[j] = ch2[rand()%4];
      ((char *)in)[j] = 0x0;
      
      printf("n = %d\n", n);
      
      mode = 2;
      
      cb_compress(in, pack, n, mode);
      stime = clock()/CLOCKS_PER_SEC;
      times(&mytime);
      users = (float) mytime.tms_utime;
      syss  = (float) mytime.tms_stime;
      
      /* amino_translate_ascii 5*REP times */
      for(k=0; k<5*REP; k++)
        cb_amino_translate_ascii(pack, n, a0, a1, a2);

      times(&mytime);
      userf = (float) mytime.tms_utime;
      sysf  = (float) mytime.tms_stime;
      ftime = clock()/CLOCKS_PER_SEC;
	
      printf("mode = %d user = %7.1f  sys = %7.1f  secs (approx) = %d\n", 
	      mode, userf-users, sysf-syss, ftime-stime);
	
	
      utot += userf-users;
      stot += sysf-syss;
      totsec += ftime-stime;
	
      free(in);
      free(pack);
      free(a0);
      free(a1);
      free(a2);
    }
    printf("total clock ticks:user = %9.1f sys = %9.1f\n", utot, stot);
    printf("total approx secs = %d for amino_translate_ascii\n\n", totsec);
    
    /* =================================================================== */
    /* countn ascii */
    
    printf("times() in clock ticks to do %d countn_asciis:\n", 3*REP);
  
    utot = stot = 0.0;
    totsec = 0;
  
    for(i=0; i<18; i++)
    {
      n = 0x00000100 << i;
      in = (long *) calloc(n/(sizeof(long))+1, sizeof(long));
      if(in==NULL)
      {
        printf("test_amino_translate_ascii: calloc error, exiting...\n");
        exit(0);
      }
      ((char *)in)[n] = 0x0;
      
      /* create random string */
      for(j=0; j<n; j++)
        ((char *)in)[j] = ch2[rand()%4];
      ((char *)in)[j] = 0x0;
      
      printf("n = %d\n", n);
      
      stime = clock()/CLOCKS_PER_SEC;
      times(&mytime);
      users = (float) mytime.tms_utime;
      syss  = (float) mytime.tms_stime;
      
      /* countn_ascii 3*REP times */
      for(k=0; k<3*REP; k++)
        cb_countn_ascii(in, n, res);
	
      times(&mytime);
      userf = (float) mytime.tms_utime;
      sysf  = (float) mytime.tms_stime;
      ftime = clock()/CLOCKS_PER_SEC;
	
      printf("user = %7.1f  sys = %7.1f  secs (approx) = %d\n", 
	      userf-users, sysf-syss, ftime-stime);
	
	
      utot += userf-users;
      stot += sysf-syss;
      totsec += ftime-stime;
	
      free(in);
    }
    printf("total clock ticks:user = %9.1f sys = %9.1f\n", utot, stot);
    printf("total approx secs = %d for countn_ascii\n\n", totsec);
    
  }
  cb_version(version);
  printf("CBL version: %s\n\n", version);
  
  return 0;
}

char codon_trans(char *three, int index) /* input a codon string, return amino acid */
{
  switch(toupper(three[index]))
  {
    case 'A': switch(toupper(three[index+1]))
              {
                case 'A': switch(toupper(three[index+2]))
                          {
			    case 'A': return 'K';
			    case 'C': return 'N';
			    case 'G': return 'K';
			    case 'T': return 'N';
			    case 'U': return 'N';
			    default: return 'Z';
			  }
		case 'C': switch(toupper(three[index+2]))
                          {
			    case 'A': return 'T';
			    case 'C': return 'T';
			    case 'G': return 'T';
			    case 'T': return 'T';
			    case 'U': return 'T';
			    default: return 'Z';
			  }
		case 'G': switch(toupper(three[index+2]))
                          {
			    case 'A': return 'R';
			    case 'C': return 'S';
			    case 'G': return 'R';
			    case 'T': return 'S';
			    case 'U': return 'S';
			    default: return 'Z';
			  }
		case 'T': switch(toupper(three[index+2]))
                          {
			    case 'A': return 'I';
			    case 'C': return 'I';
			    case 'G': return 'M';
			    case 'T': return 'I';
			    case 'U': return 'I';
			    default: return 'Z';
			  }
		case 'U': switch(toupper(three[index+2]))
                          {
			    case 'A': return 'I';
			    case 'C': return 'I';
			    case 'G': return 'M';
			    case 'T': return 'I';
			    case 'U': return 'I';
			    default: return 'Z';
			  }
		default: return 'Z';
	      }
	      
    case 'C': switch(toupper(three[index+1]))
              {
                case 'A': switch(toupper(three[index+2]))
                          {
			    case 'A': return 'Q';
			    case 'C': return 'H';
			    case 'G': return 'Q';
			    case 'T': return 'H';
			    case 'U': return 'H';
			    default: return 'Z';
			  }
		case 'C': switch(toupper(three[index+2]))
                          {
			    case 'A': return 'P';
			    case 'C': return 'P';
			    case 'G': return 'P';
			    case 'T': return 'P';
			    case 'U': return 'P';
			    default: return 'Z';
			  }
		case 'G': switch(toupper(three[index+2]))
                          {
			    case 'A': return 'R';
			    case 'C': return 'R';
			    case 'G': return 'R';
			    case 'T': return 'R';
			    case 'U': return 'R';
			    default: return 'Z';
			  }
		case 'T': switch(toupper(three[index+2]))
                          {
			    case 'A': return 'L';
			    case 'C': return 'L';
			    case 'G': return 'L';
			    case 'T': return 'L';
			    case 'U': return 'L';
			    default: return 'Z';
			  }
		case 'U': switch(toupper(three[index+2]))
                          {
			    case 'A': return 'L';
			    case 'C': return 'L';
			    case 'G': return 'L';
			    case 'T': return 'L';
			    case 'U': return 'L';
			    default: return 'Z';
			  }
		default: return 'Z';
	      }
	      
    case 'G': switch(toupper(three[index+1]))
              {
                case 'A': switch(toupper(three[index+2]))
                          {
			    case 'A': return 'E';
			    case 'C': return 'D';
			    case 'G': return 'E';
			    case 'T': return 'D';
			    case 'U': return 'D';
			    default: return 'Z';
			  }
		case 'C': switch(toupper(three[index+2]))
                          {
			    case 'A': return 'A';
			    case 'C': return 'A';
			    case 'G': return 'A';
			    case 'T': return 'A';
			    case 'U': return 'A';
			    default: return 'Z';
			  }
		case 'G': switch(toupper(three[index+2]))
                          {
			    case 'A': return 'G';
			    case 'C': return 'G';
			    case 'G': return 'G';
			    case 'T': return 'G';
			    case 'U': return 'G';
			    default: return 'Z';
			  }
		case 'T': switch(toupper(three[index+2]))
                          {
			    case 'A': return 'V';
			    case 'C': return 'V';
			    case 'G': return 'V';
			    case 'T': return 'V';
			    case 'U': return 'V';
			    default: return 'Z';
			  }
		case 'U': switch(toupper(three[index+2]))
                          {
			    case 'A': return 'V';
			    case 'C': return 'V';
			    case 'G': return 'V';
			    case 'T': return 'V';
			    case 'U': return 'V';
			    default: return 'Z';
			  }
		default: return 'Z';
	      }
	      
    case 'T': switch(toupper(three[index+1]))
              {
                case 'A': switch(toupper(three[index+2]))
                          {
			    case 'A': return 'O';
			    case 'C': return 'Y';
			    case 'G': return 'O';
			    case 'T': return 'Y';
			    case 'U': return 'Y';
			    default: return 'Z';
			  }
		case 'C': switch(toupper(three[index+2]))
                          {
			    case 'A': return 'S';
			    case 'C': return 'S';
			    case 'G': return 'S';
			    case 'T': return 'S';
			    case 'U': return 'S';
			    default: return 'Z';
			  }
		case 'G': switch(toupper(three[index+2]))
                          {
			    case 'A': return 'O';
			    case 'C': return 'C';
			    case 'G': return 'W';
			    case 'T': return 'C';
			    case 'U': return 'C';
			    default: return 'Z';
			  }
		case 'T': switch(toupper(three[index+2]))
                          {
			    case 'A': return 'L';
			    case 'C': return 'F';
			    case 'G': return 'L';
			    case 'T': return 'F';
			    case 'U': return 'F';
			    default: return 'Z';
			  }
		case 'U': switch(toupper(three[index+2]))
                          {
			    case 'A': return 'L';
			    case 'C': return 'F';
			    case 'G': return 'L';
			    case 'T': return 'F';
			    case 'U': return 'F';
			    default: return 'Z';
			  }
		default: return 'Z';
	      }
	      
    case 'U': switch(toupper(three[index+1]))
              {
                case 'A': switch(toupper(three[index+2]))
                          {
			    case 'A': return 'O';
			    case 'C': return 'Y';
			    case 'G': return 'O';
			    case 'T': return 'Y';
			    case 'U': return 'Y';
			    default: return 'Z';
			  }
		case 'C': switch(toupper(three[index+2]))
                          {
			    case 'A': return 'S';
			    case 'C': return 'S';
			    case 'G': return 'S';
			    case 'T': return 'S';
			    case 'U': return 'S';
			    default: return 'Z';
			  }
		case 'G': switch(toupper(three[index+2]))
                          {
			    case 'A': return 'O';
			    case 'C': return 'C';
			    case 'G': return 'W';
			    case 'T': return 'C';
			    case 'U': return 'C';
			    default: return 'Z';
			  }
		case 'T': switch(toupper(three[index+2]))
                          {
			    case 'A': return 'L';
			    case 'C': return 'F';
			    case 'G': return 'L';
			    case 'T': return 'F';
			    case 'U': return 'F';
			    default: return 'Z';
			  }
		case 'U': switch(toupper(three[index+2]))
                          {
			    case 'A': return 'L';
			    case 'C': return 'F';
			    case 'G': return 'L';
			    case 'T': return 'F';
			    case 'U': return 'F';
			    default: return 'Z';
			  }
		default: return 'Z';
	      }
    
    default: return 'Z';
  }
}

