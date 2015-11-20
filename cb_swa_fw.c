/* cb_swa_fw.c                                 http://cbl.sourceforge.net
 *
 * See original man page at the bottom.
 *
 * Copyright (C) 2004 University of Alaska Fairbanks
 * Institute of Arctic Biology (IAB)
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
 * $Id: cb_swa_fw.c,v 1.6 2005/01/13 22:51:01 jlong777 Exp $
 */

#include "cbl.h"
#include "cb_macro.h"
#include <stdio.h>
#include <stdlib.h>

void gotoh_score(long *, long *, long *, long, long, long, long, 
                 long *, long *, long *, long *, long *, long *);

void gotoh_align(long *, long *, long, long *, long, long, long **, long **, 
                 long **, long *, long *, long *, long *, long, long);

void cb_swa_fw(long *dbl, long dbllen, long *dbs, long dbslen, long eg, 
               long og, long *sslookup, long *smax, long **algl, long **algm, 
               long **algs, long *alglen, long *algstl, long *algsts, 
	       long *errno)
{
  /* note - this version has different internal routines than the one 
   *        described in the original Cray version, see the man page
   *        for the portable cb_swa_fw
   */
  
  long *swtab, *gaph, *gapv, eog, sav_i, sav_j, *tb;
  
  *errno = 0;
  
  tb = (long *) malloc((dbllen+1)*(dbslen/4+8));
  if(tb==0x0)
  {
    printf("cb_swa_fw: memory allocation error, returning...\n");
    *errno = 1;
    return;
  }
#ifdef LONG32
  gaph = (long *) calloc(dbslen+1, 4);
  if(gaph==0x0)
  {
    printf("cb_swa_fw: memory allocation error, returning...\n");
    *errno = 2;
    return;
  }
  gapv = (long *) calloc(dbslen+1, 4);
  if(gapv==0x0)
  {
    printf("cb_swa_fw: memory allocation error, returning...\n");
    *errno = 3;
    return;
  }
#endif
#ifdef LONG64
  gaph = (long *) calloc(dbslen+1, 8);
  if(gaph==0x0)
  {
    printf("cb_swa_fw: memory allocation error, returning...\n");
    *errno = 2;
    return;
  }
  gapv = (long *) calloc(dbslen+1, 8);
  if(gapv==0x0)
  {
    printf("cb_swa_fw: memory allocation error, returning...\n");
    *errno = 3;
    return;
  }
#endif

  eog = eg + og;

  gotoh_score(tb, gaph, gapv, dbllen, dbslen, eg, eog, smax, sslookup, dbl, dbs,
              &sav_i, &sav_j);
  gotoh_align(tb, dbl, dbllen, dbs, dbslen, *smax, algl, algm, algs, 
              alglen, algstl, algsts, errno, sav_i, sav_j);

  free(gaph);
  free(gapv);
  free(tb);

/* This commented-out code is the version that follows the Cray documentation
   for this routine. The code above for the portable version does not create
   an swtab matrix, only a traceback matrix using 2 bits per cell.

  long *swtab, *gaph, *gapv, eog;
#ifdef LONG32
  swtab = (long *) calloc((dbllen+1)*(dbslen+1), 4);
  if(swtab == 0x0)
  {
    printf("swtab allocation error, returning...\n");
    *smax = 0;
    return;
  }
  gaph  = (long *) calloc(dbslen+1, 4);
  if(gaph == 0x0)
  {
    printf("gapv allocation error, returning...\n");
    *smax = 0;
    free(swtab);
    return;
  }
  gapv  = (long *) calloc(dbslen+1, 4);
  if(gapv == 0x0)
  {
    printf("gapv allocation error, returning...\n");
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
    printf("swtab allocation error, returning...\n");
    *smax = 0;
    return;
  }
  gaph  = (long *) calloc(dbslen+1, 8);
  if(gaph == 0x0)
  {
    printf("gapv allocation error, returning...\n");
    *smax = 0;
    free(swtab);
    return;
  }
  gapv  = (long *) calloc(dbslen+1, 8);
  if(gapv == 0x0)
  {
    printf("gaph allocation error, returning...\n");
    *smax = 0;
    free(swtab);
    free(gaph);
    return;
  }
#endif
  eog = eg + og;
  
  cb_swa_fw_init(dbl, dbllen, dbs, dbslen, sslookup, swtab, gaph, gapv);
  cb_sw_fw_score(swtab, gaph, gapv, dbllen, dbslen, eg, eog, smax);
  cb_swa_fw_align(swtab, dbl, dbllen, dbs, dbslen, *smax, algl, algm, algs, 
             alglen, algstl, algsts, errno);
  free(swtab);
  free(gaph);
  free(gapv);
*/
}

void gotoh_score(long *tb, long *mgaph, long *mgapv, long dbllen, long dbslen, 
		 long eg, long eog, long *smax, long *sslookup, long *dbl, 
		 long *dbs, long *sav_i, long *sav_j)
{
  long i, j;
  register long max, max1, max2, max3, rmax=0, ri=0, rj=0, r;
  register unsigned long t, trace;
  register long e, g;

  unsigned long *temp = (unsigned long *) calloc(dbslen+1, sizeof(long));

  e = eg;
  g = eog;
  t = 0x0;
#ifdef LONG32
  r = dbslen/16 +1;
#endif
#ifdef LONG64 
  r = dbslen/32 +1;
#endif

  /* init */
  for(i=0; i<dbslen+1; i++)
  { 
    mgapv[i] = mgaph[i] = -g - (i-1)*e;
    temp[i] = 0x0;
  }
#ifdef LONG32
  for(i=0; i<dbslen/16+1; i++)
#endif
#ifdef LONG64
  for(i=0; i<dbslen/32+1; i++)
#endif
    tb[i] = 0x0;

  for(i=1; i<dbllen+1; i++)
  {
    max2 = 0 - g;

    for(j=1; j<dbslen+1; j++)
    {
      max =   0x0;
      trace = 0x0;
      max1 = temp[j-1] + 
             sslookup[(32*(0x1F & ((char *)dbs)[j-1]))+ (0x1F & ((char *)dbl)[i-1])];
      temp[j-1] = max2 + g;
      max3 = temp[j] - g;
      mgaph[j] = (max2 > mgaph[j-1]-e)? max2: mgaph[j-1]-e;
      mgapv[j] = (max3 > mgapv[j]-e)? max3: mgapv[j]-e;

      if(max1>max) max = max1; 
      if(mgaph[j]>max) max = mgaph[j]; 
      if(mgapv[j]>max) max = mgapv[j]; 
      if(max>=rmax)
      {
        rmax = max;
        ri = i;
        rj = j;
      }

      if(max)
      {
        if(max==mgaph[j]) trace = 0x2; /* vert (fortran), horiz (c) */
        if(max==mgapv[j]) trace = 0x1; /* horiz (fortran), vert (c) */
        if(max==max1)     trace = 0x3; /* came from corner */
#ifdef DEBUG_INFO
#ifdef LONG32
        if(max > 0x0FFFFFFF) 
	  printf("cb_sw_fw_score: cell score exceeds maximum...\n");
#endif
#ifdef LONG64
        if(max > 0x0FFFFFFFFFFFFFFF) 
          printf("cb_sw_fw_score: cell score exceeds maximum...\n");
#endif
#endif
      }
#ifdef LONG32 
      t |= (trace << (30-2*(j%16)));
      if((j%16) == 15)
      {
        tb[r*i + j/16] = t; /* tb is the 2-bit traceback matrix */
	t = 0x0;
      }
#endif
#ifdef LONG64
      t |= (trace << (62-2*(j%32))); /* 0x1F&j is j%32 */
      if((j%32) == 31)
      {
        tb[r*i + j/32] = t; /* (unsigned long)j>>5 is j/32 */
	t = 0x0;
      }
#endif
      max2 = max - g;
    }
#ifdef LONG32 
    if(dbslen%16 != 15)
      tb[r*i + dbslen/16] = t;
#endif
#ifdef LONG64
    if(dbslen%32 != 31)
      tb[r*i + dbslen/32] = t;
#endif
    t = 0x0;
    temp[dbslen] = max;
  }
  *smax  = rmax;
  *sav_i = ri;
  *sav_j = rj;
  free(temp);
}

void backfill(long *tb, long tb_i, long tb_j, long *dbl, long *dbs, long index,
              long **algl, long **algm, long **algs, long dbslen)
{
  char *algl_c = (char *)(*algl);
  char *algm_c = (char *)(*algm);
  char *algs_c = (char *)(*algs);

  if(tb_i==0 || tb_j==0 || index<0) return;

#ifdef LONG32
  switch((unsigned long)(tb[(dbslen/16 +1)*tb_i + tb_j/16] << 2*(tb_j%16)) >> 30)
#endif
#ifdef LONG64
  switch((unsigned long)(tb[(dbslen/32 +1)*tb_i + tb_j/32] << 2*(tb_j%32)) >> 62)
#endif
  {
    case 0x0:  /* end  */
               tb_i = 0;
	       tb_j = 0;
	       break;
    case 0x1:  /* horizontal (fortran), vertical (c) */
               algl_c[index] = ((char *)dbl)[tb_i-1];
               algm_c[index] = ' ';
               algs_c[index] = '-';
	       tb_i--;
               break;
    case 0x2:  /* vertical (fortran), horizontal (c) */
               algl_c[index] = '-';
               algm_c[index] = ' ';
               algs_c[index] = ((char *)dbs)[tb_j-1];
	       tb_j--;
	       break;
    default: /* corner */
             algl_c[index] = ((char *)dbl)[tb_i-1];
             algs_c[index] = ((char *)dbs)[tb_j-1];
	     if(algl_c[index] == algs_c[index])
               algm_c[index] = ':';
	     else
	       algm_c[index] = ' ';
	     tb_i--;
	     tb_j--;
             break;
  }

  index--;

  backfill(tb, tb_i, tb_j, dbl, dbs, index, algl, algm, algs, dbslen);
}

void gotoh_align(long *tb, long *dbl, long dbllen, long *dbs, long dbslen,
                 long smax, long **algl, long **algm, long **algs, long *alglen,
                 long *algstl, long *algsts, long *errno, long sav_i, long sav_j)
{
  long i, index=0, j;
  register unsigned long blah, r;

#ifdef LONG32
  r = dbslen/16 +1;
#endif
#ifdef LONG64 
  r = dbslen/32 +1;
#endif

  i = sav_i;
  j = sav_j;
  *alglen = 0;
#ifdef LONG32
  while(blah=((unsigned long)(tb[r*i + j/16] << 2*(j%16)) >> 30))
#endif
#ifdef LONG64 
  while(blah=((unsigned long)(tb[r*i + j/32] << 2*(j%32)) >> 62))
#endif
    switch(blah)
    {
      case 0x1:  i--;
	         (*alglen)++;
		 /**algstl--;*/
		 break;
      case 0x2:  j--;
                 (*alglen)++;
		 /**algsts--;*/
                 break;
      default:   i--;
                 j--;
		 /**algstl--;
		 *algsts--;*/
                 (*alglen)++;
	         break;
    }
/*printf("alglen = %d  smax = %d\n", *alglen, smax);*/

  /**algstl = sav_i - *alglen;
  *algsts = sav_j - *alglen;*/
  *algstl = i+1;
  *algsts = j+1;

  *algl = (long *) malloc(*alglen+8);
  *algm = (long *) malloc(*alglen+8);
  *algs = (long *) malloc(*alglen+8);
  if(*algl==0 || *algm==0 || *algs==0)
  {
    printf("cb_swa_fw: malloc error in align, returning...\n");
    return;
  }

  index = *alglen-1;

  /* recursively trace back through scoring matrix for alignment */
  backfill(tb, sav_i, sav_j, dbl, dbs, index, algl, algm, algs, dbslen);
}


/*
cb_swa_fw(3B)                                           Last changed: 09-26-03

NAME

        cb_swa_fw,  cb_swa_fw_init,  cb_swa_fw_align, cb_sw_fw_score,
        cb_swn_fw,  cb_swn_fw_init,  cb_swn_fw_align,
        cb_swn4_fw, cb_swn4_fw_init, cb_swn4_fw_align   - compute 
        Smith-Waterman alignment with full word scoring matrix

SYNOPSIS

        C/C++:

        #include <cbl.h>
        void cb_swa_fw( long *dbl, long dbllen, long *dbs, long dbslen, long eg,long og,
                        long *sslookup, long *smax,long **algl, long **algm, long **algs,
                        long *alglen,long *algstl,long *algsts,long *errno);

        void cb_swa_fw_align( long *swtab, long *dbl, long dbllen, long *dbs,long dbslen,
                              long smax, long **algl, long **algm,long **algs,long *alglen,
                              long *algstl, long *algsts, long *errno);

        void cb_swa_fw_init( long *dbl, long dbllen, long *dbs, long dbslen, 
                             long *sslookup,long *swtab, long *gaph,long *gapv);

        void cb_swn_fw( long *dbl, long dbllen, long *dbs, long dbslen, long eg,long og,
                        long *sslookup, long *smax,long **algl, long **algm, long **algs,
                        long *alglen,long *algstl,long *algsts,long *errno);

        void cb_swn_fw_align( long *swtab, long *dbl, long dbllen, long *dbs,long dbslen,
                              long smax, long **algl, long **algm,long **algs,long *alglen,
                              long *algstl, long *algsts, long *errno);

        void cb_swn_fw_init( long *dbl, long dbllen, long *dbs, long dbslen, 
                             long *sslookup,long *swtab, long *gaph,long *gapv);

        void cb_swn4_fw( long *dbl, long dbllen, long *dbs, long dbslen, long eg,long og,
                        long *sslookup, long *smax,long **algl, long **algm, long **algs,
                        long *alglen,long *algstl,long *algsts,long *errno);

        void cb_swn4_fw_align( long *swtab, long *dbl, long dbllen, long *dbs,long dbslen,
                              long smax, long **algl, long **algm,long **algs,long *alglen,
                              long *algstl, long *algsts, long *errno);

        void cb_swn4_fw_init( long *dbl, long dbllen, long *dbs, long dbslen, 
                             long *sslookup,long *swtab, long *gaph,long *gapv);

        void cb_sw_fw_score( long *swtab, long *gaph, long *gapv, long dbllen, long dbslen,
                             long eg,long eog,long *smax);


        Fortran:

        use cb_sw
        call cb_swa_fw(dbl, dbllen, dbs, dbslen, eg, og, sslookup, smax, algl, algm, algs,
                       alglen, algstl, algsts, errno)

        call cb_swa_fw_align(swtab, dbl, dbllen, dbs, dbslen, smax, algl, algm, algs,
                             alglen, algstl, algsts, errno)

        call cb_swa_fw_init(dbl, dbllen, dbs, dbslen, sslookup, swtab, gaph, gapv)

        call cb_swn_fw(dbl, dbllen, dbs, dbslen, eg, og, sslookup, smax, algl, algm, algs,
                       alglen, algstl, algsts, errno)

        call cb_swn_fw_align(swtab, dbl, dbllen, dbs, dbslen, smax, algl, algm, algs,
                             alglen, algstl, algsts, errno)

        call cb_swn_fw_init(dbl, dbllen, dbs, dbslen, sslookup, swtab, gaph, gapv)

        call cb_swn4_fw(dbl, dbllen, dbs, dbslen, eg, og, sslookup, smax, algl, algm, algs,
                       alglen, algstl, algsts, errno)

        call cb_swn4_fw_align(swtab, dbl, dbllen, dbs, dbslen, smax, algl, algm, algs,
                             alglen, algstl, algsts, errno)

        call cb_swn4_fw_init(dbl, dbllen, dbs, dbslen, sslookup, swtab, gaph, gapv)

        call cb_sw_fw_score(swtab, gaph, gapv, dbllen, dbslen, eg, eog, smax)


IMPLEMENTATION

        UNICOS/mp and Cray SV1 series UNICOS systems

DESCRIPTION

        cb_swa_fw accepts two strings of ASCII characters, dbl and dbs,
        which nominally represent amino acids though the actual meaning
        is determined by the contents of the substitution score lookup
        table sslookup. The substitution scores and the extend gap penalty,
        eg, and the open gap penalty, og, are used to compute the Smith-Waterman
        aligment of the two strings.  This is performed in three steps,
        one step each for cb_swa_fw_init, cb_sw_fw_score, and cb_swa_fw_align.
        cb_swa_fw is a wrapper routine that calls these three routines and
        handles internal allocation and deallocation of work space.

        cb_swa_fw_init accepts two strings, dbl and dbs, and a substitution
        score lookup table, sslookup, and fills the table swtab with the 
        substitution scores in rows and columns corresponding to the entries
        in the dbs and dbl strings.  The arrays gapv and gaph are initialized
        to zero. Row and column zero of the swtab are initialized to zero.

        cb_sw_fw_score accepts the swtab, gaph, and gapv arrays created by
        cb_swa_fw_init, along with open gap and extend gap penalty values, and
        replaces the contents of swtab with the Smith-Waterman scores and 
        traceback information for each cell. It also returns the value in the 
        largest score in smax.

        cb_swa_fw_align accepts a completed Smith-Waterman score table as created
        by cb_sw_fw_score and a score, smax, contained in the table. It returns an 
        alignment that ends on a cell with the score value supplied. In cb_swa_fw 
        this score is the maximum score returned by cb_sw_fw_score.

        The cb_swn_* routines are the same as the corresponding cb_swa_* routines
        except that the input strings, dbl and dbs, contain 2-bit nucleotide
        data. See cb_compress, mode=2, for details on the 2-bit format for
        packed nucleotide data.

        The cb_swn4_* routines are the same as the corresponding cb_swa_* routines
        except that the input strings, dbl and dbs, contain 4-bit nucleotide
        data. See cb_compress, mode=4, for details on the 4-bit format for
        packed nucleotide data.


        dbl     input string characters representing the longer of the two 
                sequences to be compared. For the cb_swa_* routines the
                input is ASCII characters, packed 8 per word. For the cb_swn_*
                routines the input is 2-bit nucleotide characters, packed
                32 per word. For the cb_swn4_* routines the input is 4-bit
                nucleotide characters, packed 16 per word. In Fortran, dbl 
                should be an INTEGER(8) array.

        dbllen  input number of characters in dbl. In Fortran, dbllen should be 
                an INTEGER(8) variable, constant, or expression.

        dbs     input string characters representing the shorter of the two 
                sequences to be compared. For the cb_swa_* routines the
                input is ASCII characters, packed 8 per word. For the cb_swn_*
                routines the input is 2-bit nucleotide characters, packed
                32 per word. For the cb_swn4_* routines the input is 4-bit
                nucleotide characters, packed 16 per word. In Fortran, dbs 
                should be an INTEGER(8) array.

        dbslen  input number of characters in dbs. In Fortran, dbslen should be 
                an INTEGER(8) variable, constant, or expression.

        eg      input score penalty for extending a gap. For cb_sw_fw_score,
                eg must be zero or positive.  In Fortran, eg should be
                an INTEGER(8) variable, constant, or expression.

        og      input score penalty for opening a gap.  If eg and og 
                are both non-zero, they should have the same sign.
                In Fortran, og should be an INTEGER(8) variable, constant, or
                expression.

        eog     input combined gap open and gap extend penalty as a zero or
                positive integer. In Fortran, eog should be an INTEGER(8) variable,
                constant, or expression.

        sslookup  input array containing the substitution matrix for 
                the type of data represented by the input strings.
                sslookup(i,j) (Fortran) or sslookup[j][i] (C) is the substitution 
                score for data entities represented by the i'th and j'th letters 
                of the alphabet or packed nucleotide codes equal to i and j.
                The declaration for the array should be

                Fortran                                 C/C++

                For cb_swa_fw and cb_swa_fw_init:

                integer(8) :: sslookup(0:31, 0:26)      long sslookup[27][32];

                For cb_swn_fw and cb_swn_fw_init:

                integer(8) :: sslookup(0:3, 0:3)        long sslookup[4][4];

                For cb_swn4_fw and cb_swn4_fw_init:

                integer(8) :: sslookup(0:15, 0:15)      long sslookup[16][16];

                Note that in the cb_swa_* case the leading size is 32 for 
                performance reasons, even though only 26 entries are used. 
                The entries in column and row zero are not used. Entries must 
                be defined for each letter combination that might appear in 
                the input strings.  The letter 'A' corresponds to a subscript
                value of 1.

                For the cb_swn_* and cb_swn4_* cases, the encoded 2-bit and
                4-bit values correspond to the subscripts used in sslookup.

        smax    the largest cell score computed by cb_sw_fw_score, In Fortran
                smax should be a INTEGER(8) variable.

        algl    an output array containing letters from dbl and "-" characters to 
                indicate gaps, corresponding to the alignment that ends at the
                location of the score smax in the score table. If the same 
                score appears more than once in the table, the one lowest and 
                rightmost in the table is used. In Fortran, algl should be specified
                as an allocatable array: INTEGER(8),DIMENSION(:),ALLOCATABLE. The
                memory for algl is allocated in cb_sw*_fw_align. Algl should not
                be allocated before calling either cb_sw*_fw_align or cb_sw*_fw.
                In C, deallocation of the memory associated with algl must be
                done with cb_free().

        algs    an output array similar to algl, but containing the alignment
                characters from the dbs string.

        algm    an output array similar to algl. algm contains a ":" character 
                in every location where the corresponding characters in algl
                and algs are the same (independent of case), and a " " (space)
                character where the corresponding characters in algl and algs 
                are different.

        alglen  output length, in characters, of the strings in algl, algs, and
                algm. In Fortran, alglen should be an INTEGER(8) variable.

        algstl  output location in dbl corresponding to the first character in
                the algl array. The positions in dbl are numbered starting at 1.
                In Fortran, algstl should be an INTEGER(8) variable.

        algsts  output location in dbs corresponding to the first character in
                the algs array. The positions in dbs are numbered starting at 1.
                In Fortran, algsts should be an INTEGER(8) variable.

        errno   output error number.  If errno is 0, no error occurred.  If
                errno is non-zero, an error in memory allocation occurred.

        swtab   output (cb_sw*_fw_init) or input (cb_sw_fw_score) table of 
                substitution scores as described above.
                output (cb_sw_fw_score) or input (cb_sw*_fw_align)  table of 
                Smith-Waterman cell scores. These are zero or positive by 
                construction, and are stored in the lower 60 bits (bits number 59-0) 
                of each word. The upper bits of the word are set as follows:

                bit 60 = 1 if score was computed as a horizontal gap extension
        
                bit 61 = 1 if score was computed as a vertical gap extension

                bit 62 = 1 if score was computed as a sequence extension from
                           the cell diagonally up and to the let.

                It is possible for more than one of the three bits to be set 
                if an equal score could be computed by more than one path.
                If the score is zero, bits 60-62 are all set to zero.

                Fortran declaration : INTEGER(8) :: swtab(0:dbslen,0:dbllen)
                C/C++ declarattion  : long swtab[dbllen+1][dbslen+1]

                The memory for swtab must be allocated before calling cb_sw*_fw_init.
                swtab is allocated internally in cb_sw*_fw.

        gaph    work array of (dbslen+1) elements used to hold accumulated gap 
                penalties in the horizontal direction. In Fortran, gaph should
                be an INTEGER(8) array. The memory for gaph must be allocated before
                calling cb_sw*_fw_init.  gaph is allocated internaly in cb_sw*_fw.

        gapv    work array of (dbslen+1) elements used to hold accumulated gap 
                penalties in the vertical direction. In Fortran, gapv should
                be an INTEGER(8) array. The memory for gapv must be allocated before
                calling cb_sw*_fw_init.  gapv is allocated internaly in cb_sw*_fw.


NOTES       

        For the cb_swa_fw and cb_swa_fw_init (ASCII input strings) routines,
        the expected format of the substitution matrix is alphabetical.  If the entries
        represent amino acids with the usual naming conventions, then the entries
        involving O and U are not used and can be set to any value.  The "*" character
        maps to the letter "J", which is otherwise not used. Entries corresponding to 
        "*" should be placed in the J row and column.

        cb_sw*_fw, cb_sw*_fw_init, cb_sw_fw_score, and cb_sw*_fw_align are 
        single-threaded (i.e. not tasked) and may be called from within a parallel 
        region.
               
        Some of the Smith-Waterman routines dynamically allocate memory for the
        user supplied arguments.  If this memory is to be deallocated later, the
        deallocation must be done properly. In Fortran, the variables should be
        deallocated with the deallocate statement.  In C, the variables must be
        deallocated using cb_free() and not the C library free function.

        cb_sw*_fw and cb_sw*_fw_align relpace the contents of the bmm register.

SEE ALSO

        cb_compress(3B), cb_free(3B), INTRO_LIBCBL(3B)

        This man page is available only online.
*/
