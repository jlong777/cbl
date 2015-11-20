/* cb_fasta_convert.c                            http://cbl.sourceforge.net
 *
 * restructure the memory image of a FASTA format file.
 * see original man page at the bottom.
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
 * $Id: cb_fasta_convert.c,v 1.4 2003/09/25 00:53:27 jlong777 Exp $
 */

#include "cbl.h"
#include "cb_macro.h"

void cb_fasta_convert(long *db, long dblen, long *dat, long datlen,
                      long *hdr, long hdrlen, ptrinfo *ptrs,
                      long ptrslen, long *nsegs, long *errno )
{
  int count, count2, i=0, j=0, k=0, pi=0;
  char *ch, *dch, *hch;
  
  ch  = (char *)db;
  dch = (char *)dat;
  hch = (char *)hdr;
  *nsegs = 0;
  
  while(i < dblen)
  {
    while((ch[i] != '>') && (i < dblen)) i++; /* go to ">" */
    
    if(ch[i] == '>')
    {
      if(pi>=ptrslen)
      {
        *errno = 3;
        return;
      }
      ptrs[pi].hdrstart = j+1;
      (*nsegs)++;
    }
    else if(i>=dblen)
    {
      if(pi>=ptrslen)
      {
        *errno = 4;
        return;
      }
      ptrs[pi].hdrstart = j+1;
      *errno = 0;
      return;
    }
    
    /* header record */
    while((ch[i]!='\n') && (i<dblen))
    {
      if(j>=hdrlen) 
      {
        *errno = 2;
        if(pi<ptrslen)
          ptrs[pi].hdrstart = j+1;
        return;
      }
      hch[j++] = ch[i++];
    }
    
    if(ch[i] == '\n') i++; /* go to data */
    else if(i>=dblen)
    {
      if(pi>=ptrslen)
      {
        *errno = 4;
        return;
      }
      ptrs[pi].hdrstart = j+1;
      *errno = 0;
      return;
    }

#ifdef LONG32
    ptrs[pi].datstart = k/4+1;
#endif
#ifdef LONG64
    ptrs[pi].datstart = k/8+1;
#endif

    count = 0;
    
    /* data record */
    while((ch[i]!='>') && (i<dblen))
    {
      if(k>datlen)
      {
        *errno = 1;
        if(pi<ptrslen)
          ptrs[pi].hdrstart = j+1;
        return;
      }
      if(ch[i]=='\n') i++;
      else
      {
        dch[k++] = ch[i++];
        count++;
      }
    }
    
    if(count)
    {
      ptrs[pi++].datleng = count;
      
      /* pad data array w/0's so it is a multiple of 4 words in length */
      count2 = count;
#ifdef LONG32
      while(count2%16)
      {
        if(k>=datlen)
        {
          *errno = 1;
          if(pi<ptrslen)
            ptrs[pi].hdrstart = j+1;
          return;
        }
        dch[k++] = 0x0;
        count2++;
      }
#endif
#ifdef LONG64
      while(count2%32)
      {
        if(k>=datlen)
        {
          *errno = 1;
          if(pi<ptrslen)
            ptrs[pi].hdrstart = j+1;
          return;
        }
        dch[k++] = 0x0;
        count2++;
      }
#endif
      if(k>datlen)
      {
        *errno = 1;
        if(pi<ptrslen)
          ptrs[pi].hdrstart = j+1;
        return;
      }
    }
  }
    
  if(pi>=ptrslen)
  {
    *errno = 4;
    return;
  }
  ptrs[pi].hdrstart = j+1;
  
  *errno = 0;
}

/*
cb_fasta_convert(3B)                                           Last changed: 01-31-03

NAME

        cb_fasta_convert - restructure the memory image of a FASTA format file
SYNOPSIS

        C/C++:

        #include <cbl.h>
        void cb_fasta_convert( long *db, long dblen, long *dat, long datlen,
                               long *hdr, long hdrlen, ptrinfo *ptrs,
                               long ptrslen, long *nsegs, long *errno );

        Fortran:

        use cb_fasta
        call cb_fasta_convert( db, dblen, dat, datlen, hdr, hdrlen,  &
                               ptrs, ptrslen, nsegs, errno)


IMPLEMENTATION

        Cray SV1 series UNICOS systems

DESCRIPTION

        cb_fasta_convert extracts and organizes data contained in a memory image
        of a FASTA style data file.  The image is assumed to be composed of lines
        consisting of ASCII characters terminated be a newline character with the
        decimal value of 10. The image contains two types of lines, header lines 
        and data lines. Header lines are distinguished by having ">" as the first 
        character of the line. This memory image may be created by transferring 
        an entire external file as a single unformatted binary record from disk 
        to memory. Thus, the usual character based read operations are replaced by 
        a combination of a binary read and formatting done by cb_fasta_convert. 
        This two step process is significantly faster, especially for large files.

        Each header line defines the beginning of a new block consisting of that
        header line and the subsequent data line(s) up to the next header line or
        the end of the image. The header line text is transferred to the hdr array
        and text from the associated data line(s) is transferred to the dat array.
        Pointers to the beginning of the text for each block, and the data text 
        length are stored in the ptrs array.

        db      input memory image to be reformatted. In Fortran, db should 
                be an INTEGER(8) array of length (dblen+7)/8.

        dblen   input number of characters in the db array. In Fortran, dblen
                should be an INTEGER(8) variable, constant, or expression.

        dat     output array holding the text of the data lines, excluding the
                trailing newline characters. Each block starts on a word
                boundary within dat. Each block is padded with zero bits if 
                necessary to fill a multiple of 4 words in dat. The padding allows
                the entire dat array to be compressed using either two or four
                bit compression and have the resulting blocks still word aligned.
                The starting location in dat of the data text for each block is
                specified by the datstart component of the ptrs structure for that
                block.  The number of characters in the text, excluding any pad bits,
                is specified by the datleng component of the same ptrs structure.
                The memory for dat must be allocated before calling cb_fasta_convert.
                In Fortran, dat should be an INTEGER(8) array of length (datlen+7)/8.

        datlen  input allocated size of the dat array in units of characters.
                In Fortran, datlen should be an INTEGER(8) variable, constant,
                or expression.

        hdr     output array holding the text of the header lines, including the
                ">" characters but excluding the trailing newline characters, 
                concatenated together. The starting location in hdr of the header 
                text for each block is specified by the hdrstart component of the 
                ptrs structure for that block. The memory for hdr must be allocated
                before calling cb_fasta_convert. In Fortran, hdr should be declared 
                as an INTEGER(8) array of length (hdrlen+7)/8.

        hdrlen  input allocated size of the hdr array in units of characters.
                In Fortran, hdrlen should be an INTEGER(8) variable, constant,
                or expression.

        ptrs    output array of structures that describe the locations in the
                dat and hdr arrays of the text associated with each block. The 
                memory for ptrs must be allocated before calling cb_fasta_convert.
                The structure contains three values:

                hdrstart - character number in hdr where the header text
                           starts. The first location is hdrstart=1.

                datstart - word number in dat where the data text starts.
                           The first location is datstart=1.

                datleng  - number of characters in the data segment

                An additional structure is written at the end of the ptrs
                array.  This has the hdrstart element set to the location
                in the hdr array one past the end of the final header. This 
                provides a simple way for users to find the end of the
                final header.

                In Fortran, ptrs should be a TYPE(PTRINFO) array of length 
                ptrslen. The type definition for PTRINFO is included in the
                cb_fasta module. For C, the typedef for prtinfo is supplied
                in the cbl.h header file.

        ptrslen input indicating that the allocated size of the ptrs array is
                3*ptrslen long words, or 24*ptrslen bytes. In Fortran ptrslen
                should be an INTEGER(8) variable, constant, or expression that
                specifies the allocated length of the prts array. This value
                needs to be at least one greater than the total number of data
                blocks in the db array.

        nsegs   output number of actual entries written to the ptrs array, equal
                to the number of header lines, and hence blocks, in the original 
                file. In Fortran nsegs should be an INTEGER(8) variable.

        errno   output error status set to one of these values:

                0 - no error
                1 - not enough space in dat array
                2 - not enough space in hdr array
                3 - not enough space in ptrs array
                4 - not enough space for the final ptrs array element

                If errno is returned with a non-zero value, the data in db was
                not completely processed. In Fortran errno should be an INTEGER(8) 
                variable.

NOTES       
        cb_fasta_convert is single-threaded (i.e. not tasked) and may be called 
        from within a parallel region.

        A wrapper routine, cb_read_fasta, is available that combines the file I/O, 
        memory allocations, and reformatting using cb_fasta_convert.

SEE ALSO

        cb_copy_bits(3B), cb_compress(3B), cb_read_fasta(3B), fread(3C),
        INTRO_LIBCBL(3B)


*/
