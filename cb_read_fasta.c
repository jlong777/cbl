/* cb_read_fasta.c                               http://cbl.sourceforge.net
 *
 * loads data from a FASTA file into memory arrays
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
 * $Id: cb_read_fasta.c,v 1.5 2004/03/09 21:10:46 jlong777 Exp $
 */
 
#include "cbl.h"
#include "cb_macro.h"

#include <stdio.h>
#include <stdlib.h>
#include <sys/stat.h>
#include <sys/types.h>
#include <unistd.h>

void cb_read_fasta(char *file, long lun, long **dat, long **hdr, 
                   ptrinfo **ptrs, long *nsegs, long *errno)
{
  long i, datsz=0, hdrsz=0;
  char *temp;
  long *db;
  struct stat sbuf;
  FILE *infile;
  
  /* get file size */
  if(stat(file, &sbuf) == -1)
  {
    fprintf(stderr, "cb_read_fasta: error getting file stats for %s, returning...\n", file);
    *errno = 1;
    return;
  }
#ifdef LONG32
  db = (long *) malloc(((sbuf.st_size+3)/4)*4);
#endif
#ifdef LONG64
  db = (long *) malloc(((sbuf.st_size+7)/8)*8);
#endif

  if(db==NULL)
  {
    fprintf(stderr, "cb_read_fasta: malloc error, returning...\n");
    *errno = 1;
    return;
  }
  
  temp = (char *)db;
  
  infile = fopen(file, "r");
  if(infile==NULL)
  {
    fprintf(stderr, "cb_read_fasta: error opening file %s, returning...\n", file);
    *errno = 1;
    return;
  }
  *errno = fread(temp, 1, sbuf.st_size, infile);
  if(*errno != sbuf.st_size)
  {
    fprintf(stderr, "cb_read_fasta: error reading file %s, returning...\n", file);
    *errno = 1;
    return;
  }
  fclose(infile);
  
  /* count the number of header lines and sizes needed for dat & hdr */
  i=0;
  *nsegs=0;
  while((temp[i]!='>') && (i<sbuf.st_size)) i++; /* start at first ">" */
  
  while(i < sbuf.st_size)
  {
    while((temp[i] != '\n') && (i < sbuf.st_size))
    {
      hdrsz++;
      i++;
    }
    
    (*nsegs)++; /* use nsegs to track number of headers seen */
    i++;

    while((temp[i] != '>') && (i < sbuf.st_size))
    {
      if(temp[i] == '\n') i++;
      else
      {
        datsz++;
        i++;
      }
    }
    /* make sure datsz (bytes) is always a multipe of 4 words */
#ifdef LONG32
    while(datsz%16) datsz++;
#endif
#ifdef LONG64
    while(datsz%32) datsz++;
#endif
  }

  /* allocate memory for the dat, hdr, and ptrs arrays */
  *dat = (long *) malloc(datsz);
#ifdef LONG32
  *hdr = (long *) malloc(((hdrsz+3)/4)*4);
#endif
#ifdef LONG64
  *hdr = (long *) malloc(((hdrsz+7)/8)*8);
#endif
  *ptrs = (ptrinfo *) calloc(*nsegs+1, sizeof(ptrinfo));
  
  if(*dat==NULL || *hdr==NULL || *ptrs==NULL)
  {
    fprintf(stderr, "cb_read_fasta: malloc error, returning...\n");
    *errno = 1;
    return;
  }
  
  /* now transfer the contents of temp to dat & hdr */
#ifdef LONG32
  cb_fasta_convert(db,    sbuf.st_size,
                  *dat,   datsz, 
                  *hdr, ((hdrsz+3)/4)*4, 
                  *ptrs, *nsegs+1, nsegs, errno);
#endif
#ifdef LONG64
  cb_fasta_convert(db,    sbuf.st_size,
                  *dat,   datsz, 
                  *hdr, ((hdrsz+7)/8)*8,
                  *ptrs, *nsegs+1, nsegs, errno);
#endif

  free(db);
  
  if(*errno)
  {
    fprintf(stderr, "cb_read_fasta: cb_fasta_convert error\n");
    *errno = 1;
  }
}

 
/*
cb_read_fasta(3B)                                        Last changed: 02-01-03

NAME

        cb_read_fasta - loads data from a FASTA file into memory arrays

SYNOPSIS

        C/C++:

        #include <cbl.h>
        void cb_read_fasta( char *file, long lun, long **seq, long **hdr, 
                            ptrinfo **ptrs, long *nsegs, long *error);

        Fortran:

        use cb_fasta
        call cb_read_fasta( file, lun, dat, hdr, ptrs, nsegs, error )


IMPLEMENTATION

        Cray SV1 series UNICOS systems

DESCRIPTION

        cb_read_fasta processes data from a file containing header lines and
        data lines of ASCII letters representing nucleotides or amino acids
        in the FASTA format. See the man page for cb_fasta_convert for more
        detail on the assumed file format.  cb_read_fasta transfers the file
        to memory as a single record, counts the number of header lines,
        allocates memory for the dat, hdr, and ptrs arrays, and moves the data
        and header information from the file image to the output arrays. The
        last step is done by a call to cb_fasta_convert.

        Each header line defines the beginning of a new block consisting of that
        header line and the subsequent data line(s) up to the next header line or
        the end of the image. The header line text is transferred to the hdr array
        and text from the associated data line(s) is transferred to the dat array.
        Pointers to the beginning of the text for each block, and the data text 
        length are stored in the ptrs array.

        file    input name of the file to be read. In Fortran, file should be  
                specified as a CHARACTER variable, constant, or expression in 
                the caller. In C, file should be a null terminated string with
                fewer than 1024 characters.

        lun     input value to be used as the Fortran unit number for the file 
                open and read operations.  This number should be different from 
                the unit numbers for any other currently open files.  If the call 
                to cb_read_fasta is in a parallel region, each instance of the call 
                should specify a different value for lun.  In Fortran, lun should 
                be an INTEGER(8) variable, consant, or expression.  In C, values of
                lun between 110 and 2047 should avoid any internal library conflicts.

        dat     output array holding the text of the data lines from the file. 
                This array should be unallocated when cb_read_fasta is called. The 
                memory for dat is allocated inside the routine. Each block of text 
                data starts on a word boundary within dat. Each block is padded with 
                zero bits if necessary to fill a multiple of 4 words in dat. The 
                padding allows the entire dat array to be compressed using either 
                two or four bit compression and have the resulting blocks still word 
                aligned. The starting location in dat of the data text for each block 
                is specified by the datstart component of the ptrs structure for that
                block.  The number of characters in the text, excluding any pad bits,
                is specified by the datleng component of the same ptrs structure.
                In Fortran, dat should be declared as a rank 1 allocatable array:
                INTEGER(8),DIMENSION(:),ALLOCATABLE. In C, deallocating the memory 
                for dat must be done with cb_free. See note below.

        hdr     output array holding the text of the header lines, including the ">" 
                characters, from the file. This array should be unallocated when 
                cb_read_fasta is called. The memory for hdr is allocated inside the 
                routine. If there is more than one block in the file, the header 
                lines are concatenated together in hdr. The starting location in 
                hdr of the header text for each block is specified by the hdrstart 
                component of the ptrs structure for that block. 
                In Fortran, hdr should be declared as a rank 1 allocatable array:
                INTEGER(8),DIMENSION(:),ALLOCATABLE. In C, deallocating the memory 
                for hdr must be done with cb_free. See note below.

        ptrs    output array of structures that describe the locations in the dat 
                and hdr arrays of the text associated with each block. The array should 
                be unallocated when cb_read_fasta is called. The memory for ptrs is 
                allocated inside the routine.

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

                In Fortran, ptrs should be declared as a rank 1 allocatable array:
                TYPE(PTRINFO),DIMENSION(:),ALLOCATABLE. The type definition for 
                PTRINFO is included in the cb_fasta module. For C, the typedef for
                ptrinfo is included in cbl.h. In C, deallocating the memory for ptrs 
                must be done with cb_free. See note below.

        nsegs   output number of entries written to the ptrs array, equal to the 
                number of header lines, and hence blocks, in the original file. 
                In Fortran, nsegs should be an INTEGER(8) variable.

        error   output error indicator. In Fortran, error should be a LOGICAL(8)
                variable. It is set to .true. if an error occurred in
                cb_read_fasta, and set to .false. if no error occurred. In C,
                error is set to 1 if an error occurred, and 0 otherwise.

NOTES       
        cb_read_fasta is single-threaded (i.e. not tasked) and may be called 
        from within a parallel region. Each parallel caller should specify a
        unique value for lun.

        cb_read_fasta allocates the output arrays dat, hdr, and ptrs, using
        the Fortan allocate statement.  This provides automatic memory block
        alignment for large arrays. Memory allocated in this manner cannot be
        deallocated with the C free() function.  The cb_free() function is supplied
        in this library to allow these arrays to be properly deallocated in a
        C function. See the man page for cb_free for details.

SEE ALSO

        cb_fasta_convert(3B), cb_free(3B), INTRO_LIBCBL(3B)

*/
