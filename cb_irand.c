/* cb_irand.c                                    http://cbl.sourceforge.net
 *
 * generates an array of words with random bits
 * see original man page at the bottom. This code is the generic version
 * with provision for 32- and 64-bit long ints.
 *
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
 * $Id: cb_irand.c,v 1.4 2003/09/25 00:53:27 jlong777 Exp $
 */

#include "cb_macro.h"
#include <stdlib.h>

void cb_irand(long *array, long n)
{
  long i;
  
  for(i=0; i<n; i++)
  {
#ifdef LONG32
    array[i] = ((rand() & 0x003FFFC0) << 10) | ((rand() & 0x001FFFE0) >> 5);
#endif

#ifdef LONG64
    array[i] = ((unsigned long)(rand() & 0x003FFFC0) << 42) | 
               ((unsigned long)(rand() & 0x001FFFE0) << 27) |
	       ((unsigned long)(rand() & 0x000FFFF0) << 12) |  
	       ((unsigned long)(rand() & 0x0007FFF8) >>  3);
#endif
  }
}

/*
cb_irand(3B)                                           Last changed: 09-17-02

NAME
        cb_irand - generates a list of random bits

SYNOPSIS

        C/C++:

        #include <cbl.h>
        void cb_irand( long *array, long n );

        Fortran:

        use cb_bits
        call cb_irand( array, n)


IMPLEMENTATION

        Cray SV1 series UNICOS systems

DESCRIPTION

        cb_irand generates a list of 64-bit words with random bit patterns
        in each word. For each word of output, two random real values are
        obtained from the system random number generator. From the first
        real value bits 20-51 are extracted and used for the upper 32 bits 
        of the result. From the second real value bits 23-54 are extracted
        and used for the lower 32 bits of the result.

        array   (output) returned array of 64-bit word of random bits.
                For Fortran, array should be an INTEGER(8) array of 
                size at least n words. Space for array must be allocated
                in the caller.

         n      (input) the number of elements in array to be filled
                with random bits.

NOTES       
        Internally cb_irand uses the Fortran intrinsic random_number. The
        seed value, and hence the sequence of values, can be modified with
        the random_seed intrinsic.

        cb_irand replaces the contents of the bmm register.

SEE ALSO

        RANDOM_NUMBER(3I), RANDOM_SEED(3I), INTRO_LIBCBL(3B)


*/
