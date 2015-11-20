/* cb_malloc.c                                   http://cbl.sourceforge.net
 *
 * Included only for portability with Cray SV1 series UNICOS systems.
 * See original man page at the bottom, allocates block aligned memory region.
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
 * $Id: cb_malloc.c,v 1.4 2003/09/25 00:53:27 jlong777 Exp $
 */
 
#include <stdlib.h>

void * cb_malloc(long nbytes)
{
  return malloc(nbytes);
}


/*
cb_malloc(3B)                                           Last changed: 01-31-03

NAME
        cb_malloc - allocate block aligned memory region.

SYNOPSIS

        C/C++:

        #include <cbl.h>
        void *cb_malloc(long nbytes);

IMPLEMENTATION

        Cray SV1 series UNICOS systems

DESCRIPTION

        cb_malloc allocates memory using the Fortran allocate statement. On
        the SV1 this memory is aligned on 512 word address boundaries if the
        amount of memory requested is large. cb_malloc can be used to 
        allocate memory appropriate for the cb_ssd_* routines. Memory allocated
        using cb_malloc must be deallocated using cb_free(), and not free().
        cb_malloc returns a void pointer that should be cast to the appropriate
        data type, just as with malloc().

        nbytes  input number of bytes of memory to be allocated.

NOTES       
        cb_malloc is single-threaded (i.e. not tasked) but contains updates of
        internal data structures that are not thread safe. If cb_malloc is called
        from within a parallel region it should be done from within a protected
        region.

        cb_malloc is intended for C programmers only. In a Fortran program, arrays
        aligned arrays can be allocated using the allocate statement.

SEE ALSO

        cb_free(3B), INTRO_LIBCBL(3B)

*/
