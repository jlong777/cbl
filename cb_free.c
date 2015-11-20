/* cb_free.c                                  http://cbl.sourceforge.net
 *
 * Included only for portability with Cray SV1 series UNICOS systems.
 * See original man page at the bottom, frees memory allocated with cb_malloc.
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
 * $Id: cb_free.c,v 1.4 2003/09/25 00:53:27 jlong777 Exp $
 */
 
#include <stdlib.h>

long cb_free(void *ptr)
{
  free(ptr);
  return 0;
}

/*
cb_free(3B)                                           Last changed: 01-31-03

NAME
        cb_free - frees memory allocated with cb_malloc

SYNOPSIS

        C/C++:

        #include <cbl.h>
        long cb_free(void *ptr);

IMPLEMENTATION

        Cray SV1 series UNICOS systems

DESCRIPTION

        cb_free releases memory that was allocated by cb_malloc, or allocated 
        by Fortran allocate statements inside other cbl routines. Memory 
        allocated with malloc() must not be deallocated with cb_free(), and 
        memory allocated with cb_malloc() must not be deallocated with free().
        Unlike the C library free() function, cb_free returns a status value
        of 0 is there was no error, and non-zero if corruption of the internal
        data structures was detected.

        ptr     a pointer to a block of memory previously allocated with
                cb_malloc or returned as an internally allocated argument
                in one of the cbl library routines.

NOTES       
        cb_free is single-threaded (i.e. not tasked) but contains updates of
        internal data structures that are not thread safe. If cb_free is called
        from within a parallel region it should be done from within a protected
        region.

        cb_free is intended for C programmers only. In a Fortran program, arrays
        allocated inside cbl library routines can be deallocated using the 
        deallocate statement.

SEE ALSO

        cb_malloc(3B), INTRO_LIBCBL(3B)
*/
