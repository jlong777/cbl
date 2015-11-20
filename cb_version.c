/* cb_version.c                                  http://cbl.sourceforge.net
 *
 * returns the version number of libcbl
 * see original man page at the bottom. This code is the generic version.
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
 * $Id: cb_version.c,v 1.12 2005/01/13 23:21:15 jlong777 Exp $
 */
 
void cb_version(char *string)
{
  char * version = "1.1c";
  
  while(*version) *string++ = *version++;
  *string = 0x0;
}
 
/*
cb_version(3B)                                           Last changed: 09-17-02

NAME
        cb_version - returns the version number of libcbl

SYNOPSIS

        C/C++:

        #include <cbl.h>
        void cb_version(char *string);

        Fortran:

        call cb_version(string)


IMPLEMENTATION

        Cray SV1 series UNICOS systems

DESCRIPTION

        cb_version returns the version number of the library used to link
        the executing program. The value is a character string of up to
        5 characters, such as "1.0".

        string  (output) For C, string should be a character 
                pointer pointing to a block of memory with at least
                6 bytes of space.  The returned string is null terminated.
                For Fortran, string should be declared CHARACTER(5).
                If the version text is shorter than 5 characters, string is
                padded with blanks at the end.
 
NOTES       

        (none)

EXAMPLES

        C/C++:

        main(){
          char *p ;
          p = (char *) malloc(6);
          cb_version(p);
          printf("%s\n", p);
        }


        Fortran:

        character(len=5) ::  string
        call cb_version(string)
        print *, string
        end
               
SEE ALSO

        INTRO_LIBCBL(3B)


*/ 
