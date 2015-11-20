/* $Id: README,v 1.5 2004/04/27 20:02:19 jlong777 Exp $ */

Portable Cray Bioinformatics Library
====================================

Copyright (C) 2003 University of Alaska Fairbanks
Arctic Region Supercomputing Center (ARSC)
http://www.arsc.edu

  This library is free software; you can redistribute it and/or modify
  it under the terms of the GNU Lesser General Public License as
  published by the Free Software Foundation; either version 2.1 of the
  License, or (at your option) any later version.

  This library is distributed in the hope that it will be useful, but
  WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
  Lesser General Public License for more details.

  You should have received a copy of the GNU Lesser General Public
  License along with this library; if not, write to the
  Free Software Foundation, Inc.
  59 Temple Place, Suite 330
  Boston, MA  02111-1307 USA


FILES in version 1.1
====================

LICESE
README
cb_amino_translate_ascii.c
cb_compress.c
cb_copy_bits.c
cb_countn_ascii.c
cb_fasta_convert.c
cb_free.c
cb_irand.c
cb_malloc.c
cb_read_fasta.c
cb_repeatn.c
cb_revcompl.c
cb_searchn.c
cb_sw_fw_score.c
cb_swa_fw.c
cb_swa_fw_align.c
cb_swa_fw_init.c
cb_uncompress.c
cb_version.c
cbl.in
configure
html - html man pages directory
makefile.in
pam2cbl.c - see bullet 7) below
test_copy.c
test_repeat.c
test_search.c
test_suite.c


BUILDING, TESTING, and INSTALLING
=================================

1) edit the "configure" file to set the following:  

TOP_DIR:    top directory for lib and include, default is "/usr/local"
VBOX:       if you have a vector machine like a CRAY SV1, SX6 etc.
            altivec and mmx don't count here (at least not yet...)
            allowable values are "VECTOR_BOX" or "NOT_VECTOR_BOX" (default)
VLEN:       the vector length of a VECTOR_BOX, or 128 for NOT_VECTOR_BOX
            SV1="64", SX6="256", default "128" is for non-vector box,
	    results vary w/64 or 256, but must be >= number of bits in a word.
BMARK:      set to "1" to run benchhmark during "make test", default "0"
DEBUG_INFO: set to "1" to print out diagnostics if any "make test" fails
CC:         your C compiler, set for 64 bit mode if 64-bit machine
CCFLAGS:    set for the highest optimization possible
AR:         archiver, whatever makes static libraries on your system
            "ar" is the default
ARFLAGS:    default is "-rs", IBM at least needs "-rsX64" for 64 bit

optionally, set DEBUG and PROF if you need these.

2) "./configure" will create the makefile you need for the build

3) "make" will build libcbl.a

4) optionally, "make test" will run the test programs (takes awhile) and, if
   BMARK is set, run a little benchmark for each, removing the executables
   when finished
   
5) "make install" puts libcbl in TOP_DIR/lib and cbl.h in TOP_DIR/include

6) "make clean" removes *o and libcbl.a

7) pam2cbl.c is code to generate an sslookup definition from a pam scoring 
   matrix. Cut and paste the output into the sslookup definition of your 
   driver for the cb_sw routine, see the source code for details.
