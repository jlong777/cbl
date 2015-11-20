/* cb_swa_fw_init.c                             http://cbl.sourceforge.net
 *
 * called by cb_swa_fw.
 *
 * Copyright (C) 2004 University of Alaska Fairbanks
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
 * $Id: cb_swa_fw_init.c,v 1.2 2004/03/10 00:22:26 jlong777 Exp $
 */
 
void cb_swa_fw_init(long *dbl, long dbllen, long *dbs, long dbslen, 
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

        
