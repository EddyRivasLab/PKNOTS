/* fillmtx.c 
 *
 * fill the no-hole matrices: VX, WX, WBX.         Dimension  (len x len)
 * and the hole matrices:     VHX, ZHX, YHX, WHX.  Dimension  (len x len x len x len)
 *
 */
                         
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>

#include "cfg.h"
#include "proto.h"
#include "squid.h"
                                                 
/* Function: FillMtx_nested()
 * 
 * Purpose:  The dynamic programming algorithm for folding
 *           an RNA using a stochastic context-free grammar.
 *           No pseudoknots
 *           
 * Arguments: s     - sequence (0..len-1) converted to integers
 *                    representing array indices of symbols
 *            len   - length of iseq
 *            icfg  - context-free grammar state transitions, integer log form
 *             wx   - DP matrix, already alloc'ed 
 *            wbx   - DP matrix, already alloc'ed 
 *             vx   - DP matrix, already alloc'ed 
 *            
 * Return:    (void)           
 *            vx, wx, wbx are filled (in that order).
 */       
void
FillMtx_nested(int *s, int len, int **icfg, int **wx, int **wbx, int **vx, int *vp, int allow_coaxials)
{
  int j;			/* row index (sequence end position)               */
  int d;		        /* column indices ("distance from diagonal")       */
  int min, mind;                /* minimum leght of hairpin loops allowed          */
  int min_kn;                   /* minimum leght of knots  allowed                 */

  for (j = 0; j < len; j++)
    {
                            
      /* only allows hairpin loops of length > lng 
       * only allows knots of length > lng_kn  
       */
      
      min    = (j - lng > 0) ? lng : j;
      min_kn = (j - lng_kn > 0) ? lng_kn : j;
      
      for (d = 0; d <= min; d++)
	{
	  wx[j][d]  =  0;
	  wbx[j][d] = -BIGINT;
	  vx[j][d]  = -BIGINT;
	}
      
      /* Initialize on d == 1 for lng == 0
       */
      if (j > 0 && lng == 0)
	{
	  wx[j][1]  =  0;
	  wbx[j][1] = -BIGINT;
	  vp[0]     = -BIGINT;
	  vp[1]     = -BIGINT;
	  vx[j][1]  = -BIGINT;    /*  hairpins of d=1 not allowed  */
	}
      
      /* Recursively calculate the rest of this row
       */
      mind = (lng > 0)? min+1 : min+2;

      for (d = mind; d <= j; d++)        
	{
	  /*  filling of VX(j)(d) without pseudoknots
	   *                      
	   */
	  FillVP(s, len, icfg, vx, vp, j, d);
	  FillVX_nested(s, len, icfg, wbx, vx, vp, j, d, allow_coaxials);

	  /* filling  of  WX(j)(d) == No hole between (i,j). Hole penalty (P11).
	   * filling  of WBX(j)(d) == No hole between (i,j) for multiloops loops. 
	   *                          It has dangling penalties (P6)
	   *                          and base pairs penalties (P10). Hole penalty (P13).
	   */
	  FillWX_nested(s, len, icfg, wx, vx, j, d);
	  FillWBX_nested(s, len, icfg, wbx, vx, j, d);


	}   /* while filling d */
    }   /* while filling j */
}
 

/* Function: FillMtx()
 * 
 * Purpose:  The dynamic programming algorithm for folding
 *           an RNA using a stochastic context-free grammar.
 *           Including pseudoknots.
 *           
 * Arguments: s     - sequence (0..len-1) converted to integers
 *                    representing array indices of symbols
 *            len   - length of iseq
 *            icfg  - context-free grammar state transitions, integer log form
 *            whx   - DP matrix, already alloc'ed 
 *            vhx   - DP matrix, already alloc'ed 
 *            zhx   - DP matrix, already alloc'ed 
 *            yhx   - DP matrix, already alloc'ed 
 *             wx   - DP matrix, already alloc'ed 
 *            wbx   - DP matrix, already alloc'ed 
 *             vx   - DP matrix, already alloc'ed 
 *            
 * Return:    (void)           
 *            vx, wx, wbx, vhx, zhx, yhx, whx are filled (in that order).
 */       
void
FillMtx(int *s, int len, int **icfg, int **wx, int **wbx, int **vx, int *vp, 
	int ****whx, int ****vhx, int ****zhx, int ****yhx, int allow_coaxials, int approx)
{
  int j;			/* row index (sequence end position)               */
  int d, d1, d2;		/* column indices ("distance from diagonal")       */
  int min, mind;                /* minimum leght of hairpin loops allowed          */
  int min_kn;                   /* minimum leght of knots  allowed                 */

  for (j = 0; j < len; j++)
    {
      /* only allows hairpin loops of length > lng 
       * only allows knots of length > lng_kn  
       */
      
      min    = (j - lng > 0) ? lng : j;
      min_kn = (j - lng_kn > 0) ? lng_kn : j;
      
      for (d = 0; d <= min; d++)
	{
	  wx[j][d]  = -BIGINT;
	  wbx[j][d] = -BIGINT;
	  vx[j][d]  = -BIGINT;
	  
	  for (d1 = 0; d1 <= d; d1++)
	    for (d2 = 0; d2 <= (d - d1); d2++)
	      {
		whx[j][d][d1][d2] = -BIGINT; 
		vhx[j][d][d1][d2] = -BIGINT; 
		zhx[j][d][d1][d2] = -BIGINT; 
		yhx[j][d][d1][d2] = -BIGINT; 
	      }
	}
      
      /* Initialize on d == 1 for lng == 0
       */
      if (j > 0 && lng == 0)
	{
	  vp[0]     = -BIGINT;
	  vp[1]     = -BIGINT;

	  vx[j][1]  = -BIGINT;    /*  hairpins of d=1 not allowed  */
	  wx[j][1]  = -BIGINT;
	  wbx[j][1] = -BIGINT;
	  
	  whx[j][1][0][0] = -BIGINT;
	  whx[j][1][0][1] = -BIGINT;
	  whx[j][1][1][0] = -BIGINT;
	  
	  vhx[j][1][0][0] = -BIGINT;
	  vhx[j][1][0][1] = -BIGINT;
	  vhx[j][1][1][0] = -BIGINT;
	  
	  zhx[j][1][0][0] = -BIGINT;
	  zhx[j][1][0][1] = -BIGINT;
	  zhx[j][1][1][0] = -BIGINT;
	  
	  yhx[j][1][0][0] = -BIGINT;
	  yhx[j][1][0][1] = -BIGINT;
	  yhx[j][1][1][0] = -BIGINT;
	}
      
      /* Recursively calculate the rest of this row
       */
      mind = (lng > 0)? min+1 : min+2;

      for (d = mind; d <= j; d++)        
	{
	  /* Initialize on d1 == 0 and d2 == 0 */
	  whx[j][d][0][0] = -BIGINT;
	  vhx[j][d][0][0] = -BIGINT;
	  zhx[j][d][0][0] = -BIGINT;
	  yhx[j][d][0][0] = -BIGINT;
	  
	  /*  filling of VX(j)(d) == no hole between (i,j) base-paired. 
	   *                       Multiloop (IS > 2) penalty (P5). Hole penalty (P12).
	   */
	  FillVP(s, len, icfg, vx, vp, j, d);
	  FillVX(s, len, icfg, wbx, vx, vp, whx, zhx, yhx, j, d, allow_coaxials, approx);
	  
	  /* filling  of  WX(j)(d) == No hole between (i,j). Hole penalty (P11).
	   * filling  of WBX(j)(d) == No hole between (i,j) for multiloops loops. 
	   *                          It has dangling penalties (P6)
	   *                          and base pairs penalties (P10). Hole penalty (P13).
	   */
	  FillWX (s, len, icfg, wx, vx, whx, yhx, j, d);
	  FillWBX(s, len, icfg, wbx, vx, whx, yhx, j, d, approx);
	  
	  
	  /*  filling of the ONE-HOLE matrices.
	   */
	  if (d <= min_kn)
	    {
	      for (d1 = 0; d1 <= d; d1++)
		for (d2 = 0; d2 <= (d-d1); d2++)
		  {
		    whx[j][d][d1][d2] = -BIGINT;
   		    vhx[j][d][d1][d2] = -BIGINT;
		    zhx[j][d][d1][d2] = -BIGINT;
		    yhx[j][d][d1][d2] = -BIGINT;
		  }
	      continue;
	    }
	  
	  for (d1 = 0; d1 <= d; d1++)
	    {

	      /* if no hole (k=l)][ (d2 = d - d1) use wbx, vx, or prohibit */
	      whx[j][d][d1][d-d1] = wbx[j][d];  
	      zhx[j][d][d1][d-d1] = vx[j][d];  
	      
	      vhx[j][d][d1][d-d1] = -BIGINT;  
	      yhx[j][d][d1][d-d1] = -BIGINT;  
	      
	      /* if  k < l */	     
	      for (d2 = 0; d2 < (d - d1); d2++) 
		{
		  if (d1 + d2 == 0) continue;  /* already assigned  */
		  
		  /* filling  of VHX(j)(d)(d1)(d2), (i,k)--(l,j). [i-j base-paired, k-l base-paired]
		   */	 
		  FillVHX(s, len, icfg, whx, vhx, j, d, d1, d2, min_kn);

		  /* filling  of ZHX(j)(d)(d1)(d2), (i,k)--(l,j). [i-j base-paired]
		   */
		  FillZHX(s, len, icfg, wbx, vx, whx, vhx, zhx, j, d, d1, d2);

		  /* filling  of YHX(j)(d)(d1)(d2), (i,k)--(l,j). [k-l base-paired]
		   */
		  FillYHX(s, len, icfg, wbx, vx, whx, vhx, yhx, j, d, d1, d2, min_kn);

		  /* filling  of WHX(j)(d)(d1)(d2), (i,k)--(l,j).
		   */
		  FillWHX(s, len, icfg, wbx, vx, whx, vhx, zhx, yhx, j, d, d1, d2);

		}   /* while filling d2 */
	    }   /* while filling d1 */
	}   /* while filling d */
    }   /* while filling j */
}
 

 

