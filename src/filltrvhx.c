/* filltrvhx.c 
 *
 * includes functions FillVHX,  TraceVHX
 * that fill and traceback the hole matix:
 *      vhx[j][d][d1][d2]  (len x len x len x len)
 *                         [(j - d, j)           are based paired]
 *                         [(j - d + d1, j - d2) are based paired]
 *
 */
                         
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>

#include "cfg.h"
#include "proto.h"
#include "protovhx.h"
#include "squid.h"

/* Function: FillVHX()
 * 
 * Purpose:  fill hole matrix VHX. Using Knots-IS(2) method.
 *           
 * Arguments: s     - sequence (0..len-1) converted to integers
 *                    representing array indices of symbols
 *            len   - length of iseq
 *           icfg   - context-free grammar state transitions, integer log form
 *            whx   - DP matrix, already alloc'ed 
 *            vhx   - DP matrix, already alloc'ed 
 *            zhx   - DP matrix, already alloc'ed 
 *            yhx   - DP matrix, already alloc'ed 
 *             wx   - DP matrix, already alloc'ed 
 *            wbx   - DP matrix, already alloc'ed 
 *             vx   - DP matrix, already alloc'ed 
 *            
 * Return:    (void)           
 *            vhx is filled.
 */       
void
FillVHX(int *s, int len, struct rnapar_2 *rnapar, 
	int **icfg, int ****whx, int ****vhx, 
	int j, int d, int d1, int d2, int min_kn)
{
  int      i,k,l;
  int mid1, mid2;  /* used for midpoint of a bifurc   */
  int         sc;  /* temporary score to check        */
  int     bestsc;  /* best score so far               */

  i = j - d;
  k = i + d1;
  l = j - d2;

  if (d1 > 0 && d2 > 0 && d - d1 - d2  > min_kn && 
      icfg[idxS][idxP(s[i],s[j])] == 1*INTSCALE
      && icfg[idxS][idxP(s[k],s[l])] == 1*INTSCALE)
    { 
      
      /* (VH1)__/ IRREDUCIBLE SURFACES of  O(2). 
       */
      bestsc = F2(s, len, rnapar, icfg, j, d, d1, d2);
    
      /* (VH2)__/ VHX + VHX */
      for (mid1 = 1; mid1 < d1; mid1++)
	for (mid2 = 1; mid2 < d2; mid2++)
	  if ((sc = vhx[j][d][d1-mid1][d2-mid2]
	       + vhx[l+mid2][d-d1-d2+mid1+mid2][mid1][mid2]) > bestsc)
	    bestsc = sc;
      
      /* (VH3)__/ REST of VHX (multiloops).  
       */
      
      /* (VH3.1) */
      if (d1 > 1 && d2 > 1 &&
	  (sc = 2*rnapar->P10P + rnapar->P5P 
	   + whx[j-1][d-2][d1-2][d2-2]) > bestsc)
	bestsc = sc; 
      
      /*  (VH3.2) */
      if (d1 > 2 && d2 > 1 &&
	  (sc = 2*rnapar->P10P + rnapar->P5P + rnapar->P6P 
	   + wkn*icfg[idxR(s[i+1])][idxP(s[j],s[i])]
	   + whx[j-1][d-3][d1-3][d2-2]) > bestsc)
	bestsc = sc; 
      
      /*  (VH3.3) */
      if (d1 > 1 && d2 > 2 &&
	  (sc = 2*rnapar->P10P + rnapar->P5P + rnapar->P6P 
	   + wkn*icfg[idxL(s[j-1])][idxP(s[j],s[i])]
	   + whx[j-2][d-3][d1-2][d2-3]) > bestsc)
	bestsc = sc; 
      
      /*  (VH3.4) */
      if (d1 > 2 && d2 > 1 &&
	  (sc = 2*rnapar->P10P + rnapar->P5P + rnapar->P6P 
	   + wkn*icfg[idxL(s[k-1])][idxP(s[l],s[k])]
	   + whx[j-1][d-2][d1-3][d2-2]) > bestsc)
	bestsc = sc; 
      
      /*  (VH3.5) */
      if (d1 > 1 && d2 > 2 &&
	  (sc = 2*rnapar->P10P + rnapar->P5P + rnapar->P6P 
	   + wkn*icfg[idxR(s[l+1])][idxP(s[l],s[k])] 
	   + whx[j-1][d-2][d1-2][d2-3]) > bestsc)
	bestsc = sc; 
      
      /*  (VH3.6) */
      if (d1 > 3 && d2 > 1 &&
	  (sc = 2*rnapar->P10P + rnapar->P5P + 2*rnapar->P6P 
	   + wkn*icfg[idxR(s[i+1])][idxP(s[j],s[i])] 
	   + wkn*icfg[idxL(s[k-1])][idxP(s[l],s[k])] 
	   + whx[j-1][d-3][d1-4][d2-2]) > bestsc)
	bestsc = sc; 
      
      /*  (VH3.7) */
      if (d1 > 2 && d2 > 2 &&
	  (sc = 2*rnapar->P10P + rnapar->P5P + 2*rnapar->P6P 
	   + wkn*icfg[idxR(s[i+1])][idxP(s[j],s[i])]
	   + wkn*icfg[idxR(s[l+1])][idxP(s[l],s[k])]
	   + whx[j-1][d-3][d1-3][d2-3]) > bestsc)
	bestsc = sc; 
      
      /*  (VH3.8) */
      if (d1 > 2 && d2 > 2 &&
	  (sc =2*rnapar->P10P + rnapar->P5P + 2*rnapar->P6P 
	   + wkn*icfg[idxR(s[i+1])][idxP(s[j],s[i])]
	   + wkn*icfg[idxL(s[j-1])][idxP(s[j],s[i])]
	   + whx[j-2][d-4][d1-3][d2-3] ) > bestsc)
	bestsc = sc; 
      
      /*  (VH3.9) */
      if (d1 > 2 && d2 > 2 &&
	  (sc = 2*rnapar->P10P + rnapar->P5P + 2*rnapar->P6P 
	   + wkn*icfg[idxL(s[k-1])][idxP(s[l],s[k])] 
	   + wkn*icfg[idxR(s[l+1])][idxP(s[l],s[k])]
	   + whx[j-1][d-2][d1-3][d2-3]) > bestsc)
	bestsc = sc; 
      
      /*  (VH3.10) */
      if (d1 > 2 && d2 > 2 &&
	  (sc = 2*rnapar->P10P + rnapar->P5P + 2*rnapar->P6P 
	   + wkn*icfg[idxL(s[k-1])][idxP(s[l],s[k])] 
	   + wkn*icfg[idxL(s[j-1])][idxP(s[j],s[i])]
	   + whx[j-2][d-3][d1-3][d2-3]) > bestsc)
	bestsc = sc; 
      
      /*  (VH3.11) */
      if (d1 > 1 && d2 > 3 &&
	  (sc = 2*rnapar->P10P + rnapar->P5P + 2*rnapar->P6P 
	   + wkn*icfg[idxR(s[l+1])][idxP(s[l],s[k])]
	   + wkn*icfg[idxL(s[j-1])][idxP(s[j],s[i])]
	   + whx[j-2][d-3][d1-2][d2-4]) > bestsc)
	bestsc = sc; 
      
      /*  (VH3.12) */
      if (d1 > 2 && d2 > 3 &&
	  (sc = 2*rnapar->P10P + rnapar->P5P + 3*rnapar->P6P 
	   + wkn*icfg[idxR(s[i+1])][idxP(s[j],s[i])] 
	   + wkn*icfg[idxR(s[l+1])][idxP(s[l],s[k])] 
	   + wkn*icfg[idxL(s[j-1])][idxP(s[j],s[i])]
	   + whx[j-2][d-4][d1-3][d2-4]) > bestsc)
	bestsc = sc; 
      
      /*  (VH3.13) */
      if (d1 > 3 && d2 > 2 &&
	  (sc = 2*rnapar->P10P + rnapar->P5P + 3*rnapar->P6P 
	   + wkn*icfg[idxR(s[i+1])][idxP(s[j],s[i])] 
	   + wkn*icfg[idxL(s[k-1])][idxP(s[l],s[k])] 
	   + wkn*icfg[idxR(s[l+1])][idxP(s[l],s[k])] 
	   + whx[j-1][d-3][d1-4][d2-3]) > bestsc)
	bestsc = sc; 
      
      /*  (VH3.14) */
      if (d1 > 3 && d2 > 2 &&
	  (sc = 2*rnapar->P10P + rnapar->P5P + 3*rnapar->P6P 
	   + wkn*icfg[idxR(s[i+1])][idxP(s[j],s[i])] 
	   + wkn*icfg[idxL(s[k-1])][idxP(s[l],s[k])]
	   + wkn*icfg[idxL(s[j-1])][idxP(s[j],s[i])]
	   + whx[j-2][d-4][d1-4][d2-3]) > bestsc)
	bestsc = sc; 
      
      /*  (VH3.15) */
      if (d1 > 2 && d2 > 3 &&
	  (sc =2*rnapar->P10P + rnapar->P5P + 3*rnapar->P6P 
	   + wkn*icfg[idxL(s[k-1])][idxP(s[l],s[k])]
	   + wkn*icfg[idxR(s[l+1])][idxP(s[l],s[k])] 
	   + wkn*icfg[idxL(s[j-1])][idxP(s[j],s[i])] 
	   + whx[j-2][d-3][d1-3][d2-4] ) > bestsc)
	bestsc = sc; 
      
      /*  (VH3.16) */
      if (d1 > 3 && d2 > 3 &&
	  (sc = 2*rnapar->P10P + rnapar->P5P + 4*rnapar->P6P 
	   + wkn*icfg[idxR(s[i+1])][idxP(s[j],s[i])] 
	   + wkn*icfg[idxL(s[k-1])][idxP(s[l],s[k])] 
	   + wkn*icfg[idxR(s[l+1])][idxP(s[l],s[k])] 
	   + wkn*icfg[idxL(s[j-1])][idxP(s[j],s[i])]
	   + whx[j-2][d-4][d1-4][d2-4]) > bestsc)
	bestsc = sc; 
      
      vhx[j][d][d1][d2] = bestsc;
      
    }  /* if (i,j) and (k,l) are base-paired */
  
  else vhx[j][d][d1][d2] = -BIGINT;  
}

/* Function: TraceVHX()
 * 
 * Purpose:  traceback hole matrix VHX. Using Knots-IS(2) method.
 *           
 * Arguments: s     - sequence (0..len-1) converted to integers
 *                    representing array indices of symbols
 *            len   - length of iseq
 *           icfg   - context-free grammar state transitions, integer log form
 *            whx   - DP matrix, already alloc'ed
 *            vhx   - DP matrix, already alloc'ed 
 *            zhx   - DP matrix, already alloc'ed 
 *            yhx   - DP matrix, already alloc'ed 
 *             wx   - DP matrix, already alloc'ed 
 *            wbx   - DP matrix, already alloc'ed 
 *             vx   - DP matrix, already alloc'ed 
 *            
 * Return:    (void)           
 */ 
void
TraceVHX(FILE *outf, int *s, int len, struct rnapar_2 *rnapar, 
	 int **icfg, int ****whx, int ****vhx, 
	 int j, int d, int d1, int d2, 
	 struct tracekn_s *curr_tr, struct traceknstack_s *dolist, int traceback)
{
  int mid1, mid2;  /* used for midpoint of a bifurc   */

  /* (VH1)__/
   */
  if (vhx[j][d][d1][d2] == VH1(s, len, rnapar, icfg, j, d, d1, d2)){
    trace_VH1(outf, vhx, j, d, d1, d2, curr_tr, dolist, traceback);
    return;
  }
  
  
  /* (VH2)__/ VHX + VHX */
  for (mid1 = d1 - 1; mid1 > 0; mid1--)
    for (mid2 = d2 - 1; mid2 > 0; mid2--)
      if (vhx[j][d][d1][d2] == VH2(s, len, rnapar, icfg, vhx, j, d, d1, d2, mid1, mid2)){
	trace_VH2(outf, vhx, j, d, d1, d2, mid1, mid2, curr_tr, dolist, traceback);
	return;
      }
  
  /* (VH3)__/  REST of VHX (multiloops) */
  
  /* (VH3.1) */
  if (vhx[j][d][d1][d2] == VH3_1(s, len, rnapar, icfg, whx, j, d, d1, d2)){
    trace_VH3_1(outf, whx, j, d, d1, d2, curr_tr, dolist, traceback);
    return;
  }
  
  /* (VH3.2) */
  else  if (vhx[j][d][d1][d2] == VH3_2(s, len, rnapar, icfg, whx, j, d, d1, d2)){
    trace_VH3_2(outf, whx, j, d, d1, d2, curr_tr, dolist, traceback);
    return;
  }
  
  /* (VH3.3) */
  else  if (vhx[j][d][d1][d2] == VH3_3(s, len, rnapar, icfg, whx, j, d, d1, d2)){
    trace_VH3_3(outf, whx, j, d, d1, d2, curr_tr, dolist, traceback);
    return;
      }
  
  /* (VH3.4) */
  else  if (vhx[j][d][d1][d2] == VH3_4(s, len, rnapar, icfg, whx, j, d, d1, d2)){
    trace_VH3_4(outf, whx, j, d, d1, d2, curr_tr, dolist, traceback);
    return;
  }
  
  /* (VH3.5) */
  else  if (vhx[j][d][d1][d2] == VH3_5(s, len, rnapar, icfg, whx, j, d, d1, d2)){
    trace_VH3_5(outf, whx, j, d, d1, d2, curr_tr, dolist, traceback);
    return;
  }
 
  /* (VH3.6) */
  else  if (vhx[j][d][d1][d2] == VH3_6(s, len, rnapar, icfg, whx, j, d, d1, d2)){
    trace_VH3_6(outf, whx, j, d, d1, d2, curr_tr, dolist, traceback);
    return;
  }
  
  /* (VH3.7) */
  else  if (vhx[j][d][d1][d2] == VH3_7(s, len, rnapar, icfg, whx, j, d, d1, d2)){
    trace_VH3_7(outf, whx, j, d, d1, d2, curr_tr, dolist, traceback);
    return;
  }
  
  /* (VH3.8) */
  else  if (vhx[j][d][d1][d2] == VH3_8(s, len, rnapar, icfg, whx, j, d, d1, d2)){
    trace_VH3_8(outf, whx, j, d, d1, d2, curr_tr, dolist, traceback);
    return;
  }
  
  /* (VH3.9) */
  else  if (vhx[j][d][d1][d2] == VH3_9(s, len, rnapar, icfg, whx, j, d, d1, d2)){
    trace_VH3_9(outf, whx, j, d, d1, d2, curr_tr, dolist, traceback);
    return;
  }
  
  /* (VH3.10) */
  else  if (vhx[j][d][d1][d2] == VH3_10(s, len, rnapar, icfg, whx, j, d, d1, d2)){
    trace_VH3_10(outf, whx, j, d, d1, d2, curr_tr, dolist, traceback);
    return;
  }
  
  /* (VH3.11) */
  else  if (vhx[j][d][d1][d2] == VH3_11(s, len, rnapar, icfg, whx, j, d, d1, d2)){
    trace_VH3_11(outf, whx, j, d, d1, d2, curr_tr, dolist, traceback);
    return;
  }
  
  /* (VH3.12) */
 else  if (vhx[j][d][d1][d2] == VH3_12(s, len, rnapar, icfg, whx, j, d, d1, d2)){
   trace_VH3_12(outf, whx, j, d, d1, d2, curr_tr, dolist, traceback);
   return;
 }
  
  /* (VH3.13) */
  else  if (vhx[j][d][d1][d2] == VH3_13(s, len, rnapar, icfg, whx, j, d, d1, d2)){
    trace_VH3_13(outf, whx, j, d, d1, d2, curr_tr, dolist, traceback);
    return;
  }
  
  /* (VH3.14) */
  else  if (vhx[j][d][d1][d2] == VH3_14(s, len, rnapar, icfg, whx, j, d, d1, d2)){
    trace_VH3_14(outf, whx, j, d, d1, d2, curr_tr, dolist, traceback);
    return;
  }
  
  /* (VH3.15) */
  else  if (vhx[j][d][d1][d2] == VH3_15(s, len, rnapar, icfg, whx, j, d, d1, d2)){
    trace_VH3_15(outf, whx, j, d, d1, d2, curr_tr, dolist, traceback);
    return;
  }
  
  /* (VH3.16) */
  else  if (vhx[j][d][d1][d2] == VH3_16(s, len, rnapar, icfg, whx, j, d, d1, d2)){
    trace_VH3_16(outf, whx, j, d, d1, d2, curr_tr, dolist, traceback);
    return;
  }
}

