/* filltrzhxyhx.c 
 *
 * includes functions FillZHX, FillYHX, TraceZHX, TraceYHX,
 * that fill and traceback the hole matices:
 *      zhx[j][d][d1][d2]  (len x len x len x len)
 *                         [(j - d, j)           are based paired]
 *                         [(j - d + d1, j - d2) unknown relation]
 *      yhx[j][d][d1][d2]  (len x len x len x len)
 *                         [(j - d, j)           unknown relation ]
 *                         [(j - d + d1, j - d2) are based paired]
 *
 */
                         
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>

#include <easel.h>

#include "pknots.h"
#include "pk_filltrzhxyhx.h"
#include "pk_irredsurf.h"
#include "pk_zhxyhxgraphs.h"
#include "pk_util.h"



/* Function: FillZHX()
 * 
 * Purpose:  fill hole matrix ZHX. Using Knots-IS(2) method.
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
 *            zhx is filled.
 */ 
void
FillZHX(ESL_DSQ *s, int len, struct rnapar_2 *rnapar, int **icfg, int **wbx, int **vx, 
	 int ****whx, int ****vhx, int ****zhx, 
	 int j, int d, int d1, int d2)
{
  int            i, k, l;
  int                mid;  /* midpoint of bifurcation         */
  int         mid1, mid2;  /* used for midpoint of a bifurc   */
  int                 sc;  /* temporary score to check        */
  int             bestsc;  /* best score so far               */

  i = j -d;
  k = i + d1;
  l = j - d2;

  if (d1 > 0 && d2 > 0 && d1 + d2 < d - 1 &&
      icfg[idxS][idxP(s[i],s[j])] == 1*INTSCALE)
    { 
      
      /* (ZH1)__/ 
       */
      bestsc = rnapar->P10P + vhx[j][d][d1][d2];
      
      /* (ZH2)__/ 
       */
      if (d1 > 0 && d2 > 0 &&
	  (sc = rnapar->P10P + 2*rnapar->P6P 
	   + wkn*icfg[idxR(s[k])][idxP(s[k-1],s[l+1])]
	   + wkn*icfg[idxL(s[l])][idxP(s[k-1],s[l+1])]
	   + vhx[j][d][d1-1][d2-1]) > bestsc)
	bestsc = sc;
      
      /* (ZH3)__/ 
       */
      if (d1 > 0 && 
	  (sc = rnapar->P10P + rnapar->P6P 
	   + wkn*icfg[idxR(s[k])][idxP(s[k-1],s[l])]
	   + vhx[j][d][d1-1][d2]) > bestsc)
	bestsc = sc; 
      
      /* (ZH4)
       */
      if (d2 > 0 && 
	  (sc = rnapar->P10P + rnapar->P6P 
	   + wkn*icfg[idxL(s[l])][idxP(s[k],s[l+1])]
	   + vhx[j][d][d1][d2-1]) > bestsc)
	bestsc = sc; 
      
      /* (ZH5)
       */
      if (d1 > 0 && 
	  (sc = rnapar->P6P + zhx[j][d][d1-1][d2]) > bestsc)
	bestsc = sc;
      
      /* (ZH6)
       */
      if (d2 > 0 && 
	  (sc = rnapar->P6P + zhx[j][d][d1][d2-1]) > bestsc)
	bestsc = sc; 
      
     for (mid = 0; mid < d1; mid++)
	{ 
	  /* (ZH7)__/ 1 WBX 1 ZHX
	   */                                       
	  if ((sc = wbx[k][mid] + zhx[j][d][d1-mid-1][d2]) > bestsc)
	    bestsc = sc;
	  
	  /* (ZH8)__/ 1 VX 1 VHX
	   */
	  /* (ZH8.1) */
	  if ((sc = 2*rnapar->P10P 
	       + wkn*(icfg[idxPS(s[k-mid-1],s[l])]
		      [idxPS(s[k-mid],s[k])] + INTSCALE)
	       + vx[k][mid] + vhx[j][d][d1-mid-1][d2]) > bestsc)
	    bestsc = sc;
	  
	  /* (ZH8.2) */
	  if (mid > 0 &&
	      (sc = 2*rnapar->P10P + rnapar->P6P 
	       + wkn*icfg[idxR(s[k])][idxP(s[k-mid],s[k-1])]
	       + wkn*icfg[idxPS(s[k-mid-1],s[l])][idxPS(s[k-mid],s[k-1])]
	       + vx[k-1][mid-1] + vhx[j][d][d1-mid-1][d2]) > bestsc)
	    bestsc = sc;
	  
	  /* (ZH8.3) */
	  if (d2 > 0 && 
	      (sc = 2*rnapar->P10P + rnapar->P6P 
	       + wkn*icfg[idxL(s[l])][idxP(s[l+1],s[k-mid-1])]
	       + wkn*icfg[idxPS(s[k-mid-1],s[l+1])][idxPS(s[k-mid],s[k])]
	       + vx[k][mid] + vhx[j][d][d1-mid-1][d2-1]) > bestsc)
	    bestsc = sc;
	  
	  /* (ZH8.4) */
	  if (mid > 0 && d2 > 0 && 
	      (sc = 2*rnapar->P10P + 2*rnapar->P6P 
	       + wkn*icfg[idxR(s[k])][idxP(s[k-mid],s[k-1])]
	       + wkn*icfg[idxL(s[l])][idxP(s[l+1],s[k-mid-1])]
	       + wkn*icfg[idxPS(s[k-mid-1],s[l+1])][idxPS(s[k-mid],s[k-1])]
	       + vx[k-1][mid-1] + vhx[j][d][d1-mid-1][d2-1]) > bestsc)
	    bestsc = sc;
	}
      
      for (mid = 0; mid < d2; mid++)
	{ 
	  /* (ZH9)__/ 1 WBX 1 ZHX
	   */  
	  if ((sc = wbx[l+mid][mid] + zhx[j][d][d1][d2-mid-1]) > bestsc)
	    bestsc = sc;
	  
	  /* (ZH10)__/ 1 VX 1 VHX 
	   */ 
	  /* (ZH10.1) */
	  if ((sc = 2*rnapar->P10P 
	       + wkn*(icfg[idxPS(s[l+mid],s[l])]
		      [idxPS(s[l+mid+1],s[k])] + INTSCALE)
	       + vx[l+mid][mid] + vhx[j][d][d1][d2-mid-1]) > bestsc)
	    bestsc = sc;
	  
	  /* (ZH10.2) */
	  if (mid > 0 &&  
	      (sc =2*rnapar->P10P + rnapar->P6P 
	       + wkn*icfg[idxL(s[l])][idxP(s[l+1],s[l+mid])]
	       + wkn*icfg[idxPS(s[l+mid],s[l+1])][idxPS(s[l+mid+1],s[k])]
	       + vx[l+mid][mid-1] + vhx[j][d][d1][d2-mid-1] ) > bestsc)
	    bestsc = sc;
	  
	  /* (ZH10.3) */
	  if (d1 > 0 && 
	      (sc = 2*rnapar->P10P + rnapar->P6P 
	       + wkn*icfg[idxR(s[k])][idxP(s[l+mid+1],s[k-1])]
	       + wkn*icfg[idxPS(s[l+mid],s[l])][idxPS(s[l+mid+1],s[k-1])]
	       + vx[l+mid][mid] + vhx[j][d][d1-1][d2-mid-1]) > bestsc)
	    bestsc = sc;
	  
	  /* (ZH10.4) */
	  if (mid > 0 && d1 > 0 &&  
	      (sc = 2*rnapar->P10P + 2*rnapar->P6P 
	       + wkn*icfg[idxR(s[k])][idxP(s[l+mid+1],s[k-1])]
	       + wkn*icfg[idxL(s[l])][idxP(s[l+1],s[l+mid])]
	       + wkn*icfg[idxPS(s[l+mid],s[l+1])][idxPS(s[l+mid+1],s[k-1])]
	       + vx[l+mid][mid-1] + vhx[j][d][d1-1][d2-mid-1]) > bestsc)
	    bestsc = sc;
	  
	}
      
      /* (ZH11)__/ 1 ZHX 1 WHX
       */
      for (mid1 = 1; mid1 <= d1; mid1++)                     
	for (mid2 = 1; mid2 <= d2; mid2++)
	  if ((sc = vhx[j][d][mid1][mid2]
	       + zhx[j-mid2][d-mid1-mid2][d1-mid1][d2-mid2]) > bestsc)
	    bestsc = sc; 
      
      /* (ZH12)__/ REST of ZHX (multiloops).  
       */
      
      /* (ZH12.1) */
      if (d1 > 0 && d2 > 0 &&
	  (sc = rnapar->P10P + rnapar->P5P 
	   + whx[j-1][d-2][d1-1][d2-1]) > bestsc)
	bestsc = sc; 
      
      /*  (ZH12.2) */
      if (d1 > 1 && d2 > 0 &&
	  (sc =  rnapar->P10P + rnapar->P5P + rnapar->P6P 
	   + wkn*icfg[idxR(s[i+1])][idxP(s[j],s[i])]
	   + whx[j-1][d-3][d1-2][d2-1]) > bestsc)
	bestsc = sc; 
      
      /*  (ZH12.3) */
      if (d1 > 0 && d2 > 1 &&
	  (sc = rnapar->P10P + rnapar->P5P + rnapar->P6P 
	   + wkn*icfg[idxL(s[j-1])][idxP(s[j],s[i])]
	   + whx[j-2][d-3][d1-1][d2-2]) > bestsc)
	bestsc = sc; 
      
      /*  (ZH12.4) */
      if (d1 > 1 && d2 > 1 &&
	  (sc = rnapar->P10P + rnapar->P5P + 2*rnapar->P6P 
	   + wkn*icfg[idxR(s[i+1])][idxP(s[j],s[i])]
	   + wkn*icfg[idxL(s[j-1])][idxP(s[j],s[i])]
	   + whx[j-2][d-4][d1-2][d2-2]) > bestsc)
	bestsc = sc; 
      
      zhx[j][d][d1][d2] = bestsc;
      
    } /* if (i,j) are base paired and hole*/
  
  else if (d1 + d2 >= d - 1 &&
	   icfg[idxS][idxP(s[i],s[j])] == 1*INTSCALE)
    zhx[j][d][d1][d2] = vx[j][d];  /* if (i,j) are base paired and no hole*/
  
  else zhx[j][d][d1][d2] = -BIGINT;
}


/* Function: TraceZHX()
 * 
 * Purpose:  traceback hole matrix ZHX. Using Knots-IS(2) method.
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
TraceZHX(FILE *outf, ESL_DSQ *s, int len, struct rnapar_2 *rnapar, int **icfg, int **wbx, int **vx, 
	 int ****whx, int ****vhx, int ****zhx, 
	 int j, int d, int d1, int d2, 
	 struct tracekn_s *curr_tr, struct traceknstack_s *dolist, int traceback)
{
  int mid, mid1, mid2;    /* used for midpoint of a bifurc   */

  /* (ZH1)__/
   */
  if (zhx[j][d][d1][d2] == ZH1(s, len, rnapar, icfg, vhx, j, d, d1, d2)){
    trace_ZH1(outf, vhx, j, d, d1, d2, curr_tr, dolist, traceback);
    return;
  }   
  
  /* (ZH2)__/
   */
  else if (zhx[j][d][d1][d2] == ZH2(s, len, rnapar, icfg, vhx, j, d, d1, d2)){
    trace_ZH2(outf, vhx, j, d, d1, d2, curr_tr, dolist, traceback);
    return;
  }    
  
  /* (ZH3)__/
   */
  else if (zhx[j][d][d1][d2] == ZH3(s, len, rnapar, icfg, vhx, j, d, d1, d2)){
    trace_ZH3(outf, vhx, j, d, d1, d2, curr_tr, dolist, traceback);
    return;
  } 
  
  /* (ZH4)__/ 
   */
  else if (zhx[j][d][d1][d2] ==  ZH4(s, len, rnapar, icfg, vhx, j, d, d1, d2)){
    trace_ZH4(outf, vhx, j, d, d1, d2, curr_tr, dolist, traceback);
    return;
  }
  
  /* (ZH5)__/ 
   */
  else if (zhx[j][d][d1][d2] == ZH5(s, len, rnapar, icfg, zhx, j, d, d1, d2)){
    trace_ZH5(outf, zhx, j, d, d1, d2, curr_tr, dolist, traceback);
    return;
  }
  
  /* (ZH6)__/
   */
  else if (zhx[j][d][d1][d2] == ZH6(s, len, rnapar, icfg, zhx, j, d, d1, d2)){
    trace_ZH6(outf, zhx, j, d, d1, d2, curr_tr, dolist, traceback);
    return;
  }     
  
  for (mid = d1-1; mid > 0; mid--) {
    
    /* (ZH7)__/ 1 WBX 1 ZHX
     */                                       
    if (zhx[j][d][d1][d2] == ZH7(s, len, rnapar, icfg, wbx, zhx, j, d, d1, d2, mid)){
      trace_ZH7(outf, wbx, zhx, j, d, d1, d2, mid, curr_tr, dolist, traceback);
      return;
    }

    /* (ZH8)__/ 1 VX 1 ZHX 
     */
    /* (ZH8.1) */
    else if (zhx[j][d][d1][d2] == ZH8_1(s, len, rnapar, icfg, vx, vhx, j, d, d1, d2, mid)){
      trace_ZH8_1(outf, vx, vhx, j, d, d1, d2, mid, curr_tr, dolist, traceback);
      return;
    }

    /* (ZH8.2) */
    else if (zhx[j][d][d1][d2] == ZH8_2(s, len, rnapar, icfg, vx, vhx, j, d, d1, d2, mid)){
      trace_ZH8_2(outf, vx, vhx, j, d, d1, d2, mid, curr_tr, dolist, traceback);
      return;
    }

    /* (ZH8.3) */
    else if (zhx[j][d][d1][d2] == ZH8_3(s, len, rnapar, icfg, vx, vhx, j, d, d1, d2, mid)){
      trace_ZH8_3(outf, vx, vhx, j, d, d1, d2, mid, curr_tr, dolist, traceback);
      return;
    }
      
    /* (ZH8.4) */
    else if (zhx[j][d][d1][d2] == ZH8_4(s, len, rnapar, icfg, vx, vhx, j, d, d1, d2, mid)){
      trace_ZH8_4(outf, vx, vhx, j, d, d1, d2, mid, curr_tr, dolist, traceback);
      return;
    }
  }
  
  for (mid = d2 - 1; mid > 0; mid--) {

    /* (ZH9)__/ 1 WBX 1 ZHX
     */  
    if  (zhx[j][d][d1][d2] == ZH9(s, len, rnapar, icfg, wbx, zhx, j, d, d1, d2, mid)){
      trace_ZH9(outf, wbx, zhx, j, d, d1, d2, mid, curr_tr, dolist, traceback);
      return;
    }	
    /* (ZH10)__/ 1 VX 1 ZHX 
     */ 
    /* (ZH10.1) */
    else if (zhx[j][d][d1][d2] == ZH10_1(s, len, rnapar, icfg, vx, vhx, j, d, d1, d2, mid)){
      trace_ZH10_1(outf, vx, vhx, j, d, d1, d2, mid, curr_tr, dolist, traceback);
      return;
    }
	    
    /* (ZH10.2) */
    else if (zhx[j][d][d1][d2] == ZH10_2(s, len, rnapar, icfg, vx, vhx, j, d, d1, d2, mid)){
      trace_ZH10_2(outf, vx, vhx, j, d, d1, d2, mid, curr_tr, dolist, traceback);
      return;
    }
	    
    /* (ZH10.3) */
    else if (zhx[j][d][d1][d2] == ZH10_3(s, len, rnapar, icfg, vx, vhx, j, d, d1, d2, mid)){
      trace_ZH10_3(outf, vx, vhx, j, d, d1, d2, mid, curr_tr, dolist, traceback);
      return;
    }
	    
    /* (ZH10.4) */
    else if (zhx[j][d][d1][d2] == ZH10_4(s, len, rnapar, icfg, vx, vhx, j, d, d1, d2, mid)){
      trace_ZH10_4(outf, vx, vhx, j, d, d1, d2, mid, curr_tr, dolist, traceback);
      return;
    }
  }
  
  /* (ZH11)__/ 1 VHX 1 ZHX
   */
  for (mid1 = d1; mid1 > 0; mid1--)
    for (mid2 = d2; mid2 > 0; mid2--)
      if (zhx[j][d][d1][d2] == ZH11(s, len, rnapar, icfg, zhx, vhx, j, d, d1, d2, mid1, mid2)){
      trace_ZH11(outf, zhx, vhx, j, d, d1, d2, mid1, mid2, curr_tr, dolist, traceback);
      return;
    }

  /* (ZH12)__/  REST of ZHX (multiloops) */
  
  /* (ZH12.1) */
  if (zhx[j][d][d1][d2] == ZH12_1(s, len, rnapar, icfg, whx, j, d, d1, d2)){
    trace_ZH12_1(outf, whx, j, d, d1, d2, curr_tr, dolist, traceback);
    return;
  }
  
  /* (ZH12.2) */
  else  if (zhx[j][d][d1][d2] == ZH12_2(s, len, rnapar, icfg, whx, j, d, d1, d2)){
    trace_ZH12_2(outf, whx, j, d, d1, d2, curr_tr, dolist, traceback);
    return;
  }
  
  /* (ZH12.3) */
  else  if (zhx[j][d][d1][d2] == ZH12_3(s, len, rnapar, icfg, whx, j, d, d1, d2)){
    trace_ZH12_3(outf, whx, j, d, d1, d2, curr_tr, dolist, traceback);
    return;
      }
  
  /* (ZH12.4) */
  else  if (zhx[j][d][d1][d2] == ZH12_4(s, len, rnapar, icfg, whx, j, d, d1, d2)){
    trace_ZH12_4(outf, whx, j, d, d1, d2, curr_tr, dolist, traceback);
    return;
  }
}

/* Function: FillYHX()
 * 
 * Purpose:  fill hole matrix YHX. Using Knots-IS(2) method.
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
 *            yhx is filled.
 */ 
void
FillYHX(ESL_DSQ *s, int len, struct rnapar_2 *rnapar, int **icfg, int **wbx, int **vx, 
	 int ****whx, int ****vhx, int ****yhx,
	 int j, int d, int d1, int d2, int min_kn)
{
  int              i,k,l;
  int                mid;  /* midpoint of bifurcation         */
  int         mid1, mid2;  /* used for midpoint of a bifurc   */
  int                 sc;  /* temporary score to check        */
  int             bestsc;  /* best score so far               */
  
  i = j - d;
  k = i + d1;
  l = j - d2;

 if (d1 > 0 && d2 > 0 && d - d1 - d2 > min_kn && 
      icfg[idxS][idxP(s[k],s[l])] == 1*INTSCALE)
    { 
      
      /* (YH1)__/ 
       */
      bestsc = rnapar->P10P + vhx[j][d][d1][d2];
      
      /* (YH2)__/  
       */
      if (d1 > 0 && d2 > 0 &&
	  (sc = rnapar->P10P + 2*rnapar->P6P 
	   + wkn*icfg[idxL(s[i])][idxP(s[i+1],s[j-1])]
	   + wkn*icfg[idxR(s[j])][idxP(s[i+1],s[j-1])]
	   + vhx[j-1][d-2][d1-1][d2-1]) > bestsc)
	bestsc = sc;
      
      /* (YH3)__/ 
       */
      if (d1 > 0 && 
	  (sc = rnapar->P10P + rnapar->P6P 
	   + wkn*icfg[idxL(s[i])][idxP(s[i+1],s[j])] 
	   + vhx[j][d-1][d1-1][d2]) > bestsc)
	bestsc = sc; 
      
      /* (YH4)__/ */
      if (d2 > 0 && 
	  (sc = rnapar->P10P + rnapar->P6P 
	   + wkn*icfg[idxR(s[j])][idxP(s[i],s[j-1])]
	   + vhx[j-1][d-1][d1][d2-1]) > bestsc)
	bestsc = sc; 
      
      /* (YH5)__/ 
       */
      if (d1 > 0 && 
	  (sc = rnapar->P6P + yhx[j][d-1][d1-1][d2]) > bestsc)
	bestsc = sc;
      
      /* (YH6)__/  */
      if (d2 > 0 && 
	  (sc = rnapar->P6P + yhx[j-1][d-1][d1][d2-1]) > bestsc)
	bestsc = sc; 

      for (mid = 0; mid < d1; mid++)
	{ 
	  /* (YH7)__/ 1 WBX 1 YHX
	   */                                       
	  if ((sc = wbx[i+mid][mid] + yhx[j][d-mid-1][d1-mid-1][d2]) > bestsc)
	    bestsc = sc;
	  
	  /* (YH8)__/ 1 VX 1 VHX 
	   */                                       
	  /* (YH8.1) */
	  if ((sc =2*rnapar->P10P 
	       + wkn*(icfg[idxPS(s[i+mid],s[i])][idxPS(s[i+mid+1],s[j])] + INTSCALE)
	       + vx[i+mid][mid] + vhx[j][d-mid-1][d1-mid-1][d2]) > bestsc)
	    bestsc = sc;
	  
	  /* (YH8.2) */
	  if (mid > 0 &&
	      (sc =2*rnapar->P10P + rnapar->P6P
	       + wkn*icfg[idxL(s[i])][idxP(s[i+1],s[i+mid])]
	       + wkn*icfg[idxPS(s[i+mid],s[i+1])][idxPS(s[i+mid+1],s[j])]
	       + vx[i+mid][mid-1] + vhx[j][d-mid-1][d1-mid-1][d2] ) > bestsc)
	    bestsc = sc;
	  
	  /* (YH8.3) */
	  if (d2 > 0 &&
	      (sc = 2*rnapar->P10P + rnapar->P6P
	       + wkn*icfg[idxR(s[j])][idxP(s[i+mid+1],s[j-1])]
	       + wkn*icfg[idxPS(s[i+mid],s[i])][idxPS(s[i+mid+1],s[j-1])]
	       + vx[i+mid][mid] + vhx[j-1][d-mid-2][d1-mid-1][d2-1]) > bestsc)
	    bestsc = sc;
	  
	  /* (YH8.4) */
	  if (mid > 0 && d2 > 0 &&
	      (sc = 2*rnapar->P10P + 2*rnapar->P6P
	       + wkn*icfg[idxL(s[i])][idxP(s[i+1],s[i+mid])]
	       + wkn*icfg[idxR(s[j])][idxP(s[i+mid+1],s[j-1])]
	       + wkn*icfg[idxPS(s[i+mid],s[i+1])][idxPS(s[i+mid+1],s[j-1])]
	       + vx[i+mid][mid-1] + vhx[j-1][d-mid-2][d1-mid-1][d2-1]) > bestsc)
	    bestsc = sc;
	}
      
      for (mid = 0; mid < d2; mid++)
	{ 
	  /* (YH9)__/ 1 WBX 1 YHX
	   */                                       
	  if ((sc = wbx[j][mid] + yhx[j-mid-1][d-mid-1][d1][d2-mid-1]) > bestsc)
	    bestsc = sc;
	  
	  /* (YH10)__/ 1 VX 1 VHX 
	   */                                       
	  /* (YH10.1) */
	  if ((sc = 2*rnapar->P10P 
	       + wkn*(icfg[idxPS(s[j-mid-1],s[i])]
		      [idxPS(s[j-mid],s[j])] + INTSCALE)
	       + vx[j][mid] + vhx[j-mid-1][d-mid-1][d1][d2-mid-1] ) > bestsc)
	    bestsc = sc;
	  
	  /* (YH10.2) */
	  if (d1 > 0 &&
	      (sc = 2*rnapar->P10P + rnapar->P6P
	       + wkn*icfg[idxL(s[i])][idxP(s[i+1],s[j-mid-1])]
	       + wkn*icfg[idxPS(s[j-mid-1],s[i+1])][idxPS(s[j-mid],s[j])]
	       + vx[j][mid] + vhx[j-mid-1][d-mid-2][d1-1][d2-mid-1]) > bestsc)
	    bestsc = sc;
	  
	  /* (YH10.3) */
	  if (mid > 0 &&
	      (sc = 2*rnapar->P10P + rnapar->P6P
	       + wkn*icfg[idxR(s[j])][idxP(s[j-mid],s[j-1])]
	       + wkn*icfg[idxPS(s[j-mid-1],s[i])][idxPS(s[j-mid],s[j-1])]
	       + vx[j-1][mid-1] + vhx[j-mid-1][d-mid-1][d1][d2-mid-1]) > bestsc)
	    bestsc = sc;
	  
	  /* (YH10.4) */
	  if (d1 > 0 && mid > 0 &&
	      (sc = 2*rnapar->P10P + 2*rnapar->P6P
	       + wkn*icfg[idxL(s[i])][idxP(s[i+1],s[j-mid-1])]
	       + wkn*icfg[idxR(s[j])][idxP(s[j-mid],s[j-1])]
	       + wkn*icfg[idxPS(s[j-mid-1],s[i+1])][idxPS(s[j-mid],s[j-1])]
	       + vx[j-1][mid-1] + vhx[j-mid-1][d-mid-2][d1-1][d2-mid-1]) > bestsc)
	    bestsc = sc;
	  
	}
      
      /* (YH11)__/ 1 VHX 1 YHX
       */
      for (mid1 = 1; mid1 <= d1; mid1++)                     
	for (mid2 = 1; mid2 <= d2; mid2++)
	  if ((sc = vhx[l+mid2][d-d1-d2+mid1+mid2][mid1][mid2]
	       + yhx[j][d][d1-mid1][d2-mid2]) > bestsc)
	    bestsc = sc;

      /* (YH12)__/ REST of YHX (multiloops).  
       */
      
      /* (YH12.1) */
      if (d1 > 0 && d2 > 0 &&
	  (sc =  rnapar->P10P + rnapar->P5P 
	   + whx[j][d][d1-1][d2-1]) > bestsc)
	bestsc = sc; 
      
      /*  (YH12.2) */
      if (d1 > 1 && d2 > 0 &&
	  (sc = rnapar->P10P + rnapar->P5P + rnapar->P6P 
	   + wkn*icfg[idxL(s[k-1])][idxP(s[l],s[k])]
	   + whx[j][d][d1-2][d2-1]) > bestsc)
	bestsc = sc; 
      
      /*  (YH12.3) */
      if (d1 > 0 && d2 > 1 &&
	  (sc = rnapar->P10P + rnapar->P5P + rnapar->P6P 
	   + wkn*icfg[idxR(s[l+1])][idxP(s[l],s[k])] 
	   + whx[j][d][d1-1][d2-2]) > bestsc)
	bestsc = sc; 
      
      /*  (YH12.4) */
      if (d1 > 1 && d2 > 1 &&
	  (sc = rnapar->P10P + rnapar->P5P + 2*rnapar->P6P 
	   + wkn*icfg[idxL(s[k-1])][idxP(s[l],s[k])]
	   + wkn*icfg[idxR(s[l+1])][idxP(s[l],s[k])] 
	   + whx[j][d][d1-2][d2-2]) > bestsc)
	bestsc = sc; 

      yhx[j][d][d1][d2] = bestsc;
      
    }  /* if (k,l) are base paired */
  
  else yhx[j][d][d1][d2] = -BIGINT;
}

/* Function: TraceYHX()
 * 
 * Purpose:  traceback hole matrix YHX. Using Knots-IS(2) method.
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
TraceYHX(FILE *outf, ESL_DSQ *s, int len, struct rnapar_2 *rnapar, int **icfg, int **wbx, int **vx, 
	 int ****whx, int ****vhx, int ****yhx,
	 int j, int d, int d1, int d2, 
	 struct tracekn_s *curr_tr, struct traceknstack_s *dolist, int traceback)
{
  int mid, mid1, mid2;  /* used for midpoint of a bifurc   */

  /* (YH1)__/  
   */
  if (yhx[j][d][d1][d2] == YH1(s, len, rnapar, icfg, vhx, j, d, d1, d2)){
    trace_YH1(outf, vhx, j, d, d1, d2, curr_tr, dolist, traceback);
    return;
  }   
  
  /* (YH2)__/ 
   */
  else if (yhx[j][d][d1][d2] == YH2(s, len, rnapar, icfg, vhx, j, d, d1, d2)){
    trace_YH2(outf, vhx, j, d, d1, d2, curr_tr, dolist, traceback);
    return;
  }
  
  /* (YH3)__/
   */
  else if (yhx[j][d][d1][d2] == YH3(s, len, rnapar, icfg, vhx, j, d, d1, d2)){
    trace_YH3(outf, vhx, j, d, d1, d2, curr_tr, dolist, traceback);
    return;
  } 

  /* (YH4)__/
   */
  else if (yhx[j][d][d1][d2] == YH4(s, len, rnapar, icfg, vhx, j, d, d1, d2)){
    trace_YH4(outf, vhx, j, d, d1, d2, curr_tr, dolist, traceback);
    return;
  }
  
  /* (YH5)__/ 
   */
  else if (yhx[j][d][d1][d2] == YH5(s, len, rnapar, icfg, yhx, j, d, d1, d2)){
    trace_YH5(outf, yhx, j, d, d1, d2, curr_tr, dolist, traceback);
    return;
  } 

  /* (YH6)__/ 
   */
  else if (yhx[j][d][d1][d2] ==  YH6(s, len, rnapar, icfg, yhx, j, d, d1, d2)){
    trace_YH6(outf, yhx, j, d, d1, d2, curr_tr, dolist, traceback);
    return;
  }
  
  for (mid = d1-1; mid > 0; mid--) {
    
    /* (YH7)__/ 1 WBX 1 YHX
     */                                       
    if (yhx[j][d][d1][d2] == YH7(s, len, rnapar, icfg, wbx, yhx, j, d, d1, d2, mid)){
      trace_YH7(outf, wbx, yhx, j, d, d1, d2, mid, curr_tr, dolist, traceback);
      return;
    }

    /* (YH8)__/ 1 VX 1 ZHX 
     */                                       
    /* (YH8.1) */
    else if (yhx[j][d][d1][d2] == YH8_1(s, len, rnapar, icfg, vx, vhx, j, d, d1, d2, mid)){
      trace_YH8_1(outf, vx, vhx, j, d, d1, d2, mid, curr_tr, dolist, traceback);
      return;
    }
  
    /* (YH8.2) */
    else if (yhx[j][d][d1][d2] == YH8_2(s, len, rnapar, icfg, vx, vhx, j, d, d1, d2, mid)){
      trace_YH8_2(outf, vx, vhx, j, d, d1, d2, mid, curr_tr, dolist, traceback);
      return;
    }

    /* (YH8.3) */
    else if (yhx[j][d][d1][d2] == YH8_3(s, len, rnapar, icfg, vx, vhx, j, d, d1, d2, mid)){
      trace_YH8_3(outf, vx, vhx, j, d, d1, d2, mid, curr_tr, dolist, traceback);
      return;
    } 

    /* (YH8.4) */
    else if (yhx[j][d][d1][d2] ==  YH8_4(s, len, rnapar, icfg, vx, vhx, j, d, d1, d2, mid)){
      trace_YH8_4(outf, vx, vhx, j, d, d1, d2, mid, curr_tr, dolist, traceback);
      return;
    } 
  }
  
  for (mid = d2 - 1; mid > 0; mid--) {

    /* (YH9)__/ 1 WBX 1 YHX
     */                                       
    if (yhx[j][d][d1][d2] == YH9(s, len, rnapar, icfg, wbx, yhx, j, d, d1, d2, mid)){
      trace_YH9(outf, wbx, yhx, j, d, d1, d2, mid, curr_tr, dolist, traceback);
      return;
    }
      
    /* (YH10)__/ 1 VX 1 YHX 
     */                                       
    
    /* (YH10.1) */
    else if (yhx[j][d][d1][d2] == YH10_1(s, len, rnapar, icfg, vx, vhx, j, d, d1, d2, mid)){
      trace_YH10_1(outf, vx, vhx, j, d, d1, d2, mid, curr_tr, dolist, traceback);
      return;
    }
     
    /* (YH10.2) */
    else if (yhx[j][d][d1][d2] == YH10_2(s, len, rnapar, icfg, vx, vhx, j, d, d1, d2, mid)){
      trace_YH10_2(outf, vx, vhx, j, d, d1, d2, mid, curr_tr, dolist, traceback);
      return;
    }
      
    /* (YH10.3) */
    else if (yhx[j][d][d1][d2] == YH10_3(s, len, rnapar, icfg, vx, vhx, j, d, d1, d2, mid)){
      trace_YH10_3(outf, vx, vhx, j, d, d1, d2, mid, curr_tr, dolist, traceback);
      return;
    }	

    /* (YH10.4) */
    else if (yhx[j][d][d1][d2] == YH10_4(s, len, rnapar, icfg, vx, vhx, j, d, d1, d2, mid)){
      trace_YH10_4(outf, vx, vhx, j, d, d1, d2, mid, curr_tr, dolist, traceback);
      return;
    }	
  }      
  
  /* (YH11)__/ 1 VHX 1 YHX
     */
  
  for (mid1 = d1; mid1 > 0; mid1--)
    for (mid2 = d2; mid2 > 0; mid2--)
      if (yhx[j][d][d1][d2] == YH11(s, len, rnapar, icfg, yhx, vhx, j, d, d1, d2, mid1, mid2)){
      trace_YH11(outf, yhx, vhx, j, d, d1, d2, mid1, mid2, curr_tr, dolist, traceback);
      return;
    }	

  /* (YH12)__/  REST of YHX (multiloops) */
  
  /* (YH12.1) */
  if (yhx[j][d][d1][d2] == YH12_1(s, len, rnapar, icfg, whx, j, d, d1, d2)){
    trace_YH12_1(outf, whx, j, d, d1, d2, curr_tr, dolist, traceback);
    return;
  }
  
  /* (YH12.2) */
  else  if (yhx[j][d][d1][d2] == YH12_2(s, len, rnapar, icfg, whx, j, d, d1, d2)){
    trace_YH12_2(outf, whx, j, d, d1, d2, curr_tr, dolist, traceback);
    return;
  }
  
  /* (YH12.3) */
  else  if (yhx[j][d][d1][d2] == YH12_3(s, len, rnapar, icfg, whx, j, d, d1, d2)){
    trace_YH12_3(outf, whx, j, d, d1, d2, curr_tr, dolist, traceback);
    return;
      }
  
  /* (YH12.4) */
  else  if (yhx[j][d][d1][d2] == YH12_4(s, len, rnapar, icfg, whx, j, d, d1, d2)){
    trace_YH12_4(outf, whx, j, d, d1, d2, curr_tr, dolist, traceback);
    return;
  }
}

