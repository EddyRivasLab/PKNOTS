/* filltrvx.c 
 *
 * includes functions FillVX and TraceVX
 * that fill and traceback the no-hole matix:
 *      vx[j][d] (len x len)
 *               [(j-d,j) are base pared]
 *
 */
                         
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>
#include <math.h>

#include <easel.h>

#include "pknots.h"
#include "pk_filltrvx.h"
#include "pk_irredsurf.h"
#include "pk_vxgraphs.h"
#include "pk_util.h"

static void FillVX_pks(ESL_DSQ *s, int len, struct rnapar_2 *rnapar, 
		       int **icfg, int **wbx, int **vx, 
		       int ****whx, int ****zhx, int ****yhx, int j, int d);
static void TraceVX_pks(FILE *outf, ESL_DSQ *s, int len, struct rnapar_2 *rnapar, 
			int **icfg, int **wbx, int **vx, 
			int ****whx, int ****zhx, int ****yhx, int j, int d, int *flag,
			struct tracekn_s *curr_tr, struct traceknstack_s *dolist, int traceback);

/* Function: FillVX() 
 * 
 * Purpose: fill no-hole matrix VX. 
 *           
 * Arguments: s     - sequence (0..len-1) converted to integers
 *                    representing array indices of symbols
 *            len   - length of iseq
 *           icfg   - context-free grammar state transitions, integer log form
 *            whx   - DP matrix, already alloc'ed 
 *            vhx   - DP matrix, already alloc'ed 
 *            zhx   - DP matrix, already alloc'ed 
 *            yhx   - DP matrix, already alloc'ed 
 *            wbx   - DP matrix, already alloc'ed 
 *             vx   - DP matrix, already alloc'ed 
 *            
 * Return:    (void)           
 *            vx is filled.
 */       
void
FillVX(ESL_DSQ *s, int len, struct rnapar_2 *rnapar, 
       int **icfg, int **wbx, int **vx, int *vp, 
       int ****whx, int ****zhx, int ****yhx, int j, int d, int allow_coaxials, int approx)
{
  /* no internal pseudoknots approximation, read it from the zuker algorithm.
   */
  FillVX_nested(s, len, rnapar, icfg, wbx, vx, vp, j, d, allow_coaxials);

  /* if approx = FALSE include internal pseudoknot  diagrams V7-V10. 
  */
  if (!approx) 
    FillVX_pks(s, len, rnapar, icfg, wbx, vx, whx, zhx, yhx, j, d); 

}

/* Function: TraceVX()
 * 
 * Purpose:  Trace back of VX diagrams with pseudoknots.
 *           
 * Arguments: s     - sequence (0..len-1) converted to integers
 *                    representing array indices of symbols
 *            len   - length of iseq
 *           icfg   - context-free grammar state transitions, integer log form
 *            whx   - DP matrix, already alloc'ed 
 *            vhx   - DP matrix, already alloc'ed 
 *            zhx   - DP matrix, already alloc'ed 
 *            yhx   - DP matrix, already alloc'ed 
 *            wbx   - DP matrix, already alloc'ed 
 *             vx   - DP matrix, already alloc'ed 
 *            
 * Return:    (void)           
 *            vx is filled.
 */       
void
TraceVX(FILE *outf, ESL_DSQ *s, int len, struct rnapar_2 *rnapar, 
	int **icfg, int **wbx, int **vx,  
	int ****whx, int ****zhx, int ****yhx, int j, int d, int approx,
	struct tracekn_s *curr_tr, struct traceknstack_s *dolist, int traceback)
{
  int fl, *flag;

  fl   = FALSE;
  flag =   &fl;

  /* no-internal-pseudoknots approximation, read it from the zuker algorithm (V1-V6).
   */
  TraceVX_nested(outf, s, len, rnapar, icfg, wbx, vx, j, d, flag, curr_tr, dolist, traceback);

  /* if approx = FALSE include internal pseudoknot  diagrams V7-V10. 
   */ 
  if (!approx && *flag == FALSE)
    TraceVX_pks(outf, s, len, rnapar, icfg, wbx, vx, whx, zhx, yhx, 
  j, d, flag, curr_tr, dolist, traceback); 
 
  if(*flag == FALSE)
    pk_fatal("something went wrong in the traceback of VX");
}

/* Function: FillVP() 
 * 
 * Purpose: calculate the best vp, for internal loops 
 *           
 * Arguments: s     - sequence (0..len-1) converted to integers
 *                    representing array indices of symbols
 *            cfg   - context-free grammar state transitions, float log2 form
 *            ntc   - context-free grammar nucleotide composition, float log2 form
 *             vx   - DP matrix, already alloc'ed 
 *            j,d   - coordenates of the matrix element
 *
 *            
 * Return:    (void)           
 *            vpbest is filled.
 */       
void 
FillVP(ESL_DSQ *s, int len, struct rnapar_2 *rnapar, 
       int **icfg, int **vx, int *vp, int j, int d)
{
  int mid;  /* midpoint of bifurcation         */
  int k;
  int sc;

  k = j - d + 2;
  
  if (d < 4) {
    vp[d] = -BIGINT;
    return;
  }
  else if (d > 4 &&
	   icfg[idxS][idxP(s[j-d],s[j])] == 1*INTSCALE && 
	   icfg[idxS][idxP(s[k],s[j-2])] == 1*INTSCALE &&
	   icfg[idxS][idxP(s[j-d+1],s[j-1])] != INTSCALE)
    vp[d-4] = rnapar->P3
      + wsf*icfg[idxPS(s[j-2],s[k])][idxPS(s[j-1],s[k-1])]
      + vx[j-2][d-4];
  else 
    vp[d-4] = -BIGINT;

  /* vp[d-5] is an special case because of the change from idcPS to idxPI 
   * it cannot be read from vp[d-4]
   */
  vp[d-5] = -BIGINT;
  /* two possibilities 
   */
  if (icfg[idxS][idxP(s[k+1],s[j-2])] == 1*INTSCALE &&
      (sc =  rnapar->P3
       + wsf*icfg[idxPI(s[j-2],s[k+1])][idxPI(s[j-1],s[k])]
       + vx[j-2][d-5]) > vp[d-5]) vp[d-5] = sc;
  /* or 
   */
  if (icfg[idxS][idxP(s[k],s[j-3])] == 1*INTSCALE &&
      (sc = rnapar->P3
       + wsf*icfg[idxPI(s[j-3],s[k])][idxPI(s[j-2],s[k-1])]
       + vx[j-3][d-5]) > vp[d-5]) vp[d-5] = sc;
  
  /* Now the rest of mid's
   */
  for (mid = 0; mid < (d-5); mid++) 
    if (icfg[idxS][idxP(s[k],s[k+mid])] == 1*INTSCALE &&
	(sc = rnapar->P3
	 + wsf*icfg[idxPI(s[k+mid],s[k])][idxPI(s[k+mid+1],s[k-1])]
	 + vx[k+mid][mid]) > vp[mid])
      vp[mid] = sc;
}

/* Function: FillVX_nested() 
 * 
 * Purpose: fill no-hole matrix VX without pseudoknots (zuker) 
 *           
 * Arguments: s     - sequence (0..len-1) converted to integers
 *                    representing array indices of symbols
 *            len   - length of iseq
 *           icfg   - context-free grammar state transitions, integer log form
 *            wbx   - DP matrix, already alloc'ed 
 *             vx   - DP matrix, already alloc'ed 
 *            
 * Return:    (void)           
 *            diagrams v1, v2, v3, v4, v5, v6 of vx are calculated.
 */       
void
FillVX_nested(ESL_DSQ *s, int len, struct rnapar_2 *rnapar,
	      int **icfg, int **wbx, int **vx, int *vp, 
	      int j, int d, int allow_coaxials)
{
  int                  i;  /* coords in mtx's                 */
  int                mid;  /* midpoint of bifurcation         */
  int         mid1, mid2;  /* used for midpoint of a bifurc   */
  int                 sc;  /* temporary score to check        */
  int             bestsc;  /* best score so far               */

  i = j - d;

  if (icfg[idxS][idxP(s[i],s[j])] == 1*INTSCALE)
    { 
      /* (V1)__/ IRREDUCIBLE SURFACES of  O(1)
       */
      bestsc =  F1(s, len, rnapar, icfg, j, d);	

      /* (V2)__/ IRREDUCIBLE SURFACES of  O(2)
       */
      /* STEMS
       */
      if ((sc = F2(s, len, rnapar, icfg, j, d, 1, 1)
	   + vx[j-1][d-2]) > bestsc)
	bestsc = sc;
      
      /* BULGES 
       */
      for (mid = 2; mid < d; mid++) {
	if ((sc = F2(s, len, rnapar, icfg, j, d, mid, 1)
	     + vx[j-1][d-1-mid]) > bestsc)
	  bestsc = sc;
	if ((sc = F2(s, len, rnapar, icfg, j, d, 1, mid)
	     + vx[j-mid][d-1-mid]) > bestsc)
	  bestsc = sc;
      }      

      if (1) {
	/* INTERNAL LOOPS (N^3)
	 */
	for (mid = 0; mid <= (d-4); mid++)     
	  if (icfg[idxS][idxP(s[j-d],s[j])] == 1*INTSCALE)
	    { 
	      if (d-mid > 32 &&
		  (sc = wsf*icfg[idxPI(s[i],s[j])][idxPI(s[i+1],s[j-1])] 
		   + wsf*(int)(PRELOG*log((d-mid-2)/30.))
		   + wsf*icfg[idxS][idxW(30)]
		   + vp[mid]) > bestsc)
		bestsc = sc;
	      if (d-mid == 4 && 
		  icfg[idxS][idxP(s[i+1],s[j-1])] != INTSCALE && 
		  (sc = wsf*icfg[idxPS(s[i],s[j])][idxPS(s[i+1],s[j-1])] 
		   + vp[mid]) > bestsc)
		bestsc = sc;
	      if (d-mid > 4 && d-mid < 33 &&
		  (sc = wsf*icfg[idxPI(s[i],s[j])][idxPI(s[i+1],s[j-1])] 
		   + wsf*icfg[idxS][idxW(d-mid-2)]
		   + vp[mid]) > bestsc)
		bestsc = sc;
	    }
      }
      else {
	/* INTERNAL LOOPS (N^4)
	 */
	for (mid1 = 2; mid1 < d; mid1++)
	  for (mid2 = 2; mid2 < (d-mid1); mid2++)
	    if ((sc = F2(s, len, rnapar, icfg, j, d, mid1, mid2)
		 + vx[j-mid2][d-mid1-mid2]) > bestsc)
	      bestsc = sc;
      }
      
      
      /* REST of VX (MULTILOOPS).  
       */
      
      for (mid = 1; mid < d-3; mid++)
	{                      
	  /* (V3)__/ 2 WBX
	   */
	  
	  /* (V3_1) */
	  if (mid < d-3 &&
	      (sc = rnapar->P10 + rnapar->P5 
	       + wbx[i+1+mid][mid] + wbx[j-1][d-mid-3]) > bestsc)
	    bestsc = sc;
	  
	  /* (V3_2) */
	  if (mid > 1 &&
	      (sc =  rnapar->P10 + rnapar->P5 + rnapar->P6 
	       + wsf*icfg[idxR(s[i+1])][idxP(s[j],s[i])]
	       + wbx[i+1+mid][mid-1] + wbx[j-1][d-mid-3]) > bestsc)
	    bestsc = sc;
	  
	  /* (V3_3) */
	  if (mid < d-4 &&
	      (sc = rnapar->P10 + rnapar->P5 + rnapar->P6 
	       + wsf*icfg[idxL(s[j-1])][idxP(s[j],s[i])]
	       + wbx[i+1+mid][mid] + wbx[j-2][d-mid-4]) > bestsc)
	    bestsc = sc;
	  
	  /* (V3_4) */
 	  if (mid < d-4 &&
	      (sc = rnapar->P10 + rnapar->P5 + 2*rnapar->P6 
	       + wsf*icfg[idxR(s[i+1])][idxP(s[j],s[i])] 
	       + wsf*icfg[idxL(s[j-1])][idxP(s[j],s[i])]
	       + wbx[i+1+mid][mid-1] + wbx[j-2][d-mid-4]) > bestsc)
	    bestsc = sc;
	  
	  /* (V4)__/ 1VX 1 WBX 
	   */ 

	  /* (V4_1) */
	  if (mid < d-3 &&
	      (sc = 2*rnapar->P10 + rnapar->P5 
	       + wsf*(icfg[idxPS(s[i],s[j])][idxPS(s[i+1],s[i+1+mid])] + INTSCALE)
	       + vx[i+1+mid][mid] + wbx[j-1][d-mid-3]) > bestsc)
	    bestsc = sc;
	  
	  /* (V4_2) */
	  if (mid < d-4 &&
	      (sc = 2*rnapar->P10 + rnapar->P5 + rnapar->P6 
	       + wsf*icfg[idxL(s[j-1])][idxP(s[j],s[i])]
	       + wsf*icfg[idxPS(s[i],s[j])][idxPS(s[i+1],s[i+1+mid])] 
	       + vx[i+1+mid][mid] + wbx[j-2][d-mid-4]) > bestsc)
	    bestsc = sc;

	  /* (V4_3) */
	  if (mid < d-4 &&
	      (sc = 2*rnapar->P10 + rnapar->P5 + rnapar->P6
	       + wsf*icfg[idxR(s[i+mid+2])][idxP(s[i+1],s[i+mid+1])]
	       + wsf*icfg[idxPS(s[i],s[j])][idxPS(s[i+1],s[i+1+mid])] 
	       + vx[i+1+mid][mid] + wbx[j-1][d-mid-4]) > bestsc)
	    bestsc = sc;
	  
	  /* (V4_4) */
	  if (mid < d-5 &&
	      (sc = 2*rnapar->P10 + rnapar->P5 + 2*rnapar->P6 
	       + wsf*icfg[idxL(s[j-1])][idxP(s[j],s[i])]
	       + wsf*icfg[idxR(s[i+mid+2])][idxP(s[i+1],s[i+mid+1])]
	       + wsf*icfg[idxPS(s[i],s[j])][idxPS(s[i+1],s[i+1+mid])] 
	       + vx[i+1+mid][mid] + wbx[j-2][d-mid-5]) > bestsc)
	    bestsc = sc;
	  
	  /* (V4_5) */
	  if (mid < d-4 &&
	      (sc = 2*rnapar->P10 + rnapar->P5 + rnapar->P6
	       + wsf*icfg[idxPS(s[i],s[j])][idxPS(s[i+2],s[i+2+mid])]
	       + vx[i+2+mid][mid] + wbx[j-1][d-mid-4]) > bestsc)
	    bestsc = sc;
	  
	  /* (V4_6) */
	  if (mid < d-5 &&
	      (sc = 2*rnapar->P10 + rnapar->P5 + 2*rnapar->P6 
	       + wsf*icfg[idxL(s[j-1])][idxP(s[j],s[i])]
	       + wsf*icfg[idxPS(s[i],s[j])][idxPS(s[i+2],s[i+2+mid])] 
	       + vx[i+2+mid][mid] + wbx[j-2][d-mid-5]) > bestsc)
	    bestsc = sc;

	  /* (V4_7) */
	  if (mid < d-4 &&
	      (sc = 2*rnapar->P10 + rnapar->P5 + 2*rnapar->P6
	       + wsf*icfg[idxR(s[i+mid+3])][idxP(s[i+2],s[i+mid+2])]
	       + wsf*icfg[idxPS(s[i],s[j])][idxPS(s[i+2],s[i+2+mid])] 
	       + vx[i+2+mid][mid] + wbx[j-1][d-mid-5]) > bestsc)
	    bestsc = sc;

	  /* (V4_8) */
	  if (mid < d-6 &&
	      (sc = 2*rnapar->P10 + rnapar->P5 + 3*rnapar->P6 
	       + wsf*icfg[idxL(s[j-1])][idxP(s[j],s[i])]
	       + wsf*icfg[idxR(s[i+mid+3])][idxP(s[i+2],s[i+mid+2])]
	       + wsf*icfg[idxPS(s[i],s[j])][idxPS(s[i+2],s[i+2+mid])] 
	       + vx[i+2+mid][mid] + wbx[j-2][d-mid-6]) > bestsc)
	    bestsc = sc;
	  
	  /* (V5)__/ 1 WBX 1VX 
	   */ 
	  
	  /* (V5_1) */
	  if (mid < d-3 &&
	      (sc = 2*rnapar->P10 + rnapar->P5 
	       + wsf*(icfg[idxPS(s[i],s[j])][idxPS(s[i+2+mid],s[j-1])] + INTSCALE)
	       + wbx[i+1+mid][mid] + vx[j-1][d-mid-3]) > bestsc)
	    bestsc = sc;
	  
	  /* (V5_2)  */
	  if (mid > 1 && mid < d-3 &&
	      (sc = 2*rnapar->P10 + rnapar->P5 + rnapar->P6 
	       + wsf*icfg[idxR(s[i+1])][idxP(s[j],s[i])]
	       + wsf*icfg[idxPS(s[i],s[j])][idxPS(s[i+2+mid],s[j-1])]
	       + wbx[i+1+mid][mid-1] + vx[j-1][d-mid-3]) > bestsc)
	    bestsc = sc;
	  
	  /* (V5_3)  */
	  if (mid > 1 && mid < d-3 &&
	      (sc = 2*rnapar->P10 + rnapar->P5 + rnapar->P6  
	       + wsf*icfg[idxL(s[i+mid+1])][idxP(s[i+2+mid],s[j-1])]
	       + wsf*icfg[idxPS(s[i],s[j])][idxPS(s[i+2+mid],s[j-1])]
	       + wbx[i+mid][mid-1] + vx[j-1][d-mid-3]) > bestsc)
	    bestsc = sc;
	  
	  /* (V5_4) */
	  if (mid > 1 && mid < d-3 &&
	      (sc = 2*rnapar->P10 + rnapar->P5 + 2*rnapar->P6 
	       + wsf*icfg[idxR(s[i+1])][idxP(s[j],s[i])]
	       + wsf*icfg[idxL(s[i+mid+1])][idxP(s[i+2+mid],s[j-1])]
	       + wsf*icfg[idxPS(s[i],s[j])][idxPS(s[i+2+mid],s[j-1])]
	       + wbx[i+mid][mid-2] + vx[j-1][d-mid-3]) > bestsc)
	    bestsc = sc;
	  
	  /* (V5_5) */
	  if (mid < d-4 &&
	      (sc = 2*rnapar->P10 + rnapar->P5 + rnapar->P6
	       + wsf*icfg[idxPS(s[i],s[j])][idxPS(s[i+2+mid],s[j-2])]
	       + wbx[i+1+mid][mid] + vx[j-2][d-mid-4]) > bestsc)
	    bestsc = sc;
	  
	  /* (V5_6) */
	  if (mid > 1 && mid < d-4 &&
	      (sc = 2*rnapar->P10 + rnapar->P5 + 2*rnapar->P6 
	       + wsf*icfg[idxR(s[i+1])][idxP(s[j],s[i])]
	       + wsf*icfg[idxPS(s[i],s[j])][idxPS(s[i+2+mid],s[j-2])]
	       + wbx[i+1+mid][mid-1] + vx[j-2][d-mid-4]) > bestsc)
	    bestsc = sc;
	  
	  /* (V5_7)  */
	  if (mid > 1 && mid < d-4 &&
	      (sc = 2*rnapar->P10 + rnapar->P5 + 2*rnapar->P6  
	       + wsf*icfg[idxL(s[i+mid+1])][idxP(s[i+2+mid],s[j-2])]
	       + wsf*icfg[idxPS(s[i],s[j])][idxPS(s[i+2+mid],s[j-2])]
	       + wbx[i+mid][mid-1] + vx[j-2][d-mid-4]) > bestsc)
	    bestsc = sc;
	  
	  /* (V5_8)  */
	  if (mid > 2 && mid < d-4 &&
	      (sc = 2*rnapar->P10 + rnapar->P5 + 3*rnapar->P6 
	       + wsf*icfg[idxR(s[i+1])][idxP(s[j],s[i])]
	       + wsf*icfg[idxL(s[i+mid+1])][idxP(s[i+2+mid],s[j-2])]
	       + wsf*icfg[idxPS(s[i],s[j])][idxPS(s[i+2+mid],s[j-2])]
	       + wbx[i+mid][mid-2] + vx[j-2][d-mid-4]) > bestsc)
	  bestsc = sc;
	  
	  /* (V6)__/ 2 VX 
	   */
	  if (allow_coaxials) 
	    for (mid1 = 1; mid1 < d-mid-2; mid1++)
	      for (mid2 = 1; mid2 < d-mid-mid1-1; mid2++)
		{
		  /* (V6_1) */
		  if (mid < d-3 && mid1 == 1 && mid2 == 1 &&
		      (sc = 3*rnapar->P10 + rnapar->P5
		       + wsf*(icfg[idxPS(s[i+mid+mid1],s[i+mid1])]
			      [idxPS(s[i+mid+mid1+1],s[j-mid2])] + INTSCALE)
		       + vx[i+mid+mid1][mid] + vx[j-mid2][d-mid-mid1-mid2-1]) > bestsc)
		    bestsc = sc;
		  
		  /* (V6_2a) */
		  if (mid < d-3 && mid1 == 2 && mid2 == 1 &&
		      (sc = 3*rnapar->P10 + rnapar->P5 + rnapar->P6*(mid1-1) 
		       + wsf*icfg[idxR(s[i+1])][idxP(s[j],s[i])]
		       + wsf*icfg[idxPS(s[i+mid+mid1],s[i+mid1])]
		       [idxPS(s[i+mid+mid1+1],s[j-mid2])]
		       + vx[i+mid+mid1][mid] + vx[j-mid2][d-mid-mid1-mid2-1]) > bestsc)
		    bestsc = sc;
		
		  /* (V6_2b)  */
		  if (mid < d-3 && mid1 == 2 && mid2 == 1 &&
		      (sc = 3*rnapar->P10 + rnapar->P5 + rnapar->P6*(mid1-1) 
		       + wsf*icfg[idxL(s[i+1])][idxP(s[i+mid1],s[i+mid1+mid])]
		       + wsf*icfg[idxPS(s[i+mid+mid1],s[i+mid1])]
		       [idxPS(s[i+mid+mid1+1],s[j-mid2])]
		       + vx[i+mid+mid1][mid] + vx[j-mid2][d-mid-mid1-mid2-1]) > bestsc)
		  bestsc = sc;
		  
		  /* (V6_2c) */
		  if (mid < d-3 && mid1 > 2 && mid2 == 1 &&
		      (sc = 3*rnapar->P10 + rnapar->P5 + rnapar->P6*(mid1-1) 
		       + wsf*icfg[idxL(s[i+1])][idxP(s[i+mid1],s[i+mid1+mid])]
		       + wsf*icfg[idxL(s[i+mid1-1])][idxP(s[i+mid1],s[i+mid1+mid])]
		       + wsf*icfg[idxPS(s[i+mid+mid1],s[i+mid1])]
		       [idxPS(s[i+mid+mid1+1],s[j-mid2])]
		       + vx[i+mid+mid1][mid] + vx[j-mid2][d-mid-mid1-mid2-1]) > bestsc)
		  bestsc = sc;
		  
		  /* (V6_3a) */
		  if (mid < d-3 && mid1 == 1 && mid2 == 2 &&
		      (sc =3*rnapar->P10 + rnapar->P5 + rnapar->P6*(mid2-1) 
		       + wsf*icfg[idxL(s[j-1])][idxP(s[j],s[i])]
		     + wsf*icfg[idxPS(s[i+mid+mid1],s[i+mid1])]
		       [idxPS(s[i+mid+mid1+1],s[j-mid2])]
		       + vx[i+mid+mid1][mid] + vx[j-mid2][d-mid-mid1-mid2-1] ) > bestsc)
		    bestsc = sc;
		  
		  /* (V6_3b) */
		  if (mid < d-3 && mid1 == 1 && mid2 == 2 &&
		      (sc = 3*rnapar->P10 + rnapar->P5 + rnapar->P6*(mid2-1) 
		       + wsf*icfg[idxR(s[j-1])][idxP(s[i+mid+mid1+1],s[j-mid2])]
		       + wsf*icfg[idxPS(s[i+mid+mid1],s[i+mid1])]
		       [idxPS(s[i+mid+mid1+1],s[j-mid2])]
		     + vx[i+mid+mid1][mid] + vx[j-mid2][d-mid-mid1-mid2-1]) > bestsc)
		    bestsc = sc;
		  
		  /* (V6_3c) */
		  if (mid < d-3 && mid1 == 1 && mid2 > 2 &&
		    (sc = 3*rnapar->P10 + rnapar->P5 + rnapar->P6*(mid2-1) 
		     + wsf*icfg[idxL(s[j-1])][idxP(s[j],s[i])]
		     + wsf*icfg[idxR(s[j-mid2+1])][idxP(s[i+mid+mid1+1],s[j-mid2])]
		     + wsf*icfg[idxPS(s[i+mid+mid1],s[i+mid1])]
		     [idxPS(s[i+mid+mid1+1],s[j-mid2])]
		     + vx[i+mid+mid1][mid] + vx[j-mid2][d-mid-mid1-mid2-1]) > bestsc)
		    bestsc = sc;
		  
		/* (V6_4a) */
		  if (mid < d-3 && mid1 == 2 && mid2 == 2 &&
		      (sc = 3*rnapar->P10 + rnapar->P5 + rnapar->P6*(mid1+mid2-2) 
		       + wsf*icfg[idxR(s[i+1])][idxP(s[j],s[i])] 
		       + wsf*icfg[idxL(s[j-1])][idxP(s[j],s[i])]
		       + wsf*icfg[idxPS(s[i+mid+mid1],s[i+mid1])]
		               [idxPS(s[i+mid+mid1+1],s[j-mid2])]
		       + vx[i+mid+mid1][mid] + vx[j-mid2][d-mid-mid1-mid2-1]) > bestsc)
		    bestsc = sc;
		  
		  /* (V6_4b) */
		  if (mid < d-3 && mid1 == 2 && mid2 == 2 &&
		      (sc = 3*rnapar->P10 + rnapar->P5 + rnapar->P6*(mid1+mid2-2) 
		       + wsf*icfg[idxR(s[i+1])][idxP(s[j],s[i])]
		       + wsf*icfg[idxR(s[j-1])][idxP(s[i+mid+mid1+1],s[j-mid2])]
		     + wsf*icfg[idxPS(s[i+mid+mid1],s[i+mid1])]
		       [idxPS(s[i+mid+mid1+1],s[j-mid2])]
		       + vx[i+mid+mid1][mid] + vx[j-mid2][d-mid-mid1-mid2-1]) > bestsc)
		    bestsc = sc;
		  
		  /* (V6_4c) */
		  if (mid < d-3 && mid1 == 2 && mid2 == 2 &&
		    (sc =3*rnapar->P10 + rnapar->P5 + rnapar->P6*(mid1+mid2-2) 
		     + wsf*icfg[idxL(s[j-1])][idxP(s[j],s[i])]
		     + wsf*icfg[idxL(s[i+1])][idxP(s[i+mid1],s[i+mid1+mid])]
		     + wsf*icfg[idxPS(s[i+mid+mid1],s[i+mid1])]
		     [idxPS(s[i+mid+mid1+1],s[j-mid2])]
		     + vx[i+mid+mid1][mid] + vx[j-mid2][d-mid-mid1-mid2-1] ) > bestsc)
		    bestsc = sc;
		  
		/* (V6_4d) */
		  if (mid < d-3 && mid1 == 2 && mid2 == 2 &&
		      (sc = 3*rnapar->P10 + rnapar->P5 + rnapar->P6*(mid1+mid2-2) 
		       + wsf*icfg[idxL(s[i+1])][idxP(s[i+mid1],s[i+mid1+mid])]
		       + wsf*icfg[idxR(s[j-1])][idxP(s[i+mid+mid1+1],s[j-mid2])]
		       + wsf*icfg[idxPS(s[i+mid+mid1],s[i+mid1])]
		       [idxPS(s[i+mid+mid1+1],s[j-mid2])]
		     + vx[i+mid+mid1][mid] + vx[j-mid2][d-mid-mid1-mid2-1]) > bestsc)
		    bestsc = sc;
		  
		  /* (V6_4e) */
		  if (mid < d-3 && mid1 > 2 && mid2 == 2 &&
		      (sc = 3*rnapar->P10 + rnapar->P5 + rnapar->P6*(mid1+mid2-2) 
		     + wsf*icfg[idxR(s[i+1])][idxP(s[j],s[i])] 
		       + wsf*icfg[idxL(s[j-1])][idxP(s[j],s[i])]
		       + wsf*icfg[idxL(s[i+mid1-1])][idxP(s[i+mid1],s[i+mid1+mid])]
		       + wsf*icfg[idxPS(s[i+mid+mid1],s[i+mid1])]
		       [idxPS(s[i+mid+mid1+1],s[j-mid2])]
		       + vx[i+mid+mid1][mid] + vx[j-mid2][d-mid-mid1-mid2-1]) > bestsc)
		    bestsc = sc;
		  
		  /* (V6_4f) */
		  if (mid < d-3 && mid1 > 2 && mid2 == 2 &&
		      (sc = 3*rnapar->P10 + rnapar->P5 + rnapar->P6*(mid1+mid2-2) 
		       + wsf*icfg[idxR(s[i+1])][idxP(s[j],s[i])]
		       + wsf*icfg[idxL(s[i+mid1-1])][idxP(s[i+mid1],s[i+mid1+mid])]
		       + wsf*icfg[idxR(s[j-1])][idxP(s[i+mid+mid1+1],s[j-mid2])]
		       + wsf*icfg[idxPS(s[i+mid+mid1],s[i+mid1])]
		       [idxPS(s[i+mid+mid1+1],s[j-mid2])]
		     + vx[i+mid+mid1][mid] + vx[j-mid2][d-mid-mid1-mid2-1]) > bestsc)
		    bestsc = sc;
		  
		  /* (V6_4g) */
		  if (mid < d-3 && mid1 == 2 && mid2 > 2 &&
		      (sc = 3*rnapar->P10 + rnapar->P5 + rnapar->P6*(mid1+mid2-2) 
		       + wsf*icfg[idxR(s[i+1])][idxP(s[j],s[i])]
		       + wsf*icfg[idxL(s[j-1])][idxP(s[j],s[i])]
		       + wsf*icfg[idxR(s[j-mid2+1])][idxP(s[i+mid+mid1+1],s[j-mid2])]
		     + wsf*icfg[idxPS(s[i+mid+mid1],s[i+mid1])]
		       [idxPS(s[i+mid+mid1+1],s[j-mid2])]
		       + vx[i+mid+mid1][mid] + vx[j-mid2][d-mid-mid1-mid2-1]) > bestsc)
		    bestsc = sc;
		  
		  /* (V6_4h) */
		  if (mid < d-3 && mid1 == 2 && mid2 > 2 &&
		    (sc = 3*rnapar->P10 + rnapar->P5 + rnapar->P6*(mid1+mid2-2) 
		     + wsf*icfg[idxL(s[j-1])][idxP(s[j],s[i])]
		     + wsf*icfg[idxL(s[i+1])][idxP(s[i+mid1],s[i+mid1+mid])]
		     + wsf*icfg[idxR(s[j-mid2+1])][idxP(s[i+mid+mid1+1],s[j-mid2])]
		     + wsf*icfg[idxPS(s[i+mid+mid1],s[i+mid1])]
		     [idxPS(s[i+mid+mid1+1],s[j-mid2])]
		     + vx[i+mid+mid1][mid] + vx[j-mid2][d-mid-mid1-mid2-1]) > bestsc)
		    bestsc = sc;
		  
		  /* (V6_4i) */
		  if (mid < d-3 && mid1 > 2 && mid2 > 2 &&
		      (sc = 3*rnapar->P10 + rnapar->P5 + rnapar->P6*(mid1+mid2-2) 
		       + wsf*icfg[idxR(s[i+1])][idxP(s[j],s[i])]
		       + wsf*icfg[idxL(s[j-1])][idxP(s[j],s[i])]
		       + wsf*icfg[idxL(s[i+mid1-1])][idxP(s[i+mid1],s[i+mid1+mid])]
		       + wsf*icfg[idxR(s[j-mid2+1])][idxP(s[i+mid+mid1+1],s[j-mid2])]
		       + wsf*icfg[idxPS(s[i+mid+mid1],s[i+mid1])]
		       [idxPS(s[i+mid+mid1+1],s[j-mid2])]
		       + vx[i+mid+mid1][mid] + vx[j-mid2][d-mid-mid1-mid2-1]) > bestsc)
		    bestsc = sc;
		  
		  /* (V6_5) */
		  if (mid < d-4 && mid1 == 1 && mid2 == 1 &&
		      (sc = 3*rnapar->P10 + rnapar->P5 + rnapar->P6
		       + wsf*icfg[idxPS(s[i+mid+mid1],s[i+mid1])]
		       [idxPS(s[i+mid+mid1+2],s[j-mid2])] 
		       + vx[i+mid+mid1][mid] + vx[j-mid2][d-mid-mid1-mid2-2]) > bestsc)
		    bestsc = sc;
		  
		  /* (V6_6a) */
		  if (mid < d-4 && mid1 == 2 && mid2 == 1 &&
		      (sc = 3*rnapar->P10 + rnapar->P5 + rnapar->P6*mid1 
		       + wsf*icfg[idxR(s[i+1])][idxP(s[j],s[i])]
		       + wsf*icfg[idxPS(s[i+mid+mid1],s[i+mid1])]
		       [idxPS(s[i+mid+mid1+2],s[j-mid2])]
		       + vx[i+mid+mid1][mid] + vx[j-mid2][d-mid-mid1-mid2-2]) > bestsc)
		    bestsc = sc;
		  
		  /* (V6_6b)  */
		  if (mid < d-4 && mid1 == 2 && mid2 == 1 &&
		      (sc = 3*rnapar->P10 + rnapar->P5 + rnapar->P6*mid1 
		       + wsf*icfg[idxL(s[i+1])][idxP(s[i+mid1],s[i+mid1+mid])]
		       + wsf*icfg[idxPS(s[i+mid+mid1],s[i+mid1])]
		       [idxPS(s[i+mid+mid1+2],s[j-mid2])]
		       + vx[i+mid+mid1][mid] + vx[j-mid2][d-mid-mid1-mid2-2]) > bestsc)
		    bestsc = sc;
		
		  /* (V6_6c) */
		  if (mid < d-4 && mid1 > 2 && mid2 == 1 &&
		      (sc = 3*rnapar->P10 + rnapar->P5 + rnapar->P6*mid1 
		       + wsf*icfg[idxL(s[i+1])][idxP(s[i+mid1],s[i+mid1+mid])]
		       + wsf*icfg[idxL(s[i+mid1-1])][idxP(s[i+mid1],s[i+mid1+mid])]
		       + wsf*icfg[idxPS(s[i+mid+mid1],s[i+mid1])]
		       [idxPS(s[i+mid+mid1+2],s[j-mid2])]
		       + vx[i+mid+mid1][mid] + vx[j-mid2][d-mid-mid1-mid2-2]) > bestsc)
		    bestsc = sc;
		  
		  /* (V6_7a) */
		  if (mid < d-4 && mid1 == 1 && mid2 == 2 &&
		      (sc = 3*rnapar->P10 + rnapar->P5 + rnapar->P6*mid2
		       + wsf*icfg[idxL(s[j-1])][idxP(s[j],s[i])]
		       + wsf*icfg[idxPS(s[i+mid+mid1],s[i+mid1])]
		       [idxPS(s[i+mid+mid1+2],s[j-mid2])]
		       + vx[i+mid+mid1][mid] + vx[j-mid2][d-mid-mid1-mid2-2]) > bestsc)
		    bestsc = sc;
		  
		  /* (V6_7b) */
		  if (mid < d-4 && mid1 == 1 && mid2 == 2 &&
		      (sc = 3*rnapar->P10 + rnapar->P5 + rnapar->P6*mid2 
		       + wsf*icfg[idxR(s[j-1])][idxP(s[i+mid+mid1+1],s[j-mid2])]
		       + wsf*icfg[idxPS(s[i+mid+mid1],s[i+mid1])]
		       [idxPS(s[i+mid+mid1+2],s[j-mid2])]
		       + vx[i+mid+mid1][mid] + vx[j-mid2][d-mid-mid1-mid2-2]) > bestsc)
		    bestsc = sc;
		  
		  /* (V6_7c) */
		  if (mid < d-4 && mid1 == 1 && mid2 > 2 &&
		      (sc = 3*rnapar->P10 + rnapar->P5 + rnapar->P6*mid2 
		       + wsf*icfg[idxL(s[j-1])][idxP(s[j],s[i])]
		       + wsf*icfg[idxR(s[j-mid2+1])][idxP(s[i+mid+mid1+1],s[j-mid2])]
		       + wsf*icfg[idxPS(s[i+mid+mid1],s[i+mid1])]
		       [idxPS(s[i+mid+mid1+2],s[j-mid2])]
		       + vx[i+mid+mid1][mid] + vx[j-mid2][d-mid-mid1-mid2-2]) > bestsc)
		    bestsc = sc;
		  
		  /* (V6_8a) */
		  if (mid < d-4 && mid1 == 2 && mid2 == 2 &&
		      (sc = 3*rnapar->P10 + rnapar->P5 + rnapar->P6*(mid1+mid2-1) 
		       + wsf*icfg[idxR(s[i+1])][idxP(s[j],s[i])] 
		       + wsf*icfg[idxL(s[j-1])][idxP(s[j],s[i])]
		       + wsf*icfg[idxPS(s[i+mid+mid1],s[i+mid1])]
		       [idxPS(s[i+mid+mid1+2],s[j-mid2])]
		       + vx[i+mid+mid1][mid] + vx[j-mid2][d-mid-mid1-mid2-2]) > bestsc)
		    bestsc = sc;
		  
		  /* (V6_8b) */
		  if (mid < d-4 && mid1 == 2 && mid2 == 2 &&
		      (sc = 3*rnapar->P10 + rnapar->P5 + rnapar->P6*(mid1+mid2-1) 
		       + wsf*icfg[idxR(s[i+1])][idxP(s[j],s[i])]
		       + wsf*icfg[idxR(s[j-1])][idxP(s[i+mid+mid1+1],s[j-mid2])]
		       + wsf*icfg[idxPS(s[i+mid+mid1],s[i+mid1])]
		       [idxPS(s[i+mid+mid1+2],s[j-mid2])]
		       + vx[i+mid+mid1][mid] + vx[j-mid2][d-mid-mid1-mid2-2]) > bestsc)
		    bestsc = sc;
		  
		  /* (V6_8c) */
		  if (mid < d-4 && mid1 == 2 && mid2 == 2 &&
		      (sc = 3*rnapar->P10 + rnapar->P5 + rnapar->P6*(mid1+mid2-1) 
		       + wsf*icfg[idxL(s[j-1])][idxP(s[j],s[i])]
		       + wsf*icfg[idxL(s[i+1])][idxP(s[i+mid1],s[i+mid1+mid])]
		       + wsf*icfg[idxPS(s[i+mid+mid1],s[i+mid1])]
		       [idxPS(s[i+mid+mid1+2],s[j-mid2])]
		       + vx[i+mid+mid1][mid] + vx[j-mid2][d-mid-mid1-mid2-2]) > bestsc)
		    bestsc = sc;
		  
		  /* (V6_8d) */
		  if (mid < d-4 && mid1 == 2 && mid2 == 2 &&
		      (sc = 3*rnapar->P10 + rnapar->P5 + rnapar->P6*(mid1+mid2-1) 
		       + wsf*icfg[idxL(s[i+1])][idxP(s[i+mid1],s[i+mid1+mid])]
		       + wsf*icfg[idxR(s[j-1])][idxP(s[i+mid+mid1+1],s[j-mid2])]
		       + wsf*icfg[idxPS(s[i+mid+mid1],s[i+mid1])]
		       [idxPS(s[i+mid+mid1+2],s[j-mid2])]
		       + vx[i+mid+mid1][mid] + vx[j-mid2][d-mid-mid1-mid2-2]) > bestsc)
		    bestsc = sc;
		  
		  /* (V6_8e) */
		  if (mid < d-4 && mid1 > 2 && mid2 == 2 &&
		      (sc =3*rnapar->P10 + rnapar->P5 + rnapar->P6*(mid1+mid2-1) 
		       + wsf*icfg[idxR(s[i+1])][idxP(s[j],s[i])] 
		       + wsf*icfg[idxL(s[j-1])][idxP(s[j],s[i])]
		       + wsf*icfg[idxL(s[i+mid1-1])][idxP(s[i+mid1],s[i+mid1+mid])]
		       + wsf*icfg[idxPS(s[i+mid+mid1],s[i+mid1])]
		       [idxPS(s[i+mid+mid1+2],s[j-mid2])]
		       + vx[i+mid+mid1][mid] + vx[j-mid2][d-mid-mid1-mid2-2] ) > bestsc)
		    bestsc = sc;
		  
		  /* (V6_8f) */
		  if (mid < d-4 && mid1 > 2 && mid2 == 2 &&
		      (sc = 3*rnapar->P10 + rnapar->P5 + rnapar->P6*(mid1+mid2-1) 
		       + wsf*icfg[idxR(s[i+1])][idxP(s[j],s[i])]
		       + wsf*icfg[idxL(s[i+mid1-1])][idxP(s[i+mid1],s[i+mid1+mid])]
		       + wsf*icfg[idxR(s[j-1])][idxP(s[i+mid+mid1+1],s[j-mid2])]
		       + wsf*icfg[idxPS(s[i+mid+mid1],s[i+mid1])]
		       [idxPS(s[i+mid+mid1+2],s[j-mid2])]
		       + vx[i+mid+mid1][mid] + vx[j-mid2][d-mid-mid1-mid2-2]) > bestsc)
		    bestsc = sc;
		  
		  /* (V6_8g) */
		  if (mid < d-4 && mid1 == 2 && mid2 > 2 &&
		      (sc = 3*rnapar->P10 + rnapar->P5 + rnapar->P6*(mid1+mid2-1) 
		       + wsf*icfg[idxR(s[i+1])][idxP(s[j],s[i])]
		       + wsf*icfg[idxL(s[j-1])][idxP(s[j],s[i])]
		       + wsf*icfg[idxR(s[j-mid2+1])][idxP(s[i+mid+mid1+1],s[j-mid2])]
		       + wsf*icfg[idxPS(s[i+mid+mid1],s[i+mid1])]
		       [idxPS(s[i+mid+mid1+2],s[j-mid2])]
		       + vx[i+mid+mid1][mid] + vx[j-mid2][d-mid-mid1-mid2-2]) > bestsc)
		    bestsc = sc;
		  
		/* (V6_8h) */
		  if (mid < d-4 && mid1 == 2 && mid2 > 2 &&
		      (sc = 3*rnapar->P10 + rnapar->P5 + rnapar->P6*(mid1+mid2-1) 
		       + wsf*icfg[idxL(s[j-1])][idxP(s[j],s[i])]
		       + wsf*icfg[idxL(s[i+1])][idxP(s[i+mid1],s[i+mid1+mid])]
		       + wsf*icfg[idxR(s[j-mid2+1])][idxP(s[i+mid+mid1+1],s[j-mid2])]
		       + wsf*icfg[idxPS(s[i+mid+mid1],s[i+mid1])]
		       [idxPS(s[i+mid+mid1+2],s[j-mid2])]
		       + vx[i+mid+mid1][mid] + vx[j-mid2][d-mid-mid1-mid2-2]) > bestsc)
		    bestsc = sc;
		  
		  /* (V6_8i) */
		  if (mid < d-4 && mid1 > 2 && mid2 > 2 &&
		      (sc = 3*rnapar->P10 + rnapar->P5 + rnapar->P6*(mid1+mid2-1) 
		       + wsf*icfg[idxR(s[i+1])][idxP(s[j],s[i])]
		       + wsf*icfg[idxL(s[j-1])][idxP(s[j],s[i])]
		       + wsf*icfg[idxL(s[i+mid1-1])][idxP(s[i+mid1],s[i+mid1+mid])]
		       + wsf*icfg[idxR(s[j-mid2+1])][idxP(s[i+mid+mid1+1],s[j-mid2])]
		       + wsf*icfg[idxPS(s[i+mid+mid1],s[i+mid1])]
		       [idxPS(s[i+mid+mid1+2],s[j-mid2])]
		       + vx[i+mid+mid1][mid] + vx[j-mid2][d-mid-mid1-mid2-2]) > bestsc)
		    bestsc = sc;
		}
	}
      vx[j][d] = bestsc;
    }  /* if (i,j) are base-paired, else prohibited  */
  
  else vx[j][d] = -BIGINT;
}

/* Function: TraceVX_nested()
 * 
 * Purpose:  Trace back of VX diagrams without pseudoknots.
 *           
 * Arguments: s     - sequence (0..len-1) converted to integers
 *                    representing array indices of symbols
 *            len   - length of iseq
 *           icfg   - context-free grammar state transitions, integer log form
 *            wbx   - DP matrix, already alloc'ed 
 *             vx   - DP matrix, already alloc'ed 
 *            
 * Return:    (void)           
 *            vx is filled.
 */       
void
TraceVX_nested(FILE *outf, ESL_DSQ *s, int len, struct rnapar_2 *rnapar,
	       int **icfg, int **wbx, int **vx, int j, int d,  
	       int *flag, struct tracekn_s *curr_tr, struct traceknstack_s *dolist, int traceback)
{
  int                  i;  /* coords in mtx's                 */
  int    mid, mid1, mid2;  /* used for midpoint of a bifurc   */
  
  i = j - d;

  /* (V1)__/  IRREDUCIBLE SURFACES O(1) 
   */
  if (vx[j][d] == V1(s, len, rnapar, icfg, j, d)){
    trace_V1(outf, vx, j, d, curr_tr, dolist, traceback);
    *flag = TRUE;
    return;
  }
  
  /* (V2)__/  IRREDUCIBLE SURFACES O(2) 
   */
  for (mid1 = 1; mid1 < d-1; mid1++)
    for (mid2 = 1; mid2 < d-mid1; mid2++)
      if (vx[j][d] == V2(s, len, rnapar, icfg, vx, j, d, mid1, mid2)) {
	trace_V2(outf, vx, s, len, rnapar, icfg, j, d, mid1, mid2, curr_tr, dolist, traceback);
	*flag = TRUE;
	return;
      }
  
  /* REST OF VX (MULTILOOPS)
   */
  
  /* 2 NO-HOLE STRUCTURES (four types of diagrams: V3-V6) */
  for (mid = 1; mid < d-3; mid++)
    {
      /* (V3)__/  2 WBX 
       */
      
      /* (V3_1) */
      if (vx[j][d] == V3_1(s, len, rnapar, icfg, wbx, j, d, mid)){
	trace_V3_1(outf, wbx, j, d, mid, curr_tr, dolist, traceback);
	*flag = TRUE;
	return;
      }
      
      /* (V3_2) */
      else if (vx[j][d] == V3_2(s, len, rnapar, icfg, wbx, j, d, mid)){
	trace_V3_2(outf, wbx, j, d, mid, curr_tr, dolist, traceback);
	*flag = TRUE;
	return;
      }
      
      /* (V3_3) */
      else if (vx[j][d] == V3_3(s, len, rnapar, icfg, wbx, j, d, mid)){
	trace_V3_3(outf, wbx, j, d, mid, curr_tr, dolist, traceback);
	*flag = TRUE;
	return;
      }
      
      /* (V3_4) */
      else if (vx[j][d] == V3_4(s, len, rnapar, icfg, wbx, j, d, mid)){
	trace_V3_4(outf, wbx, j, d, mid, curr_tr, dolist, traceback);
	*flag = TRUE;
	return;
      }
      
      /*  (V4)__/ 1VX 1 WBX 
       */ 
      
      /* (V4_1) */
      else if (vx[j][d] == V4_1(s, len, rnapar, icfg, wbx, vx, j, d, mid)){
	trace_V4_1(outf, wbx, vx, j, d, mid, curr_tr, dolist, traceback);
	*flag = TRUE;
	return;
      }
      
      /* (V4_2) */
      else if (vx[j][d] == V4_2(s, len, rnapar, icfg, wbx, vx, j, d, mid)){
	trace_V4_2(outf, wbx, vx, j, d, mid, curr_tr, dolist, traceback);
	*flag = TRUE;
	return;
      }
      
      /* (V4_3) */
      else if (vx[j][d] == V4_3(s, len, rnapar, icfg, wbx, vx, j, d, mid)){
	trace_V4_3(outf, wbx, vx, j, d, mid, curr_tr, dolist, traceback);
	*flag = TRUE;
	return;
      }
      
      /* (V4_4) */
      else if (vx[j][d] == V4_4(s, len, rnapar, icfg, wbx, vx, j, d, mid)){
	trace_V4_4(outf, wbx, vx, j, d, mid, curr_tr, dolist, traceback);
	*flag = TRUE;
	return;
      }
      
      /* (V4_5) */
      else if (vx[j][d] == V4_5(s, len, rnapar, icfg, wbx, vx, j, d, mid)){
	trace_V4_5(outf, wbx, vx, j, d, mid, curr_tr, dolist, traceback);
	*flag = TRUE;
	return;
      }
      
      /* (V4_6) */
      else if (vx[j][d] == V4_6(s, len, rnapar, icfg, wbx, vx, j, d, mid)){
	trace_V4_6(outf, wbx, vx, j, d, mid, curr_tr, dolist, traceback);
	*flag = TRUE;
	return;
      }
      
      /* (V4_7) */
      else if (vx[j][d] == V4_7(s, len, rnapar, icfg, wbx, vx, j, d, mid)){
	trace_V4_7(outf, wbx, vx, j, d, mid, curr_tr, dolist, traceback);
	*flag = TRUE;
	return;
      }
      
      /* (V4_8) */
      else if (vx[j][d] == V4_8(s, len, rnapar, icfg, wbx, vx, j, d, mid)){
	trace_V4_8(outf, wbx, vx, j, d, mid, curr_tr, dolist, traceback);
	*flag = TRUE;
	return;
      }
      
      /* (V5)__/ 1 WBX 1VX 
       */ 
      
      /* (V5_1) */
      else if (vx[j][d] == V5_1(s, len, rnapar, icfg, wbx, vx, j, d, mid)){
	trace_V5_1(outf, wbx, vx, j, d, mid, curr_tr, dolist, traceback);
	*flag = TRUE;
	return;
      }
      
      
      /* (V5_2)  */
      else if (vx[j][d] == V5_2(s, len, rnapar, icfg, wbx, vx, j, d, mid)){
	trace_V5_2(outf, wbx, vx, j, d, mid, curr_tr, dolist, traceback);
	*flag = TRUE;
	return;
      }
      
      /* (V5_3)  */
      else if (vx[j][d] == V5_3(s, len, rnapar, icfg, wbx, vx, j, d, mid)){
	trace_V5_3(outf, wbx, vx, j, d, mid, curr_tr, dolist, traceback);
	*flag = TRUE;
	return;
      }
      
      /* (V5_4) */
      else if (vx[j][d] == V5_4(s, len, rnapar, icfg, wbx, vx, j, d, mid)){
	trace_V5_4(outf, wbx, vx, j, d, mid, curr_tr, dolist, traceback);
	*flag = TRUE;
	return;
      }
      
      /* (V5_5) */
      else if (vx[j][d] == V5_5(s, len, rnapar, icfg, wbx, vx, j, d, mid)){
	trace_V5_5(outf, wbx, vx, j, d, mid, curr_tr, dolist, traceback);
	*flag = TRUE;
	return;
      }
      
      
      /* (V5_6)  */
      else if (vx[j][d] == V5_6(s, len, rnapar, icfg, wbx, vx, j, d, mid)){
	trace_V5_6(outf, wbx, vx, j, d, mid, curr_tr, dolist, traceback);
	*flag = TRUE;
	return;
      }
      
      /* (V5_7)  */
      else if (vx[j][d] == V5_7(s, len, rnapar, icfg, wbx, vx, j, d, mid)){
	trace_V5_7(outf, wbx, vx, j, d, mid, curr_tr, dolist, traceback);
	*flag = TRUE;
	return;
      }
      
      /* (V5_8) */
      else if (vx[j][d] == V5_8(s, len, rnapar, icfg, wbx, vx, j, d, mid)){
	trace_V5_8(outf, wbx, vx, j, d, mid, curr_tr, dolist, traceback);
	*flag = TRUE;
	return;
      }
      
      /* (V6)__/  2 VX 
       */
      for (mid1 = 1; mid1 < d-mid-2; mid1++)
	for (mid2 = 1; mid2 < d-mid-mid1-1; mid2++) {
	  
	  /* (V6_1) */
	  if (vx[j][d] == V6_1(s, len, rnapar, icfg, vx, j, d, mid, mid1, mid2)){
	    trace_V6_CC(outf, vx, j, d, mid, mid1, mid2, curr_tr, dolist, traceback);
	    *flag = TRUE;
	    return;
	  }
	  
	  /* (V6_2a) */
	  else if (vx[j][d] == V6_2a(s, len, rnapar, icfg, vx, j, d, mid, mid1, mid2)){
	    trace_V6_CC(outf, vx, j, d, mid, mid1, mid2, curr_tr, dolist, traceback);
	    *flag = TRUE;
	    return;
	  }
	  
	  /* (V6_2b) */
	  else if (vx[j][d] == V6_2b(s, len, rnapar, icfg, vx, j, d, mid, mid1, mid2)){
	    trace_V6_CC(outf, vx, j, d, mid, mid1, mid2, curr_tr, dolist, traceback);
	    *flag = TRUE;
	    return;
	  }
	  
	  /* (V6_2c) */
	  else if (vx[j][d] == V6_2c(s, len, rnapar, icfg, vx, j, d, mid, mid1, mid2)){
	    trace_V6_CC(outf, vx, j, d, mid, mid1, mid2, curr_tr, dolist, traceback);
	    *flag = TRUE;
	    return;
	  }
	  
	  /* (V6_3a) */
	  else if (vx[j][d] == V6_3a(s, len, rnapar, icfg, vx, j, d, mid, mid1, mid2)){
	    trace_V6_CC(outf, vx, j, d, mid, mid1, mid2, curr_tr, dolist, traceback);
	    *flag = TRUE;
	    return;
	  }
	  
	  /* (V6_3b) */
	  else if (vx[j][d] == V6_3b(s, len, rnapar, icfg, vx, j, d, mid, mid1, mid2)){
	    trace_V6_CC(outf, vx, j, d, mid, mid1, mid2, curr_tr, dolist, traceback);
	    *flag = TRUE;
	    return;
	  }
	  
	  /* (V6_3c) */
	  else if (vx[j][d] == V6_3c(s, len, rnapar, icfg, vx, j, d, mid, mid1, mid2)){
	    trace_V6_CC(outf, vx, j, d, mid, mid1, mid2, curr_tr, dolist, traceback);
	    *flag = TRUE;
	    return;
	  }
	  
	  /* (V6_4a) */
 	  else if (vx[j][d] == V6_4a(s, len, rnapar, icfg, vx, j, d, mid, mid1, mid2)){
	    trace_V6_CC(outf, vx, j, d, mid, mid1, mid2, curr_tr, dolist, traceback);
	    *flag = TRUE;
	    return;
	  }
	  
	  /* (V6_4b) */
 	  else if (vx[j][d] == V6_4b(s, len, rnapar, icfg, vx, j, d, mid, mid1, mid2)){
	    trace_V6_CC(outf, vx, j, d, mid, mid1, mid2, curr_tr, dolist, traceback);
	    *flag = TRUE;
	    return;
	  }
	  
	  /* (V6_4c) */
	  else if (vx[j][d] == V6_4c(s, len, rnapar, icfg, vx, j, d, mid, mid1, mid2)){
	    trace_V6_CC(outf, vx, j, d, mid, mid1, mid2, curr_tr, dolist, traceback);
	    *flag = TRUE;
	    return;
	  }
	  
	  /* (V6_4d) */
	  else if (vx[j][d] == V6_4d(s, len, rnapar, icfg, vx, j, d, mid, mid1, mid2)){
	    trace_V6_CC(outf, vx, j, d, mid, mid1, mid2, curr_tr, dolist, traceback);
	    *flag = TRUE;
	    return;
	  }
	  
	  /* (V6_4e) */
 	  else if (vx[j][d] == V6_4e(s, len, rnapar, icfg, vx, j, d, mid, mid1, mid2)){
	    trace_V6_CC(outf, vx, j, d, mid, mid1, mid2, curr_tr, dolist, traceback);
	    *flag = TRUE;
	    return;
	  }
	  
	  /* (V6_4f) */
	  else if (vx[j][d] == V6_4f(s, len, rnapar, icfg, vx, j, d, mid, mid1, mid2)){
	    trace_V6_CC(outf, vx, j, d, mid, mid1, mid2, curr_tr, dolist, traceback);
	    *flag = TRUE;
	    return;
	  }
	  
	  /* (V6_4g) */
 	  else if (vx[j][d] == V6_4g(s, len, rnapar, icfg, vx, j, d, mid, mid1, mid2)){
	    trace_V6_CC(outf, vx, j, d, mid, mid1, mid2, curr_tr, dolist, traceback);
	    *flag = TRUE;
	    return;
	  }
	  
	  /* (V6_4h) */
 	  else if (vx[j][d] == V6_4h(s, len, rnapar, icfg, vx, j, d, mid, mid1, mid2)){
	    trace_V6_CC(outf, vx, j, d, mid, mid1, mid2, curr_tr, dolist, traceback);
	    *flag = TRUE;
	    return;
	  }
	  
	  /* (V6_4i) */
 	  else if (vx[j][d] == V6_4i(s, len, rnapar, icfg, vx, j, d,  mid, mid1, mid2)){
	    trace_V6_CC(outf, vx, j, d, mid, mid1, mid2, curr_tr, dolist, traceback);
	    *flag = TRUE;
	    return;
	  }
	  
	  /* (V6_5) */
	  if (vx[j][d] == V6_5(s, len, rnapar, icfg, vx, j, d, mid, mid1, mid2)){
	    trace_V6_NCC(outf, vx, j, d, mid, mid1, mid2, curr_tr, dolist, traceback);
	    *flag = TRUE;
	    return;
	  }
	  
	  /* (V6_6a) */
	  else if (vx[j][d] == V6_6a(s, len, rnapar, icfg, vx, j, d, mid, mid1, mid2)){
	    trace_V6_NCC(outf, vx, j, d, mid, mid1, mid2, curr_tr, dolist, traceback);
	    *flag = TRUE;
	    return;
	  }
	  
	  /* (V6_6b) */
	  else if (vx[j][d] == V6_6b(s, len, rnapar, icfg, vx, j, d, mid, mid1, mid2)){
	    trace_V6_NCC(outf, vx, j, d, mid, mid1, mid2, curr_tr, dolist, traceback);
	    *flag = TRUE;
	    return;
	  }
	  
	  /* (V6_6c) */
	  else if (vx[j][d] == V6_6c(s, len, rnapar, icfg, vx, j, d, mid, mid1, mid2)){
	    trace_V6_NCC(outf, vx, j, d, mid, mid1, mid2, curr_tr, dolist, traceback);
	    *flag = TRUE;
	    return;
	  }
	  
	  /* (V6_7a) */
	  else if (vx[j][d] == V6_7a(s, len, rnapar, icfg, vx, j, d, mid, mid1, mid2)){
	    trace_V6_NCC(outf, vx, j, d, mid, mid1, mid2, curr_tr, dolist, traceback);
	    *flag = TRUE;
	    return;
	  }
	  
	  /* (V6_7b) */
	  else if (vx[j][d] == V6_7b(s, len, rnapar, icfg, vx, j, d, mid, mid1, mid2)){
	    trace_V6_NCC(outf, vx, j, d, mid, mid1, mid2, curr_tr, dolist, traceback);
	    *flag = TRUE;
	    return;
	  }
	  
	  /* (V6_7c) */
	  else if (vx[j][d] == V6_7c(s, len, rnapar, icfg, vx, j, d, mid, mid1, mid2)){
	    trace_V6_NCC(outf, vx, j, d, mid, mid1, mid2, curr_tr, dolist, traceback);
	    *flag = TRUE;
	    return;
	  }
	  
	  /* (V6_8a) */
 	  else if (vx[j][d] == V6_8a(s, len, rnapar, icfg, vx, j, d, mid, mid1, mid2)){
	    trace_V6_NCC(outf, vx, j, d, mid, mid1, mid2, curr_tr, dolist, traceback);
	    *flag = TRUE;
	    return;
	  }
	  
	  /* (V6_8b) */
 	  else if (vx[j][d] == V6_8b(s, len, rnapar, icfg, vx, j, d, mid, mid1, mid2)){
	    trace_V6_NCC(outf, vx, j, d, mid, mid1, mid2, curr_tr, dolist, traceback);
	    *flag = TRUE;
	    return;
	  }
	  
	  /* (V6_8c) */
	  else if (vx[j][d] == V6_8c(s, len, rnapar, icfg, vx, j, d, mid, mid1, mid2)){
	    trace_V6_NCC(outf, vx, j, d, mid, mid1, mid2, curr_tr, dolist, traceback);
	    *flag = TRUE;
	    return;
	  }
	  
	  /* (V6_8d) */
	  else if (vx[j][d] == V6_8d(s, len, rnapar, icfg, vx, j, d, mid, mid1, mid2)){
	    trace_V6_NCC(outf, vx, j, d, mid, mid1, mid2, curr_tr, dolist, traceback);
	    *flag = TRUE;
	    return;
	  }
	  
	  /* (V6_8e) */
 	  else if (vx[j][d] == V6_8e(s, len, rnapar, icfg, vx, j, d, mid, mid1, mid2)){
	    trace_V6_NCC(outf, vx, j, d, mid, mid1, mid2, curr_tr, dolist, traceback);
	    *flag = TRUE;
	    return;
	  }
	  
	  /* (V6_8f) */
	  else if (vx[j][d] == V6_8f(s, len, rnapar, icfg, vx, j, d, mid, mid1, mid2)){
	    trace_V6_NCC(outf, vx, j, d, mid, mid1, mid2, curr_tr, dolist, traceback);
	    *flag = TRUE;
	    return;
	  }
	  
	  /* (V6_8g) */
 	  else if (vx[j][d] == V6_8g(s, len, rnapar, icfg, vx, j, d, mid, mid1, mid2)){
	    trace_V6_NCC(outf, vx, j, d, mid, mid1, mid2, curr_tr, dolist, traceback);
	    *flag = TRUE;
	    return;
	  }
	  
	  /* (V6_8h) */
 	  else if (vx[j][d] == V6_8h(s, len, rnapar, icfg, vx, j, d, mid, mid1, mid2)){
	    trace_V6_NCC(outf, vx, j, d, mid, mid1, mid2, curr_tr, dolist, traceback);
	    *flag = TRUE;
	    return;
	  }
	  
	  /* (V6_8i) */
 	  else if (vx[j][d] == V6_8i(s, len, rnapar, icfg, vx, j, d,  mid, mid1, mid2)){
	    trace_V6_NCC(outf, vx, j, d, mid, mid1, mid2, curr_tr, dolist, traceback);
	    *flag = TRUE;
	    return;
	  }
	} /* mid1,mid2 loops */
    } /* mid loop */
}

/* Function: FillVX_pks()
 * 
 * Purpose: fill internal pseudoknots in vx_ Using Knots-IS(2) method.
 *
 * 2 ONE-HOLE STRUCTURES (4 types of diagrams V7-V10)
 *
 * (add P13 penalty for using a hole structure)
 * (weigth contributions with wkn instead of wsf as before)
 * (parameters P11, P12, P13 have no tuner counterparts)
 *
 *
 * one connects (i+1, i+1+mid1) to (i+1+mid-mid2, i+1+mid) 
 *              (i+1+mid, mid, mid1, mid2)
 * one connects (i+1+mid1+1, i+mid-mid2) to (i+2+mid,j-1)
 *              (j-1, d-mid1-3, mid-mid1-mid2-2, d-mid-3)
 *           
 * Arguments: s     - sequence (0..len-1) converted to integers
 *                    representing array indices of symbols
 *            len   - length of iseq
 *           icfg   - context-free grammar state transitions, integer log form
 *            whx   - DP matrix, already alloc'ed 
 *            vhx   - DP matrix, already alloc'ed 
 *            zhx   - DP matrix, already alloc'ed 
 *            yhx   - DP matrix, already alloc'ed 
 *            wbx   - DP matrix, already alloc'ed 
 *             vx   - DP matrix, already alloc'ed 
 *            
 * Return:    (void)           
 *            structures V7, V8, V9, V10 of vx are filled.
 */     
void 
FillVX_pks(ESL_DSQ *s, int len, struct rnapar_2 *rnapar,
	   int **icfg, int **wbx, int **vx, 
	   int ****whx, int ****zhx, int ****yhx, int j, int d)
{
  int               i;
  int              sc;
  int          bestsc;
  int mid, mid1, mid2;
  
  i = j - d;
  
  if (icfg[idxS][idxP(s[i],s[j])] == 1*INTSCALE){ 

   bestsc = vx[j][d];

    for (mid = 2; mid < d-4; mid++)
      for (mid1 = 0; mid1 < mid-2; mid1++)
	for (mid2 = 0; mid2 < mid-mid1-2; mid2++)
	  {
	    /* (V7)__/ 2 WHX
	     */
	    
	    /* (V7.1) */
	    if ((sc = rnapar->P10P + rnapar->P5P + rnapar->P13 + 3*rnapar->P6P 
		 + whx[i+1+mid][mid][mid1][mid2] 
		 + whx[j-1][d-mid1-4][mid-mid1-mid2-3][d-mid-5]) > bestsc)
	      bestsc = sc;
	    
	    /* (V7.2) */
	    if (mid1 > 0 &&
		(sc = rnapar->P10P + rnapar->P5P + rnapar->P13 + 4*rnapar->P6P 
		 + wkn*icfg[idxR(s[i+1])][idxP(s[j],s[i])]
		 + whx[i+1+mid][mid-1][mid1-1][mid2] 
		 + whx[j-1][d-mid1-4][mid-mid1-mid2-3][d-mid-5]) > bestsc)
	      bestsc = sc;
	    
	    /* (V7.3)  */
	    if (d-mid > 5 &&
		(sc = rnapar->P10P + rnapar->P5P + rnapar->P13 + 4*rnapar->P6P 
		 + wkn*icfg[idxL(s[j-1])][idxP(s[j],s[i])]
		 + whx[i+1+mid][mid][mid1][mid2] 
		 + whx[j-2][d-mid1-5][mid-mid1-mid2-3][d-mid-6]) > bestsc)
	      bestsc = sc;
	    
	    /* (V7.4) */
	    if (mid1 > 0 && d-mid > 5 &&
		(sc = rnapar->P10P + rnapar->P5P + rnapar->P13 + 5*rnapar->P6P 
		 + wkn*icfg[idxR(s[i+1])][idxP(s[j],s[i])]
		 + wkn*icfg[idxL(s[j-1])][idxP(s[j],s[i])]
		 + whx[i+1+mid][mid-1][mid1-1][mid2] 
		 + whx[j-2][d-mid1-5][mid-mid1-mid2-3][d-mid-6]) > bestsc)
	      bestsc = sc;
	    
	    /* (V8)__/ 1 ZHX 1 WHX
	     */
	    
	    /* (V8.1) */
	    if ((sc = 2*rnapar->P10P + rnapar->P5P + rnapar->P13 + 3*rnapar->P6P 
		 + wkn*icfg[idxPS(s[i],s[j])][idxPS(s[i+1],s[i+mid+1])] 
		 + zhx[i+1+mid][mid][mid1][mid2] 
		 + whx[j-1][d-mid1-4][mid-mid1-mid2-3][d-mid-5]) > bestsc)
	      bestsc = sc;
	    
	    /* (V8.2)  */
	    if (d-mid > 5 &&
		(sc = 2*rnapar->P10P + rnapar->P5P + rnapar->P13 + 4*rnapar->P6P 
		 + wkn*icfg[idxL(s[j-1])][idxP(s[j],s[i])]
		 + wkn*icfg[idxPS(s[i],s[j])][idxPS(s[i+1],s[i+mid+1])]
		 + zhx[i+1+mid][mid][mid1][mid2] 
		 + whx[j-2][d-mid1-5][mid-mid1-mid2-3][d-mid-6]) > bestsc)
	      bestsc = sc;
	    
	    /* (V8.3) */
	    if (d-mid > 5 &&
		(sc = 2*rnapar->P10P + rnapar->P5P + rnapar->P13 + 4*rnapar->P6P 
		 + wkn*icfg[idxPS(s[i],s[j])][idxPS(s[i+2],s[i+mid+2])] 
		 + zhx[i+2+mid][mid][mid1][mid2] 
		 + whx[j-1][d-mid1-5][mid-mid1-mid2-3][d-mid-6]) > bestsc)
	      bestsc = sc;
	    
	    /* (V8.4)  */
	    if (d-mid > 6 &&
		(sc = 2*rnapar->P10P + rnapar->P5P + rnapar->P13 + 5*rnapar->P6P 
		 + wkn*icfg[idxL(s[j-1])][idxP(s[j],s[i])]
		 + wkn*icfg[idxPS(s[i],s[j])][idxPS(s[i+2],s[i+mid+2])]
		 + zhx[i+2+mid][mid][mid1][mid2] 
		 + whx[j-2][d-mid1-6][mid-mid1-mid2-3][d-mid-7]) > bestsc)
	      bestsc = sc;
	    
	    /* (V9)__/ 1 WHX 1 ZHX 
	     */
	    
	    /* (V9.1) */
	    if ((sc = 2*rnapar->P10P + rnapar->P5P + rnapar->P13 + 3*rnapar->P6P
		 + wkn*icfg[idxPS(s[j-1],s[i+3+mid1])][idxPS(s[j],s[i])]
		 + whx[i+1+mid][mid][mid1][mid2] 
		 + zhx[j-1][d-mid1-4][mid-mid1-mid2-3][d-mid-5]) > bestsc)
	      bestsc = sc;
	    
	    /* (V9.2) */
	    if (mid1 > 0 &&
		(sc = 2*rnapar->P10P + rnapar->P5P + rnapar->P13 + 4*rnapar->P6P 
		 + wkn*icfg[idxR(s[i+1])][idxP(s[j],s[i])]
		 + wkn*icfg[idxPS(s[j-1],s[i+3+mid1])][idxPS(s[j],s[i])]
		 + whx[i+1+mid][mid-1][mid1-1][mid2] 
		 + zhx[j-1][d-mid1-4][mid-mid1-mid2-3][d-mid-5] ) > bestsc)
	      bestsc = sc;

	    /* (V9.3) */
	    if (d-mid > 5 &&
		(sc = 2*rnapar->P10P + rnapar->P5P + rnapar->P13 + 4*rnapar->P6P
		 + wkn*icfg[idxPS(s[j-2],s[i+3+mid1])][idxPS(s[j],s[i])]
		 + whx[i+1+mid][mid][mid1][mid2] 
		 + zhx[j-2][d-mid1-5][mid-mid1-mid2-3][d-mid-6]) > bestsc)
	      bestsc = sc;
	    
	    /* (V9.4) */
	    if (mid1 > 0 && d-mid > 5 &&
		(sc = 2*rnapar->P10P + rnapar->P5P + rnapar->P13 + 5*rnapar->P6P 
		 + wkn*icfg[idxR(s[i+1])][idxP(s[j],s[i])]
		 + wkn*icfg[idxPS(s[j-2],s[i+3+mid1])][idxPS(s[j],s[i])]
		 + whx[i+1+mid][mid-1][mid1-1][mid2] 
		 + zhx[j-2][d-mid1-5][mid-mid1-mid2-3][d-mid-6]) > bestsc)
	      bestsc = sc;
	    
	    /* (V10)__/ 2 YHX 
	     */
	    
	    /* (V10.1) */
	    if ((sc = 3*rnapar->P10P + rnapar->P5P + rnapar->P13 + 3*rnapar->P6P
		 + wkn*icfg[idxR(s[i+mid1+2])][idxP(s[i+1+mid-mid2],s[i+1+mid1])]
		 + wkn*icfg[idxL(s[i+mid+3])][idxP(s[i+mid+4],s[i+mid-mid2])]
		 + wkn*icfg[idxPS(s[i+mid-mid2],s[i+mid+4])]
		           [idxPS(s[i+1+mid-mid2],s[i+1+mid1])] 
		 + yhx[i+1+mid][mid][mid1][mid2]
		 + yhx[j-1][d-mid1-4][mid-mid1-mid2-3][d-mid-5]) > bestsc)
	      bestsc = sc;
	    
	    /* (V10.2) */
	    if (mid > 2 && mid1 > 0 && 
		(sc = 3*rnapar->P10P + rnapar->P5P + rnapar->P13 + 4*rnapar->P6P 
		 + wkn*icfg[idxR(s[i+1])][idxP(s[j],s[i])]
		 + wkn*icfg[idxR(s[i+mid1+2])][idxP(s[i+1+mid-mid2],s[i+1+mid1])]
		 + wkn*icfg[idxL(s[i+mid+3])][idxP(s[i+mid+4],s[i+mid-mid2])]
		 + wkn*icfg[idxPS(s[i+mid-mid2],s[i+mid+4])]
		           [idxPS(s[i+1+mid-mid2],s[i+1+mid1])]
		 + yhx[i+1+mid][mid-1][mid1-1][mid2] 
		 + yhx[j-1][d-mid1-4][mid-mid1-mid2-3][d-mid-5]) > bestsc)
	      bestsc = sc;
	    
	    /* (V10.3) */
	    if (mid > 1 && d-mid > 5 &&
		(sc = 3*rnapar->P10P + rnapar->P5P + rnapar->P13 + 4*rnapar->P6P 
		 + wkn*icfg[idxL(s[j-1])][idxP(s[j],s[i])]
		 + wkn*icfg[idxR(s[i+mid1+2])][idxP(s[i+1+mid-mid2],s[i+1+mid1])]
		 + wkn*icfg[idxL(s[i+mid+3])][idxP(s[i+mid+4],s[i+mid-mid2])]
		 + wkn*icfg[idxPS(s[i+mid-mid2],s[i+mid+4])]
		           [idxPS(s[i+1+mid-mid2],s[i+1+mid1])] 
		 + yhx[i+1+mid][mid][mid1][mid2] 
		 + yhx[j-2][d-mid1-5][mid-mid1-mid2-3][d-mid-6]) > bestsc)
	      bestsc = sc;
	    
	    /* (V10.4) */
	    if (mid > 2 && mid1 > 0 && d-mid > 5 &&
		(sc = 3*rnapar->P10P + rnapar->P5P + rnapar->P13 + 5*rnapar->P6P 
		 + wkn*icfg[idxR(s[i+1])][idxP(s[j],s[i])]
		 + wkn*icfg[idxL(s[j-1])][idxP(s[j],s[i])]
		 + wkn*icfg[idxR(s[i+mid1+2])][idxP(s[i+1+mid-mid2],s[i+1+mid1])]
		 + wkn*icfg[idxL(s[i+mid+3])][idxP(s[i+mid+4],s[i+mid-mid2])]
		 + wkn*icfg[idxPS(s[i+mid-mid2],s[i+mid+4])]
		           [idxPS(s[i+1+mid-mid2],s[i+1+mid1])]
		 + yhx[i+1+mid][mid-1][mid1-1][mid2]
		 + yhx[j-2][d-mid1-5][mid-mid1-mid2-3][d-mid-6]) > bestsc)
	      bestsc = sc;

	    /* (V10.5) */
	    if (mid-mid1-mid2 > 3 &&
		(sc = 3*rnapar->P10P + rnapar->P5P + rnapar->P13 + 4*rnapar->P6P
		 + wkn*icfg[idxR(s[i+mid1+2])][idxP(s[i+1+mid-mid2],s[i+1+mid1])]
		 + wkn*icfg[idxL(s[i+mid+3])][idxP(s[i+mid+4],s[i+mid-mid2-1])]
		 + wkn*icfg[idxPS(s[i+mid-mid2-1],s[i+mid+4])]
		           [idxPS(s[i+1+mid-mid2],s[i+1+mid1])] 
		 + yhx[i+1+mid][mid][mid1][mid2]
		 + yhx[j-1][d-mid1-4][mid-mid1-mid2-4][d-mid-5] ) > bestsc)
	      bestsc = sc;
	    
	    /* (V10.6) */
	    if (mid > 2 && mid1 > 0 && mid-mid1-mid2 > 3 &&
		(sc = 3*rnapar->P10P + rnapar->P5P + rnapar->P13 + 5*rnapar->P6P 
		 + wkn*icfg[idxR(s[i+1])][idxP(s[j],s[i])]
		 + wkn*icfg[idxR(s[i+mid1+2])][idxP(s[i+1+mid-mid2],s[i+1+mid1])]
		 + wkn*icfg[idxL(s[i+mid+3])][idxP(s[i+mid+4],s[i+mid-mid2-1])]
		 + wkn*icfg[idxPS(s[i+mid-mid2-1],s[i+mid+4])]
		           [idxPS(s[i+1+mid-mid2],s[i+1+mid1])]
		 + yhx[i+1+mid][mid-1][mid1-1][mid2] 
		 + yhx[j-1][d-mid1-4][mid-mid1-mid2-4][d-mid-5]) > bestsc)
	      bestsc = sc;
	    
	    /* (V10.7) */
	    if (mid > 1 && d-mid > 5 && mid-mid1-mid2 > 3 &&
		(sc = 3*rnapar->P10P + rnapar->P5P + rnapar->P13 + 5*rnapar->P6P 
		 + wkn*icfg[idxL(s[j-1])][idxP(s[j],s[i])]
		 + wkn*icfg[idxR(s[i+mid1+2])][idxP(s[i+1+mid-mid2],s[i+1+mid1])]
		 + wkn*icfg[idxL(s[i+mid+3])][idxP(s[i+mid+4],s[i+mid-mid2-1])]
		 + wkn*icfg[idxPS(s[i+mid-mid2-1],s[i+mid+4])]
		           [idxPS(s[i+1+mid-mid2],s[i+1+mid1])] 
		 + yhx[i+1+mid][mid][mid1][mid2] 
		 + yhx[j-2][d-mid1-5][mid-mid1-mid2-4][d-mid-6]) > bestsc)
	      bestsc = sc;
	    
	    /* (V10.8) */
	    if (mid > 2 && mid1 > 0 && d-mid > 5 && mid-mid1-mid2 > 3 &&
		(sc = 3*rnapar->P10P + rnapar->P5P + rnapar->P13 + 6*rnapar->P6P 
		 + wkn*icfg[idxR(s[i+1])][idxP(s[j],s[i])]
		 + wkn*icfg[idxL(s[j-1])][idxP(s[j],s[i])]
		 + wkn*icfg[idxR(s[i+mid1+2])][idxP(s[i+1+mid-mid2],s[i+1+mid1])]
		 + wkn*icfg[idxL(s[i+mid+3])][idxP(s[i+mid+4],s[i+mid-mid2-1])]
		 + wkn*icfg[idxPS(s[i+mid-mid2-1],s[i+mid+4])]
		           [idxPS(s[i+1+mid-mid2],s[i+1+mid1])]
		 + yhx[i+1+mid][mid-1][mid1-1][mid2]
		 + yhx[j-2][d-mid1-5][mid-mid1-mid2-4][d-mid-6]) > bestsc)
	      bestsc = sc;
	  }
    vx[j][d] = bestsc;
 }  /* if (i,j) are base-paired, else prohibited  */
  
  else vx[j][d] = -BIGINT;
}

/* Function: TraceVX_pks()
 * 
 * Purpose: traceback internal pseudoknots in vx. Using Knots-IS(2) method.
 *
 * 2 ONE-HOLE STRUCTURES (4 types of diagrams V7-V10)
 *
 * (add P13 penalty for using a hole structure)
 * (weigth contributions with wkn instead of wsf as before)
 * (parameters P11, P12, P13 have no tuner counterparts)
 *
 *
 * one connects (i+1, i+1+mid1) to (i+1+mid-mid2, i+1+mid) 
 *              (i+1+mid, mid, mid1, mid2)
 * one connects (i+1+mid1+1, i+mid-mid2) to (i+2+mid,j-1)
 *              (j-1, d-mid1-3, mid-mid1-mid2-2, d-mid-3)
 *           
 * Arguments: s     - sequence (0..len-1) converted to integers
 *                    representing array indices of symbols
 *            len   - length of iseq
 *           icfg   - context-free grammar state transitions, integer log form
 *            whx   - DP matrix, already alloc'ed 
 *            vhx   - DP matrix, already alloc'ed 
 *            zhx   - DP matrix, already alloc'ed 
 *            yhx   - DP matrix, already alloc'ed 
 *            wbx   - DP matrix, already alloc'ed 
 *             vx   - DP matrix, already alloc'ed 
 *            
 * Return:    (void)           
 *            structures V7, V8, V9, V10 of vx are traced back.
 */      
void 
TraceVX_pks(FILE *outf, ESL_DSQ *s, int len, struct rnapar_2 *rnapar,
	    int **icfg, int **wbx, int **vx, 
	    int ****whx, int ****zhx, int ****yhx, int j, int d, int *flag,
	    struct tracekn_s *curr_tr, struct traceknstack_s *dolist, int traceback)
{
  int               i;
  int mid, mid1, mid2;
  
  i = j - d;
  
    for (mid = 2; mid < d-4; mid++)
      for (mid1 = 0; mid1 < mid-2; mid1++)
	for (mid2 = 0; mid2 < mid-mid1-2; mid2++){

	  /* (V7)__/2 WHX 
	   */
	  
	  /* (V7.1) */
	  if (vx[j][d] == V7_1(s, len, rnapar, icfg, whx, j, d, mid, mid1, mid2)){
	    trace_V7_1(outf, whx, j, d, mid, mid1, mid2, curr_tr, dolist, traceback);
	    *flag = TRUE;
	    return;
	  }
	  
	  /* (V7.2) */
	  else if (vx[j][d] == V7_2(s, len, rnapar, icfg, whx, j, d, mid, mid1, mid2)){
	    trace_V7_2(outf, whx, j, d, mid, mid1, mid2, curr_tr, dolist, traceback);
	    *flag = TRUE;
	    return;
	  }
	  
	  /* (V7.3) */
	  else if (vx[j][d] == V7_3(s, len, rnapar, icfg, whx, j, d, mid, mid1, mid2)){
	    trace_V7_3(outf, whx, j, d, mid, mid1, mid2, curr_tr, dolist, traceback);
	    *flag = TRUE;
	    return;
	  }
	  
	  /* (V7.4) */
	  else if (vx[j][d] == V7_4(s, len, rnapar, icfg, whx, j, d, mid, mid1, mid2)){
	    trace_V7_4(outf, whx, j, d, mid, mid1, mid2, curr_tr, dolist, traceback);
	    *flag = TRUE;
	    return;
	  }
	  
	  /* (V8)__/  1 ZHX 1 WHX
	   */
	  
	  /* (V8.1) */
	  else if (vx[j][d] == V8_1(s, len, rnapar, icfg, whx, zhx,  j, d, mid, mid1, mid2)){
	    trace_V8_1(outf, whx, zhx, j, d, mid, mid1, mid2, curr_tr, dolist, traceback);
	    *flag = TRUE;
	    return;
	  }
	  
	  /* (V8.2) */
	  else if (vx[j][d] == V8_2(s, len, rnapar, icfg, whx, zhx, j, d, mid, mid1, mid2)){
	    trace_V8_2(outf, whx, zhx, j, d, mid, mid1, mid2, curr_tr, dolist, traceback);
	    *flag = TRUE;
	    return;
	  }

	  /* (V8.3) */
	  else if (vx[j][d] == V8_3(s, len, rnapar, icfg, whx, zhx,  j, d, mid, mid1, mid2)){
	    trace_V8_3(outf, whx, zhx, j, d, mid, mid1, mid2, curr_tr, dolist, traceback);
	    *flag = TRUE;
	    return;
	  }
	  
	  /* (V8.4) */
	  else if (vx[j][d] == V8_4(s, len, rnapar, icfg, whx, zhx, j, d, mid, mid1, mid2)){
	    trace_V8_4(outf, whx, zhx, j, d, mid, mid1, mid2, curr_tr, dolist, traceback);
	    *flag = TRUE;
	    return;
	  }
	  
	  /* (V9)__/ 1 WHX 1 ZHX 
	   */
	  
	  /* (V9.1) */
	  else if (vx[j][d] == V9_1(s, len, rnapar, icfg, whx, zhx, j, d, mid, mid1, mid2)){
	    trace_V9_1(outf, whx, zhx, j, d, mid, mid1, mid2, curr_tr, dolist, traceback);
	    *flag = TRUE;
	    return;
	  }
	  
	  /* (V9.2) */
	  else if (vx[j][d] == V9_2(s, len, rnapar, icfg, whx, zhx, j, d, mid, mid1, mid2)){
	    trace_V9_2(outf, whx, zhx, j, d, mid, mid1, mid2, curr_tr, dolist, traceback);
	    *flag = TRUE;
	    return;
	  }

	  /* (V9.3) */
	  else if (vx[j][d] == V9_3(s, len, rnapar, icfg, whx, zhx, j, d, mid, mid1, mid2)){
	    trace_V9_3(outf, whx, zhx, j, d, mid, mid1, mid2, curr_tr, dolist, traceback);
	    *flag = TRUE;
	    return;
	  }
	  
	  /* (V9.4) */
	  else if (vx[j][d] == V9_4(s, len, rnapar, icfg, whx, zhx, j, d, mid, mid1, mid2)){
	    trace_V9_4(outf, whx, zhx, j, d, mid, mid1, mid2, curr_tr, dolist, traceback);
	    *flag = TRUE;
	    return;
	  }
	  
	  /* (V10)__/ 2 YHX 
	   */
	  
	  /* (V10.1) */
	  else if (vx[j][d] == V10_1(s, len, rnapar, icfg, yhx, j, d, mid, mid1, mid2)){
	    trace_V10_1(outf, yhx, j, d, mid, mid1, mid2, curr_tr, dolist, traceback);
	    *flag = TRUE;
	    return;
	  }
	  
	  /* (V10.2) */
	  else if (vx[j][d] ==  V10_2(s, len, rnapar, icfg, yhx, j, d, mid, mid1, mid2)){
	    trace_V10_2(outf, yhx, j, d, mid, mid1, mid2, curr_tr, dolist, traceback);
	    *flag = TRUE;
	    return;
	  }
	  
	  /* (V10.3) */
	  else if (vx[j][d] == V10_3(s, len, rnapar, icfg, yhx, j, d, mid, mid1, mid2)){
	    trace_V10_3(outf, yhx, j, d, mid, mid1, mid2, curr_tr, dolist, traceback);
	    *flag = TRUE;
	    return;
	  }
	  
	  /* (V10.4) */
	  else if (vx[j][d] ==  V10_4(s, len, rnapar, icfg, yhx, j, d, mid, mid1, mid2)){
	    trace_V10_4(outf, yhx, j, d, mid, mid1, mid2, curr_tr, dolist, traceback);
	    *flag = TRUE;
	    return;
	  }

	  /* (V10.5) */
	  else if (vx[j][d] == V10_5(s, len, rnapar, icfg, yhx, j, d, mid, mid1, mid2)){
	    trace_V10_5(outf, yhx, j, d, mid, mid1, mid2, curr_tr, dolist, traceback);
	    *flag = TRUE;
	    return;
	  }
	  
	  /* (V10.6) */
	  else if (vx[j][d] ==  V10_6(s, len, rnapar, icfg, yhx, j, d, mid, mid1, mid2)){
	    trace_V10_6(outf, yhx, j, d, mid, mid1, mid2, curr_tr, dolist, traceback);
	    *flag = TRUE;
	    return;
	  }
	  
	  /* (V10.7) */
	  else if (vx[j][d] == V10_7(s, len, rnapar, icfg, yhx, j, d, mid, mid1, mid2)){
	    trace_V10_7(outf, yhx, j, d, mid, mid1, mid2, curr_tr, dolist, traceback);
	    *flag = TRUE;
	    return;
	  }
	  
	  /* (V10.8) */
	  else if (vx[j][d] ==  V10_8(s, len, rnapar, icfg, yhx, j, d, mid, mid1, mid2)){
	    trace_V10_8(outf, yhx, j, d, mid, mid1, mid2, curr_tr, dolist, traceback);
	    *flag = TRUE;
	    return;
	  }
	} 
} 




