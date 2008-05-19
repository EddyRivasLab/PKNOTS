/* fillvxscan.c 
 *
 * includes functions FillVX 
 * fills the no-hole matix:
 *      vx[j][d] (len x len)
 *               [(j-d,j) are base pared]
 *
 */
                         
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>

#include "cfg.h"
#include "proto.h"
#include "squid.h"
                                                 
#ifdef MEMDEBUG
#include "dbmalloc.h"
#endif

/* Function: FillVPScan() 
 * 
 * Purpose: fill VP for internal loops. 
 *           
 * Arguments: s     - sequence (0..len-1) converted to integers
 *                    representing array indices of symbols
 *            len   - length of iseq
 *            cfg   - context-free grammar state transitions, float log2 form
 *             vx   - DP matrix, already alloc'ed 
 *            
 * Return:    (void)           
 *            vp is calculates
 */       
void
FillVPScan(int *s, int len, int win, int **icfg, int **vx, int *vp, int j, int jmod, int d)
{
  int mid;  /* midpoint of bifurcation         */
  int k;
  int kmod;
  int sc;

  k = j - d + 2;
  kmod = jmod - d + 2;
  
  if (d < 4) {
    vp[d] = -BIGINT;
    return;
  }
  else if (d > 4 &&
	   icfg[idxS][idxP(s[j-d],s[j])] == 1*INTSCALE && 
	   icfg[idxS][idxP(s[k],s[j-2])] == 1*INTSCALE &&
	   icfg[idxS][idxP(s[j-d+1],s[j-1])] != INTSCALE)
    vp[d-4] = wsf*EPARAM3
      + wsf*icfg[idxPS(s[j-2],s[k])][idxPS(s[j-1],s[k-1])]
      + vx[(jmod-2<0)? jmod-2+win:jmod-2][d-4];
  else 
    vp[d-4] = -BIGINT;

  /* vp[d-5] is an special case because of the change from idcPS to idxPI 
   * it cannot be read from vp[d-4]
   */
  vp[d-5] = -BIGINT;
  /* two possibilities 
   */
  if (icfg[idxS][idxP(s[k+1],s[j-2])] == 1*INTSCALE &&
      (sc =  wsf*EPARAM3
       + wsf*icfg[idxPI(s[j-2],s[k+1])][idxPI(s[j-1],s[k])]
       + vx[(jmod-2<0)? jmod-2+win:jmod-2][d-5]) > vp[d-5]) vp[d-5] = sc;
  /* or 
   */
  if (icfg[idxS][idxP(s[k],s[j-3])] == 1*INTSCALE &&
      (sc = wsf*EPARAM3
       + wsf*icfg[idxPI(s[j-3],s[k])][idxPI(s[j-2],s[k-1])]
       + vx[(jmod-3<0)? jmod-3+win:jmod-3][d-5]) > vp[d-5]) vp[d-5] = sc;
  
  /* Now the rest of mid's
   */
  for (mid = 0; mid < (d-5); mid++) 
    if (icfg[idxS][idxP(s[k],s[k+mid])] == 1*INTSCALE &&
	(sc = wsf*EPARAM3
	 + wsf*icfg[idxPI(s[k+mid],s[k])][idxPI(s[k+mid+1],s[k-1])]
	 + vx[(kmod+mid<0)? kmod+mid+win:kmod+mid][mid]) > vp[mid])
      vp[mid] = sc;
}


/* Function: FillVX_nestedScan() 
 * 
 * Purpose: fill no-hole matrix VX without pseudoknots. 
 *           
 * Arguments: s     - sequence (0..len-1) converted to integers
 *                    representing array indices of symbols
 *            len   - length of iseq
 *            cfg   - context-free grammar state transitions, float log2 form
 *            wbx   - DP matrix, already alloc'ed 
 *             vx   - DP matrix, already alloc'ed 
 *            
 * Return:    (void)           
 *            diagrams v1, v2, v3, v4, v5, v6 of vx are calculated.
 */       
void
FillVX_nestedScan(int *s, int len, int win, int **icfg,
		  int **wbx, int **vx, int *vp, int j, int jmod, int d)
{
  int                  i;  /* coords in mtx's                 */
  int                mid;  /* midpoint of bifurcation         */
  int                 sc;  /* temporary score to check        */
  int             bestsc;  /* best score so far               */

  i = j - d;

  if (icfg[idxS][idxP(s[i],s[j])] == 1*INTSCALE)
    { 
      /* (V1)__/ IRREDUCIBLE SURFACES of  O(1)
       */
      bestsc =  wsf*F1(s, len, icfg, j, d);	
      
      /* (V2)__/ IRREDUCIBLE SURFACES of  O(2)
       */
      /* STEMS
       */
      if (d > 1 &&
	  (sc = wsf*F2(s, len, icfg, j, d, 1, 1)
	   + vx[(jmod-1<0)? jmod-1+win:jmod-1][d-2]) > bestsc)
	bestsc = sc;

      /* BULGES 
       */
      for (mid = 2; mid < d; mid++) {
	if ((sc = wsf*F2(s, len, icfg, j, d, mid, 1)
	     + vx[(jmod-1<0)? jmod-1+win:jmod-1][d-1-mid]) > bestsc)
	  bestsc = sc;
	if ((sc = wsf*F2(s, len, icfg, j, d, 1, mid)
	     + vx[(jmod-mid<0)? jmod-mid+win:jmod-mid][d-1-mid]) > bestsc)
	  bestsc = sc;
      }
      
      /* INTERNAL LOOPS (N^3)
       */
      for (mid = 0; mid <= (d-4); mid++) {
	if (d-mid > 32 &&
	    (sc = wsf*icfg[idxPI(s[i],s[j])][idxPI(s[i+1],s[j-1])] 
	     + wsf*(int)(PRELOG*log((d-mid-2)/30))
	     + wsf*icfg[idxS][idxW(30)]
	     + vp[mid]) > bestsc)
	  bestsc = sc;
	else if (d-mid == 4 && 
		 icfg[idxS][idxP(s[i+1],s[j-1])] != INTSCALE && 
		 (sc = wsf*icfg[idxPS(s[i],s[j])][idxPS(s[i+1],s[j-1])] 
		  + vp[mid]) > bestsc)
	  bestsc = sc;
	else if (d-mid > 4 && d-mid < 33 &&
		 (sc = wsf*icfg[idxPI(s[i],s[j])][idxPI(s[i+1],s[j-1])] 
		  + wsf*icfg[idxS][idxW(d-mid-2)]
		  + vp[mid]) > bestsc)
	  bestsc = sc;
      }
      /* INTERNAL LOOPS (N^4)
	 for (mid1 = 2; mid1 < d; mid1++)
	 for (mid2 = 2; mid2 < (d-mid1); mid2++)
	 if ((sc = wsf*F2(s, len, icfg, j, d, mid1, mid2)
	 + vx[j-mid2][d-mid1-mid2]) > bestsc)
	 bestsc = sc;
      */      
      
      /* REST of VX (MULTILOOPS).  
       */
      
      for (mid = 1; mid < d-2; mid++)
	{                      
	  /* (V3)__/ 2 WBX
	   */
	  
	  /* (V3_1) */
	  if (mid < d-3 &&
	      (sc = P10 + P5 
	       + wbx[(jmod-d+1+mid<0)? jmod-d+1+mid+win:jmod-d+1+mid][mid] 
	       + wbx[(jmod-1<0)? jmod-1+win:jmod-1][d-mid-3]) > bestsc)

	    bestsc = sc;
	  
	  /* (V3_2) */
	  if (mid > 1 &&
	      (sc =  P10 + P5 + P6 
	       + wsf*icfg[idxR(s[i+1])][idxP(s[j],s[i])]
	       + wbx[(jmod-d+1+mid<0)? jmod-d+1+mid+win:jmod-d+1+mid][mid-1] 
	       + wbx[(jmod-1<0)? jmod-1+win:jmod-1][d-mid-3]) > bestsc)
	    bestsc = sc;
	  
	  /* (V3_3) */
	  if (mid < d-4 &&
	      (sc = P10 + P5 + P6 
	       + wsf*icfg[idxL(s[j-1])][idxP(s[j],s[i])]
	       + wbx[(jmod-d+1+mid<0)? jmod-d+1+mid+win:jmod-d+1+mid][mid] 
	       + wbx[(jmod-2<0)? jmod-2+win:jmod-2][d-mid-4]) > bestsc)
	    bestsc = sc;
	  
	  /* (V3_4) */
	  if (mid < d-4 &&
	      (sc = P10 + P5 + 2*P6 
	       + wsf*icfg[idxR(s[i+1])][idxP(s[j],s[i])] 
	       + wsf*icfg[idxL(s[j-1])][idxP(s[j],s[i])]
	       + wbx[(jmod-d+1+mid<0)? jmod-d+1+mid+win:jmod-d+1+mid][mid-1] 
	       + wbx[(jmod-2<0)? jmod-2+win:jmod-2][d-mid-4]) > bestsc)
	    bestsc = sc;
	  
	  /* (V4)__/ 1VX 1 WBX 
	   */ 

	  /* (V4_1) */
	  if (mid < d-3 &&
	      (sc = 2*P10 + P5 
	       + wsf*(icfg[idxPS(s[i],s[j])][idxPS(s[i+1],s[i+1+mid])] + INTSCALE)
	       + vx[(jmod-d+1+mid<0)? jmod-d+1+mid+win:jmod-d+1+mid][mid] 
	       + wbx[(jmod-1<0)? jmod-1+win:jmod-1][d-mid-3]) > bestsc)
	    bestsc = sc;
	  
	  /* (V4_2) */
	  if (mid < d-4 &&
	      (sc = 2*P10 + P5 + P6 
	       + wsf*icfg[idxL(s[j-1])][idxP(s[j],s[i])]
	       + wsf*icfg[idxPS(s[i],s[j])][idxPS(s[i+1],s[i+1+mid])] 
	       + vx[(jmod-d+1+mid<0)? jmod-d+1+mid+win:jmod-d+1+mid][mid] 
	       + wbx[(jmod-2<0)? jmod-2+win:jmod-2][d-mid-4]) > bestsc)
	    bestsc = sc;

	  /* (V4_3) */
	  if (mid < d-4 &&
	      (sc = 2*P10 + P5 + P6
	       + wsf*icfg[idxR(s[i+mid+2])][idxP(s[i+1],s[i+mid+1])]
	       + wsf*icfg[idxPS(s[i],s[j])][idxPS(s[i+1],s[i+1+mid])] 
	       + vx[(jmod-d+1+mid<0)? jmod-d+1+mid+win:jmod-d+1+mid][mid] 
	       + wbx[(jmod-1<0)? jmod-1+win:jmod-1][d-mid-4]) > bestsc)
	    bestsc = sc;
	  
	  /* (V4_4) */
	  if (mid < d-5 &&
	      (sc = 2*P10 + P5 + 2*P6 
	       + wsf*icfg[idxL(s[j-1])][idxP(s[j],s[i])]
	       + wsf*icfg[idxR(s[i+mid+2])][idxP(s[i+1],s[i+mid+1])]
	       + wsf*icfg[idxPS(s[i],s[j])][idxPS(s[i+1],s[i+1+mid])] 
	       + vx[(jmod-d+1+mid<0)? jmod-d+1+mid+win:jmod-d+1+mid][mid] 
	       + wbx[(jmod-2<0)? jmod-2+win:jmod-2][d-mid-5]) > bestsc)
	    bestsc = sc;
	  
	  /* (V4_5) */
	  if (mid < d-4 &&
	      (sc = 2*P10 + P5 + P6
	       + wsf*icfg[idxPS(s[i],s[j])][idxPS(s[i+2],s[i+2+mid])]
	       + vx[(jmod-d+2+mid<0)? jmod-d+2+mid+win:jmod-d+2+mid][mid] 
	       + wbx[(jmod-1<0)? jmod-1+win:jmod-1][d-mid-4]) > bestsc)
	    bestsc = sc;
	  
	  /* (V4_6) */
	  if (mid < d-5 &&
	      (sc = 2*P10 + P5 + 2*P6 
	       + wsf*icfg[idxL(s[j-1])][idxP(s[j],s[i])]
	       + wsf*icfg[idxPS(s[i],s[j])][idxPS(s[i+2],s[i+2+mid])] 
	       + vx[(jmod-d+2+mid<0)? jmod-d+2+mid+win:jmod-d+2+mid][mid] 
	       + wbx[(jmod-2<0)? jmod-2+win:jmod-2][d-mid-5]) > bestsc)
	    bestsc = sc;

	  /* (V4_7) */
	  if (mid < d-4 &&
	      (sc = 2*P10 + P5 + 2*P6
	       + wsf*icfg[idxR(s[i+mid+3])][idxP(s[i+2],s[i+mid+2])]
	       + wsf*icfg[idxPS(s[i],s[j])][idxPS(s[i+2],s[i+2+mid])] 
	       + vx[(jmod-d+2+mid<0)? jmod-d+2+mid+win:jmod-d+2+mid][mid] 
	       + wbx[(jmod-1<0)? jmod-1+win:jmod-1][d-mid-5]) > bestsc)
	    bestsc = sc;
	  
	  /* (V4_8) */
	  if (mid < d-6 &&
	      (sc = 2*P10 + P5 + 3*P6 
	       + wsf*icfg[idxL(s[j-1])][idxP(s[j],s[i])]
	       + wsf*icfg[idxR(s[i+mid+3])][idxP(s[i+2],s[i+mid+2])]
	       + wsf*icfg[idxPS(s[i],s[j])][idxPS(s[i+2],s[i+2+mid])] 
	       + vx[(jmod-d+1+mid<0)? jmod-d+1+mid+win:jmod-d+1+mid][mid] 
	       + wbx[(jmod-2<0)? jmod-2+win:jmod-2][d-mid-6]) > bestsc)
	    bestsc = sc;
  
	  /* (V5)__/ 1 WBX 1VX 
	   */ 
	  
	  /* (V5_1) */
	  if (mid < d-3 &&
	      (sc = 2*P10 + P5 
	       + wsf*(icfg[idxPS(s[j-1],s[i+2+mid])][idxPS(s[j],s[i])] + INTSCALE)
	       + wbx[(jmod-d+1+mid<0)? jmod-d+1+mid+win:jmod-d+1+mid][mid] 
	       + vx[(jmod-1<0)? jmod-1+win:jmod-1][d-mid-3]) > bestsc)
	    bestsc = sc;
	  
	  /* (V5_2)  */
	  if (mid > 1 && mid < d-3 &&
	      (sc = 2*P10 + P5 + P6 
	       + wsf*icfg[idxR(s[i+1])][idxP(s[j],s[i])]
	       + wsf*icfg[idxPS(s[j-1],s[i+2+mid])][idxPS(s[j],s[i])]
	       + wbx[(jmod-d+1+mid<0)? jmod-d+1+mid+win:jmod-d+1+mid][mid-1] 
	       + vx[(jmod-1<0)? jmod-1+win:jmod-1][d-mid-3]) > bestsc)
	    bestsc = sc;
	  
	  /* (V5_3)  */
	  if (mid > 1 && mid < d-3 &&
	      (sc = 2*P10 + P5 + P6  
	       + wsf*icfg[idxL(s[i+mid+1])][idxP(s[i+2+mid],s[j-1])]
	       + wsf*icfg[idxPS(s[j-1],s[i+2+mid])][idxPS(s[j],s[i])] 
	       + wbx[(jmod-d+mid<0)? jmod-d+mid+win:jmod-d+mid][mid-1] 
	       + vx[(jmod-1<0)? jmod-1+win:jmod-1][d-mid-3]) > bestsc)
	    bestsc = sc;
	  
	  /* (V5_4) */
	  if (mid > 1 && mid < d-3 &&
	      (sc = 2*P10 + P5 + 2*P6 
	       + wsf*icfg[idxR(s[i+1])][idxP(s[j],s[i])]
	       + wsf*icfg[idxL(s[i+mid+1])][idxP(s[i+2+mid],s[j-1])]
	       + wsf*icfg[idxPS(s[j-1],s[i+2+mid])][idxPS(s[j],s[i])]
	       + wbx[(jmod-d+mid<0)? jmod-d+mid+win:jmod-d+mid][mid-2] 
	       + vx[(jmod-1<0)? jmod-1+win:jmod-1][d-mid-3]) > bestsc)
	    bestsc = sc;
	  
	  /* (V5_5) */
	  if (mid < d-4 &&
	      (sc = 2*P10 + P5 + P6
	       + wsf*icfg[idxPS(s[j-2],s[i+2+mid])][idxPS(s[j],s[i])]
	       + wbx[(jmod-d+1+mid<0)? jmod-d+1+mid+win:jmod-d+1+mid][mid] 
	       + vx[(jmod-2<0)? jmod-2+win:jmod-2][d-mid-4]) > bestsc)
	    bestsc = sc;
	  
	  /* (V5_6) */
	  if (mid > 1 && mid < d-4 &&
	      (sc = 2*P10 + P5 + 2*P6 
	       + wsf*icfg[idxR(s[i+1])][idxP(s[j],s[i])]
	       + wsf*icfg[idxPS(s[j-2],s[i+2+mid])][idxPS(s[j],s[i])]
	       + wbx[(jmod-d+1+mid<0)? jmod-d+1+mid+win:jmod-d+1+mid][mid-1] 
	       + vx[(jmod-2<0)? jmod-2+win:jmod-2][d-mid-4]) > bestsc)
	    bestsc = sc;
	  
	  /* (V5_7)  */
	  if (mid > 1 && mid < d-4 &&
	      (sc = 2*P10 + P5 + 2*P6  
	       + wsf*icfg[idxL(s[i+mid+1])][idxP(s[i+2+mid],s[j-2])]
	       + wsf*icfg[idxPS(s[j-2],s[i+2+mid])][idxPS(s[j],s[i])] 
	       + wbx[(jmod-d+mid<0)? jmod-d+mid+win:jmod-d+mid][mid-1] 
	       + vx[(jmod-2<0)? jmod-2+win:jmod-2][d-mid-4]) > bestsc)
	    bestsc = sc;
	  
	  /* (V5_8)  */
	  if (mid > 2 && mid < d-4 &&
	      (sc = 2*P10 + P5 + 3*P6 
	       + wsf*icfg[idxR(s[i+1])][idxP(s[j],s[i])]
	       + wsf*icfg[idxL(s[i+mid+1])][idxP(s[i+2+mid],s[j-2])]
	       + wsf*icfg[idxPS(s[j-2],s[i+2+mid])][idxPS(s[j],s[i])]
	       + wbx[(jmod-d+mid<0)? jmod-d+mid+win:jmod-d+mid][mid-2] 
	       + vx[(jmod-2<0)? jmod-2+win:jmod-2][d-mid-4]) > bestsc)
	  bestsc = sc;
	}      
      vx[jmod][d] = bestsc;
    }
  else vx[jmod][d] = -BIGINT;
  
}




