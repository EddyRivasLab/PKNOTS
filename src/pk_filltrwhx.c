/* filltrwhx.c 
 *
 * includes functions FillWHX,  TraceWHX
 * that fill and traceback the hole matix:
 *      whx[j][d][d1][d2]  (len x len x len x len)
 *                         [(j - d, j)           unknown relation]
 *                         [(j - d + d1, j - d2) unknown relation]
 *
 */
                         
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>

#include <easel.h>

#include "pknots.h"
#include "pk_filltrwbx.h"
#include "pk_irredsurf.h"
#include "pk_whxgraphs.h"
#include "pk_util.h"


/* Function: FillWHX()
 * 
 * Purpose:  fill hole matrix WHX. Using Knots-IS(2) method.
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
 *            whx is filled.
 */       
void
FillWHX(ESL_DSQ *s, int len, struct rnapar_2 *rnapar, int **icfg, int **wbx, int **vx, 
	 int ****whx, int ****vhx, int ****zhx, int ****yhx,
	 int j, int d, int d1, int d2)
{
  int            i, k, l;  /* coords in mtx's                 */
  int                mid;  /* midpoint of bifurcation         */
  int         mid1, mid2;  /* used for midpoint of a bifurc   */
  int                 sc;  /* temporary score to check        */
  int             bestsc;  /* best score so far               */

  i = j - d;
  k = i + d1;
  l = j - d2;

  if (d1 + d2 >= d - 1 )
    {
      whx[j][d][d1][d2] = wbx[j][d];
      return;
    }
  
  /* (WH1)__/ STRUCTURES WITH ONE VHX (sixteen).
   */ 

  /* (WH1.1) */
  bestsc = 2*rnapar->P10P + vhx[j][d][d1][d2];
  
  /* (WH1.2) */
  if (d1 > 1 &&
      (sc = 2*rnapar->P10P + 2*rnapar->P6P 
       + wkn*icfg[idxL(s[i])][idxP(s[i+1],s[j])] 
       + wkn*icfg[idxR(s[k])][idxP(s[k-1],s[l])] 
       + vhx[j][d-1][d1-2][d2]) > bestsc)
    bestsc = sc;
  
  /* (WH1.3) */
  if (d2 > 1 &&
      (sc = 2*rnapar->P10P + 2*rnapar->P6P 
       + wkn*icfg[idxL(s[l])][idxP(s[k],s[l+1])]
       + wkn*icfg[idxR(s[j])][idxP(s[i],s[j-1])]
       + vhx[j-1][d-1][d1][d2-2]) > bestsc)
    bestsc = sc;
  
  /* (WH1.4) */
  if (d1 > 0 && d2 > 0 &&
      (sc = 2*rnapar->P10P + 2*rnapar->P6P 
       + wkn*icfg[idxL(s[i])][idxP(s[i+1],s[j])]
       + wkn*icfg[idxL(s[l])][idxP(s[k],s[l+1])]
       + vhx[j][d-1][d1-1][d2-1]) > bestsc)
    bestsc = sc;
  
  /* (WH1.5) */
  if (d1 > 0 && d2 > 0 &&
      (sc = 2*rnapar->P10P + 2*rnapar->P6P 
       + wkn*icfg[idxR(s[k])][idxP(s[k-1],s[l])]
       + wkn*icfg[idxR(s[j])][idxP(s[i],s[j-1])]
       + vhx[j-1][d-1][d1-1][d2-1]) > bestsc)
    bestsc = sc;

  /* (WH1.6) */
  if (d1 > 1 && d2 > 0 &&
      (sc = 2*rnapar->P10P + 3*rnapar->P6P 
       + wkn*icfg[idxL(s[i])][idxP(s[i+1],s[j])]
       + wkn*icfg[idxR(s[k])][idxP(s[k-1],s[l+1])]
       + wkn*icfg[idxL(s[l])][idxP(s[k-1],s[l+1])] 
       + vhx[j][d-1][d1-2][d2-1]) > bestsc)
    bestsc = sc;
  
  /* (WH1.7) */
  if (d1 > 1 && d2 > 0 &&
      (sc = 2*rnapar->P10P + 3*rnapar->P6P 
       + wkn*icfg[idxL(s[i])][idxP(s[i+1],s[j-1])]
       + wkn*icfg[idxR(s[j])][idxP(s[i+1],s[j-1])]
       + wkn*icfg[idxR(s[k])][idxP(s[k-1],s[l])] 
       + vhx[j-1][d-2][d1-2][d2-1]) > bestsc)
    bestsc = sc;
  
  /* (WH1.8) */
  if (d1 > 0 && d2 > 1 &&
      (sc = 2*rnapar->P10P + 3*rnapar->P6P 
       + wkn*icfg[idxL(s[i])][idxP(s[i+1],s[j-1])]
       + wkn*icfg[idxR(s[j])][idxP(s[i+1],s[j-1])] 
       + wkn*icfg[idxL(s[l])][idxP(s[k],s[l+1])]
       + vhx[j-1][d-2][d1-1][d2-2]) > bestsc)
    bestsc = sc;
  
  /* (WH1.9) */
  if (d1 > 0 && d2 > 1 &&
      (sc = 2*rnapar->P10P + 3*rnapar->P6P 
       + wkn*icfg[idxR(s[k])][idxP(s[k-1],s[l+1])]
       + wkn*icfg[idxL(s[l])][idxP(s[k-1],s[l+1])] 
       + wkn*icfg[idxR(s[j])][idxP(s[i],s[j-1])]
       + vhx[j-1][d-1][d1-1][d2-2]) > bestsc)
    bestsc = sc;
  
  /* (WH1.10) */
  if (d1 > 1 && d2 > 1 &&
      (sc = 2*rnapar->P10P + 4*rnapar->P6P 
       + wkn*icfg[idxL(s[i])][idxP(s[i+1],s[j-1])]
       + wkn*icfg[idxR(s[j])][idxP(s[i+1],s[j-1])]
       + wkn*icfg[idxR(s[k])][idxP(s[k-1],s[l+1])]
       + wkn*icfg[idxL(s[l])][idxP(s[k-1],s[l+1])] 
       + vhx[j-1][d-2][d1-2][d2-2]) > bestsc)
    bestsc = sc;
  
  /*(WH2)__/ STRUCTURES WITH ONE ZHX (four) */ 
  
  /* (WH2.1) */
  if ((sc = rnapar->P10P + zhx[j][d][d1][d2]) > bestsc)
    bestsc = sc;
  
  /* (WH2.2) */
  if (d1 > 0 && d2 > 0 &&
      (sc = rnapar->P10P + 2*rnapar->P6P 
       + wkn*icfg[idxL(s[i])][idxP(s[i+1],s[j-1])]
       + wkn*icfg[idxR(s[j])][idxP(s[i+1],s[j-1])] 
       + zhx[j-1][d-2][d1-1][d2-1]) > bestsc)
    bestsc = sc;
  
  /* (WH2.3) */
  if (d1 > 0 && 
      (sc = rnapar->P10P + rnapar->P6P 
       + wkn*icfg[idxL(s[i])][idxP(s[i+1],s[j])]
       + zhx[j][d-1][d1-1][d2]) > bestsc)
    bestsc = sc;
  
  /* (WH2.4) */
  if (d2 > 0 &&
      (sc = rnapar->P10P + rnapar->P6P 
       + wkn*icfg[idxR(s[j])][idxP(s[i],s[j-1])]
       + zhx[j-1][d-1][d1][d2-1]) > bestsc)
    bestsc = sc;
  
  /* (WH3)__/ STRUCTURES WITH ONE YHX (four) */ 
  
  /* (WH3.1) */
  if ((sc = rnapar->P10P + yhx[j][d][d1][d2]) > bestsc)
    bestsc = sc;
  
  /* (WH3.2) */
  if (d1 > 0 && d2 > 0 &&
      (sc = rnapar->P10P + 2*rnapar->P6P 
       + wkn*icfg[idxR(s[k])][idxP(s[k-1],s[l+1])]
       + wkn*icfg[idxL(s[l])][idxP(s[k-1],s[l+1])] 
       + yhx[j][d][d1-1][d2-1]) > bestsc)
    bestsc = sc;
  
  /* (WH3.3) */
  if (d2 > 0 &&
      (sc = rnapar->P10P + rnapar->P6P 
       + wkn*icfg[idxL(s[l])][idxP(s[k],s[l+1])] 
       + yhx[j][d][d1][d2-1]) > bestsc)
    bestsc = sc;
  
  /* (WH3.4) */
  if (d1 > 0 &&
      (sc = rnapar->P10P + rnapar->P6P 
       + wkn*icfg[idxR(s[k])][idxP(s[k-1],s[l])]
       + yhx[j][d][d1-1][d2]) > bestsc)
    bestsc = sc;
  
  /* (WH4)__/ STRUCTURES WITH ONE WHX (four) */ 
  
  /* (WH4.1) */
  if (d1 > 0 &&
      (sc = rnapar->P6P + whx[j][d-1][d1-1][d2] ) > bestsc)
    bestsc = sc;
  
  /* (WH4.2) */
  if (d2 > 0 &&
      (sc = rnapar->P6P + whx[j-1][d-1][d1][d2-1]) > bestsc)
    bestsc = sc;
  
  /* (WH4.3) */ 
  if (d2 > 0 &&
      (sc = rnapar->P6P + whx[j][d][d1][d2-1]) > bestsc)
    bestsc = sc;
  
  /* (WH4.4) */
  if (d1 > 0 &&
      (sc = rnapar->P6P + whx[j][d][d1-1][d2]) > bestsc)
    bestsc = sc;
  
  /* (WH5)__/ 2 WBX  
   */
  if ((sc = wbx[j][d2] + wbx[k][d1]) > bestsc)
    bestsc = sc;
  
  for (mid = 0; mid < d1; mid++)
    { 
      /* (WH6)__/  1 WBX 1WHX
       */                                       
      if ((sc = wbx[i+mid][mid] + whx[j][d-mid-1][d1-mid-1][d2]) > bestsc)
	bestsc = sc;
      
      /* (WH7)__/  1 VX 1 ZHX
       */   
      /* (WH7.1) */                                   
      if ((sc = 2*rnapar->P10P
	   + wkn*(icfg[idxPS(s[i+mid],s[i])]
		  [idxPS(s[i+mid+1],s[j])] + INTSCALE)
	   + vx[i+mid][mid] + zhx[j][d-mid-1][d1-mid-1][d2]) > bestsc)
	bestsc = sc;
      
      /* (WH7.2) */                                   
      if (mid > 0 &&
	  (sc = 2*rnapar->P10P + rnapar->P6P
	   + wkn*icfg[idxL(s[i])][idxP(s[i+1],s[i+mid])]
	   + wkn*icfg[idxPS(s[i+mid],s[i+1])][idxPS(s[i+mid+1],s[j])]
	   + vx[i+mid][mid-1] + zhx[j][d-mid-1][d1-mid-1][d2]) > bestsc)
	bestsc = sc;
      
     /* (WH7.3) */                                   
      if (d2 > 0 &&
	  (sc = 2*rnapar->P10P + rnapar->P6P
	   + wkn*icfg[idxR(s[j])][idxP(s[i+mid+1],s[j-1])]
	   + wkn*icfg[idxPS(s[i+mid],s[i])][idxPS(s[i+mid+1],s[j-1])]
	   + vx[i+mid][mid] + zhx[j-1][d-mid-2][d1-mid-1][d2-1]) > bestsc)
	bestsc = sc;
      
      /* (WH7.4) */                                   
      if (mid > 0 && d2 > 0 &&
	  (sc = 2*rnapar->P10P + 2*rnapar->P6P
	   + wkn*icfg[idxL(s[i])][idxP(s[i+1],s[i+mid])]
	   + wkn*icfg[idxR(s[j])][idxP(s[i+mid+1],s[j-1])]
	   + wkn*icfg[idxPS(s[i+mid],s[i+1])][idxPS(s[i+mid+1],s[j-1])]
	   + vx[i+mid][mid-1] + zhx[j-1][d-mid-2][d1-mid-1][d2-1]) > bestsc)
	bestsc = sc;
      
      /* (WH8)__/ 1 WBX 1 WHX
       */                                       
      if ((sc = wbx[k][mid] + whx[j][d][d1-mid-1][d2]) > bestsc)
	bestsc = sc;
      
      /* (WH9)__/ 1 VX 1 YHX
      */
      /* (WH9.1) */                                   
      if ((sc = 2*rnapar->P10P
	   + wkn*(icfg[idxPS(s[k-mid-1],s[l])]
		  [idxPS(s[k-mid],s[k])] + INTSCALE)
	   + vx[k][mid] + yhx[j][d][d1-mid-1][d2]) > bestsc)
	bestsc = sc;
      
     /* (WH9.2) */                                   
      if (mid > 0 &&
	  (sc = 2*rnapar->P10P + rnapar->P6P
	   + wkn*icfg[idxR(s[k])][idxP(s[k-mid],s[k-1])]
	   + wkn*icfg[idxPS(s[k-mid-1],s[l])][idxPS(s[k-mid],s[k-1])]
	   + vx[k-1][mid-1] + yhx[j][d][d1-mid-1][d2]) > bestsc)
	bestsc = sc;
      
     /* (WH9.3) */                                   
      if (d2 > 0 &&
	  (sc = 2*rnapar->P10P + rnapar->P6P
	   + wkn*icfg[idxL(s[l])][idxP(s[l+1],s[k-mid-1])]
	   + wkn*icfg[idxPS(s[k-mid-1],s[l+1])][idxPS(s[k-mid],s[k])]
	   + vx[k][mid] + yhx[j][d][d1-mid-1][d2-1]) > bestsc)
	bestsc = sc;
      
     /* (WH9.4) */                                   
      if (mid > 0 && d2 > 0 &&
	  (sc = 2*rnapar->P10P + 2*rnapar->P6P
	   + wkn*icfg[idxR(s[k])][idxP(s[k-mid],s[k-1])]
	   + wkn*icfg[idxL(s[l])][idxP(s[l+1],s[k-mid-1])]
	   + wkn*icfg[idxPS(s[k-mid-1],s[l+1])][idxPS(s[k-mid],s[k-1])]
	   + vx[k-1][mid-1] + yhx[j][d][d1-mid-1][d2-1]) > bestsc)
	bestsc = sc;
    }
  
  for (mid = 0; mid < d2; mid++)
    { 
      /* (WH10)__/ 1 WBX 1 WHX
       */                                       
      if ((sc = wbx[j][mid] + whx[j-mid-1][d-mid-1][d1][d2-mid-1]) > bestsc)
	bestsc = sc;
      
      /* (WH11)__/ 1 VX 1 ZHX
       */                                       
      /* (WH11.1) */                                   
      if ((sc = 2*rnapar->P10P
	   + wkn*(icfg[idxPS(s[j-mid-1],s[i])]
		  [idxPS(s[j-mid],s[j])] + INTSCALE)
	   + vx[j][mid] + zhx[j-mid-1][d-mid-1][d1][d2-mid-1]) > bestsc)
	bestsc = sc;
      
      /* (WH11.2) */                                   
      if (d1 > 0 &&
	  (sc = 2*rnapar->P10P + rnapar->P6P
	   + wkn*icfg[idxL(s[i])][idxP(s[i+1],s[j-mid-1])]
	   + wkn*icfg[idxPS(s[j-mid-1],s[i+1])][idxPS(s[j-mid],s[j])]
	   + vx[j][mid] + zhx[j-mid-1][d-mid-2][d1-1][d2-mid-1]) > bestsc)
	bestsc = sc;
      
      /* (WH11.3) */                                   
      if (mid > 0 &&
	  (sc = 2*rnapar->P10P + rnapar->P6P
	   + wkn*icfg[idxR(s[j])][idxP(s[j-mid],s[j-1])]
	   + wkn*icfg[idxPS(s[j-mid-1],s[i])][idxPS(s[j-mid],s[j-1])]
	   + vx[j-1][mid-1] + zhx[j-mid-1][d-mid-1][d1][d2-mid-1]) > bestsc)
	bestsc = sc;
      
      /* (WH11.4) */                                         
      if (d1 > 0 && mid > 0 &&
	  (sc = 2*rnapar->P10P + 2*rnapar->P6P
	   + wkn*icfg[idxL(s[i])][idxP(s[i+1],s[j-mid-1])]
	   + wkn*icfg[idxR(s[j])][idxP(s[j-mid],s[j-1])]
	   + wkn*icfg[idxPS(s[j-mid-1],s[i+1])][idxPS(s[j-mid],s[j-1])]
	   + vx[j-1][mid-1] + zhx[j-mid-1][d-mid-2][d1-1][d2-mid-1]) > bestsc)
	bestsc = sc;
      
      /* (WH12)__/ 1 WBX 1 WHX
       */                                       
      if ((sc = wbx[l+mid][mid] + whx[j][d][d1][d2-mid-1]) > bestsc)
	bestsc = sc;
      
      /* (WH13)__/ 1 VX 1 YHX
       */ 
      /* (WH13.1) */                                   
      if ((sc = 2*rnapar->P10P
	   + wkn*(icfg[idxPS(s[l+mid],s[l])]
		  [idxPS(s[l+mid+1],s[k])] + INTSCALE)
	   + vx[l+mid][mid] + yhx[j][d][d1][d2-mid-1]) > bestsc)
	bestsc = sc;
      
      /* (WH13.2) */                                   
      if (d1 > 0 &&
	  (sc = 2*rnapar->P10P + rnapar->P6P
	   + wkn*icfg[idxL(s[k])][idxP(s[l+mid+1],s[k-1])]
	   + wkn*icfg[idxPS(s[l+mid],s[l])][idxPS(s[l+mid+1],s[k-1])]
	   + vx[l+mid][mid] + yhx[j][d][d1-1][d2-mid-1]) > bestsc)
	bestsc = sc;
      
      /* (WH13.3) */                                   
      if (mid > 0 &&
	  (sc = 2*rnapar->P10P + rnapar->P6P
	   + wkn*icfg[idxR(s[l])][idxP(s[l+1],s[l+mid])]
	   + wkn*icfg[idxPS(s[l+mid],s[l+1])][idxPS(s[l+mid+1],s[k])]
	   + vx[l+mid][mid-1] + yhx[j][d][d1][d2-mid-1]) > bestsc)
	bestsc = sc;
      
      /* (WH13.4) */                                   
      if (d1 > 0 && mid > 0 &&
	  (sc = 2*rnapar->P10P + 2*rnapar->P6P
	   + wkn*icfg[idxL(s[k])][idxP(s[l+mid+1],s[k-1])]
	   + wkn*icfg[idxR(s[l])][idxP(s[l+1],s[l+mid])]
	   + wkn*icfg[idxPS(s[l+mid],s[l+1])][idxPS(s[l+mid+1],s[k-1])]
	   + vx[l+mid][mid-1] + yhx[j][d][d1-1][d2-mid-1]) > bestsc)
	bestsc = sc;
    }
  
  for (mid1 = 0; mid1 < d1; mid1++)                     
    for (mid2 = 0; mid2 < d2; mid2++)
      {
	/* (WH14)__/  STRUCTURE WITH ONE YHX AND ONE ZHX. (only one pair)
	 */ 
	if ((sc = yhx[j][d][d1-mid1][d2-mid2]
	     + zhx[l+mid2][d-d1-d2+mid1+mid2][mid1][mid2]) > bestsc)
	  bestsc = sc;
	
	/* (WH15)__/  2 WHX (inclusive bifurcation)
	 */
	if (d1-mid1 > 1 && d2-mid2 > 2 &&
	    (sc = rnapar->P5P
	     + whx[j][d][mid1][mid2] 
	     + whx[j-mid2-1][d-mid1-mid2-2][d1-mid1-1][d2-mid2-1]) > bestsc)
	  bestsc = sc;
  
	/* (WH16)__/  2 WHX (crossed bifurcation)
	 */
	if (d1-mid1 > 1 && d2-mid2 > 2 &&
	    (sc = rnapar->P12 
	     + whx[l+mid2][d-d2+mid2][mid1][mid2] 
	     + whx[j][d-mid1-2][d1-mid1-2][d2-mid2-3]) > bestsc)
	  bestsc = sc;
      }
  
  for (mid1 = 0; mid1 < d2-1; mid1++)
    for (mid2 = 0; mid2 < d2-mid1-1; mid2++) 
      {    
	/* (WH17)__/ 2 WHX */
	if ((sc = rnapar->P12 
	     + whx[j-mid2-1][d-mid2-1][d1][d2-mid1-mid2-2] 
	     + whx[j][d2][mid1][mid2]) > bestsc)
	  bestsc = sc;
	
	/* (WH18)__/  1 ZHX 1 YHX
	 */
	/* (WH18.1) */                                   
	if ((sc = rnapar->P12 + 2*rnapar->P10P 
	     + wkn*(icfg[idxPS(s[j-mid2-1],s[i])]
		    [idxPS(s[j-mid2],s[l+mid1])] + INTSCALE)
	     + zhx[j-mid2-1][d-mid2-1][d1][d2-mid1-mid2-2] 
	     + yhx[j][d2][mid1][mid2]) > bestsc)
	  bestsc = sc;
	
	/* (WH18.2) */
	if (d1 > 0 &&
	    (sc = rnapar->P12 + 2*rnapar->P10P + rnapar->P6P
	     + wkn*icfg[idxL(s[i])][idxP(s[i+1],s[j-mid2-1])]
	     + wkn*icfg[idxPS(s[j-mid2-1],s[i+1])][idxPS(s[j-mid2],s[l+mid1])]
	     + zhx[j-mid2-1][d-mid2-2][d1-1][d2-mid1-mid2-2] 
	     + yhx[j][d2][mid1][mid2]) > bestsc)
	  bestsc = sc;
	
	/* (WH18.3) */
	if (d2-mid1-mid2 > 2 &&
	    (sc = rnapar->P12 + 2*rnapar->P10P + rnapar->P6P
	     + wkn*icfg[idxR(s[l+mid1+1])][idxP(s[j-mid2],s[l+mid1])]
	     + wkn*icfg[idxPS(s[j-mid2-1],s[i])][idxPS(s[j-mid2],s[l+mid1])]
	     + zhx[j-mid2-1][d-mid2-1][d1][d2-mid1-mid2-3] 
	     + yhx[j][d2][mid1][mid2]) > bestsc)
	  bestsc = sc;
	
	/* (WH18.4) */
	if (d1 > 0 && d2-mid1-mid2 > 2 &&
	    (sc = rnapar->P12 + 2*rnapar->P10P + 2*rnapar->P6P
	     + wkn*icfg[idxL(s[i])][idxP(s[i+1],s[j-mid2-1])]
	     + wkn*icfg[idxR(s[l+mid1+1])][idxP(s[j-mid2],s[l+mid1])]
	     + wkn*icfg[idxPS(s[j-mid2-1],s[i+1])][idxPS(s[j-mid2],s[l+mid1])]
	     + zhx[j-mid2-1][d-mid2-2][d1-1][d2-mid1-mid2-3] 
	     + yhx[j][d2][mid1][mid2]) > bestsc)
	  bestsc = sc;
	
	/* (WH19)__/ 2 YHX 
	 */
	/* (WH19.1) */                                   
	if ((sc =  rnapar->P12 + 2*rnapar->P10P 
	     + wkn*(icfg[idxPS(s[l+mid1],s[j-mid2])]
		    [idxPS(s[l+mid1+1],s[k])] + INTSCALE)
	     + yhx[j-mid2-1][d-mid2-1][d1][d2-mid1-mid2-2] 
	     + yhx[j][d2][mid1][mid2]) > bestsc)
	  bestsc = sc;
	
	/* (WH19.2) */
	if (d1 > 0 &&
	    (sc = rnapar->P12 + 2*rnapar->P10P + rnapar->P6P
	     + wkn*icfg[idxR(s[k])][idxP(s[l+mid1+1],s[k-1])]
	     + wkn*icfg[idxPS(s[l+mid1],s[j-mid2])][idxPS(s[l+mid1+1],s[k-1])]
	     + yhx[j-mid2-1][d-mid2-1][d1-1][d2-mid1-mid2-2] 
	     + yhx[j][d2][mid1][mid2]) > bestsc)
	  bestsc = sc;
 
	/* (WH19.3) */
	if (d2-mid1-mid2 > 2 &&
	    (sc = rnapar->P12 + 2*rnapar->P10P + rnapar->P6P
	     + wkn*icfg[idxL(s[j-mid2-1])][idxP(s[j-mid2],s[l+mid1])]
	     + wkn*icfg[idxPS(s[l+mid1],s[j-mid2])][idxPS(s[l+mid1+1],s[k])]
	     + yhx[j-mid2-2][d-mid2-2][d1][d2-mid1-mid2-3] 
	     + yhx[j][d2][mid1][mid2]) > bestsc)
	  bestsc = sc;
	
  
	/* (WH19.4) */
	if (d1 > 0 && d2-mid1-mid2 > 2 &&
	    (sc = rnapar->P12 + 2*rnapar->P10P + 2*rnapar->P6P
	     + wkn*icfg[idxR(s[k])][idxP(s[l+mid1+1],s[k-1])]
	     + wkn*icfg[idxL(s[j-mid2-1])][idxP(s[j-mid2],s[l+mid1])]
	     + wkn*icfg[idxPS(s[l+mid1],s[j-mid2])][idxPS(s[l+mid1+1],s[k-1])]
	     + yhx[j-mid2-2][d-mid2-2][d1-1][d2-mid1-mid2-3] 
	     + yhx[j][d2][mid1][mid2]) > bestsc)
	  bestsc = sc;
      }

  for (mid1 = 0; mid1 < d1-1; mid1++)                   
    for (mid2 = 0; mid2 < d1-mid1-1; mid2++)  
      {                       
	/* (WH20)__/ 2 WHX
	 */
	if ((sc = rnapar->P12 
	     + whx[j][d-mid1-1][d1-mid1-mid2-2][d2] 
	     + whx[k][d1][mid1][mid2]) > bestsc)
	  bestsc = sc;
	
	/* (WH21)__/ 1 ZHX 1 YHX
	 */
	/* (WH21.1) */                                   
	if ((sc = rnapar->P12 + 2*rnapar->P10P 
	     + wkn*(icfg[idxPS(s[i+mid1],s[k-mid2])]
		    [idxPS(s[i+mid1+1],s[j])] + INTSCALE)
	     + zhx[j][d-mid1-1][d1-mid1-mid2-2][d2] 
	     + yhx[k][d1][mid1][mid2]) > bestsc)
	  bestsc = sc;
	
	/* (WH21.2) */
	if (d2 > 0 &&
	    (sc = rnapar->P12 + 2*rnapar->P10P + rnapar->P6P
	     + wkn*icfg[idxR(s[j])][idxP(s[i+mid1+1],s[j-1])]
	     + wkn*icfg[idxPS(s[i+mid1],s[k-mid2])][idxPS(s[i+mid1+1],s[j-1])]
	     + zhx[j-1][d-mid1-2][d1-mid1-mid2-2][d2-1] 
	     + yhx[k][d1][mid1][mid2]) > bestsc)
	  bestsc = sc;
	
	/* (WH21.3) */
	if (d1-mid1-mid2 > 2 &&
	    (sc = rnapar->P12 + 2*rnapar->P10P + rnapar->P6P
	     + wkn*icfg[idxL(s[k-mid2-1])][idxP(s[ k-mid2],s[i+mid1])]
	     + wkn*icfg[idxPS(s[i+mid1],s[k-mid2])][idxPS(s[i+mid1+1],s[j])] 
	     + zhx[j][d-mid1-1][d1-mid1-mid2-3][d2] 
	     + yhx[k][d1][mid1][mid2]) > bestsc)
	  bestsc = sc;

	/* (WH21.4) */
	if (d2 > 0 && d1-mid1-mid2 > 2 &&
	    (sc = rnapar->P12 + 2*rnapar->P10P + 2*rnapar->P6P
	     + wkn*icfg[idxR(s[j])][idxP(s[i+mid1+1],s[j-1])]
	     + wkn*icfg[idxL(s[k-mid2-1])][idxP(s[ k-mid2],s[i+mid1])]
	     + wkn*icfg[idxPS(s[i+mid1],s[k-mid2])][idxPS(s[i+mid1+1],s[j-1])]
	     + zhx[j-1][d-mid1-2][d1-mid1-mid2-3][d2-1] 
	     + yhx[k][d1][mid1][mid2]) > bestsc)
	  bestsc = sc;

	/* (WH22)__/ 2 YHX 
	 */
	/* (WH22.1) */                                   
	if ((sc = rnapar->P12 + 2*rnapar->P10P 
	     + wkn*(icfg[idxPS(s[k-mid2-1],s[l])]
		    [idxPS(s[k-mid2],s[i+mid1])] + INTSCALE)
	     + yhx[j][d-mid1-1][d1-mid1-mid2-2][d2] 
	     + yhx[k][d1][mid1][mid2]) > bestsc)
	  bestsc = sc;
	
	/* (WH22.2) */
	if (d2 > 0 &&
	    (sc = rnapar->P12 + 2*rnapar->P10P + rnapar->P6P
	     + wkn*icfg[idxL(s[l])][idxP(s[l+1],s[k-mid2-1])]
	     + wkn*icfg[idxPS(s[k-mid2-1],s[l+1])][idxPS(s[k-mid2],s[i+mid1])]
	     + yhx[j][d-mid1-1][d1-mid1-mid2-2][d2-1] 
	     + yhx[k][d1][mid1][mid2]) > bestsc)
	  bestsc = sc;
	
	/* (WH22.3) */
	if (d1-mid1-mid2 > 2 &&
	    (sc = rnapar->P12 + 2*rnapar->P10P + rnapar->P6P 
	     + wkn*icfg[idxR(s[i+mid1+1])][idxP(s[k-mid2],s[i+mid1])]
	     + wkn*icfg[idxPS(s[k-mid2-1],s[l])][idxPS(s[k-mid2],s[i+mid1])]
	     + yhx[j][d-mid1-2][d1-mid1-mid2-3][d2] 
	     + yhx[k][d1][mid1][mid2]) > bestsc)
	  bestsc = sc;

	/* (WH22.4) */
	if (d2 > 0 && d1-mid1-mid2 > 2 &&
	    (sc = rnapar->P12 + 2*rnapar->P10P + 2*rnapar->P6P
	     + wkn*icfg[idxL(s[l])][idxP(s[l+1],s[k-mid2-1])]
	     + wkn*icfg[idxR(s[i+mid1+1])][idxP(s[k-mid2],s[i+mid1])]
	     + wkn*icfg[idxPS(s[k-mid2-1],s[l+1])][idxPS(s[k-mid2],s[i+mid1])]
	     + yhx[j][d-mid1-2][d1-mid1-mid2-3][d2-1] 
	     + yhx[k][d1][mid1][mid2]) > bestsc)
	  bestsc = sc;
      }
  whx[j][d][d1][d2] = bestsc;
}

/* Function: TraceWHX()
 * 
 * Purpose:  traceback hole matrix WHX. Using Knots-IS(2) method.
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
TraceWHX(FILE *outf, ESL_DSQ *s, int len, struct rnapar_2 *rnapar, int **icfg, int **wbx, int **vx, 
	 int ****whx, int ****vhx, int ****zhx, int ****yhx,
	 int j, int d, int d1, int d2, 
	 struct tracekn_s *curr_tr, struct traceknstack_s *dolist, int traceback)
{
  int         i, k, l;   /* coords in mtx's                 */
  int mid, mid1, mid2;  /* used for midpoint of a bifurc   */

  i = j - d;
  k = i + d1;
  l = j - d2;

  /* (WH1)__/ 1 VHX
   */ 
  /* (WH1.1) */
  if (whx[j][d][d1][d2] == WH1_1(s, len, rnapar, icfg, vhx, j, d, d1, d2)){
    trace_WH1_1(outf, vhx, j, d, d1, d2, curr_tr, dolist, traceback);
    return;
  }   
  
  /* (WH1.2) */
  else if (whx[j][d][d1][d2] == WH1_2(s, len, rnapar, icfg, vhx, j, d, d1, d2)){
    trace_WH1_2(outf, vhx, j, d, d1, d2, curr_tr, dolist, traceback);
    return;
  }   
  
  /* (WH1.3) */
  else if (whx[j][d][d1][d2] == WH1_3(s, len, rnapar, icfg, vhx, j, d, d1, d2)){ 
   trace_WH1_3(outf, vhx, j, d, d1, d2, curr_tr, dolist, traceback);
   return;
  }
  
  /* (WH1.4) */
  else if (whx[j][d][d1][d2] == WH1_4(s, len, rnapar, icfg, vhx, j, d, d1, d2)){
    trace_WH1_4(outf, vhx, j, d, d1, d2, curr_tr, dolist, traceback);
    return;
  }
  
  /* (WH1.5) */
  if (whx[j][d][d1][d2] == WH1_5(s, len, rnapar, icfg, vhx, j, d, d1, d2)){    
    trace_WH1_5(outf, vhx, j, d, d1, d2, curr_tr, dolist, traceback);
    return;
  }
  
  /* (WH1.6) */
  else if (whx[j][d][d1][d2] == WH1_6(s, len, rnapar, icfg, vhx, j, d, d1, d2)){ 
    trace_WH1_6(outf, vhx, j, d, d1, d2, curr_tr, dolist, traceback);
    return;
  }
  
  /* (WH1.7) */
  else if (whx[j][d][d1][d2] == WH1_7(s, len, rnapar, icfg, vhx, j, d, d1, d2)){ 
    trace_WH1_7(outf, vhx, j, d, d1, d2, curr_tr, dolist, traceback);
    return;
  }
  
  /* (WH1.8) */
  else if (whx[j][d][d1][d2] == WH1_8(s, len, rnapar, icfg, vhx, j, d, d1, d2)){ 
    trace_WH1_8(outf, vhx, j, d, d1, d2, curr_tr, dolist, traceback);
    return;
  }
  
  /* (WH1.9) */
  else if (whx[j][d][d1][d2] == WH1_9(s, len, rnapar, icfg, vhx, j, d, d1, d2)){ 
    trace_WH1_9(outf, vhx, j, d, d1, d2, curr_tr, dolist, traceback);
    return;
  }
  
  /* (WH1.10) */
  else if (whx[j][d][d1][d2] == WH1_10(s, len, rnapar, icfg, vhx, j, d, d1, d2)){
    trace_WH1_10(outf, vhx, j, d, d1, d2, curr_tr, dolist, traceback);
    return;
  }
  
  /* (WH2)__/ 1 ZHX
   */ 
  /* (WH2.1) */
  else if (whx[j][d][d1][d2] == WH2_1(s, len, rnapar, icfg, zhx, j, d, d1, d2)){ 
    trace_WH2_1(outf, zhx, j, d, d1, d2, curr_tr, dolist, traceback);
    return;
  }
  
  /* (WH2.2) */
  else if (whx[j][d][d1][d2] == WH2_2(s, len, rnapar, icfg, zhx, j, d, d1, d2)){ 
    trace_WH2_2(outf, zhx, j, d, d1, d2, curr_tr, dolist, traceback);
    return;
  }
  
  /* (WH2.3) */
  else if (whx[j][d][d1][d2] == WH2_3(s, len, rnapar, icfg, zhx, j, d, d1, d2)){ 
    trace_WH2_3(outf, zhx, j, d, d1, d2, curr_tr, dolist, traceback);
    return;
  }
  
  /* (WH2.4) */
  else if (whx[j][d][d1][d2] == WH2_4(s, len, rnapar, icfg, zhx, j, d, d1, d2)){ 
    trace_WH2_4(outf, zhx, j, d, d1, d2, curr_tr, dolist, traceback);
    return;
  }
  
  /* (WH3)__/ 1 YHX
   */ 
  /* (WH3.1) */
  else if (whx[j][d][d1][d2] == WH3_1(s, len, rnapar, icfg, yhx, j, d, d1, d2)){ 
    trace_WH3_1(outf, yhx, j, d, d1, d2, curr_tr, dolist, traceback);
    return;
  }
  
  /* (WH3.2) */
  else if (whx[j][d][d1][d2] == WH3_2(s, len, rnapar, icfg, yhx, j, d, d1, d2)){ 
    trace_WH3_2(outf, yhx, j, d, d1, d2, curr_tr, dolist, traceback);
    return;
  }

  /* (WH3.3) */
  else if (whx[j][d][d1][d2] == WH3_3(s, len, rnapar, icfg, yhx, j, d, d1, d2)){
    trace_WH3_3(outf, yhx, j, d, d1, d2, curr_tr, dolist, traceback);
    return;
  }
  
  /* (WH3.4) */
  else if (whx[j][d][d1][d2] == WH3_4(s, len, rnapar, icfg, yhx, j, d, d1, d2)){ 
    trace_WH3_4(outf, yhx, j, d, d1, d2, curr_tr, dolist, traceback);
    return;
  }
  
  /* (WH4)__/ 1 WHX
   */ 
  /* (WH4.1) */
  else if (whx[j][d][d1][d2] == WH4_1(s, len, rnapar, icfg, whx, j, d, d1, d2)){ 
    trace_WH4_1(outf, whx, j, d, d1, d2, curr_tr, dolist, traceback);
    return;
  }
  
  /* (WH4.2) */
  else if (whx[j][d][d1][d2] == WH4_2(s, len, rnapar, icfg, whx, j, d, d1, d2)){
    trace_WH4_2(outf, whx, j, d, d1, d2, curr_tr, dolist, traceback);
    return;
  }
  
  /* (WH4.3) */ 
  else if (whx[j][d][d1][d2] == WH4_3(s, len, rnapar, icfg, whx, j, d, d1, d2)){
    trace_WH4_3(outf, whx, j, d, d1, d2, curr_tr, dolist, traceback);
    return;
  }
  
  /* (WH4.4) */
  else if (whx[j][d][d1][d2] == WH4_4(s, len, rnapar, icfg, whx, j, d, d1, d2)){ 
    trace_WH4_4(outf, whx, j, d, d1, d2, curr_tr, dolist, traceback);
    return;
  }
  
  /* (WH5)__/ 2 WBX  
   */
  else if (whx[j][d][d1][d2] == WH5(s, len, rnapar, icfg, wbx, j, d, d1, d2)){
    trace_WH5(outf, wbx, j, d, d1, d2, curr_tr, dolist, traceback);
    return;
  }
  
  for (mid = d1-1; mid > 0; mid--) {
    
    /* (WH6)__/  1 WBX, 1 WHX
     */                                       
    if (whx[j][d][d1][d2] == WH6(s, len, rnapar, icfg, wbx, whx, j, d, d1, d2, mid)){ 
      trace_WH6(outf, wbx, whx, j, d, d1, d2, mid, curr_tr, dolist, traceback);
      return;
    }
	  
    /* (WH7)__/  1 VX, 1 ZHX
     */                                       
    /* (WH7.1) */                                   
    else if (whx[j][d][d1][d2] == WH7_1(s, len, rnapar, icfg, vx, zhx, j, d, d1, d2, mid)){ 
      trace_WH7_1(outf, vx, zhx, j, d, d1, d2, mid, curr_tr, dolist, traceback);
      return;
    }
	  
    /* (WH7.2) */                                   
    else if (whx[j][d][d1][d2] == WH7_2(s, len, rnapar, icfg, vx, zhx, j, d, d1, d2, mid)){ 
      trace_WH7_2(outf, vx, zhx, j, d, d1, d2, mid, curr_tr, dolist, traceback);
      return;
    }
	  
    /* (WH7.3) */                                   
    else if (whx[j][d][d1][d2] == WH7_3(s, len, rnapar, icfg, vx, zhx, j, d, d1, d2, mid)){ 
      trace_WH7_3(outf, vx, zhx, j, d, d1, d2, mid, curr_tr, dolist, traceback);
      return;
    }
	  
    /* (WH7.4) */                                   
    else if (whx[j][d][d1][d2] == WH7_4(s, len, rnapar, icfg, vx, zhx, j, d, d1, d2, mid)){ 
      trace_WH7_4(outf, vx, zhx, j, d, d1, d2, mid, curr_tr, dolist, traceback);
      return;
    }
	  
    /* (WH8)__/ 1 WBX 1 WHX
     */                                       
    else if (whx[j][d][d1][d2] == WH8(s, len, rnapar, icfg, wbx, whx, j, d, d1, d2, mid)){  
      trace_WH8(outf, wbx, whx, j, d, d1, d2, mid, curr_tr, dolist, traceback);
      return;
    }

    /* (WH9)__/ 1 VX 1 VHX
     */
    /* (WH9.1) */                                   
    else if (whx[j][d][d1][d2] == WH9_1(s, len, rnapar, icfg, vx, vhx, j, d, d1, d2, mid)){  
      trace_WH9_1(outf, vx, whx, j, d, d1, d2, mid, curr_tr, dolist, traceback);
      return;
    }

    /* (WH9.2) */                                   
    else if (whx[j][d][d1][d2] == WH9_2(s, len, rnapar, icfg, vx, vhx, j, d, d1, d2, mid)) { 
      trace_WH9_2(outf, vx, whx, j, d, d1, d2, mid, curr_tr, dolist, traceback);
      return;
    }

    /* (WH9.3) */                                   
    else if (whx[j][d][d1][d2] == WH9_3(s, len, rnapar, icfg, vx, vhx, j, d, d1, d2, mid)){  
      trace_WH9_3(outf, vx, whx, j, d, d1, d2, mid, curr_tr, dolist, traceback);
      return;
    }

    /* (WH9.4) */                                   
    else if (whx[j][d][d1][d2] == WH9_4(s, len, rnapar, icfg, vx, vhx, j, d, d1, d2, mid)){   
      trace_WH9_4(outf, vx, whx, j, d, d1, d2, mid, curr_tr, dolist, traceback);
      return;
    }
  }

  for (mid = d2 - 1; mid > 0; mid--) {

    /* (WH10)__/ 1 WBX 1 WHX
     */                                       
    if (whx[j][d][d1][d2] == WH10(s, len, rnapar, icfg, wbx, whx, j, d, d1, d2, mid)){ 
      trace_WH10(outf, wbx, whx, j, d, d1, d2, mid, curr_tr, dolist, traceback);
      return;
    }
      
    /* (WH11)__/ 1 VX 1 ZHX
     */                                       
    /* (WH11.1) */                                   
    else if (whx[j][d][d1][d2] == WH11_1(s, len, rnapar, icfg, vx, zhx, j, d, d1, d2, mid)){ 
      trace_WH11_1(outf, vx, zhx, j, d, d1, d2, mid, curr_tr, dolist, traceback);
      return;
    }
      
    /* (WH11.2) */                                   
    else if (whx[j][d][d1][d2] == WH11_2(s, len, rnapar, icfg, vx, zhx, j, d, d1, d2, mid)){ 
      trace_WH11_2(outf, vx, zhx, j, d, d1, d2, mid, curr_tr, dolist, traceback);
      return;
    }
      
    /* (WH11.3) */                                   
    else if (whx[j][d][d1][d2] == WH11_3(s, len, rnapar, icfg, vx, zhx, j, d, d1, d2, mid)){ 
      trace_WH11_3(outf, vx, zhx, j, d, d1, d2, mid, curr_tr, dolist, traceback);
      return;
    }
      
    /* (WH11.4) */                                       
    else if (whx[j][d][d1][d2] == WH11_4(s, len, rnapar, icfg, vx, zhx, j, d, d1, d2, mid)){ 
      trace_WH11_4(outf, vx, zhx, j, d, d1, d2, mid, curr_tr, dolist, traceback);
      return;
    }
      
    /* (WH12)__/ 1 WBX 1 WHX
     */                                       
    else if (whx[j][d][d1][d2] == WH12(s, len, rnapar, icfg, wbx, whx, j, d, d1, d2, mid)){     
      trace_WH12(outf, wbx, whx, j, d, d1, d2, mid, curr_tr, dolist, traceback);
      return;
    }
      
    /* (WH13)__/ 1 VX 1 YHX
     */ 
    /* (WH13.1) */                                   
    else if (whx[j][d][d1][d2] == WH13_1(s, len, rnapar, icfg, vx, yhx, j, d, d1, d2, mid)){     
      trace_WH13_1(outf, vx, yhx, j, d, d1, d2, mid, curr_tr, dolist, traceback);
      return;
    }
      
    /* (WH13.2) */                                   
    else if (whx[j][d][d1][d2] == WH13_2(s, len, rnapar, icfg, vx, yhx, j, d, d1, d2, mid)){     
      trace_WH13_2(outf, vx, yhx, j, d, d1, d2, mid, curr_tr, dolist, traceback);
      return;
    }
      
    /* (WH13.3) */                                   
    else if (whx[j][d][d1][d2] == WH13_3(s, len, rnapar, icfg, vx, yhx, j, d, d1, d2, mid)){     
      trace_WH13_3(outf, vx, yhx, j, d, d1, d2, mid, curr_tr, dolist, traceback);
      return;
    }
      
    /* (WH13.4) */                                       
    else if (whx[j][d][d1][d2] == WH13_4(s, len, rnapar, icfg, vx, yhx, j, d, d1, d2, mid)){     
      trace_WH13_4(outf, vx, yhx, j, d, d1, d2, mid, curr_tr, dolist, traceback);
      return;
    }
  }

  for ((mid1 = d1 - 1); (mid1 >= 0); mid1--)
    for ((mid2 = d2 - 1); (mid2 >= 0); mid2--) {
      
      /* (WH14)__/ STRUCTURE WITH ONE YHX AND ONE ZHX. (only one pair)
       */ 
      if (whx[j][d][d1][d2] == WH14(s, len, rnapar, icfg, zhx, yhx, j, d, d1, d2, mid1, mid2)){    
	trace_WH14(outf, zhx, yhx, j, d, d1, d2, mid1, mid2, curr_tr, dolist, traceback);
	return;
      } 
	  
      /* (WH15)__/  2 WHX (inclusive bifurcation)
       */
      else if (whx[j][d][d1][d2] == WH15(s, len, rnapar, icfg, whx, j, d, d1, d2, mid1, mid2)){    
	trace_WH15(outf, whx, j, d, d1, d2, mid1, mid2, curr_tr, dolist, traceback);
	return;
      } 

      /* (WH16)__/  2 WHX (crossed bifurcation)
       */
      else if (whx[j][d][d1][d2] == WH16(s, len, rnapar, icfg, whx, j, d, d1, d2, mid1, mid2)){    
	trace_WH16(outf, whx, j, d, d1, d2, mid1, mid2, curr_tr, dolist, traceback);
	return;
      } 
    }

  for ((mid1 = d2 - 2); (mid1 >= 0); mid1--) 
    for (mid2 = d2 - mid1 - 2; mid2 >= 0; mid2--) {

      /* (WH17)__/  2 WHX 
       */
      if (whx[j][d][d1][d2] == WH17(s, len, rnapar, icfg, whx, j, d, d1, d2, mid1, mid2)){         
	trace_WH17(outf, whx, j, d, d1, d2, mid1, mid2, curr_tr, dolist, traceback);
	return;
      }
	  
      /* (WH18)__/  1 ZHX 1 YHX
       */  
      /* (WH18.1) */                                   
      else if (whx[j][d][d1][d2] == WH18_1(s, len, rnapar, icfg, zhx, yhx, j, d, d1, d2, mid1, mid2)){ 
	trace_WH18_1(outf, zhx, yhx, j, d, d1, d2, mid1, mid2, curr_tr, dolist, traceback);
	return;
      }
      
      /* (WH18.2) */
      else if (whx[j][d][d1][d2] == WH18_2(s, len, rnapar, icfg, zhx, yhx, j, d, d1, d2, mid1, mid2)){  
       	trace_WH18_2(outf, zhx, yhx, j, d, d1, d2, mid1, mid2, curr_tr, dolist, traceback);
	return;
      }

      /* (WH18.3) */
      else if (whx[j][d][d1][d2] == WH18_3(s, len, rnapar, icfg, zhx, yhx, j, d, d1, d2, mid1, mid2)){
	trace_WH18_3(outf, zhx, yhx, j, d, d1, d2, mid1, mid2, curr_tr, dolist, traceback);
	return;
      }
      
      /* (WH18.4) */
      else if (whx[j][d][d1][d2] == WH18_4(s, len, rnapar, icfg, zhx, yhx, j, d, d1, d2, mid1, mid2)){     
	trace_WH18_4(outf, zhx, yhx, j, d, d1, d2, mid1, mid2, curr_tr, dolist, traceback);
	return;
      }

      /* (WH19)__/ 2 YHX 
       */
      /* (WH19.1) */                                   
      else if (whx[j][d][d1][d2] == WH19_1(s, len, rnapar, icfg, yhx, j, d, d1, d2, mid1, mid2)){        
	trace_WH19_1(outf, yhx, j, d, d1, d2, mid1, mid2, curr_tr, dolist, traceback);
	return;
      }
	  
      /* (WH19.2) */
      else if (whx[j][d][d1][d2] == WH19_2(s, len, rnapar, icfg, yhx, j, d, d1, d2, mid1, mid2)){         
	trace_WH19_2(outf, yhx, j, d, d1, d2, mid1, mid2, curr_tr, dolist, traceback);
      }

      /* (WH19.3) */
      else if (whx[j][d][d1][d2] == WH19_3(s, len, rnapar, icfg, yhx, j, d, d1, d2, mid1, mid2)){         
	trace_WH19_3(outf, yhx, j, d, d1, d2, mid1, mid2, curr_tr, dolist, traceback);
	return;
      }
	  
      /* (WH19.4) */
      else if (whx[j][d][d1][d2] == WH19_4(s, len, rnapar, icfg, yhx, j, d, d1, d2, mid1, mid2)){       
	trace_WH19_4(outf, yhx, j, d, d1, d2, mid1, mid2, curr_tr, dolist, traceback);
	return;
      }
    }

  for ((mid1 = d1 - 2); (mid1 >= 0); mid1--) 
    for ((mid2 = d1 - mid1 - 2); (mid2 >= 0); mid2--) {

      /* (WH20)__/  2 WHX */
      if (whx[j][d][d1][d2] == WH20(s, len, rnapar, icfg, whx, j, d, d1, d2, mid1, mid2)){        
	trace_WH20(outf, whx, j, d, d1, d2, mid1, mid2, curr_tr, dolist, traceback);
	return;
      }
      
      /* (WH21)__/ 1 ZHX 1 YHX
       */
      /* (WH21.1) */                                   
      else if (whx[j][d][d1][d2] == WH21_1(s, len, rnapar, icfg, zhx, yhx, j, d, d1, d2, mid1, mid2)){
	trace_WH21_1(outf, zhx, yhx, j, d, d1, d2, mid1, mid2, curr_tr, dolist, traceback);
	return;
      }

      /* (WH21.2) */
      else if (whx[j][d][d1][d2] == WH21_2(s, len, rnapar, icfg, zhx, yhx, j, d, d1, d2, mid1, mid2)){
	trace_WH21_2(outf, zhx, yhx, j, d, d1, d2, mid1, mid2, curr_tr, dolist, traceback);
	return;
      }

      /* (WH21.3) */
      else if (whx[j][d][d1][d2] == WH21_3(s, len, rnapar, icfg, zhx, yhx, j, d, d1, d2, mid1, mid2)){
	trace_WH21_3(outf, zhx, yhx, j, d, d1, d2, mid1, mid2, curr_tr, dolist, traceback);
	return;
      }
      
      /* (WH21.4) */
      else if (whx[j][d][d1][d2] == WH21_4(s, len, rnapar, icfg, zhx, yhx, j, d, d1, d2, mid1, mid2)){
	trace_WH21_4(outf, zhx, yhx, j, d, d1, d2, mid1, mid2, curr_tr, dolist, traceback);
	return;
      }
      
      /* (WH22)__/ 2 YHX 
       */
      /* (WH22.1) */                                   
      else if (whx[j][d][d1][d2] == WH22_1(s, len, rnapar, icfg, yhx, j, d, d1, d2, mid1, mid2)){        
	trace_WH22_1(outf, yhx, j, d, d1, d2, mid1, mid2, curr_tr, dolist, traceback);
	return;
      }

      /* (WH22.2) */
      else if (whx[j][d][d1][d2] == WH22_2(s, len, rnapar, icfg, yhx, j, d, d1, d2, mid1, mid2)){        
	trace_WH22_2(outf, yhx, j, d, d1, d2, mid1, mid2, curr_tr, dolist, traceback);
	return;
      }

      /* (WH22.3) */
      else if (whx[j][d][d1][d2] == WH22_3(s, len, rnapar, icfg, yhx, j, d, d1, d2, mid1, mid2)){        
	trace_WH22_3(outf, yhx, j, d, d1, d2, mid1, mid2, curr_tr, dolist, traceback);
	return;
      }

      /* (WH22.4) */
      else if (whx[j][d][d1][d2] == WH22_4(s, len, rnapar, icfg, yhx, j, d, d1, d2, mid1, mid2)){        
	trace_WH22_4(outf, yhx, j, d, d1, d2, mid1, mid2, curr_tr, dolist, traceback);
	return;
      }
    }
}


