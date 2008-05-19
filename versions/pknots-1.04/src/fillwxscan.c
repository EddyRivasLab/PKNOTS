/* 
 * fillwxscan.c 
 *
 * includes functions FillWX, 
 * fills the no-hole matrix:
 *       wx[j][d]  (len x j) [(j-d,j) are not inside a base paired structure]
 *                           [(j-d,j) not necessarely based paired]
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

/* Function: FillWX_nestedScan()
 * 
 * Purpose:   fill no-hole matrix WX without pseudoknots.
 *           
 * Arguments: s     - sequence (0..len-1) converted to integers
 *                    representing array indices of symbols
 *            len   - length of iseq
 *            cfg   - context-free grammar state transitions, float log2 form
 *             wx   - DP matrix, already alloc'ed 
 *             vx   - DP matrix, already alloc'ed 
 *            
 * Return:    (void)           
 *            wx is filled.
 */       
void
FillWX_nestedScan(int *s, int len, int win, int **icfg, int **wx, int **vx, 
		  int j, int jmod, int d)
{  
  int                  i;  /* coords in mtx's                 */
  int                mid;  /* midpoint of bifurcation         */
  int                 sc;  /* temporary score to check        */
  int             bestsc;  /* temporary best score            */
 
  i = j - d;

  /* (W1)__/ PAIR (i,j).
   */
  bestsc = vx[jmod][d];	

  /* (W2)__/ PAIR (i+1,j-1).
   */
  if (d > 1 &&
      (sc = wsf*icfg[idxL(s[i])][idxP(s[i+1],s[j-1])]
       + wsf*icfg[idxR(s[j])][idxP(s[i+1],s[j-1])] 
       + vx[(jmod-1<0)? jmod-1+win:jmod-1][d-2]) > bestsc)
    bestsc = sc;
  
  /* (W3)__/  SINGLET-LEFT  
   */
   if (d > 0 &&
       (sc = wsf*icfg[idxL(s[i])][idxP(s[i+1],s[j])] 
	+ vx[jmod][d-1]) > bestsc)
     bestsc = sc; 
  
  /* (W4)__/  SINGLET-RIGHT. 
   */
  if (d > 0 &&
      (sc = wsf*icfg[idxR(s[j])][idxP(s[i],s[j-1])] 
       + vx[(jmod-1<0)? jmod-1+win:jmod-1][d-1]) > bestsc)
    bestsc = sc; 
  
  /* (W5)__/  SINGLET-LEFT. 
   */
  if (d > 1 &&
      (sc = wx[jmod][d-1]) > bestsc)
    bestsc = sc;
  
  /* (W6)__/  SINGLET-RIGHT.
   */
  if (d > 0 &&
      (sc = wx[(jmod-1<0)? jmod-1+win:jmod-1][d-1]) > bestsc)
    bestsc = sc; 
  
  for (mid = 1; mid < d; mid++) 
    {    
      /* (W7)__/  2 WX 
       */
      if ((sc = wx[(jmod-d+mid<0)? jmod-d+mid+win:jmod-d+mid][mid] 
	   + wx[jmod][d-mid-1]) > bestsc)
	bestsc = sc; 
      
      /* (W8)__/ 2 VX 
       */
      /* (W8.1) */
      if ((sc = wsf*(icfg[idxPS(s[i+mid],s[i])][idxPS(s[i+mid+1],s[j])] + INTSCALE)
	   + vx[(jmod-d+mid<0)? jmod-d+mid+win:jmod-d+mid][mid] 
	   + vx[jmod][d-mid-1]) > bestsc)
	bestsc = sc; 
      
      /* (W8.2) */
      if ((sc = wsf*icfg[idxL(s[i])][idxP(s[i+1],s[i+mid])]
	   + wsf*icfg[idxPS(s[i+mid],s[i+1])][idxPS(s[i+mid+1],s[j])]
	   + vx[(jmod-d+mid<0)? jmod-d+mid+win:jmod-d+mid][mid-1] 
	   + vx[jmod][d-mid-1]) > bestsc)
	bestsc = sc; 
      
      /* (W8.3) */
      if (d-mid > 1 &&
	  (sc = wsf*icfg[idxR(s[j])][idxP(s[i+mid+1],s[j-1])]
	   + wsf*icfg[idxPS(s[i+mid],s[i])][idxPS(s[i+mid+1],s[j-1])]
	   + vx[(jmod-d+mid<0)? jmod-d+mid+win:jmod-d+mid][mid] 
	   + vx[(jmod-1<0)? jmod-1+win:jmod-1][d-mid-2]) > bestsc)
	bestsc = sc; 
      
      /* (W8.4)__/ i dangles off i+1,i+mid / j dangles off i+mid+1,j-1  */
      if (d-mid > 1 &&
	  (sc = wsf*icfg[idxL(s[i])][idxP(s[i+1],s[i+mid])]
	   + wsf*icfg[idxR(s[j])][idxP(s[i+mid+1],s[j-1])]
	   + wsf*icfg[idxPS(s[i+mid],s[i+1])][idxPS(s[i+mid+1],s[j-1])]
	   + vx[(jmod-d+mid<0)? jmod-d+mid+win:jmod-d+mid][mid-1] 
	   + vx[(jmod-1<0)? jmod-1+win:jmod-1][d-mid-2]) > bestsc)
	bestsc = sc; 
      
      /* (W8.5) */
      if (d-mid > 1 &&
	  (sc = wsf*icfg[idxPS(s[i+mid],s[i])][idxPS(s[i+mid+2],s[j])]
	   + vx[(jmod-d+mid<0)? jmod-d+mid+win:jmod-d+mid][mid] 
	   + vx[jmod][d-mid-2]) > bestsc)
	bestsc = sc; 
      
      /* (W8.6) */
      if (d-mid > 1 &&
	  (sc = wsf*icfg[idxL(s[i])][idxP(s[i+1],s[i+mid])]
	   + wsf*icfg[idxPS(s[i+mid],s[i+1])][idxPS(s[i+mid+2],s[j])]
	   + vx[(jmod-d+mid<0)? jmod-d+mid+win:jmod-d+mid][mid-1] 
	   + vx[jmod][d-mid-2]) > bestsc)
	bestsc = sc; 
      
      /* (W8.7) */
      if (d-mid > 2 &&
	  (sc = wsf*icfg[idxR(s[j])][idxP(s[i+mid+2],s[j-1])]
	   + wsf*icfg[idxPS(s[i+mid],s[i])][idxPS(s[i+mid+2],s[j-1])]
	   + vx[(jmod-d+mid<0)? jmod-d+mid+win:jmod-d+mid][mid] 
	   + vx[(jmod-1<0)? jmod-1+win:jmod-1][d-mid-3]) > bestsc)
	bestsc = sc; 
      
      /* (W8.8)__/ i dangles off i+1,i+mid / j dangles off i+mid+1,j-1  */
      if (d-mid > 2 &&
	  (sc = wsf*icfg[idxL(s[i])][idxP(s[i+1],s[i+mid])]
	   + wsf*icfg[idxR(s[j])][idxP(s[i+mid+2],s[j-1])]
	   + wsf*icfg[idxPS(s[i+mid],s[i+1])][idxPS(s[i+mid+2],s[j-1])]
	   + vx[(jmod-d+mid<0)? jmod-d+mid+win:jmod-d+mid][mid-1] 
	   + vx[(jmod-1<0)? jmod-1+win:jmod-1][d-mid-3]) > bestsc)
	bestsc = sc; 
    }
  
  wx[jmod][d] = bestsc;
}




