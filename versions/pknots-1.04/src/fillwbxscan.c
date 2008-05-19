/* 
 * fillwbxscan.c
 *
 * includes functions FillWBX
 * fills the no-hole matrix:
 *           wbx[j][d] (len x j) [(j-d,j) are     inside a base paired structure]
 *                     [(j-d,j) not necessarely based paired]
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

/* Function: FillWBX_nestedScan()
 * 
 * Purpose:   fill no-hole matrix WBX without pseudoknots.
 *           
 * Arguments: s     - sequence (0..len-1) converted to integers
 *                    representing array indices of symbols
 *            len   - length of iseq
 *            cfg   - context-free grammar state transitions, float log2 form
 *            wbx   - DP matrix, already alloc'ed 
 *             vx   - DP matrix, already alloc'ed 
 *            
 * Return:    (void)           
 *            wbx is filled.
 */       

void
FillWBX_nestedScan(int *s, int len, int win, int **icfg, int **wbx, int **vx, 
		   int j, int jmod, int d)
{
  int                  i;  /* coords in mtx's                 */
  int                mid;  /* midpoint of bifurcation         */
  int                 sc;  /* temporary score to check        */
  int             bestsc;  /* temporary best score            */

  i = j - d;

  /* (WB1)__/ PAIR (i,j).
   */
  bestsc = P10 + vx[jmod][d];	

  /* (WB2)__/ PAIR (i+1,j-1).
   */
  if (d > 1 &&
      (sc = P10 + 2*P6 
       + wsf*icfg[idxL(s[i])][idxP(s[i+1],s[j-1])]
       + wsf*icfg[idxR(s[j])][idxP(s[i+1],s[j-1])] 
       + vx[(jmod-1<0)? jmod-1+win:jmod-1][d-2]) > bestsc)
    bestsc = sc;
  
  /* (WB3)__/  SINGLET-LEFT  
   */
   if (d > 0 &&
       (sc = P10 + P6 
	+ wsf*icfg[idxL(s[i])][idxP(s[i+1],s[j])] 
	+ vx[jmod][d-1]) > bestsc)
     bestsc = sc; 
  
  /* (WB4)__/  SINGLET-RIGHT. 
   */
  if (d > 0 &&
      (sc = P10 + P6 
       + wsf*icfg[idxR(s[j])][idxP(s[i],s[j-1])] 
       + vx[(jmod-1<0)? jmod-1+win:jmod-1][d-1]) > bestsc)
    bestsc = sc; 
  
  /* (WB5)__/  SINGLET-LEFT. 
   */
  if (d > 1 &&
      (sc = P6 + wbx[jmod][d-1]) > bestsc)
    bestsc = sc;
  
  /* (WB6)__/  SINGLET-RIGHT.
   */
  if (d > 0 &&
      (sc = P6 + wbx[(jmod-1<0)? jmod-1+win:jmod-1][d-1]) > bestsc)
    bestsc = sc; 
  
  for (mid = 1; mid < d; mid++) 
    {    
      /* (WB7)__/  2 WBX 
       */
      if ((sc = wbx[(jmod-d+mid<0)? jmod-d+mid+win:jmod-d+mid][mid] 
	   + wbx[jmod][d-mid-1]) > bestsc)
	bestsc = sc; 
      
      /* (WB8)__/ 2 VX 
       */
      /* (WB8.1) */
      if ((sc = 2*P10 
	   + wsf*(icfg[idxPS(s[i+mid],s[i])][idxPS(s[i+mid+1],s[j])] + INTSCALE)
	   + vx[(jmod-d+mid<0)? jmod-d+mid+win:jmod-d+mid][mid] 
	   + vx[jmod][d-mid-1]) > bestsc)
	bestsc = sc; 
      
      /* (WB8.2) */
      if ((sc = 2*P10 + P6    
	   + wsf*icfg[idxL(s[i])][idxP(s[i+1],s[i+mid])]
	   + wsf*icfg[idxPS(s[i+mid],s[i+1])][idxPS(s[i+mid+1],s[j])]
	   + vx[(jmod-d+mid<0)? jmod-d+mid+win:jmod-d+mid][mid-1] 
	   + vx[jmod][d-mid-1]) > bestsc)
	bestsc = sc; 
      
      /* (WB8.3) */
      if (d-mid > 1 &&
	  (sc = 2*P10 + P6  
	   + wsf*icfg[idxR(s[j])][idxP(s[i+mid+1],s[j-1])]
	   + wsf*icfg[idxPS(s[i+mid],s[i])][idxPS(s[i+mid+1],s[j-1])]
	   + vx[(jmod-d+mid<0)? jmod-d+mid+win:jmod-d+mid][mid] 
	   + vx[(jmod-1<0)? jmod-1+win:jmod-1][d-mid-2]) > bestsc)
	bestsc = sc; 
      
      /* (WB8.4) */
      if (d-mid > 1 &&
	  (sc = 2*P10 + 2*P6 
	   + wsf*icfg[idxL(s[i])][idxP(s[i+1],s[i+mid])]
	   + wsf*icfg[idxR(s[j])][idxP(s[i+mid+1],s[j-1])]
	   + wsf*icfg[idxPS(s[i+mid],s[i+1])][idxPS(s[i+mid+1],s[j-1])]
	   + vx[(jmod-d+mid<0)? jmod-d+mid+win:jmod-d+mid][mid-1] 
	   + vx[(jmod-1<0)? jmod-1+win:jmod-1][d-mid-2]) > bestsc)
	bestsc = sc; 

      /* (WB8.5) */
      if (d-mid > 1 &&
	  (sc = 2*P10 + P6 
	   + wsf*icfg[idxPS(s[i+mid],s[i])][idxPS(s[i+mid+2],s[j])]
	   + vx[(jmod-d+mid<0)? jmod-d+mid+win:jmod-d+mid][mid] 
	   + vx[jmod][d-mid-2]) > bestsc)
	bestsc = sc; 
      
      /* (WB8.6) */
      if (d-mid > 1 &&
	  (sc = 2*P10 + 2*P6
	   + wsf*icfg[idxL(s[i])][idxP(s[i+1],s[i+mid])]
	   + wsf*icfg[idxPS(s[i+mid],s[i+1])][idxPS(s[i+mid+2],s[j])]
	   + vx[(jmod-d+mid<0)? jmod-d+mid+win:jmod-d+mid][mid-1] 
	   + vx[jmod][d-mid-2]) > bestsc)
	bestsc = sc; 
      
      /* (WB8.7) */
      if (d-mid > 2 &&
	  (sc = 2*P10 + 2*P6
	   + wsf*icfg[idxR(s[j])][idxP(s[i+mid+2],s[j-1])]
	   + wsf*icfg[idxPS(s[i+mid],s[i])][idxPS(s[i+mid+2],s[j-1])]
	   + vx[(jmod-d+mid<0)? jmod-d+mid+win:jmod-d+mid][mid] 
	   + vx[(jmod-1<0)? jmod-1+win:jmod-1][d-mid-3]) > bestsc)
	bestsc = sc; 
      
      /* (WB8.8)__/ i dangles off i+1,i+mid / j dangles off i+mid+1,j-1  */
      if (d-mid > 2 &&
	  (sc = 2*P10 + 3*P6
	   + wsf*icfg[idxL(s[i])][idxP(s[i+1],s[i+mid])]
	   + wsf*icfg[idxR(s[j])][idxP(s[i+mid+2],s[j-1])]
	   + wsf*icfg[idxPS(s[i+mid],s[i+1])][idxPS(s[i+mid+2],s[j-1])]
	   + vx[(jmod-d+mid<0)? jmod-d+mid+win:jmod-d+mid][mid-1] 
	   + vx[(jmod-1<0)? jmod-1+win:jmod-1][d-mid-3]) > bestsc)
	bestsc = sc; 
    }
  
  wbx[jmod][d] = bestsc;
}



