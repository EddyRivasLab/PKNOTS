/* 
 * filltrwx.c 
 *
 * includes functions FillWX,  TraceWX.
 * that fill and traceback the no-hole matrix:
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
#include "protowx.h"
#include "squid.h"

static void FillWX_pks(int *s, int len, int **icfg, int **wx, int **vx, 
		       int ****whx, int ****yhx, int j, int d);
static void TraceWX_pks(FILE *outf, int *s, int len, int **icfg, int **wx, int **vx, 
			int ****whx, int ****yhx, int j, int d, int *flag,
			struct tracekn_s *curr_tr, struct traceknstack_s *dolist, int traceback);

/* Function: FillWX()
 * 
 * Purpose:   fill no-hole matrix WX. Using Knots-IS(2) method.
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
 *            wx is filled.
 */       
void
FillWX(int *s, int len, int **icfg, int **wx, int **vx, 
       int ****whx, int ****yhx, int j, int d)
{

  /* no pseudoknots. (includes: W1-W8).
   */
  FillWX_nested(s, len, icfg, wx, vx, j, d);

  /* pseudoknot diagrams (W9-W10).
   */
  FillWX_pks(s, len, icfg, wx, vx, whx, yhx, j, d);
}

/* Function: TraceWX()
 * 
 * Purpose:  traceback wx.
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
 *             vx   - DP matrix, already alloc'ed 
 *            
 * Return:    (void)           
 *            traceback wx.
 */       
void
TraceWX(FILE *outf, int *s, int len, int **icfg, int **wx, int **vx, 
	int ****whx, int ****yhx, int j, int d,  
	struct tracekn_s *curr_tr, struct traceknstack_s *dolist, int traceback)
{
  int fl, *flag;

  fl   = FALSE;
  flag =   &fl;

  /* no pseudoknots. (includes: W1-W8).
   */
  TraceWX_nested(outf, s, len, icfg, wx, vx, j, d, flag, curr_tr, dolist, traceback);

  /* pseudoknot diagrams (W9-W12).
   */
  if (*flag == FALSE)
    TraceWX_pks(outf, s, len, icfg, wx, vx, whx, yhx, 
		j, d, flag, curr_tr, dolist, traceback);
  

  if(*flag == FALSE)
    Die("something went wrong in the traceback of WX");
}

/* Function: FillWX_nested()
 * 
 * Purpose:   fill no-hole matrix WX without pseudoknots (zuker).
 *           
 * Arguments: s     - sequence (0..len-1) converted to integers
 *                    representing array indices of symbols
 *            len   - length of iseq
 *           icfg   - context-free grammar state transitions, integer log form
 *             wx   - DP matrix, already alloc'ed 
 *             vx   - DP matrix, already alloc'ed 
 *            
 * Return:    (void)           
 *            wx is filled.
 */       
void
FillWX_nested(int *s, int len, int **icfg, int **wx, int **vx, int j, int d)
{
  int                  i;  /* coords in mtx's                 */
  int                mid;  /* midpoint of bifurcation         */
  int                 sc;  /* temporary score to check        */
  int             bestsc;  /* best score so far               */

  i = j - d;

  /* (W1)__/ PAIR (i,j).
   */
  bestsc = vx[j][d];	

  /* (W2)__/ PAIR (i+1,j-1).
   */
  if (d > 1 &&
      (sc = wsf*icfg[idxL(s[i])][idxP(s[i+1],s[j-1])]
       + wsf*icfg[idxR(s[j])][idxP(s[i+1],s[j-1])] 
       + vx[j-1][d-2]) > bestsc)
    bestsc = sc;
  
  /* (W3)__/  SINGLET-LEFT  
   */
   if (d > 0 &&
       (sc = wsf*icfg[idxL(s[i])][idxP(s[i+1],s[j])] 
	+ vx[j][d-1]) > bestsc)
     bestsc = sc; 
  
  /* (W4)__/  SINGLET-RIGHT. 
   */
  if (d > 0 &&
      (sc = wsf*icfg[idxR(s[j])][idxP(s[i],s[j-1])] 
       + vx[j-1][d-1]) > bestsc)
    bestsc = sc; 
  
  /* (W5)__/  SINGLET-LEFT. 
   */
  if (d > 1 &&
      (sc = wx[j][d-1]) > bestsc)
    bestsc = sc;
  
  /* (W6)__/  SINGLET-RIGHT.
   */
  if (d > 0 &&
      (sc = wx[j-1][d-1]) > bestsc)
    bestsc = sc; 
  
  for (mid = 1; mid < d; mid++) 
    {    
      /* (W7)__/  2 WX 
       */
      if ((sc = wx[i+mid][mid] + wx[j][d-mid-1]) > bestsc)
	bestsc = sc; 
      
      /* (W8)__/ 2 VX 
       */
      /* (W8.1) */
      if ((sc = wsf*(icfg[idxPS(s[i+mid],s[i])][idxPS(s[i+mid+1],s[j])] + INTSCALE)
	   + vx[i+mid][mid] + vx[j][d-mid-1]) > bestsc)
	bestsc = sc; 
      
      /* (W8.2) */
      if ((sc = wsf*icfg[idxL(s[i])][idxP(s[i+1],s[i+mid])]
	   + wsf*icfg[idxPS(s[i+mid],s[i+1])][idxPS(s[i+mid+1],s[j])]
	   + vx[i+mid][mid-1] + vx[j][d-mid-1]) > bestsc)
	bestsc = sc; 
      
      /* (W8.3) */
      if (d-mid > 1 &&
	  (sc = wsf*icfg[idxR(s[j])][idxP(s[i+mid+1],s[j-1])]
	   + wsf*icfg[idxPS(s[i+mid],s[i])][idxPS(s[i+mid+1],s[j-1])]
	   + vx[i+mid][mid] + vx[j-1][d-mid-2]) > bestsc)
	bestsc = sc; 
      
      /* (W8.4)__/ i dangles off i+1,i+mid / j dangles off i+mid+1,j-1  */
      if (d-mid > 1 &&
	  (sc = wsf*icfg[idxL(s[i])][idxP(s[i+1],s[i+mid])]
	   + wsf*icfg[idxR(s[j])][idxP(s[i+mid+1],s[j-1])]
	   + wsf*icfg[idxPS(s[i+mid],s[i+1])][idxPS(s[i+mid+1],s[j-1])]
	   + vx[i+mid][mid-1] + vx[j-1][d-mid-2]) > bestsc)
	bestsc = sc; 
      
      /* (W8.5) */
      if (d-mid > 1 &&
	  (sc = wsf*icfg[idxPS(s[i+mid],s[i])][idxPS(s[i+mid+2],s[j])]
	   + vx[i+mid][mid] + vx[j][d-mid-2]) > bestsc)
	bestsc = sc; 
      
      /* (W8.6) */
      if (d-mid > 1 &&
	  (sc = wsf*icfg[idxL(s[i])][idxP(s[i+1],s[i+mid])]
	   + wsf*icfg[idxPS(s[i+mid],s[i+1])][idxPS(s[i+mid+2],s[j])]
	   + vx[i+mid][mid-1] + vx[j][d-mid-2]) > bestsc)
	bestsc = sc; 
      
      /* (W8.7) */
      if (d-mid > 2 &&
	  (sc = wsf*icfg[idxR(s[j])][idxP(s[i+mid+2],s[j-1])]
	   + wsf*icfg[idxPS(s[i+mid],s[i])][idxPS(s[i+mid+2],s[j-1])]
	   + vx[i+mid][mid] + vx[j-1][d-mid-3]) > bestsc)
	bestsc = sc; 
      
      /* (W8.8)__/ i dangles off i+1,i+mid / j dangles off i+mid+1,j-1  */
    if (d-mid > 2 &&
	  (sc = wsf*icfg[idxL(s[i])][idxP(s[i+1],s[i+mid])]
	   + wsf*icfg[idxR(s[j])][idxP(s[i+mid+2],s[j-1])]
	   + wsf*icfg[idxPS(s[i+mid],s[i+1])][idxPS(s[i+mid+2],s[j-1])]
	   + vx[i+mid][mid-1] + vx[j-1][d-mid-3]) > bestsc)
      bestsc = sc; 
    }
  wx[j][d] = bestsc;
}


/* Function: TraceWX_nested()
 * 
 * Purpose:  traceback wx without pseudoknots (zuker).
 *           
 * Arguments: s     - sequence (0..len-1) converted to integers
 *                    representing array indices of symbols
 *            len   - length of iseq
 *           icfg   - context-free grammar state transitions, integer log form
 *             wx   - DP matrix, already alloc'ed 
 *             vx   - DP matrix, already alloc'ed 
 *            
 * Return:    (void)           
 *            traceback wx.
 */       
void
TraceWX_nested(FILE *outf, int *s, int len, int **icfg, 
	       int **wx, int **vx, int j, int d, int *flag, 
	       struct tracekn_s *curr_tr, struct traceknstack_s *dolist, int traceback)
{
  int mid;  /* used for midpoint of a bifurc   */

  /* (W1) */
  if (wx[j][d] ==  W1(s, len, icfg, vx, j, d)){
    trace_W1(outf, vx, j, d, curr_tr, dolist, traceback);
    *flag = TRUE;
    return;
  } 
  
  /* (W2) */
  else if (wx[j][d] == W2(s, len, icfg, vx, j, d)){
    trace_W2(outf, vx, j, d, curr_tr, dolist, traceback);
    *flag = TRUE;
    return;
  } 
  
  /* (W3) */
  else if (wx[j][d] == W3(s, len, icfg, vx, j, d)){
    trace_W3(outf, vx, j, d, curr_tr, dolist, traceback);
    *flag = TRUE;
    return;
  } 
  
  /* (W4) */
  else if (wx[j][d] == W4(s, len, icfg, vx, j, d)){
    trace_W4(outf, vx, j, d, curr_tr, dolist, traceback);
    *flag = TRUE;
    return;
  } 
  
  /* (W5) */
  else if (wx[j][d] == W5(s, len, icfg, wx, j, d)){
    trace_W5(outf, wx, j, d, curr_tr, dolist, traceback);
    *flag = TRUE;
    return;
  } 
  
  /* (W6) */
  else if (wx[j][d] == W6(s, len, icfg, wx, j, d)){
    trace_W6(outf, wx, j, d, curr_tr, dolist, traceback);
    *flag = TRUE;
    return;
  } 
  
  for (mid = d-1; mid > 0; mid--) {
    /* (W7) */         
    if (wx[j][d] == W7(s, len, icfg, wx, j, d, mid)){
      trace_W7(outf, wx, j, d, mid, curr_tr, dolist, traceback);
      *flag = TRUE;
      return;
    } 
    
    /* (W8) */
    /* (W8.1) */
    else if (wx[j][d] == W8_1(s, len, icfg, vx, j, d, mid)){
      trace_W8_1(outf, vx, j, d, mid, curr_tr, dolist, traceback);
      *flag = TRUE;
      return;
    } 

    /* (W8.2) */
    else if (wx[j][d] == W8_2(s, len, icfg, vx, j, d, mid)){
      trace_W8_2(outf, vx, j, d, mid, curr_tr, dolist, traceback);
      *flag = TRUE;
      return;
    } 

    /* (W8.3) */
    else if (wx[j][d] == W8_3(s, len, icfg, vx, j, d, mid)){
      trace_W8_3(outf, vx, j, d, mid, curr_tr, dolist, traceback);
      *flag = TRUE;
      return;
    } 

    /* (W8.4) */
    else if (wx[j][d] == W8_4(s, len, icfg, vx, j, d, mid)){
      trace_W8_4(outf, vx, j, d, mid, curr_tr, dolist, traceback);
      *flag = TRUE;
      return;
    } 

    /* (W8.5) */
    else if (wx[j][d] == W8_5(s, len, icfg, vx, j, d, mid)){
      trace_W8_5(outf, vx, j, d, mid, curr_tr, dolist, traceback);
      *flag = TRUE;
      return;
    } 

    /* (W8.6) */
    else if (wx[j][d] == W8_6(s, len, icfg, vx, j, d, mid)){
      trace_W8_6(outf, vx, j, d, mid, curr_tr, dolist, traceback);
      *flag = TRUE;
      return;
    } 

    /* (W8.7) */
    else if (wx[j][d] == W8_7(s, len, icfg, vx, j, d, mid)){
      trace_W8_7(outf, vx, j, d, mid, curr_tr, dolist, traceback);
      *flag = TRUE;
      return;
    } 

    /* (W8.8) */
    else if (wx[j][d] == W8_8(s, len, icfg, vx, j, d, mid)){
      trace_W8_8(outf, vx, j, d, mid, curr_tr, dolist, traceback);
      *flag = TRUE;
      return;
    } 
  }
}

/* Function: FillWX_pks()
 * 
 * Purpose: fill internal pseudoknots in wx. Using Knots-IS(2) method.
 *
 * BIFURCATIONS (with hole) (2 types of diagrams W9-W10)
 * (add P13 penalty for wx)
 *
 *     one connects (i,i+mid1) and (i+mid-mid2,i+mid). 
 *                  (i+mid,mid,mid1,mid2).
 *
 * Another connects (i+mid1+1,i+mid-mid2-1) and (i+mid+1,j). 
 *                  (j,d-mid1-1,mid-mid1-mid2-2,d-mid-1).	
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
 *             vx   - DP matrix, already alloc'ed 
 *            
 * Return:    (void)           
 *            structures W9-W12 of vx are filled.
 *
 */
void
FillWX_pks(int *s, int len, int **icfg, int **wx, int **vx, 
	   int ****whx, int ****yhx, int j, int d)
{
  int i;
  int mid, mid1, mid2;
  int sc, bestsc;

  i = j - d;

  bestsc = wx[j][d];

  for (mid = 3; mid < (d-2); mid++)            
    for (mid1 = 0; mid1 < (mid-1); mid1++)                
      for (mid2 = 0; mid2 < (mid-mid1-2); mid2++)
	{    
	  /* (W9)__/ 2 WHX
	   */
	  if ((sc = P11 + 3*P6P
	       + whx[i+mid][mid][mid1][mid2] 
	       + whx[j][d-mid1-2][mid-mid1-mid2-3][d-mid-3]) > bestsc)
	    bestsc = sc; 
	  
	  /* (W10)__/ 1 YHX 1 YHX
	   */
	  /* (W10.1) */
	  if ((sc =  P11 + 2*P10P + 3*P6P
	       + wkn*icfg[idxR(s[i+mid1+1])][idxP(s[i+mid-mid2],s[i+mid1])]     
	       + wkn*icfg[idxL(s[i+mid+2])][idxP(s[i+mid+3],s[i+mid-mid2-1])]     
	       + wkn*icfg[idxPS(s[i+mid-mid2-1],s[i+mid+3])]
	                 [idxPS(s[i+mid-mid2],s[i+mid1])] 
	       + yhx[i+mid][mid][mid1][mid2] 
	       + yhx[j][d-mid1-2][mid-mid1-mid2-3][d-mid-3]) > bestsc)
	    bestsc = sc; 
	  
	  /* (W10.2) */
	  if (mid-mid1-mid2 > 3 &&
	      (sc = P11 + 2*P10P + 4*P6P
	       + wkn*icfg[idxR(s[i+mid1+1])][idxP(s[i+mid-mid2],s[i+mid1])]     
	       + wkn*icfg[idxL(s[i+mid+2])][idxP(s[i+mid+3],s[i+mid-mid2-2])]     
	       + wkn*icfg[idxPS(s[i+mid-mid2-2],s[i+mid+3])]
	                 [idxPS(s[i+mid-mid2],s[i+mid1])] 
	       + yhx[i+mid][mid][mid1][mid2] 
	       + yhx[j][d-mid1-2][mid-mid1-mid2-4][d-mid-3]) > bestsc)
	    bestsc = sc; 
	}
  wx[j][d] = bestsc;
}

/* Function: TraceWX_pks()
 * 
 * Purpose: traceback internal pseudoknots in wx. Using Knots-IS(2) method.
 *
 * BIFURCATIONS (with hole) (2 types of diagrams W9-W10)
 * (add P11 penalty for wx)
 *
 *     one connects (i,i+mid1) and (i+mid-mid2,i+mid). 
 *                  (i+mid,mid,mid1,mid2).
 *
 * Another connects (i+mid1+1,i+mid-mid2-1) and (i+mid+1,j). 
 *                  (j,d-mid1-1,mid-mid1-mid2-2,d-mid-1).	
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
 *             vx   - DP matrix, already alloc'ed 
 *            
 * Return:    (void)           
 *            structures W9-W13 of vx are traced back.
 *
 */
void
TraceWX_pks(FILE *outf, int *s, int len, int **icfg, int **wx, int **vx, 
	    int ****whx, int ****yhx, int j, int d, int *flag,
	    struct tracekn_s *curr_tr, struct traceknstack_s *dolist, int traceback)
{
  int mid, mid1, mid2;

  for (mid = d-3; mid >= 3; mid--)
    for (mid1 = mid-2; mid1 >= 0; mid1--)   
      for (mid2 = mid-mid1-3; mid2 >= 0; mid2--) {

	/* (W9)__/  2 WHX.
	 */
	if (wx[j][d] == W9(s, len, icfg, whx, j, d, mid, mid1, mid2)){
	  trace_W9(outf, whx, j, d, mid, mid1, mid2, curr_tr, dolist, traceback);
	  *flag = TRUE;
	  return;
	} 
	
	/* (W10)__/ YHX 1 YHX
	 */
	/* (W10.1) */
	else if (wx[j][d] == W10_1(s, len, icfg, yhx, j, d, mid, mid1, mid2)){
	  trace_W10_1(outf, yhx, j, d, mid, mid1, mid2, curr_tr, dolist, traceback);
	  *flag = TRUE;
	  return;
	} 

	/* (W10.2) */
	else if (wx[j][d] == W10_2(s, len, icfg, yhx, j, d, mid, mid1, mid2)){
	  trace_W10_2(outf, yhx, j, d, mid, mid1, mid2, curr_tr, dolist, traceback);
	  *flag = TRUE;
	  return;
	} 
      }
}

