/* 
 * filltrwbx.c 
 *
 * includes functions FillWBX,  Trace WBX
 * that fill and traceback the no-hole matrix:
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
#include "protowbx.h"
#include "squid.h"

static void FillWBX_pks(int *s, int len, int **icfg, int **wbx, int **vx, 
			int ****whx, int ****yhx, int j, int d);
static void TraceWBX_pks(FILE *outf, int *s, int len, int **icfg, int **wbx, int **vx, 
			 int ****whx, int ****yhx, int j, int d, int *flag, 
			 struct tracekn_s *curr_tr, struct traceknstack_s *dolist, int traceback);

/* Function: FillWBX()
 * 
 * Purpose:   fill no-hole matrix WBX. Using Knots-IS(2) method.
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
 *            wbx is filled.
 */       
void
FillWBX(int *s, int len, int **icfg, int **wbx, int **vx, 
	int ****whx,  int ****yhx, int j, int d, int approx)
{

  /* no internal pseudoknots approximation, read it from the zuker algorithm (WB1-WB8).
   */
  FillWBX_nested(s, len, icfg, wbx, vx, j, d);

  /* if approx = FALSE include internal pseudoknot  diagrams WB9-WB10. 
   */
  if (!approx) 
    FillWBX_pks(s, len, icfg, wbx, vx, whx, yhx, j, d);
}

/* Function: TraceWBX()
 * 
 * Purpose:  traceback wbx.
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
 *            traceback wbx.
 */       
void
TraceWBX(FILE *outf, int *s, int len, int **icfg, int **wbx, int **vx, 
	 int ****whx, int ****yhx, int j, int d, int approx,
	 struct tracekn_s *curr_tr, struct traceknstack_s *dolist, int traceback)
{
  int fl, *flag;

  fl   = FALSE;
  flag =   &fl;

  /* no internal pseudoknots approximation, read it from the zuker algorithm.
   */
  TraceWBX_nested(outf, s, len, icfg, wbx, vx, j, d, flag, curr_tr, dolist, traceback);

  /* if approx = FALSE include internal pseudoknot  diagrams WB9-WB10. 
   */
  if (!approx && *flag == FALSE){
    TraceWBX_pks(outf,s, len, icfg, wbx, vx, whx, yhx, 
		 j, d, flag, curr_tr, dolist, traceback);
  }
  /*  if(*flag == FALSE)
    Die("something went wrong in the traceback of WBX");*/
}

/* Function: FillWBX_nested()
 * 
 * Purpose:   fill no-hole matrix WBX without pseudoknots (zuker).
 *           
 * Arguments: s     - sequence (0..len-1) converted to integers
 *                    representing array indices of symbols
 *            len   - length of iseq
 *           icfg   - context-free grammar state transitions, integer log form
 *            wbx   - DP matrix, already alloc'ed 
 *             vx   - DP matrix, already alloc'ed 
 *            
 * Return:    (void)           
 *            wbx is filled.
 */       
void
FillWBX_nested(int *s, int len, int **icfg, int **wbx, int **vx, int j, int d)
{
  int                  i;  /* coords in mtx's                 */
  int                mid;  /* midpoint of bifurcation         */
  int                 sc;  /* temporary score to check        */
  int             bestsc;  /* best score so far               */

  i = j - d;

  /* (WB1)__/ PAIR (i,j).
   */
  bestsc = P10 + vx[j][d];	

  /* (WB2)__/ PAIR (i+1,j-1).
   */
  if (d > 1 &&
      (sc = P10 + 2*P6 
       + wsf*icfg[idxL(s[i])][idxP(s[i+1],s[j-1])]
       + wsf*icfg[idxR(s[j])][idxP(s[i+1],s[j-1])] 
       + vx[j-1][d-2]) > bestsc)
    bestsc = sc;
  
  /* (WB3)__/  SINGLET-LEFT  
   */
   if (d > 0 &&
       (sc = P10 + P6 
	+ wsf*icfg[idxL(s[i])][idxP(s[i+1],s[j])] 
	+ vx[j][d-1]) > bestsc)
     bestsc = sc; 
  
  /* (WB4)__/  SINGLET-RIGHT. 
   */
  if (d > 0 &&
      (sc = P10 + P6 
       + wsf*icfg[idxR(s[j])][idxP(s[i],s[j-1])] 
       + vx[j-1][d-1]) > bestsc)
    bestsc = sc; 
  
  /* (WB5)__/  SINGLET-LEFT. 
   */
  if (d > 1 &&
      (sc = P6 + wbx[j][d-1]) > bestsc)
    bestsc = sc;
  
  /* (WB6)__/  SINGLET-RIGHT.
   */
  if (d > 0 &&
      (sc = P6 + wbx[j-1][d-1]) > bestsc)
    bestsc = sc; 
  
  for (mid = 1; mid < d; mid++) 
    {    
      /* (WB7)__/  2 WBX 
       */
      if ((sc = wbx[i+mid][mid] + wbx[j][d-mid-1]) > bestsc)
	bestsc = sc; 
      
      /* (WB8)__/ 2 VX 
       */
      /* (WB8.1) */
      if ((sc = 2*P10 
	   + wsf*(icfg[idxPS(s[i+mid],s[i])][idxPS(s[i+mid+1],s[j])] + INTSCALE)
	   + vx[i+mid][mid] + vx[j][d-mid-1]) > bestsc)
	bestsc = sc; 
      
      /* (WB8.2) */
      if ((sc = 2*P10 + P6    
	   + wsf*icfg[idxL(s[i])][idxP(s[i+1],s[i+mid])]
	   + wsf*icfg[idxPS(s[i+mid],s[i+1])][idxPS(s[i+mid+1],s[j])]
	   + vx[i+mid][mid-1] + vx[j][d-mid-1]) > bestsc)
	bestsc = sc; 
      
      /* (WB8.3) */
      if (d-mid > 1 &&
	  (sc = 2*P10 + P6  
	   + wsf*icfg[idxR(s[j])][idxP(s[i+mid+1],s[j-1])]
	   + wsf*icfg[idxPS(s[i+mid],s[i])][idxPS(s[i+mid+1],s[j-1])]
	   + vx[i+mid][mid] + vx[j-1][d-mid-2]) > bestsc)
	bestsc = sc; 
      
      /* (WB8.4) */
      if (d-mid > 1 &&
	  (sc = 2*P10 + 2*P6 
	   + wsf*icfg[idxL(s[i])][idxP(s[i+1],s[i+mid])]
	   + wsf*icfg[idxR(s[j])][idxP(s[i+mid+1],s[j-1])]
	   + wsf*icfg[idxPS(s[i+mid],s[i+1])][idxPS(s[i+mid+1],s[j-1])]
	   + vx[i+mid][mid-1] + vx[j-1][d-mid-2]) > bestsc)
	bestsc = sc; 

      /* (WB8.5) */
      if (d-mid > 1 &&
	  (sc = 2*P10 + P6 
	   + wsf*icfg[idxPS(s[i+mid],s[i])][idxPS(s[i+mid+2],s[j])]
	   + vx[i+mid][mid] + vx[j][d-mid-2]) > bestsc)
	bestsc = sc; 
      
      /* (WB8.6) */
      if (d-mid > 1 &&
	  (sc = 2*P10 + 2*P6
	   + wsf*icfg[idxL(s[i])][idxP(s[i+1],s[i+mid])]
	   + wsf*icfg[idxPS(s[i+mid],s[i+1])][idxPS(s[i+mid+2],s[j])]
	   + vx[i+mid][mid-1] + vx[j][d-mid-2]) > bestsc)
	bestsc = sc; 
      
      /* (WB8.7) */
      if (d-mid > 2 &&
	  (sc = 2*P10 + 2*P6
	   + wsf*icfg[idxR(s[j])][idxP(s[i+mid+2],s[j-1])]
	   + wsf*icfg[idxPS(s[i+mid],s[i])][idxPS(s[i+mid+2],s[j-1])]
	   + vx[i+mid][mid] + vx[j-1][d-mid-3]) > bestsc)
	bestsc = sc; 
      
      /* (WB8.8)__/ i dangles off i+1,i+mid / j dangles off i+mid+1,j-1  */
      if (d-mid > 2 &&
	  (sc = 2*P10 + 3*P6
	   + wsf*icfg[idxL(s[i])][idxP(s[i+1],s[i+mid])]
	   + wsf*icfg[idxR(s[j])][idxP(s[i+mid+2],s[j-1])]
	   + wsf*icfg[idxPS(s[i+mid],s[i+1])][idxPS(s[i+mid+2],s[j-1])]
	   + vx[i+mid][mid-1] + vx[j-1][d-mid-3]) > bestsc)
	bestsc = sc; 
    }
  wbx[j][d] = bestsc;
}


/* Function: TraceWBX_nested()
 * 
 * Purpose:  traceback wbx without pseudoknots (zuker).
 *           
 * Arguments: s     - sequence (0..len-1) converted to integers
 *                    representing array indices of symbols
 *            len   - length of iseq
 *           icfg   - context-free grammar state transitions, integer log form
 *            wbx   - DP matrix, already alloc'ed 
 *             vx   - DP matrix, already alloc'ed 
 *            
 * Return:    (void)           
 *            traceback wbx.
 */       
void
TraceWBX_nested(FILE *outf, int *s, int len, int **icfg, 
		int **wbx, int **vx, int j, int d, int *flag, 
		struct tracekn_s *curr_tr, struct traceknstack_s *dolist, int traceback)
{
  int mid;  /* used for midpoint of a bifurc   */

  /* (WB1) */
  if (wbx[j][d] ==  WB1(s, len, icfg, vx, j, d)){
    trace_WB1(outf, vx, j, d, curr_tr, dolist, traceback);
    *flag = TRUE;
    return;
  } 
  
  /* (WB2) */
  else if (wbx[j][d] == WB2(s, len, icfg, vx, j, d)){
    trace_WB2(outf, vx, j, d, curr_tr, dolist, traceback);
    *flag = TRUE;
    return;
  } 
  
  /* (WB3) */
  else if (wbx[j][d] == WB3(s, len, icfg, vx, j, d)){
    trace_WB3(outf, vx, j, d, curr_tr, dolist, traceback);
    *flag = TRUE;
    return;
  } 
  
  /* (WB4) */
  else if (wbx[j][d] == WB4(s, len, icfg, vx, j, d)){
    trace_WB4(outf, vx, j, d, curr_tr, dolist, traceback);
    *flag = TRUE;
    return;
  } 
  
  /* (WB5) */
  else if (wbx[j][d] == WB5(s, len, icfg, wbx, j, d)){
    trace_WB5(outf, wbx, j, d, curr_tr, dolist, traceback);
    *flag = TRUE;
    return;
  } 
  
  /* (WB6) */
  else if (wbx[j][d] == WB6(s, len, icfg, wbx, j, d)){
    trace_WB6(outf, wbx, j, d, curr_tr, dolist, traceback);
    *flag = TRUE;
    return;
  } 
  
  for (mid = d-1; mid > 0; mid--) {
    /* (WB7) */         
    if (wbx[j][d] == WB7(s, len, icfg, wbx, j, d, mid)){
      trace_WB7(outf, wbx, j, d, mid, curr_tr, dolist, traceback);
      *flag = TRUE;
      return;
    } 
   
    /* (WB8) */
    /* (WB8.1) */
    else if (wbx[j][d] == WB8_1(s, len, icfg, vx, j, d, mid)){
      trace_WB8_1(outf, vx, j, d, mid, curr_tr, dolist, traceback);
      *flag = TRUE;
      return;
    } 

    /* (WB8.2) */
    else if (wbx[j][d] == WB8_2(s, len, icfg, vx, j, d, mid)){
      trace_WB8_2(outf, vx, j, d, mid, curr_tr, dolist, traceback);
      *flag = TRUE;
      return;
    } 

    /* (WB8.3) */
    else if (wbx[j][d] == WB8_3(s, len, icfg, vx, j, d, mid)){
      trace_WB8_3(outf, vx, j, d, mid, curr_tr, dolist, traceback);
      *flag = TRUE;
      return;
    } 

    /* (WB8.4) */
    else if (wbx[j][d] == WB8_4(s, len, icfg, vx, j, d, mid)){
      trace_WB8_4(outf, vx, j, d, mid, curr_tr, dolist, traceback);
      *flag = TRUE;
      return;
    } 

    /* (WB8.5) */
    else if (wbx[j][d] == WB8_5(s, len, icfg, vx, j, d, mid)){
      trace_WB8_5(outf, vx, j, d, mid, curr_tr, dolist, traceback);
      *flag = TRUE;
      return;
    } 

    /* (WB8.6) */
    else if (wbx[j][d] == WB8_6(s, len, icfg, vx, j, d, mid)){
      trace_WB8_6(outf, vx, j, d, mid, curr_tr, dolist, traceback);
      *flag = TRUE;
      return;
    } 

    /* (WB8.7) */
    else if (wbx[j][d] == WB8_7(s, len, icfg, vx, j, d, mid)){
      trace_WB8_7(outf, vx, j, d, mid, curr_tr, dolist, traceback);
      *flag = TRUE;
      return;
    } 

    /* (WB8.8) */
    else if (wbx[j][d] == WB8_8(s, len, icfg, vx, j, d, mid)){
      trace_WB8_8(outf, vx, j, d, mid, curr_tr, dolist, traceback);
      *flag = TRUE;
      return;
    } 
  }
}

/* Function: FillWBX_pks()
 * 
 * Purpose: fill internal pseudoknots in wbx. Using Knots-IS(2) method.
 *
 * BIFURCATIONS (with hole) (2 types of diagrams WB9-WB10)
 * (add P13 penalty for wbx)
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
 *            wbx   - DP matrix, already alloc'ed 
 *             vx   - DP matrix, already alloc'ed 
 *            
 * Return:    (void)           
 *            structures WB9-WB10 of vx are filled.
 *
 */
void
FillWBX_pks(int *s, int len, int **icfg, int **wbx, int **vx, 
	    int ****whx, int ****yhx, int j, int d)
{
  int i;
  int mid, mid1, mid2;
  int sc, bestsc;

  i = j - d;

  bestsc = wbx[j][d];

  for (mid = 3; mid < (d-2); mid++)            
    for (mid1 = 0; mid1 < (mid-1); mid1++)                
      for (mid2 = 0; mid2 < (mid-mid1-2); mid2++)
	{    
	  /* (WB9)__/ 2 WHX
	   */
	  if ((sc = P13 + 3*P6P
	       + whx[i+mid][mid][mid1][mid2] 
	       + whx[j][d-mid1-2][mid-mid1-mid2-3][d-mid-3]) > bestsc)
	    bestsc = sc; 
	  
	  /* (WB10)__/ 1 YHX 1 YHX
	   */
	  /* (WB10.1) */
	  if ((sc = P13 + 2*P10P + 3*P6P
	       + wkn*icfg[idxR(s[i+mid1+1])][idxP(s[i+mid-mid2],s[i+mid1])]     
	       + wkn*icfg[idxL(s[i+mid+2])][idxP(s[i+mid+3],s[i+mid-mid2-1])]     
	       + wkn*icfg[idxPS(s[i+mid-mid2-1],s[i+mid+3])]
	                 [idxPS(s[i+mid-mid2],s[i+mid1])] 
	       + yhx[i+mid][mid][mid1][mid2] 
	       + yhx[j][d-mid1-2][mid-mid1-mid2-3][d-mid-3]) > bestsc)
	    bestsc = sc; 
	  
	  /* (WB10.2) */
	  if (mid-mid1-mid2 > 3 &&
	      (sc = P13 + 2*P10P + 4*P6P
	       + wkn*icfg[idxR(s[i+mid1+1])][idxP(s[i+mid-mid2],s[i+mid1])]     
	       + wkn*icfg[idxL(s[i+mid+2])][idxP(s[i+mid+3],s[i+mid-mid2-2])]     
	       + wkn*icfg[idxPS(s[i+mid-mid2-2],s[i+mid+3])]
	       [idxPS(s[i+mid-mid2],s[i+mid1])] 
	       + yhx[i+mid][mid][mid1][mid2] 
	       + yhx[j][d-mid1-2][mid-mid1-mid2-4][d-mid-3]) > bestsc)
	    bestsc = sc; 
	}
  wbx[j][d] = bestsc;
}

/* Function: TraceWBX_pks()
 * 
 * Purpose: traceback internal pseudoknots in wbx. Using Knots-IS(2) method.
 *
 * BIFURCATIONS (with hole) (2 types of diagrams WB9-WB10)
 * (add P13 penalty for wbx)
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
 *            wbx   - DP matrix, already alloc'ed 
 *             vx   - DP matrix, already alloc'ed 
 *            
 * Return:    (void)           
 *            structures WB9-WB13 of vx are traced back.
 *
 */
void
TraceWBX_pks(FILE *outf, int *s, int len, int **icfg, int **wbx, int **vx, 
	     int ****whx, int ****yhx, int j, int d, 
	     int *flag, struct tracekn_s *curr_tr, struct traceknstack_s *dolist, int traceback)
{
  int mid, mid1, mid2;

  for (mid = d-3; mid >= 3; mid--)
    for (mid1 = mid-2; mid1 >= 0; mid1--)   
      for (mid2 = mid-mid1-3; mid2 >= 0; mid2--) {

	/* (WB9)__/  2 WHX.
	 */
	if (wbx[j][d] == WB9(s, len, icfg, whx, j, d, mid, mid1, mid2)){
	  trace_WB9(outf, whx, j, d, mid, mid1, mid2, curr_tr, dolist, traceback);
	  *flag = TRUE;
	  return;
 	} 
	
	/* (WB10)__/ YHX 1 YHX
	 */
	/* (WB10.1) */
	else if (wbx[j][d] == WB10_1(s, len, icfg, yhx, j, d, mid, mid1, mid2)){
	  trace_WB10_1(outf, yhx, j, d, mid, mid1, mid2, curr_tr, dolist, traceback);
	  *flag = TRUE;
	  return;
	} 

	/* (WB10.2) */
	else if (wbx[j][d] == WB10_2(s, len, icfg, yhx, j, d, mid, mid1, mid2)){
	  trace_WB10_2(outf, yhx, j, d, mid, mid1, mid2, curr_tr, dolist, traceback);
	  *flag = TRUE;
	  return;
	} 
      }
}

