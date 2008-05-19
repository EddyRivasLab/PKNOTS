/* tracemtx.c 
 *
 * traces back the best score for
 * the no-hole matrices:  VX, WX, WBX.         Dimension  (len x len)
 * and the hole matrices: VHX, ZHX, YHX, WHX.  Dimension  (len x len x len x len)
 *
 */
                                      
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "cfg.h"
#include "proto.h"
#include "squid.h"


/* Function: TraceMtx_nested()
 * 
 * Purpose:  Recover an RNA structure as a tree of trace_s structures,
 *           by tracing back from some point in the matrix j, d, ROOT.
 *           emitl and emitr in the trace structure are 0..N-1 positions
 *           in the sequence; type is dpcP, dpcL, etc, (state type).
 *           
 * Args:     s            - sequence, in integer form (A=0, C=1, G=2, U=3)
 *           len          - length of s, which is 0..len-1
 *           icfg         - the SCFG
 *           wx, wbx, vx  - the matrices
 *           j            - end point of the match in vertical
 *           d            - end point of the match in horizontal
 *
 *           ret_trace - RETURN: traceback tree. No pseudoknots
 *           
 * Return:   void. ret_trace is allocated here and must be free'd by caller,
 *           using FreeTracekn().
 */
void
TraceMtx_nested(FILE *outf, int *s, int len, int **icfg, 
		int **wx, int **wbx, int **vx, int j, int d,
		struct tracekn_s **ret_trace, int traceback)
{
  struct tracekn_s          *tr;  /* the traceback tree under construction  */
  struct tracekn_s     *curr_tr;  /* ptr to node of tr we're working on     */
  struct traceknstack_s *dolist;  /* pushdown stack of active tr nodes      */
  int                     i,k,l;  /* coords in mtx's                        */
  int                    d1, d2;  /* coords in mtx's                        */
  int                   typemtx;  /* type of matrix : VX, WX, WBX           */
  int                 fl, *flag;  /* flags for the traceback                */
                                     
  /* Initialize.
   * Start at j, d, d1 = (int)(d/2), d2 = d - (int)(d/2) - 1.
   */
  tr     = InitTracekn();       /* start a trace tree */
  dolist = InitTraceknstack();	/* start a stack for traversing the trace tree */

  if (len == 1) d = 0;  /* redefine d for len == 1 */
  else if (d < 0)
    Die("check your traceback assingments");

  
  i = j - d;
  k = i + (int)(d/2);
  l = k + 1;
  
  
  fl   = FALSE;
  flag =   &fl;
  if (traceback) fprintf(outf,"---------------------------------------------------\n");

 /* assigned state to (i,j): dpcS = (relation unknown, trace wx)
   *                          dpcB = (relation unknown, trace wbx)
   *                          dpcP = (paired)
   *                          dpcL = (unpaired)
   */

  curr_tr = AttachTracekn(tr, i, j, k, l, dpcS, dpcS); 
  PushTraceknstack(dolist, curr_tr);

  /* Recursion. While there's active nodes in the stack, trace from them.
   * 
   * This is cribbed from FillMtx_nested(); it's almost the exact reverse.
   * We know the best score, we just have to figure out where it came from.
   */
  while ((curr_tr = PopTraceknstack(dolist)) != NULL)
   {
     /* get some useful numbers, mostly for clarity */
      i = curr_tr->emiti;
      j = curr_tr->emitj;
      k = curr_tr->emitk;
      l = curr_tr->emitl;
      
      d  = j - i;
      
      if (d == 0 ) continue;

      d1 = (int)(d/2);
      d2 = d - d1 - 1;
      
      k = i + d1;
      l = j - d2;
      

      if (traceback) fprintf(outf,"---------------------------------------------------\n");
      if (traceback) fprintf(outf,"%d %d %d %d  \n", j, d, d1, d2); 
      
      /* Determine the matrix we have to trace down */
      
      if (curr_tr->type1 == dpcP)
	typemtx = VX;

      else if (curr_tr->type1 == dpcS) 
	typemtx = WX;

      else if (curr_tr->type1 == dpcB) 
	typemtx = WBX;

      switch (typemtx){

      /*************************************
       * TRACE VX (pair exterior, no hole) 
       *************************************/
      case VX: 
	if (traceback) fprintf(outf,"tracing VX  %d   \n", vx[j][d]);
	
	TraceVX_nested(outf, s, len, icfg, wbx, vx, j, d, flag, curr_tr, dolist, traceback);
	if (*flag == FALSE)
	  Die("something went wrong in the traceback of vx");
	break;

      /************************************* 
       * TRACE WX  
       *************************************/
      case WX: 
	if (traceback) fprintf(outf,"tracing WX %d \n", wx[j][d]);
	
	TraceWX_nested(outf, s, len, icfg, wx, vx, j, d, flag, curr_tr, dolist, traceback);
	if (*flag == FALSE)
	  Die("something went wrong in the traceback of wx");
	break;
      
      /************************************* 
       * TRACE WBX  
       *************************************/
      case WBX: 
	if (traceback) fprintf(outf,"tracing WBX %d \n", wbx[j][d]);
	
	TraceWBX_nested(outf, s, len, icfg, wbx, vx, j, d, flag, curr_tr, dolist, traceback);
	if (*flag == FALSE)
	  Die("something went wrong in the traceback of wbx");
	break;
      
     default:
	Die("invalid traceback matrix assignement");
      }
   } /* while something is in the trace stack */
  
  FreeTracekn(curr_tr);
  FreeTraceknstack(dolist);
  *ret_trace = tr;
}

/* Function: TraceMtx()
 * 
 * Purpose:  Recover an RNA structure as a tree of trace_s structures,
 *           by tracing back from some point in the matrix j, d, ROOT.
 *           emitl and emitr in the trace structure are 0..N-1 positions
 *           in the sequence; type is dpcP, dpcL, etc, (state type).
 *           
 * Args:     s         - sequence, in integer form (A=0, C=1, G=2, U=3)
 *           len       - length of s, which is 0..len-1
 *           icfg      - the SCFG
 *           whx       - the matrix
 *           j         - end point of the match in vertical
 *           d         - end point of the match in horizontal
 *           d1
 *           d2
 *           ret_trace - RETURN: traceback tree
 *           
 * Return:   void. ret_trace is allocated here and must be free'd by caller,
 *           using FreeTracekn().
 */
void
TraceMtx(FILE *outf, int *s, int len, int **icfg, 
	 int **wx, int **wbx, int **vx, 
	 int ****whx, int ****vhx, int ****zhx, int ****yhx,
	 int j, int d, int d1, int d2, int approx,  
	 struct tracekn_s **ret_trace, int traceback)
{
  struct tracekn_s          *tr;  /* the traceback tree under construction  */
  struct tracekn_s     *curr_tr;  /* ptr to node of tr we're working on     */
  struct traceknstack_s *dolist;  /* pushdown stack of active tr nodes      */
  int                   i, k, l;  /* coords in mtx's                        */
  int                   typemtx;  /* type of matrix : VHX, VX, ZHX, YHX, WX, WBX, WHX */
                                     
  /* Initialize.
   * Start at j, d, d1, d2.
   */
  tr     = InitTracekn();       /* start a trace tree */
  dolist = InitTraceknstack();	/* start a stack for traversing the trace tree */

  if (len == 1) {d = 0; d1 = 0;  d2= 0;} /* redefine d, d1, d2 for len == 1 */
  else if (d < 0 || d1 < 0 || d2 < 0)
    Die("check your traceback assingments");
  
  i = j - d;
  k = i + d1;
  l = j - d2;
      
  /* assigned state to (i,j): dpcS = (relation unknown, trace wx)
   *                          dpcB = (relation unknown, trace wbx)
   *                          dpcP = (paired)
   *                          dpcL = (unpaired)
   */

  curr_tr = AttachTracekn(tr, i, j, k, l, dpcS, dpcS); 
  PushTraceknstack(dolist, curr_tr);

  /* Recursion. While there's active nodes in the stack, trace from them.
   * 
   * This is cribbed from FillMtx(); it's almost the exact reverse.
   * We know the best score, we just have to figure out where it came from.
   */
  while ((curr_tr = PopTraceknstack(dolist)) != NULL)
   {
     /* get some useful numbers, mostly for clarity */
     i = curr_tr->emiti;
     j = curr_tr->emitj;
     k = curr_tr->emitk;
     l = curr_tr->emitl;

     d  = j - i;
     d1 = k - i;
     d2 = j - l;
     
     if (d < 0 || d1 < 0 || d2 < 0)
       Die("check your traceback assingments");
     
     if (d == 0) continue;
     
     else if (d1 + d2 >= d-1){
       d1 = (int)(d/2);
       d2 = d - d1 - 1;
       
       k = i + d1;
       l = j - d2;
     }
     if (traceback) fprintf(outf,"---------------------------------------------------\n");
     if (traceback) fprintf(outf,"%d %d %d %d  \n", j, d, d1, d2); 
     
     /* Determine the matrix we have to traceback */
     
      if (d1 == 0 && d2 == 0) continue;
      
      else if (d1+d2 >= d-1 && curr_tr->type1 == dpcP)
	typemtx = VX;
      
      else if (d1+d2 >= d-1 && curr_tr->type1 == dpcS) 
	typemtx = WX;

      else if (d1+d2 >= d-1 && curr_tr->type1 == dpcB) 
	typemtx = WBX;

      else if (d1+d2 < d-1 && 
	       curr_tr->type1 == dpcP && curr_tr->type2 == dpcP)
	typemtx = VHX;

      else if (d1+d2 < d-1 && curr_tr->type1 == dpcP) 
	typemtx = ZHX;

      else if (d1+d2 < d-1 && curr_tr->type2 == dpcP)
	typemtx = YHX;

      else 
	typemtx = WHX;

      switch (typemtx){
      /***********************
       * TRACE VHX (pair both) 
       ***********************/
      case VHX: 
	if (traceback) fprintf(outf,"tracing VHX  %d \n", vhx[j][d][d1][d2]); 
	  
	TraceVHX(outf, s, len, icfg, whx, vhx, 
		 j, d, d1, d2, curr_tr, dolist, traceback);
	break;

      /*************************************
       * TRACE VX (pair exterior, no hole) 
       *************************************/
      case VX: 
	if (zhx[j][d][d1][d2] != vx[j][d])
	  Die("zhx is not correctly assigned");
	  
	if (traceback) fprintf(outf,"tracing VX  %d   \n", vx[j][d]);
	
	TraceVX(outf, s, len, icfg, wbx, vx, whx, zhx, yhx, 
		j, d, approx, curr_tr, dolist, traceback);
	break;
      
      /*************************************
       * TRACE ZHX hole (d1+d2 < d-1) 
       *************************************/
      case ZHX: 
	if (traceback) fprintf(outf,"tracing ZHX  %d \n", zhx[j][d][d1][d2]);
	
	TraceZHX(outf, s, len, icfg, wbx, vx, whx, vhx, zhx,  
		 j, d, d1, d2, curr_tr, dolist, traceback);
	break;

      /************************************* 
       * TRACE YHX (pair interior) 
       *************************************/
      case YHX: 
	if (traceback) fprintf(outf,"tracing YHX %d \n", yhx[j][d][d1][d2]);
	
	TraceYHX(outf, s, len, icfg, wbx, vx, whx, vhx, yhx, 
		 j, d, d1, d2, curr_tr, dolist, traceback);
	break;

      /************************************* 
       * TRACE WX  
       *************************************/
      case WX: 
	if (traceback) fprintf(outf,"tracing WX %d \n", wx[j][d]);
	
	TraceWX(outf, s, len, icfg, wx, vx, whx, yhx, 
		j, d, curr_tr, dolist, traceback);
	break;
      
      /************************************* 
       * TRACE WBX  
       *************************************/
      case WBX: 
	if (traceback) fprintf(outf,"tracing WBX %d \n", wbx[j][d]);
	
	TraceWBX(outf, s, len, icfg, wbx, vx, whx, yhx, 
		 j, d, approx, curr_tr, dolist, traceback);
	break;
      
      /************************************* 
       * TRACE WHX  
       *************************************/
      case WHX: 
	if (traceback) fprintf(outf,"tracing whx %d \n", whx[j][d][d1][d2]); 
	
	TraceWHX(outf, s, len, icfg, wbx, vx, whx, vhx, zhx, yhx, 
		 j, d, d1, d2, curr_tr, dolist, traceback);
	break;

      default:
	Die("invalid traceback matrix assignement");
      }
   } /* while something is in the trace stack */
  
  FreeTracekn(curr_tr);
  FreeTraceknstack(dolist);
  *ret_trace = tr;
}

