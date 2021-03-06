/* wxgraphs.c 
 *
 * includes functions to calculate all the diagrams 
 * that fill and traceback the no-hole matrix:
 *      wx[j][d] (len x len)
 *               [(j-d,j) are base pared]
 *
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>

#include <easel.h>
#include <esl_sqio.h>

#include "pknots.h"
#include "pk_irredsurf.h"
#include "pk_wbxgraphs.h"
#include "pk_trace.h"
#include "pk_util.h"


/* (W1)__/ PAIR (i,j).
 */
int
W1(ESL_DSQ *s, int len, struct rnapar_2 *rnapar, int **icfg, int **vx, int j, int d)
{ 
  int W1;
  int i;

  i = j - d;
  
  W1 = vx[j][d];
  
  return W1;
}
void
trace_W1(FILE *outf, int **vx,  int j, int d, 
	 struct tracekn_s *curr_tr, struct traceknstack_s *dolist, int traceback)
{  
  int i;

  i = j - d;

  if (traceback) fprintf(outf," (W1) trace vx  %d %d %d\n", j, d, vx[j][d]);
  
  PushTraceknstack(dolist, AttachTracekn(curr_tr, i, j, i+(int)(d/2), 
					 i+(int)(d/2)+1, dpcP, dpcS));
}

/* (W2)__/ PAIR (i+1,j-1). Connect to (i+1,j-1) (j-1,d-2)
 */
int
W2(ESL_DSQ *s, int len, struct rnapar_2 *rnapar, int **icfg, int **vx, int j, int d)
{ 
  int W2;
  int i;

  i = j - d;
  
  if (d > 1)
    W2 = wsf*icfg[idxL(s[i])][idxP(s[i+1],s[j-1])]
      + wsf*icfg[idxR(s[j])][idxP(s[i+1],s[j-1])] 
      + vx[j-1][d-2];
  else
    W2 = -BIGINT;
  
  return W2;
}
void
trace_W2(FILE *outf, int **vx, int j, int d, 
	 struct tracekn_s *curr_tr, struct traceknstack_s *dolist, int traceback)
{
  int i;
  
  i = j - d;
  
  if (d > 1){
    if (traceback) fprintf(outf," (W2) SR SL, trace vx %d %d %d\n", 
			   j-1, d-2,  vx[j-1][d-2]);
    
    PushTraceknstack(dolist, AttachTracekn(curr_tr, i+1, j-1, i+1+(int)((d-2)/2), 
					   i+2+(int)((d-2)/2), dpcP, dpcS));
  }
}

/* (W3)__/  SINGLET-LEFT [(i+1,j) base-paired]. 
 */
int
W3(ESL_DSQ *s, int len, struct rnapar_2 *rnapar, int **icfg, int **vx, int j, int d)
{ 
  int W3;
  int i;

  i = j - d;
  
  if (d > 0)
    W3 = wsf*icfg[idxL(s[i])][idxP(s[i+1],s[j])] 
      + vx[j][d-1];
  else
    W3 = -BIGINT;
  
  return W3;
}
void
trace_W3(FILE *outf, int **vx, int j, int d, 
	 struct tracekn_s *curr_tr, struct traceknstack_s *dolist, int traceback)
{
  int i;
  
  i = j - d;

  if (d > 0){
    if (traceback) fprintf(outf," (W3) SL, trace vx  %d %d %d\n", j, d-1, vx[j][d-1]);
    
    PushTraceknstack(dolist, AttachTracekn(curr_tr, i+1, j, i+1+(int)((d-1)/2), 
					   i+2+(int)((d-1)/2), dpcP, dpcS));
  }
}

/* (W4)__/  SINGLET-RIGHT. [(i,j-1) base-paired].  
 */
int
W4(ESL_DSQ *s, int len, struct rnapar_2 *rnapar, int **icfg, int **vx, int j, int d)
{ 
  int W4;
  int i;

  i = j - d;
  
  if (d > 0)
    W4 = wsf*icfg[idxR(s[j])][idxP(s[i],s[j-1])] 
      + vx[j-1][d-1];
  else
    W4 = -BIGINT;
  
  return W4;
}
void
trace_W4(FILE *outf, int **vx, int j, int d, 
	 struct tracekn_s *curr_tr, struct traceknstack_s *dolist, int traceback)
{
  int i;
  
  i = j - d;

  if (d > 0){
    if (traceback) fprintf(outf," (W4) SR, trace vx %d %d %d\n", j-1, d-1, vx[j-1][d-1]);
    
    PushTraceknstack(dolist, AttachTracekn(curr_tr, i, j-1, i+(int)((d-1)/2), 
					   i+1+(int)((d-1)/2), dpcP, dpcS));
  }
}

/* (W5)__/  SINGLET-LEFT. Connect (i+1,j) (j,d-1) 
 */
int
W5(ESL_DSQ *s, int len, struct rnapar_2 *rnapar, int **icfg, int **wx, int j, int d)
{ 
  int W5;
  int i;

  i = j - d;
  
  if (d > 1)
    W5 = wx[j][d-1];
  else
    W5 = -BIGINT;
  
  return W5;
}
void
trace_W5(FILE *outf, int **wx, int j, int d, 
	 struct tracekn_s *curr_tr, struct traceknstack_s *dolist, int traceback)
{
  int i;
  
  i = j - d;

  if (d > 1){
    if (traceback) fprintf(outf," (W5) SL, trace wx  %d %d %d\n",  j, d-1, wx[j][d-1]);
    
    PushTraceknstack(dolist, AttachTracekn(curr_tr, i+1, j, i+1+(int)((d-1)/2), 
					   i+2+(int)((d-1)/2), dpcS, dpcS));
  }
}

/* (W6)__/  SINGLET-RIGHT. Connect (i,j-1) (j-1,d-1) 
 */
int
W6(ESL_DSQ *s, int len, struct rnapar_2 *rnapar, int **icfg, int **wx, int j, int d)
{ 
  int W6;
  int i;

  i = j - d;
  
  if (d > 0)
    W6 = wx[j-1][d-1];
  else
    W6 = -BIGINT;
  
  return W6;
}
void
trace_W6(FILE *outf, int **wx, int j, int d, 
	 struct tracekn_s *curr_tr, struct traceknstack_s *dolist, int traceback)
{
  int i;
  
  i = j - d;

  if (d > 0){
    if (traceback) fprintf(outf," (W6) SR, trace wx %d %d %d\n", j-1, d-1, wx[j-1][d-1]);
    
    PushTraceknstack(dolist, AttachTracekn(curr_tr, i, j-1, i+(int)((d-1)/2), 
					   i+1+(int)((d-1)/2), dpcS, dpcS));
  }
}

/* BIFURCATIONS (no hole) (two type of diagrams W7, W8)
 *
 *     One  connects i to i+mid.   (i+mid,mid)
 *
 * Another  connects i+mid+1 to j. (j,d-mid-1)
 */         

/* (W7)__/  2 WX 
 */
int
W7(ESL_DSQ *s, int len, struct rnapar_2 *rnapar, int **icfg, int **wx, int j, int d, int mid)
{ 
  int W7;
  int i;

  i = j - d;
  
  W7 = wx[i+mid][mid] + wx[j][d-mid-1];
  
  return W7;
}
void
trace_W7(FILE *outf, int **wx, int j, int d, int mid,
	 struct tracekn_s *curr_tr, struct traceknstack_s *dolist, int traceback)
{
  int i;
  
  i = j - d;

  if (traceback) fprintf(outf," (W7) trace wx  %d %d %d\n", i+mid, mid, wx[i+mid][mid]);
  if (traceback) fprintf(outf," (W7) trace wx  %d %d %d\n", j, d-mid-1, wx[j][d-mid-1]);
  
  PushTraceknstack(dolist, AttachTracekn(curr_tr, i, i+mid, i+(int)(mid/2), 
					 i+(int)(mid/2)+1, dpcS, dpcS));
  PushTraceknstack(dolist, AttachTracekn(curr_tr, i+mid+1, j, i+mid+1+(int)((d-mid-1)/2),
					 i+mid+2+(int)((d-mid-1)/2), dpcS, dpcS));
}

/* (W8)
 * 2 VX 
 */

/* W8_1 to W8_4
 * contiguous coaxial pairs (i,i+mid) (i+mid+1,j)
 *                  ...add stack[i+mid][i][i+mid+1][j]
 * check also for external danglings
 */

/* (W8.1)__/ no danglings / coaxial = stack + 1  */
int
W8_1(ESL_DSQ *s, int len, struct rnapar_2 *rnapar, int **icfg, int **vx, int j, int d, int mid)
{ 
  int W8_1;
  int i;

  i = j - d;
  
  W8_1 = wsf*(icfg[idxPS(s[i+mid],s[i])][idxPS(s[i+mid+1],s[j])] + INTSCALE)
    + vx[i+mid][mid] + vx[j][d-mid-1];
  
  return W8_1;
}
void
trace_W8_1(FILE *outf, int **vx, int j, int d, int mid,
	   struct tracekn_s *curr_tr, struct traceknstack_s *dolist, int traceback)
{
  int i;
  
  i = j - d;

  if (traceback) fprintf(outf," (W8.1) trace vx %d %d %d\n", i+mid, mid, vx[i+mid][mid]);
  if (traceback) fprintf(outf," (W8.1) trace vx %d %d %d\n", j, d-mid-1, vx[j][d-mid-1]);
  
  PushTraceknstack(dolist, AttachTracekn(curr_tr, i, i+mid, i+(int)(mid/2), 
					 i+(int)(mid/2)+1, dpcP, dpcS));
  PushTraceknstack(dolist, AttachTracekn(curr_tr, i+mid+1, j, i+mid+1+(int)((d-mid-1)/2),
					 i+mid+2+(int)((d-mid-1)/2), dpcP, dpcS));
}

/* (W8.2)__/ i dangles off i+1,i+mid */
int
W8_2(ESL_DSQ *s, int len, struct rnapar_2 *rnapar, int **icfg, int **vx, int j, int d, int mid)
{ 
  int W8_2;
  int i;

  i = j - d;
  
  W8_2 = wsf*icfg[idxL(s[i])][idxP(s[i+1],s[i+mid])]
    + wsf*icfg[idxPS(s[i+mid],s[i+1])][idxPS(s[i+mid+1],s[j])]
    + vx[i+mid][mid-1] + vx[j][d-mid-1];

  return W8_2;
}
void
trace_W8_2(FILE *outf, int **vx, int j, int d, int mid,
	   struct tracekn_s *curr_tr, struct traceknstack_s *dolist, int traceback)
{
  int i;
  
  i = j - d;

  if (traceback) fprintf(outf," (W8.2) dangle i, trace vx %d %d %d\n", 
			 i+mid, mid-1, vx[i+mid][mid-1]);
  if (traceback) fprintf(outf," (W8.2) trace vx %d %d %d\n", 
			 j, d-mid-1, vx[j][d-mid-1]);
  
  PushTraceknstack(dolist, AttachTracekn(curr_tr, i+1, i+mid, i+(int)(mid/2)+1, 
					 i+(int)(mid/2)+2, dpcP, dpcS));
  PushTraceknstack(dolist,AttachTracekn(curr_tr, i+mid+1, j, i+mid+1+(int)((d-mid-1)/2),
					i+mid+2+(int)((d-mid-1)/2), dpcP, dpcS));
}

/* (W8.3)__/ j dangles off i+mid+1,j-1 */
int
W8_3(ESL_DSQ *s, int len, struct rnapar_2 *rnapar, int **icfg, int **vx, int j, int d, int mid)
{ 
  int W8_3;
  int i;

  i = j - d;
  
  if (d-mid > 1)
    W8_3 = wsf*icfg[idxR(s[j])][idxP(s[i+mid+1],s[j-1])]
      + wsf*icfg[idxPS(s[i+mid],s[i])][idxPS(s[i+mid+1],s[j-1])]
      + vx[i+mid][mid] + vx[j-1][d-mid-2];
  else
    W8_3 = -BIGINT;
  
  return W8_3;
}
void
trace_W8_3(FILE *outf, int **vx, int j, int d, int mid,
	   struct tracekn_s *curr_tr, struct traceknstack_s *dolist, int traceback)
{
  int i;
  
  i = j - d;

  if (d-mid > 1){
    if (traceback) fprintf(outf," (W8.3) trace vx %d %d %d\n", 
			   i+mid, mid, vx[i+mid][mid]);
    if (traceback) fprintf(outf," (W8.3) dangle j, trace vx %d %d %d\n", 
			   j-1, d-mid-2, vx[j-1][d-mid-2]);
    
    PushTraceknstack(dolist, AttachTracekn(curr_tr, i, i+mid, i+(int)(mid/2), 
					   i+(int)(mid/2)+1, dpcP, dpcS));
    PushTraceknstack(dolist, AttachTracekn(curr_tr, i+mid+1, j-1, i+mid+1+(int)((d-mid-2)/2),
					   i+mid+2+(int)((d-mid-2)/2), dpcP, dpcS));
  }
}

/* (W8.4)__/ i dangles off i+1,i+mid / j dangles off i+mid+1,j-1 */
int
W8_4(ESL_DSQ *s, int len, struct rnapar_2 *rnapar, int **icfg, int **vx, int j, int d, int mid)
{ 
  int W8_4;
  int i;

  i = j - d;
  
  if (d-mid > 1)
    W8_4 = wsf*icfg[idxL(s[i])][idxP(s[i+1],s[i+mid])]
      + wsf*icfg[idxR(s[j])][idxP(s[i+mid+1],s[j-1])]
      + wsf*icfg[idxPS(s[i+mid],s[i+1])][idxPS(s[i+mid+1],s[j-1])]
      + vx[i+mid][mid-1] + vx[j-1][d-mid-2];
  else
    W8_4 = -BIGINT;
  
  return W8_4;
}
void
trace_W8_4(FILE *outf, int **vx, int j, int d, int mid,
	   struct tracekn_s *curr_tr, struct traceknstack_s *dolist, int traceback)
{
  int i;
  
  i = j - d;

  if (d-mid > 1){
    if (traceback) fprintf(outf," (W8.4) dangle i, trace vx %d %d %d\n", 
			   i+mid, mid-1, vx[i+mid][mid-1]);
    if (traceback) fprintf(outf," (W8.4) dangle j, trace vx %d %d %d\n", 
			   j-1, d-mid-2, vx[j-1][d-mid-2]);
    
    PushTraceknstack(dolist, AttachTracekn(curr_tr, i+1, i+mid, i+(int)(mid/2)+1, 
					   i+(int)(mid/2)+2, dpcP, dpcS));
    PushTraceknstack(dolist,
		     AttachTracekn(curr_tr, i+mid+1, j-1, i+mid+1+(int)((d-mid-2)/2),
				   i+mid+2+(int)((d-mid-2)/2), dpcP, dpcS));
  }
}
 
/* W8_5 to W8_8
 * non-contiguous coaxial pairs (i,i+mid) (i+mid+2,j)
 *                  ...add stack[i+mid][i][i+mid+2][j]
 * check also for external danglings
 */

/* (W8.5)__/ no danglings (do not add +1 becuase the stacking is not perfect) */
int
W8_5(ESL_DSQ *s, int len, struct rnapar_2 *rnapar, int **icfg, int **vx, int j, int d, int mid)
{ 
  int W8_5;
  int i;

  i = j - d;
  
  if (d-mid > 1)
    W8_5 = wsf*icfg[idxPS(s[i+mid],s[i])][idxPS(s[i+mid+2],s[j])]
      + vx[i+mid][mid] + vx[j][d-mid-2];
   else
    W8_5 = -BIGINT;
   
  return W8_5;
}
void
trace_W8_5(FILE *outf, int **vx, int j, int d, int mid,
	   struct tracekn_s *curr_tr, struct traceknstack_s *dolist, int traceback)
{
  int i;
  
  i = j - d;

  if (d-mid > 1){
    if (traceback) fprintf(outf," (W8.5) trace vx %d %d %d\n", 
			   i+mid, mid, vx[i+mid][mid]);
    if (traceback) fprintf(outf," (W8.5) trace vx %d %d %d\n", 
			   j, d-mid-2, vx[j][d-mid-2]);
    
    PushTraceknstack(dolist, AttachTracekn(curr_tr, i, i+mid, i+(int)(mid/2), 
					   i+(int)(mid/2)+1, dpcP, dpcS));
    PushTraceknstack(dolist, AttachTracekn(curr_tr, i+mid+2, j, i+mid+2+(int)((d-mid-2)/2),
					   i+mid+3+(int)((d-mid-2)/2), dpcP, dpcS));
  }
}

/* (W8.6)__/ i dangles off i+1,i+mid */
int
W8_6(ESL_DSQ *s, int len, struct rnapar_2 *rnapar, int **icfg, int **vx, int j, int d, int mid)
{ 
  int W8_6;
  int i;

  i = j - d;
  
  if (d-mid > 1)
    W8_6 = wsf*icfg[idxL(s[i])][idxP(s[i+1],s[i+mid])]
      + wsf*icfg[idxPS(s[i+mid],s[i+1])][idxPS(s[i+mid+2],s[j])]
      + vx[i+mid][mid-1] + vx[j][d-mid-2];
  else
    W8_6 = -BIGINT;

  return W8_6;
}
void
trace_W8_6(FILE *outf, int **vx, int j, int d, int mid,
	   struct tracekn_s *curr_tr, struct traceknstack_s *dolist, int traceback)
{
  int i;
  
  i = j - d;

  if (d-mid > 1){
    if (traceback) fprintf(outf," (W8.6) dangle i, trace vx %d %d %d\n", 
			   i+mid, mid-1, vx[i+mid][mid-1]);
    if (traceback) fprintf(outf," (W8.6) trace vx %d %d %d\n", 
			   j, d-mid-2, vx[j][d-mid-2]);
    
    PushTraceknstack(dolist, AttachTracekn(curr_tr, i+1, i+mid, i+(int)(mid/2)+1, 
					   i+(int)(mid/2)+2, dpcP, dpcS));
    PushTraceknstack(dolist,AttachTracekn(curr_tr, i+mid+2, j, i+mid+2+(int)((d-mid-2)/2),
					  i+mid+3+(int)((d-mid-2)/2), dpcP, dpcS));
  }
}

/* (W8.7)__/ j dangles off i+mid+1,j-1 */
int
W8_7(ESL_DSQ *s, int len, struct rnapar_2 *rnapar, int **icfg, int **vx, int j, int d, int mid)
{ 
  int W8_7;
  int i;

  i = j - d;
  
  if (d-mid > 2)
    W8_7 = wsf*icfg[idxR(s[j])][idxP(s[i+mid+2],s[j-1])]
      + wsf*icfg[idxPS(s[i+mid],s[i])][idxPS(s[i+mid+2],s[j-1])]
      + vx[i+mid][mid] + vx[j-1][d-mid-3];
  else
    W8_7 = -BIGINT;
  
  return W8_7;
}
void
trace_W8_7(FILE *outf, int **vx, int j, int d, int mid,
	   struct tracekn_s *curr_tr, struct traceknstack_s *dolist, int traceback)
{
  int i;
  
  i = j - d;

  if (d-mid > 2){
    if (traceback) fprintf(outf," (W8.7) trace vx %d %d %d\n", 
			   i+mid, mid, vx[i+mid][mid]);
    if (traceback) fprintf(outf," (W8.7) dangle j, trace vx %d %d %d\n", 
			   j-1, d-mid-3, vx[j-1][d-mid-3]);
    
    PushTraceknstack(dolist, AttachTracekn(curr_tr, i, i+mid, i+(int)(mid/2), 
					   i+(int)(mid/2)+1, dpcP, dpcS));
    PushTraceknstack(dolist, AttachTracekn(curr_tr, i+mid+2, j-1, i+mid+2+(int)((d-mid-3)/2),
					   i+mid+3+(int)((d-mid-3)/2), dpcP, dpcS));
  }
}

/* (W8.8)__/ i dangles off i+1,i+mid / j dangles off i+mid+2,j-1 */
int
W8_8(ESL_DSQ *s, int len, struct rnapar_2 *rnapar, int **icfg, int **vx, int j, int d, int mid)
{ 
  int W8_8;
  int i;

  i = j - d;
  
  if (d-mid > 2)
    W8_8 = wsf*icfg[idxL(s[i])][idxP(s[i+1],s[i+mid])]
      + wsf*icfg[idxR(s[j])][idxP(s[i+mid+2],s[j-1])]
      + wsf*icfg[idxPS(s[i+mid],s[i+1])][idxPS(s[i+mid+2],s[j-1])]
      + vx[i+mid][mid-1] + vx[j-1][d-mid-3];
  else
    W8_8 = -BIGINT;
  
  return W8_8;
}
void
trace_W8_8(FILE *outf, int **vx, int j, int d, int mid,
	   struct tracekn_s *curr_tr, struct traceknstack_s *dolist, int traceback)
{
  int i;
  
  i = j - d;

  if (d-mid > 1){
    if (traceback) fprintf(outf," (W8.8) dangle i, trace vx %d %d %d\n", 
			   i+mid, mid-1, vx[i+mid][mid-1]);
    if (traceback) fprintf(outf," (W8.8) dangle j, trace vx %d %d %d\n", 
			   j-1, d-mid-3, vx[j-1][d-mid-3]);
    
    PushTraceknstack(dolist, AttachTracekn(curr_tr, i+1, i+mid, i+(int)(mid/2)+1, 
					   i+(int)(mid/2)+2, dpcP, dpcS));
    PushTraceknstack(dolist,
		     AttachTracekn(curr_tr, i+mid+2, j-1, i+mid+2+(int)((d-mid-3)/2),
				   i+mid+3+(int)((d-mid-3)/2), dpcP, dpcS));
  }
}

/* BIFURCATIONS (with hole) (2 types of diagrams W9-W10)
 * (add P11 penalty for wx)
 *
 *     one connects (i,i+mid1) and (i+mid-mid2,i+mid). 
 *                  (i+mid,mid,mid1,mid2).
 *
 * Another connects (i+mid1+2,i+mid-mid2-1) and (i+mid+3,j). 
 *                  (j,d-mid1-2,mid-mid1-mid2-3,d-mid-3).	
 *
 *  [notice that we have left 1 base on the left and 2 on the right already unparied
 *   for pks topological constrains.]
 */
/* (W9)__/ 2 WHX
 */
int
W9(ESL_DSQ *s, int len, struct rnapar_2 *rnapar, int **icfg, int ****whx, int j, int d, int mid, int mid1, int mid2)
{ 
  int W9;
  int i;

  i = j - d;
 
  W9 =  rnapar->P11 + 3*rnapar->P6P
    + whx[i+mid][mid][mid1][mid2] 
    + whx[j][d-mid1-2][mid-mid1-mid2-3][d-mid-3];
  
  return W9;
}
void
trace_W9(FILE *outf, int ****whx, int j, int d, int mid, int mid1, int mid2,
	 struct tracekn_s *curr_tr, struct traceknstack_s *dolist, int traceback)
{
  int i;
  
  i = j - d;

  if (traceback) fprintf(outf," (W9) trace whx %d %d %d %d %d\n", 
			 i+mid, mid, mid1, mid2,
			 whx[i+mid][mid][mid1][mid2]);
  if (traceback) fprintf(outf," (W9) trace whx %d %d %d %d %d\n", 
			 j, d-mid1-2, mid-mid1-mid2-3, d-mid-3,
			 whx[j][d-mid1-2][mid-mid1-mid2-3][d-mid-3]);
  
  PushTraceknstack(dolist, AttachTracekn(curr_tr, i, i+mid, i+mid1, 
					 i+mid-mid2, dpcS, dpcS));
  PushTraceknstack(dolist, AttachTracekn(curr_tr, i+mid1+2, j, i+mid-mid2-1,
					 i+mid+3, dpcS, dpcS));
}
 
/* (W10)__/ 1 YHX 1 YHX 
 * coaxial pairs (i+mid1,i+mid-mid2) (i+mid-mid2-1,i+mid+3)
 * ....add stack[i+mid-mid2-1][i+mid+3][i+mid-mid2][i+mid1]
 */
/* (W10.1)__/ contiguous coaxial: i+mid-mid2-1, i+mid-mid2     /
            / i+mid1+1 dangles right off i+mid-mid2,i+mid1     /
            / i+mid+2  dangles left  off i+mid+3,i+mid-mid2-1  */
int
W10_1(ESL_DSQ *s, int len, struct rnapar_2 *rnapar, int **icfg, int ****yhx, 
      int j, int d, int mid, int mid1, int mid2)
{ 
  int W10_1;
  int i;

  i = j - d;
  
  W10_1 = rnapar->P11 + 2*rnapar->P10P + 3*rnapar->P6P
    + wkn*icfg[idxR(s[i+mid1+1])][idxP(s[i+mid-mid2],s[i+mid1])]     
    + wkn*icfg[idxL(s[i+mid+2])][idxP(s[i+mid+3],s[i+mid-mid2-1])]     
    + wkn*icfg[idxPS(s[i+mid-mid2-1],s[i+mid+3])]
              [idxPS(s[i+mid-mid2],s[i+mid1])] 
    + yhx[i+mid][mid][mid1][mid2] 
    + yhx[j][d-mid1-2][mid-mid1-mid2-3][d-mid-3];
  
  return W10_1;
}
void
trace_W10_1(FILE *outf, int ****yhx, int j, int d, int mid, int mid1, int mid2,
	    struct tracekn_s *curr_tr, struct traceknstack_s *dolist, int traceback)
{
  int i;
  
  i = j - d;

  if (traceback) fprintf(outf," (W10.1) trace yhx %d %d %d %d %d\n", 
			 i+mid, mid, mid1, mid2,
			 yhx[i+mid][mid][mid1][mid2]);
  if (traceback) fprintf(outf," (W10.1) trace yhx %d %d %d %d %d\n", 
			 j, d-mid1-2, mid-mid1-mid2-3, d-mid-3, 
			 yhx[j][d-mid1-2][mid-mid1-mid2-3][d-mid-3]);
  
  PushTraceknstack(dolist, AttachTracekn(curr_tr, i, i+mid, 
					 i+mid1, i+mid-mid2, dpcS, dpcP));
  PushTraceknstack(dolist, AttachTracekn(curr_tr, i+mid1+2, j, 
					 i+mid-mid2-1, i+mid+3, dpcS, dpcP));
}

/* (W10.2)__/ non-contiguous coaxial: i+mid-mid2-2, i+mid-mid2 /
            / i+mid1+1 dangles right off i+mid-mid2,i+mid1     /
            / i+mid+2  dangles left  off i+mid+3,i+mid-mid2-1  */
int
W10_2(ESL_DSQ *s, int len, struct rnapar_2 *rnapar, int **icfg, int ****yhx, 
      int j, int d, int mid, int mid1, int mid2)
{ 
  int W10_2;
  int i;

  i = j - d;
  
  if (mid-mid1-mid2 > 3)
    W10_2 = rnapar->P11 + 2*rnapar->P10P + 4*rnapar->P6P
      + wkn*icfg[idxR(s[i+mid1+1])][idxP(s[i+mid-mid2],s[i+mid1])]     
      + wkn*icfg[idxL(s[i+mid+2])][idxP(s[i+mid+3],s[i+mid-mid2-2])]     
      + wkn*icfg[idxPS(s[i+mid-mid2-2],s[i+mid+3])]
                [idxPS(s[i+mid-mid2],s[i+mid1])] 
      + yhx[i+mid][mid][mid1][mid2] 
      + yhx[j][d-mid1-2][mid-mid1-mid2-4][d-mid-3];
  else
    W10_2 = -BIGINT;

  return W10_2;
}
void
trace_W10_2(FILE *outf, int ****yhx, int j, int d, int mid, int mid1, int mid2,
	    struct tracekn_s *curr_tr, struct traceknstack_s *dolist, int traceback)
{
  int i;
  
  i = j - d;

  if (mid-mid1-mid2 > 3){
    if (traceback) fprintf(outf," (W10.2) trace yhx %d %d %d %d %d\n", 
			   i+mid, mid, mid1, mid2,
			   yhx[i+mid][mid][mid1][mid2]);
    if (traceback) fprintf(outf," (W10.2) trace yhx %d %d %d %d %d\n", 
			   j, d-mid1-2, mid-mid1-mid2-4, d-mid-3, 
			   yhx[j][d-mid1-2][mid-mid1-mid2-4][d-mid-3]);
    
    PushTraceknstack(dolist, AttachTracekn(curr_tr, i, i+mid, 
					   i+mid1, i+mid-mid2, dpcS, dpcP));
    PushTraceknstack(dolist, AttachTracekn(curr_tr, i+mid1+2, j, 
					   i+mid-mid2-2, i+mid+3, dpcS, dpcP));
  }
}

