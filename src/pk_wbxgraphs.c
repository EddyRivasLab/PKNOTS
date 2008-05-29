/* wbxgraphs.c 
 *
 * includes functions to calculate all the diagrams 
 * that fill and traceback the no-hole matrix:
 *      wbx[j][d] (len x len)
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

/* (WB1)__/ PAIR (i,j).
 */
int
WB1(ESL_DSQ *s, int len, struct rnapar_2 *rnapar, int **icfg, int **vx, int j, int d)
{ 
  int WB1;
  int i;

  i = j - d;
  
  WB1 = rnapar->P10 + vx[j][d];
  
  return WB1;
}
void
trace_WB1(FILE *outf, int **vx,  int j, int d, 
	  struct tracekn_s *curr_tr, struct traceknstack_s *dolist, int traceback)
{  
  int i;

  i = j - d;

  if (traceback) fprintf(outf," (WB1) trace vx  %d %d %d\n", j, d, vx[j][d]);
  
  PushTraceknstack(dolist, AttachTracekn(curr_tr, i, j, i+(int)(d/2), 
					 i+(int)(d/2)+1, dpcP, dpcS));
}

/* (WB2)__/ PAIR (i+1,j-1). Connect to (i+1,j-1) (j-1,d-2)
 */
int
WB2(ESL_DSQ *s, int len, struct rnapar_2 *rnapar, int **icfg, int **vx, int j, int d)
{ 
  int WB2;
  int i;

  i = j - d;
  
  if (d > 1)
    WB2 = rnapar->P10 + 2*rnapar->P6 
      + wsf*icfg[idxL(s[i])][idxP(s[i+1],s[j-1])]
      + wsf*icfg[idxR(s[j])][idxP(s[i+1],s[j-1])] 
      + vx[j-1][d-2];
  else
    WB2 = -BIGINT;
  
  return WB2;
}
void
trace_WB2(FILE *outf, int **vx, int j, int d, 
	  struct tracekn_s *curr_tr, struct traceknstack_s *dolist, int traceback)
{
  int i;
  
  i = j - d;
  
  if (d > 1){
    if (traceback) fprintf(outf," (WB2) SR SL, trace vx %d %d %d\n", j-1, d-2,  vx[j-1][d-2]);
    
    PushTraceknstack(dolist, AttachTracekn(curr_tr, i+1, j-1, i+1+(int)((d-2)/2), 
					   i+2+(int)((d-2)/2), dpcP, dpcS));
  }
}

/* (WB3)__/  SINGLET-LEFT [(i+1,j) base-paired]. 
 */
int
WB3(ESL_DSQ *s, int len, struct rnapar_2 *rnapar, int **icfg, int **vx, int j, int d)
{ 
  int WB3;
  int i;

  i = j - d;
  
  if (d > 0)
    WB3 = rnapar->P10 + rnapar->P6 
      + wsf*icfg[idxL(s[i])][idxP(s[i+1],s[j])] 
      + vx[j][d-1];
  else
    WB3 = -BIGINT;
  
  return WB3;
}
void
trace_WB3(FILE *outf, int **vx, int j, int d, 
	  struct tracekn_s *curr_tr, struct traceknstack_s *dolist, int traceback)
{
  int i;
  
  i = j - d;

  if (d > 0){
    if (traceback) fprintf(outf," (WB3) SL, trace vx  %d %d %d\n", j, d-1, vx[j][d-1]);
    
    PushTraceknstack(dolist, AttachTracekn(curr_tr, i+1, j, i+1+(int)((d-1)/2), 
					   i+2+(int)((d-1)/2), dpcP, dpcS));
  }
}

/* (WB4)__/  SINGLET-RIGHT. [(i,j-1) base-paired].  
 */
int
WB4(ESL_DSQ *s, int len, struct rnapar_2 *rnapar, int **icfg, int **vx, int j, int d)
{ 
  int WB4;
  int i;

  i = j - d;
  
  if (d > 0)
    WB4 = rnapar->P10 + rnapar->P6 
      + wsf*icfg[idxR(s[j])][idxP(s[i],s[j-1])] 
      + vx[j-1][d-1];
  else
    WB4 = -BIGINT;
  
  return WB4;
}
void
trace_WB4(FILE *outf, int **vx, int j, int d, 
	  struct tracekn_s *curr_tr, struct traceknstack_s *dolist, int traceback)
{
  int i;
  
  i = j - d;

  if (d > 0){
    if (traceback) fprintf(outf," (WB4) SR, trace vx %d %d %d\n", j-1, d-1, vx[j-1][d-1]);
    
    PushTraceknstack(dolist, AttachTracekn(curr_tr, i, j-1, i+(int)((d-1)/2), 
					   i+1+(int)((d-1)/2), dpcP, dpcS));
  }
}
 
/* (WB5)__/  SINGLET-LEFT. Connect (i+1,j) (j,d-1) 
 */
int
WB5(ESL_DSQ *s, int len, struct rnapar_2 *rnapar, int **icfg, int **wbx, int j, int d)
{ 
  int WB5;
  int i;

  i = j - d;
  
  if (d > 1)
    WB5 = rnapar->P6 + wbx[j][d-1];
  else
    WB5 = -BIGINT;
  
  return WB5;
}
void
trace_WB5(FILE *outf, int **wbx, int j, int d, 
	  struct tracekn_s *curr_tr, struct traceknstack_s *dolist, int traceback)
{
  int i;
  
  i = j - d;

  if (d > 1){
    if (traceback) fprintf(outf," (WB5) SL, trace wbx  %d %d %d\n",  j, d-1, wbx[j][d-1]);
    
    PushTraceknstack(dolist, AttachTracekn(curr_tr, i+1, j, i+1+(int)((d-1)/2), 
					   i+2+(int)((d-1)/2), dpcB, dpcS));
  }
}

/* (WB6)__/  SINGLET-RIGHT. Connect (i,j-1) (j-1,d-1) 
 */
int
WB6(ESL_DSQ *s, int len, struct rnapar_2 *rnapar, int **icfg, int **wbx, int j, int d)
{ 
  int WB6;
  int i;

  i = j - d;
  
  if (d > 0)
    WB6 = rnapar->P6 + wbx[j-1][d-1];
  else
    WB6 = -BIGINT;
  
  return WB6;
}
void
trace_WB6(FILE *outf, int **wbx, int j, int d, 
	  struct tracekn_s *curr_tr, struct traceknstack_s *dolist, int traceback)
{
  int i;
  
  i = j - d;

  if (d > 0){
    if (traceback) fprintf(outf," (WB6) SR, trace wbx %d %d %d\n", j-1, d-1, wbx[j-1][d-1]);
    
    PushTraceknstack(dolist, AttachTracekn(curr_tr, i, j-1, i+(int)((d-1)/2), 
					   i+1+(int)((d-1)/2), dpcB, dpcS));
  }
}


/* BIFURCATIONS (no hole) (two type of diagrams WB7, WB8)
 *
 *     One  connects i to i+mid.   (i+mid,mid)
 *
 * Another  connects i+mid+1 to j. (j,d-mid-1)
 */         

/* (WB7)__/  2 WBX 
 */
int
WB7(ESL_DSQ *s, int len, struct rnapar_2 *rnapar, int **icfg, int **wbx, int j, int d, int mid)
{ 
  int WB7;
  int i;

  i = j - d;
  
  WB7 = wbx[i+mid][mid] + wbx[j][d-mid-1];
  
  return WB7;
}
void
trace_WB7(FILE *outf, int **wbx, int j, int d, int mid,
	  struct tracekn_s *curr_tr, struct traceknstack_s *dolist, int traceback)
{
  int i;
  
  i = j - d;

  if (traceback) fprintf(outf," (WB7) trace wbx  %d %d %d\n", i+mid, mid, wbx[i+mid][mid]);
  if (traceback) fprintf(outf," (WB7) trace wbx  %d\n", wbx[j][d-mid-1]);
  
  PushTraceknstack(dolist, AttachTracekn(curr_tr, i, i+mid, i+(int)(mid/2), 
					 i+(int)(mid/2)+1, dpcB, dpcS));
  PushTraceknstack(dolist, AttachTracekn(curr_tr, i+mid+1, j, i+mid+1+(int)((d-mid-1)/2),
					 i+mid+2+(int)((d-mid-1)/2), dpcB, dpcS));
}

/* (WB8)
 * 2 VX 
 */

/* WB8_1 to WB8_4
 * contiguous coaxial pairs (i,i+mid) (i+mid+1,j)
 * ...add stack[i+mid][i][i+mid+1][j]
 * check also for external danglings
 */

/* (WB8.1)__/ no danglings / coaxial = stack + 1  */
int
WB8_1(ESL_DSQ *s, int len, struct rnapar_2 *rnapar, int **icfg, int **vx, int j, int d, int mid)
{ 
  int WB8_1;
  int i;

  i = j - d;
  
  WB8_1 = 2*rnapar->P10 
    + wsf*(icfg[idxPS(s[i+mid],s[i])][idxPS(s[i+mid+1],s[j])] + INTSCALE)
    + vx[i+mid][mid] + vx[j][d-mid-1];
  
  return WB8_1;
}
void
trace_WB8_1(FILE *outf, int **vx, int j, int d, int mid,
	    struct tracekn_s *curr_tr, struct traceknstack_s *dolist, int traceback)
{
  int i;
  
  i = j - d;

  if (traceback) fprintf(outf," (WB8.1) trace vx %d %d %d\n", 
			 i+mid, mid, vx[i+mid][mid]);
  if (traceback) fprintf(outf," (WB8.1) trace vx %d %d %d\n", 
			 j, d-mid-1, vx[j][d-mid-1]);
  
  PushTraceknstack(dolist, AttachTracekn(curr_tr, i, i+mid, i+(int)(mid/2), 
					 i+(int)(mid/2)+1, dpcP, dpcS));
  PushTraceknstack(dolist, AttachTracekn(curr_tr, i+mid+1, j, i+mid+1+(int)((d-mid-1)/2),
					 i+mid+2+(int)((d-mid-1)/2), dpcP, dpcS));
}

/* (WB8.2)__/ i dangles off i+1,i+mid */
int
WB8_2(ESL_DSQ *s, int len, struct rnapar_2 *rnapar, int **icfg, int **vx, int j, int d, int mid)
{ 
  int WB8_2;
  int i;

  i = j - d;
  
  WB8_2 = 2*rnapar->P10 + rnapar->P6    
    + wsf*icfg[idxL(s[i])][idxP(s[i+1],s[i+mid])]
    + wsf*icfg[idxPS(s[i+mid],s[i+1])][idxPS(s[i+mid+1],s[j])]
    + vx[i+mid][mid-1] + vx[j][d-mid-1];

  return WB8_2;
}
void
trace_WB8_2(FILE *outf, int **vx, int j, int d, int mid,
	    struct tracekn_s *curr_tr, struct traceknstack_s *dolist, int traceback)
{
  int i;
  
  i = j - d;

  if (traceback) fprintf(outf," (WB8.2) dangle i, trace vx %d %d %d\n", 
			 i+mid, mid-1, vx[i+mid][mid-1]);
  if (traceback) fprintf(outf," (WB8.2) trace vx %d %d %d\n", 
			 j, d-mid-1, vx[j][d-mid-1]);
  
  PushTraceknstack(dolist, AttachTracekn(curr_tr, i+1, i+mid, i+(int)(mid/2)+1, 
					 i+(int)(mid/2)+2, dpcP, dpcS));
  PushTraceknstack(dolist,AttachTracekn(curr_tr, i+mid+1, j, i+mid+1+(int)((d-mid-1)/2),
					i+mid+2+(int)((d-mid-1)/2), dpcP, dpcS));
}

/* (WB8.3)__/ j dangles off i+mid+1,j-1 */
int
WB8_3(ESL_DSQ *s, int len, struct rnapar_2 *rnapar, int **icfg, int **vx, int j, int d, int mid)
{ 
  int WB8_3;
  int i;

  i = j - d;
  
  if (d-mid > 1)
    WB8_3 = 2*rnapar->P10 + rnapar->P6  
      + wsf*icfg[idxR(s[j])][idxP(s[i+mid+1],s[j-1])]
      + wsf*icfg[idxPS(s[i+mid],s[i])][idxPS(s[i+mid+1],s[j-1])]
      + vx[i+mid][mid] + vx[j-1][d-mid-2];
  else
    WB8_3 = -BIGINT;
  
  return WB8_3;
}
void
trace_WB8_3(FILE *outf, int **vx, int j, int d, int mid,
	    struct tracekn_s *curr_tr, struct traceknstack_s *dolist, int traceback)
{
  int i;
  
  i = j - d;

  if (d-mid > 1){
    if (traceback) fprintf(outf," (WB8.3) trace vx %d %d %d\n", 
			   i+mid, mid, vx[i+mid][mid]);
    if (traceback) fprintf(outf," (WB8.3) dangle j, trace vx %d %d %d\n", 
			   j-1, d-mid-2, vx[j-1][d-mid-2]);
    
    PushTraceknstack(dolist, AttachTracekn(curr_tr, i, i+mid, i+(int)(mid/2), 
					   i+(int)(mid/2)+1, dpcP, dpcS));
    PushTraceknstack(dolist, AttachTracekn(curr_tr, i+mid+1, j-1, i+mid+1+(int)((d-mid-2)/2),
					   i+mid+2+(int)((d-mid-2)/2), dpcP, dpcS));
  }
}

/* (WB8.4)__/ i dangles off i+1,i+mid / j dangles off i+mid+1,j-1 */
int
WB8_4(ESL_DSQ *s, int len, struct rnapar_2 *rnapar, int **icfg, int **vx, int j, int d, int mid)
{ 
  int WB8_4;
  int i;

  i = j - d;
  
  if (d-mid > 1)
    WB8_4 = 2*rnapar->P10 + 2*rnapar->P6 
      + wsf*icfg[idxL(s[i])][idxP(s[i+1],s[i+mid])]
      + wsf*icfg[idxR(s[j])][idxP(s[i+mid+1],s[j-1])]
      + wsf*icfg[idxPS(s[i+mid],s[i+1])][idxPS(s[i+mid+1],s[j-1])]
      + vx[i+mid][mid-1] + vx[j-1][d-mid-2];
  else
    WB8_4 = -BIGINT;
  
  return WB8_4;
}
void
trace_WB8_4(FILE *outf, int **vx, int j, int d, int mid,
	    struct tracekn_s *curr_tr, struct traceknstack_s *dolist, int traceback)
{
  int i;
  
  i = j - d;

  if (d-mid > 1){
    if (traceback) fprintf(outf," (WB8.4) dangle i, trace vx %d %d %d\n", 
			   i+mid, mid-1, vx[i+mid][mid-1]);
    if (traceback) fprintf(outf," (WB8.4) dangle j, trace vx %d %d %d\n", 
			   j-1, d-mid-2, vx[j-1][d-mid-2]);
    
    PushTraceknstack(dolist, AttachTracekn(curr_tr, i+1, i+mid, i+(int)(mid/2)+1, 
					   i+(int)(mid/2)+2, dpcP, dpcS));
    PushTraceknstack(dolist,
		     AttachTracekn(curr_tr, i+mid+1, j-1, i+mid+1+(int)((d-mid-2)/2),
				   i+mid+2+(int)((d-mid-2)/2), dpcP, dpcS));
  }
}

/* WB8_5 to WB8_8
 * non-contiguous coaxial pairs (i,i+mid) (i+mid+2,j)
 *                  ...add stack[i+mid][i][i+mid+2][j]
 * check also for external danglings
 */

/* (WB8.5)__/ no danglings (do not add +1 becuase the stacking is not perfect) */
int
WB8_5(ESL_DSQ *s, int len, struct rnapar_2 *rnapar, int **icfg, int **vx, int j, int d, int mid)
{ 
  int WB8_5;
  int i;

  i = j - d;
  
  if (d-mid > 1)
    WB8_5 = 2*rnapar->P10 + rnapar->P6 
      + wsf*icfg[idxPS(s[i+mid],s[i])][idxPS(s[i+mid+2],s[j])]
      + vx[i+mid][mid] + vx[j][d-mid-2];
   else
    WB8_5 = -BIGINT;
   
  return WB8_5;
}
void
trace_WB8_5(FILE *outf, int **vx, int j, int d, int mid,
	    struct tracekn_s *curr_tr, struct traceknstack_s *dolist, int traceback)
{
  int i;
  
  i = j - d;

  if (d-mid > 1){
    if (traceback) fprintf(outf," (WB8.5) trace vx %d %d %d\n", 
			   i+mid, mid, vx[i+mid][mid]);
    if (traceback) fprintf(outf," (WB8.5) trace vx %d %d %d\n", 
			   j, d-mid-2, vx[j][d-mid-2]);
    
    PushTraceknstack(dolist, AttachTracekn(curr_tr, i, i+mid, i+(int)(mid/2), 
					   i+(int)(mid/2)+1, dpcP, dpcS));
    PushTraceknstack(dolist, AttachTracekn(curr_tr, i+mid+2, j, i+mid+2+(int)((d-mid-2)/2),
					   i+mid+3+(int)((d-mid-2)/2), dpcP, dpcS));
  }
}

/* (WB8.6)__/ i dangles off i+1,i+mid */
int
WB8_6(ESL_DSQ *s, int len, struct rnapar_2 *rnapar, int **icfg, int **vx, int j, int d, int mid)
{ 
  int WB8_6;
  int i;

  i = j - d;
  
  if (d-mid > 1)
    WB8_6 = 2*rnapar->P10 + 2*rnapar->P6
      + wsf*icfg[idxL(s[i])][idxP(s[i+1],s[i+mid])]
      + wsf*icfg[idxPS(s[i+mid],s[i+1])][idxPS(s[i+mid+2],s[j])]
      + vx[i+mid][mid-1] + vx[j][d-mid-2];
  else
    WB8_6 = -BIGINT;

  return WB8_6;
}
void
trace_WB8_6(FILE *outf, int **vx, int j, int d, int mid,
	    struct tracekn_s *curr_tr, struct traceknstack_s *dolist, int traceback)
{
  int i;
  
  i = j - d;

  if (d-mid > 1){
    if (traceback) fprintf(outf," (WB8.6) dangle i, trace vx %d %d %d\n", 
			   i+mid, mid-1, vx[i+mid][mid-1]);
    if (traceback) fprintf(outf," (WB8.6) trace vx %d %d %d\n", 
			   j, d-mid-2, vx[j][d-mid-2]);
    
    PushTraceknstack(dolist, AttachTracekn(curr_tr, i+1, i+mid, i+(int)(mid/2)+1, 
					   i+(int)(mid/2)+2, dpcP, dpcS));
    PushTraceknstack(dolist,AttachTracekn(curr_tr, i+mid+2, j, i+mid+2+(int)((d-mid-2)/2),
					  i+mid+3+(int)((d-mid-2)/2), dpcP, dpcS));
  }
}

/* (WB8.7)__/ j dangles off i+mid+1,j-1 */
int
WB8_7(ESL_DSQ *s, int len, struct rnapar_2 *rnapar, int **icfg, int **vx, int j, int d, int mid)
{ 
  int WB8_7;
  int i;

  i = j - d;
  
  if (d-mid > 2)
    WB8_7 = 2*rnapar->P10 + 2*rnapar->P6
      + wsf*icfg[idxR(s[j])][idxP(s[i+mid+2],s[j-1])]
      + wsf*icfg[idxPS(s[i+mid],s[i])][idxPS(s[i+mid+2],s[j-1])]
      + vx[i+mid][mid] + vx[j-1][d-mid-3];
  else
    WB8_7 = -BIGINT;
  
  return WB8_7;
}
void
trace_WB8_7(FILE *outf, int **vx, int j, int d, int mid,
	    struct tracekn_s *curr_tr, struct traceknstack_s *dolist, int traceback)
{
  int i;
  
  i = j - d;

  if (d-mid > 2){
    if (traceback) fprintf(outf," (WB8.7) trace vx %d %d %d\n", 
			   i+mid, mid, vx[i+mid][mid]);
    if (traceback) fprintf(outf," (WB8.7) dangle j, trace vx %d %d %d\n", 
			   j-1, d-mid-3, vx[j-1][d-mid-3]);
    
    PushTraceknstack(dolist, AttachTracekn(curr_tr, i, i+mid, i+(int)(mid/2), 
					   i+(int)(mid/2)+1, dpcP, dpcS));
    PushTraceknstack(dolist, AttachTracekn(curr_tr, i+mid+2, j-1, i+mid+2+(int)((d-mid-3)/2),
					   i+mid+3+(int)((d-mid-3)/2), dpcP, dpcS));
  }
}

/* (WB8.8)__/ i dangles off i+1,i+mid / j dangles off i+mid+2,j-1 */
int
WB8_8(ESL_DSQ *s, int len, struct rnapar_2 *rnapar, int **icfg, int **vx, int j, int d, int mid)
{ 
  int WB8_8;
  int i;

  i = j - d;
  
  if (d-mid > 2)
    WB8_8 = 2*rnapar->P10 + 3*rnapar->P6
      + wsf*icfg[idxL(s[i])][idxP(s[i+1],s[i+mid])]
      + wsf*icfg[idxR(s[j])][idxP(s[i+mid+2],s[j-1])]
      + wsf*icfg[idxPS(s[i+mid],s[i+1])][idxPS(s[i+mid+2],s[j-1])]
      + vx[i+mid][mid-1] + vx[j-1][d-mid-3];
  else
    WB8_8 = -BIGINT;
  
  return WB8_8;
}
void
trace_WB8_8(FILE *outf, int **vx, int j, int d, int mid,
	    struct tracekn_s *curr_tr, struct traceknstack_s *dolist, int traceback)
{
  int i;
  
  i = j - d;

  if (d-mid > 1){
    if (traceback) fprintf(outf," (WB8.8) dangle i, trace vx %d %d %d\n", 
			   i+mid, mid-1, vx[i+mid][mid-1]);
    if (traceback) fprintf(outf," (WB8.8) dangle j, trace vx %d %d %d\n", 
			   j-1, d-mid-3, vx[j-1][d-mid-3]);
    
    PushTraceknstack(dolist, AttachTracekn(curr_tr, i+1, i+mid, i+(int)(mid/2)+1, 
					   i+(int)(mid/2)+2, dpcP, dpcS));
    PushTraceknstack(dolist,
		     AttachTracekn(curr_tr, i+mid+2, j-1, i+mid+2+(int)((d-mid-3)/2),
				   i+mid+3+(int)((d-mid-3)/2), dpcP, dpcS));
  }
}
 
/* BIFURCATIONS (with hole) (2 types of diagrams WB9-WB10)
 * (add P13 penalty for wbx)
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
/* (WB9)__/ 2 WHX
 */
int
WB9(ESL_DSQ *s, int len, struct rnapar_2 *rnapar, int **icfg, int ****whx, int j, int d, int mid, int mid1, int mid2)
{ 
  int WB9;
  int i;

  i = j - d;
  
  WB9 = rnapar->P13 + 3*rnapar->P6P
    + whx[i+mid][mid][mid1][mid2] 
    + whx[j][d-mid1-2][mid-mid1-mid2-3][d-mid-3];
  
  return WB9;
}
void
trace_WB9(FILE *outf, int ****whx, int j, int d, int mid, int mid1, int mid2,
	 struct tracekn_s *curr_tr, struct traceknstack_s *dolist, int traceback)
{
  int i;
  
  i = j - d;

  if (traceback) fprintf(outf," (WB9) trace whx %d %d %d %d %d\n", i+mid, mid, mid1, mid2,
			 whx[i+mid][mid][mid1][mid2]);
  if (traceback) fprintf(outf," (WB9) trace whx %d %d %d %d %d\n", 
			 j, d-mid1-2, mid-mid1-mid2-3, d-mid-3,
			 whx[j][d-mid1-2][mid-mid1-mid2-3][d-mid-3]);
  
  PushTraceknstack(dolist, AttachTracekn(curr_tr, i, i+mid, i+mid1, 
					 i+mid-mid2, dpcS, dpcS));
  PushTraceknstack(dolist, AttachTracekn(curr_tr, i+mid1+2, j, i+mid-mid2-1,
					 i+mid+3, dpcS, dpcS));
}
 
/* (WB10)__/ 1 YHX 1 YHX 
 * coaxial pairs (i+mid1,i+mid-mid2) (i+mid-mid2-1,i+mid+3)
 * ....add stack[i+mid-mid2-1][i+mid+3][i+mid-mid2][i+mid1]
 */
/* (WB10.1)__/ contiguous coaxial: i+mid-mid2-1, i+mid-mid2     /
             / i+mid1+1 dangles right off i+mid-mid2,i+mid1     /
             / i+mid+2  dangles left  off i+mid+3,i+mid-mid2-1  */
int
WB10_1(ESL_DSQ *s, int len, struct rnapar_2 *rnapar, int **icfg, int ****yhx, 
       int j, int d, int mid, int mid1, int mid2)
{ 
  int WB10_1;
  int i;

  i = j - d;
  
  WB10_1 = rnapar->P13 + 2*rnapar->P10P + 3*rnapar->P6P
    + wkn*icfg[idxR(s[i+mid1+1])][idxP(s[i+mid-mid2],s[i+mid1])]     
    + wkn*icfg[idxL(s[i+mid+2])][idxP(s[i+mid+3],s[i+mid-mid2-1])]     
    + wkn*icfg[idxPS(s[i+mid-mid2-1],s[i+mid+3])]
    [idxPS(s[i+mid-mid2],s[i+mid1])] 
    + yhx[i+mid][mid][mid1][mid2] 
    + yhx[j][d-mid1-2][mid-mid1-mid2-3][d-mid-3];
  
  return WB10_1;
}
void
trace_WB10_1(FILE *outf, int ****yhx, int j, int d, int mid, int mid1, int mid2,
	     struct tracekn_s *curr_tr, struct traceknstack_s *dolist, int traceback)
{
  int i;
  
  i = j - d;

  if (traceback) fprintf(outf," (WB10.1) trace yhx %d %d %d %d %d\n", i+mid, mid, mid1, mid2,
			 yhx[i+mid][mid][mid1][mid2]);
  if (traceback) fprintf(outf," (WB10.1) trace yhx %d %d %d %d %d\n", 
			 j, d-mid1-2, mid-mid1-mid2-3, d-mid-3, 
			 yhx[j][d-mid1-2][mid-mid1-mid2-3][d-mid-3]);
  
  PushTraceknstack(dolist, AttachTracekn(curr_tr, i, i+mid, 
					 i+mid1, i+mid-mid2, dpcS, dpcP));
  PushTraceknstack(dolist, AttachTracekn(curr_tr, i+mid1+2, j, 
					 i+mid-mid2-1, i+mid+3, dpcS, dpcP));
}

/* (WB10.2)__/ non-contiguous coaxial: i+mid-mid2-2, i+mid-mid2 /
            / i+mid1+1 dangles right off i+mid-mid2,i+mid1     /
            / i+mid+2  dangles left  off i+mid+3,i+mid-mid2-1  */
int
WB10_2(ESL_DSQ *s, int len, struct rnapar_2 *rnapar, int **icfg, int ****yhx, 
       int j, int d, int mid, int mid1, int mid2)
{ 
  int WB10_2;
  int i;

  i = j - d;
  
  if (mid-mid1-mid2 > 3)
    WB10_2 = rnapar->P13 + 2*rnapar->P10P + 4*rnapar->P6P
      + wkn*icfg[idxR(s[i+mid1+1])][idxP(s[i+mid-mid2],s[i+mid1])]     
      + wkn*icfg[idxL(s[i+mid+2])][idxP(s[i+mid+3],s[i+mid-mid2-2])]     
      + wkn*icfg[idxPS(s[i+mid-mid2-2],s[i+mid+3])]
                [idxPS(s[i+mid-mid2],s[i+mid1])] 
      + yhx[i+mid][mid][mid1][mid2] 
      + yhx[j][d-mid1-2][mid-mid1-mid2-4][d-mid-3];
  else
    WB10_2 = -BIGINT;
  
  return WB10_2;
}
void
trace_WB10_2(FILE *outf, int ****yhx, int j, int d, int mid, int mid1, int mid2,
	     struct tracekn_s *curr_tr, struct traceknstack_s *dolist, int traceback)
{
  int i;
  
  i = j - d;

  if (mid-mid1-mid2 > 3){
    if (traceback) fprintf(outf," (WB10.2) trace yhx %d %d %d %d %d\n", i+mid, mid, mid1, mid2,
			   yhx[i+mid][mid][mid1][mid2]);
    if (traceback) fprintf(outf," (WB10.2) trace yhx %d %d %d %d %d\n", 
			   j, d-mid1-2, mid-mid1-mid2-4, d-mid-3, 
			   yhx[j][d-mid1-2][mid-mid1-mid2-4][d-mid-3]);
    
    PushTraceknstack(dolist, AttachTracekn(curr_tr, i, i+mid, 
					   i+mid1, i+mid-mid2, dpcS, dpcP));
    PushTraceknstack(dolist, AttachTracekn(curr_tr, i+mid1+2, j, 
					   i+mid-mid2-2, i+mid+3, dpcS, dpcP));
  }
}

