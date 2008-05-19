/* vxgraphs.c 
 *
 * includes functions to calculate all the diagrams 
 * that fill and traceback the no-hole matrix:
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
#include "protovx.h"
#include "squid.h"

/*  FUNCTION: V1
 *   
 *  Purpose: calculate IRREDUCIBLE SURFACES of  O(1) (hairpin loops) 
 *   
 *  Return: V1           
 */
int 
V1(int *s, int len, struct rnapar_2 *rnapar, int **icfg, int j, int d)
{ 
  int V1;
  
  V1 =  F1(s, len, rnapar, icfg, j, d);
  
  return V1;
}
void
trace_V1(FILE *outf, int **vx, int j, int d, 
	 struct tracekn_s *curr_tr, struct traceknstack_s *dolist, int traceback)
{
     if (traceback) fprintf(outf," (v1) hairpin %d, %d \n", d, vx[j][d]);
}

/*  FUNCTION: V2
 *   
 *  Purpose: calculate IRREDUCIBLE SURFACES of  O(2) 
 *   
 *  Return: V2           
 */
int 
V2(int *s, int len, struct rnapar_2 *rnapar,
   int **icfg, int **vx, int j, int d, int mid1, int mid2)
{ 
  int V2;
  
    V2 = F2(s, len, rnapar, icfg, j, d, mid1, mid2)
      + vx[j-mid2][d-mid1-mid2];
  return V2;
}

void
trace_V2(FILE *outf, int **vx, int *s, int len, struct rnapar_2 *rnapar, int **icfg, 
	 int j, int d, int mid1, int mid2, 
	 struct tracekn_s *curr_tr, struct traceknstack_s *dolist, int traceback)
{  
  int i;

  i = j - d;

  if (traceback) fprintf(outf," (v2) IS2 (%d,%d) %d\n", 
			 mid1, mid2, F2(s, len, rnapar, icfg, j, d, mid1, mid2));
  if (traceback) fprintf(outf," trace vx %d %d %d\n", 
			 j-mid2, d-mid1-mid2, vx[j-mid2][d-mid1-mid2]);
  
  curr_tr->type1 = dpcP;
  
  PushTraceknstack(dolist, AttachTracekn(curr_tr, i+mid1, j-mid2,
					 i+mid1+(int)((d-mid1-mid2)/2),
					 i+mid1+(int)((d-mid1-mid2)/2)+1,
					 dpcP, dpcS));
}

/* 2 NO-HOLE STRUCTURES (four types of diagrams: V3-V6)
 * 
 * one connects i+1,k   (i+1+mid, mid)
 * one connects k+1,j-1 (j-1, d-mid-3)
 */

/*  FUNCTIONS: V3 (V3_1, V3_2, V3_3, V3_4)
 *   
 *  Purpose: calculate diagrams with 2 WBX 
 *  
 *  Return: int V3_X           
 */
 
/* (V3_1)__/no danglings */
int 
V3_1(int *s, int len, struct rnapar_2 *rnapar, int **icfg, int **wbx, int j, int d, int mid)
{
  int V3_1;
  int    i;

  i = j - d;
  
  if (mid < d-3)
    V3_1 = rnapar->P10 + rnapar->P5 + wbx[i+1+mid][mid] + wbx[j-1][d-mid-3];
  else
    V3_1 = -BIGINT;

  return V3_1;
}
void
trace_V3_1(FILE *outf, int **wbx, int j, int d, int mid, 
	   struct tracekn_s *curr_tr, struct traceknstack_s *dolist, int traceback)
{  
  int i;

  i = j - d;

  if (mid < d-3){
    if (traceback) fprintf(outf," (v3_1) multiloop wbx %d %d %d\n", 
			   i+1+mid, mid, wbx[i+1+mid][mid]);
    if (traceback) fprintf(outf," (v3_1) multiloop wbx %d %d %d\n", 
			   j-1, d-mid-3, wbx[j-1][d-mid-3]);
    
    PushTraceknstack(dolist, AttachTracekn(curr_tr, i+1, i+1+mid,
					   i+1+(int)(mid/2),
					   i+2+(int)(mid/2), dpcB, dpcS));
    PushTraceknstack(dolist, AttachTracekn(curr_tr, i+2+mid, j-1,
					   i+2+mid+(int)((d-mid-3)/2),
					   i+3+mid+(int)((d-mid-3)/2), dpcB, dpcS));  
  }
}

/*  (V3_2)__/ i+1 dangles off i,j */           
int 
V3_2(int *s, int len, struct rnapar_2 *rnapar, int **icfg, int **wbx, int j, int d, int mid)
{
  int V3_2;
  int    i;

  i = j - d;

  if (mid > 1)
    V3_2 = rnapar->P10 + rnapar->P5 + rnapar->P6 
      + wsf*icfg[idxR(s[i+1])][idxP(s[j],s[i])]
      + wbx[i+1+mid][mid-1] + wbx[j-1][d-mid-3];
  else
    V3_2 = -BIGINT;

  return V3_2;
}
void
trace_V3_2(FILE *outf, int **wbx, int j, int d, int mid, 
	   struct tracekn_s *curr_tr, struct traceknstack_s *dolist, int traceback)
{  
  int i;

  i = j - d;

  if (mid > 1){
    if (traceback) fprintf(outf," (v3_2) multiloop wbx %d %d %d\n", 
			   i+1+mid, mid-1, wbx[i+1+mid][mid-1]);
    if (traceback) fprintf(outf," (v3_2) multiloop wbx %d %d %d\n", 
			   j-1, d-mid-3, wbx[j-1][d-mid-3]);
    
    PushTraceknstack(dolist, AttachTracekn(curr_tr, i+2, i+1+mid,
					   i+2+(int)((mid-1)/2),
					   i+3+(int)((mid-1)/2), dpcB, dpcS));
    PushTraceknstack(dolist, AttachTracekn(curr_tr, i+2+mid, j-1,
					   i+2+mid+(int)((d-mid-3)/2),
					   i+3+mid+(int)((d-mid-3)/2), dpcB, dpcS));  
  }
}

/*  (V3_3)__/ j-1 dangles off i,j */
int 
V3_3(int *s, int len, struct rnapar_2 *rnapar, int **icfg, int **wbx, int j, int d, int mid)
{
  int V3_3;
  int    i;

  i = j - d;

  if (mid < d-4)
    V3_3 = rnapar->P10 + rnapar->P5 + rnapar->P6 
      + wsf*icfg[idxL(s[j-1])][idxP(s[j],s[i])]
      + wbx[i+1+mid][mid] + wbx[j-2][d-mid-4];
  else
    V3_3 = -BIGINT;

  return V3_3;
}
void
trace_V3_3(FILE *outf, int **wbx, int j, int d, int mid, 
	   struct tracekn_s *curr_tr, struct traceknstack_s *dolist, int traceback)
{  
  int i;

  i = j - d;

  if (mid < d-4){
    if (traceback) fprintf(outf," (v3_3) multiloop wbx %d %d %d\n", 
			   i+1+mid, mid, wbx[i+1+mid][mid]);
    if (traceback) fprintf(outf," (v3_3) multiloop wbx %d %d %d\n", 
			   j-2, d-mid-4, wbx[j-2][d-mid-4]);
    
    PushTraceknstack(dolist, AttachTracekn(curr_tr, i+1, i+1+mid,
					   i+1+(int)((mid)/2),
					   i+2+(int)((mid)/2), dpcB, dpcS));
    PushTraceknstack(dolist, AttachTracekn(curr_tr, i+2+mid, j-2,
					   i+2+mid+(int)((d-mid-4)/2),
					   i+3+mid+(int)((d-mid-4)/2), dpcB, dpcS)); 
  }
}

/*  (V3_4)__/ i+1 and j-1 dangle off i,j */
int 
V3_4(int *s, int len, struct rnapar_2 *rnapar, int **icfg, int **wbx, int j, int d, int mid)
{
  int V3_4;
  int    i;

  i = j - d;

  if (mid > 1 && mid < d-4)
    V3_4 = rnapar->P10 + rnapar->P5 + 2*rnapar->P6 
      + wsf*icfg[idxR(s[i+1])][idxP(s[j],s[i])] 
      + wsf*icfg[idxL(s[j-1])][idxP(s[j],s[i])]
      + wbx[i+1+mid][mid-1] + wbx[j-2][d-mid-4];
  else
    V3_4 = -BIGINT;
  
  return V3_4;
}
void
trace_V3_4(FILE *outf, int **wbx, int j, int d, int mid, 
	   struct tracekn_s *curr_tr, struct traceknstack_s *dolist, int traceback)
{  
  int i;

  i = j - d;

  if (mid > 1 && mid < d-4){
    if (traceback) fprintf(outf," (v3_4) multiloop wbx %d %d %d\n", 
			   i+1+mid, mid-1, wbx[i+1+mid][mid-1]);
    if (traceback) fprintf(outf," (v3_4) multiloop wbx %d %d %d\n",
			   j-2, d-mid-4, wbx[j-2][d-mid-4]);
    
    PushTraceknstack(dolist, AttachTracekn(curr_tr, i+2, i+1+mid,
					   i+2+(int)((mid-1)/2),
					   i+3+(int)((mid-1)/2), dpcB, dpcS));
    PushTraceknstack(dolist, AttachTracekn(curr_tr, i+2+mid, j-2,
					   i+2+mid+(int)((d-mid-4)/2),		
					   i+3+mid+(int)((d-mid-4)/2), dpcB, dpcS));  
  }
}


/*  FUNCTIONS: V4__/ 1 VX 1 WBX
 */

/* V4_1, V4_2, V4_3, V4_4
 *   
 * contiguous coaxial pairs: (i,j) (i+1,i+1+mid)
 *                   ...add stack[i][j][i+1][i+1+mid] (+ 1 if no danglings)
 *   
 *  Return: int V4_X          
 */

/*  (V4_1)__/ no danglings / coaxial = stack + 1 */
int 
V4_1(int *s, int len, struct rnapar_2 *rnapar, int **icfg, int **wbx, int **vx, int j, int d, int mid)
{
  int V4_1;
  int    i;

  i = j - d;

  if (mid < d-3)
    V4_1 = 2*rnapar->P10 + rnapar->P5 
      + wsf*(icfg[idxPS(s[i],s[j])][idxPS(s[i+1],s[i+1+mid])] + INTSCALE)
      + vx[i+1+mid][mid] + wbx[j-1][d-mid-3];
  else
    V4_1 = -BIGINT;

  return V4_1;
}
void
trace_V4_1(FILE *outf, int **wbx, int **vx, int j, int d, int mid, 
	   struct tracekn_s *curr_tr, struct traceknstack_s *dolist, int traceback)
{  
  int i;

  i = j - d;

  if (mid < d-3){
    if (traceback) fprintf(outf," (v4_1) multiloop, trace vx  %d %d %d\n", 
			   i+1+mid, mid, vx[i+1+mid][mid]);
    if (traceback) fprintf(outf," (v4_1) multiloop, trace wbx %d %d %d\n", 
			   j-1, d-mid-3, wbx[j-1][d-mid-3]);
    
    PushTraceknstack(dolist, AttachTracekn(curr_tr, i+1, i+1+mid, i+1+(int)((mid)/2),
					   i+2+(int)((mid)/2), dpcP, dpcS));
    PushTraceknstack(dolist, AttachTracekn(curr_tr, i+2+mid, j-1, i+2+mid+(int)((d-mid-3)/2),
					   i+3+mid+(int)((d-mid-3)/2), dpcB, dpcS));  
  }
}

/* (V4_2)__/ j-1 dangles off i,j */
int 
V4_2(int *s, int len, struct rnapar_2 *rnapar, int **icfg, int **wbx, int **vx, int j, int d, int mid)
{
  int V4_2;
  int    i;

  i = j - d;

  if (mid < d-4)
    V4_2 = 2*rnapar->P10 + rnapar->P5 + rnapar->P6 
      + wsf*icfg[idxL(s[j-1])][idxP(s[j],s[i])]
      + wsf*icfg[idxPS(s[i],s[j])][idxPS(s[i+1],s[i+1+mid])] 
      + vx[i+1+mid][mid] + wbx[j-2][d-mid-4];
  else
    V4_2 = -BIGINT;

  return V4_2;
}
void
trace_V4_2(FILE *outf, int **wbx, int **vx, int j, int d, int mid, 
	   struct tracekn_s *curr_tr, struct traceknstack_s *dolist, int traceback)
{  
  int i;

  i = j - d;

  if (mid < d-4){
    if (traceback) fprintf(outf," (v4_2) multiloop, trace vx  %d %d %d\n", 
			   i+1+mid, mid, vx[i+1+mid][mid]);
    if (traceback) fprintf(outf," (v4_2) multiloop, trace wbx %d %d %d\n", 
			   j-2, d-mid-4, wbx[j-2][d-mid-4]);
    
    PushTraceknstack(dolist, AttachTracekn(curr_tr, i+1, i+1+mid, i+1+(int)((mid)/2),
					   i+2+(int)((mid)/2), dpcP, dpcS));
    PushTraceknstack(dolist, AttachTracekn(curr_tr, i+2+mid, j-2, i+2+mid+(int)((d-mid-4)/2),
					   i+3+mid+(int)((d-mid-4)/2), dpcB, dpcS));  
  }
}

/*  (V4_3)__/ i+mid+2 dangles off i+1,i+mid+1 */
int
V4_3(int *s, int len, struct rnapar_2 *rnapar, int **icfg, int **wbx, int **vx, int j, int d, int mid)
{
  int V4_3;
  int    i;

  i = j - d;

  if (mid < d-4)
    V4_3 = 2*rnapar->P10 + rnapar->P5 + rnapar->P6
      + wsf*icfg[idxR(s[i+mid+2])][idxP(s[i+1],s[i+mid+1])]
      + wsf*icfg[idxPS(s[i],s[j])][idxPS(s[i+1],s[i+1+mid])] 
      + vx[i+1+mid][mid] + wbx[j-1][d-mid-4];
  else
    V4_3 = -BIGINT;
  
  return V4_3;
}
void
trace_V4_3(FILE *outf, int **wbx, int **vx, int j, int d, int mid, 
	   struct tracekn_s *curr_tr, struct traceknstack_s *dolist, int traceback)
{  
  int i;

  i = j - d;

  if (mid < d-4){
    if (traceback) fprintf(outf," (v4_3) multiloop, trace vx  %d %d %d\n", 
			   i+1+mid, mid, vx[i+1+mid][mid]);
    if (traceback) fprintf(outf," (v4_3) multiloop, trace wbx %d %d %d\n", 
			   j-1, d-mid-4, wbx[j-1][d-mid-4]);
    
    PushTraceknstack(dolist, AttachTracekn(curr_tr, i+1, i+1+mid, i+1+(int)((mid)/2),
					   i+2+(int)((mid)/2), dpcP, dpcS));
    PushTraceknstack(dolist, AttachTracekn(curr_tr, i+3+mid, j-1, i+3+mid+(int)((d-mid-4)/2),
					   i+4+mid+(int)((d-mid-4)/2), dpcB, dpcS));  
  }
}

/*  (V4_4)__/ j-1 dangles off i,j / i+mid+2 dangles off i+1,i+mid+1 */
int 
V4_4(int *s, int len, struct rnapar_2 *rnapar, int **icfg, int **wbx, int **vx, int j, int d, int mid)
{
  int V4_4;
  int    i;

  i = j - d;

  if (mid < d-5)
    V4_4 = 2*rnapar->P10 + rnapar->P5 + 2*rnapar->P6 
      + wsf*icfg[idxL(s[j-1])][idxP(s[j],s[i])]
      + wsf*icfg[idxR(s[i+mid+2])][idxP(s[i+1],s[i+mid+1])]
      + wsf*icfg[idxPS(s[i],s[j])][idxPS(s[i+1],s[i+1+mid])] 
      + vx[i+1+mid][mid] + wbx[j-2][d-mid-5];
  else
    V4_4 = -BIGINT;
  
  return V4_4;
}
void
trace_V4_4(FILE *outf, int **wbx, int **vx, int j, int d, int mid, 
	   struct tracekn_s *curr_tr, struct traceknstack_s *dolist, int traceback)
{  
  int i;

  i = j - d;

  if (mid < d-5){
    if (traceback) fprintf(outf," (v4_4) multiloop, trace vx  %d %d %d\n", 
			   i+1+mid, mid, vx[i+1+mid][mid]);
    if (traceback) fprintf(outf," (v4_4) multiloop, trace wbx %d %d %d\n", 
			   j-2, d-mid-5, wbx[j-2][d-mid-5]);
    
    PushTraceknstack(dolist, AttachTracekn(curr_tr, i+1, i+1+mid, i+1+(int)((mid)/2),
					 i+2+(int)((mid)/2), dpcP, dpcS));
    PushTraceknstack(dolist, AttachTracekn(curr_tr, i+3+mid, j-2, i+3+mid+(int)((d-mid-5)/2),
					   i+4+mid+(int)((d-mid-5)/2), dpcB, dpcS));  
  }
}

/* V4_5, V4_6, V4_7, V4_8
 *   
 * non-contiguous coaxial pairs: (i,j) (i+2,i+2+mid)
 *                        ...add stack[i][j][i+1][i+2+mid] 
 *   
 *  Return: int V4_X          
 */

/*  (V4_5)__/ no danglings  */
int 
V4_5(int *s, int len, struct rnapar_2 *rnapar, int **icfg, int **wbx, int **vx, int j, int d, int mid)
{
  int V4_5;
  int    i;

  i = j - d;

  if (mid < d-4)
    V4_5 = 2*rnapar->P10 + rnapar->P5 + rnapar->P6
      + wsf*icfg[idxPS(s[i],s[j])][idxPS(s[i+2],s[i+2+mid])]
      + vx[i+2+mid][mid] + wbx[j-1][d-mid-4];
  else
    V4_5 = -BIGINT;

  return V4_5;
}
void
trace_V4_5(FILE *outf, int **wbx, int **vx, int j, int d, int mid, 
	   struct tracekn_s *curr_tr, struct traceknstack_s *dolist, int traceback)
{  
  int i;

  i = j - d;

  if (mid < d-4){
    if (traceback) fprintf(outf," (v4_5) multiloop, trace vx  %d %d %d\n", 
			   i+2+mid, mid, vx[i+2+mid][mid]);
    if (traceback) fprintf(outf," (v4_5) multiloop, trace wbx %d %d %d\n", 
			   j-1, d-mid-4, wbx[j-1][d-mid-4]);
    
    PushTraceknstack(dolist, AttachTracekn(curr_tr, i+2, i+2+mid, i+2+(int)((mid)/2),
					   i+3+(int)((mid)/2), dpcP, dpcS));
    PushTraceknstack(dolist, AttachTracekn(curr_tr, i+3+mid, j-1, i+3+mid+(int)((d-mid-4)/2),
					   i+4+mid+(int)((d-mid-4)/2), dpcB, dpcS));  
  }
}

/* (V4_6)__/ j-1 dangles off i,j */
int 
V4_6(int *s, int len, struct rnapar_2 *rnapar, int **icfg, int **wbx, int **vx, int j, int d, int mid)
{
  int V4_6;
  int    i;

  i = j - d;

  if (mid < d-5)
    V4_6 = 2*rnapar->P10 + rnapar->P5 + 2*rnapar->P6 
      + wsf*icfg[idxL(s[j-1])][idxP(s[j],s[i])]
      + wsf*icfg[idxPS(s[i],s[j])][idxPS(s[i+2],s[i+2+mid])] 
      + vx[i+2+mid][mid] + wbx[j-2][d-mid-5];
  else
    V4_6 = -BIGINT;

  return V4_6;
}
void
trace_V4_6(FILE *outf, int **wbx, int **vx, int j, int d, int mid, 
	   struct tracekn_s *curr_tr, struct traceknstack_s *dolist, int traceback)
{  
  int i;

  i = j - d;

  if (mid < d-5){
    if (traceback) fprintf(outf," (v4_6) multiloop, trace vx  %d %d %d\n", 
			   i+2+mid, mid, vx[i+2+mid][mid]);
    if (traceback) fprintf(outf," (v4_6) multiloop, trace wbx %d %d %d\n", 
			   j-2, d-mid-5, wbx[j-2][d-mid-5]);
    
    PushTraceknstack(dolist, AttachTracekn(curr_tr, i+2, i+2+mid, i+2+(int)((mid)/2),
					   i+3+(int)((mid)/2), dpcP, dpcS));
    PushTraceknstack(dolist, AttachTracekn(curr_tr, i+3+mid, j-2, i+3+mid+(int)((d-mid-5)/2),
					   i+4+mid+(int)((d-mid-5)/2), dpcB, dpcS));  
  }
}

/*  (V4_7)__/ i+mid+3 dangles off i+2,i+mid+2 */
int
V4_7(int *s, int len, struct rnapar_2 *rnapar, int **icfg, int **wbx, int **vx, int j, int d, int mid)
{
  int V4_7;
  int    i;

  i = j - d;

  if (mid < d-4)
    V4_7 = 2*rnapar->P10 + rnapar->P5 + 2*rnapar->P6
      + wsf*icfg[idxR(s[i+mid+3])][idxP(s[i+2],s[i+mid+2])]
      + wsf*icfg[idxPS(s[i],s[j])][idxPS(s[i+2],s[i+2+mid])] 
      + vx[i+2+mid][mid] + wbx[j-1][d-mid-5];
  else
    V4_7 = -BIGINT;
  
  return V4_7;
}
void
trace_V4_7(FILE *outf, int **wbx, int **vx, int j, int d, int mid, 
	   struct tracekn_s *curr_tr, struct traceknstack_s *dolist, int traceback)
{  
  int i;

  i = j - d;

  if (mid < d-5){
    if (traceback) fprintf(outf," (v4_7) multiloop, trace vx  %d %d %d\n", 
			   i+2+mid, mid, vx[i+2+mid][mid]);
    if (traceback) fprintf(outf," (v4_7) multiloop, trace wbx %d %d %d\n", 
			   j-1, d-mid-5, wbx[j-1][d-mid-5]);
    
    PushTraceknstack(dolist, AttachTracekn(curr_tr, i+2, i+2+mid, i+2+(int)((mid)/2),
					   i+3+(int)((mid)/2), dpcP, dpcS));
    PushTraceknstack(dolist, AttachTracekn(curr_tr, i+4+mid, j-1, i+4+mid+(int)((d-mid-5)/2),
					   i+5+mid+(int)((d-mid-5)/2), dpcB, dpcS));  
  }
}

/*  (V4_8)__/ j-1 dangles off i,j / i+mid+3 dangles off i+2,i+mid+2 */
int 
V4_8(int *s, int len, struct rnapar_2 *rnapar, int **icfg, int **wbx, int **vx, int j, int d, int mid)
{
  int V4_8;
  int    i;

  i = j - d;

  if (mid < d-6)
    V4_8 = 2*rnapar->P10 + rnapar->P5 + 3*rnapar->P6 
      + wsf*icfg[idxL(s[j-1])][idxP(s[j],s[i])]
      + wsf*icfg[idxR(s[i+mid+3])][idxP(s[i+2],s[i+mid+2])]
      + wsf*icfg[idxPS(s[i],s[j])][idxPS(s[i+2],s[i+2+mid])] 
      + vx[i+1+mid][mid] + wbx[j-2][d-mid-6];
  else
    V4_8 = -BIGINT;
  
  return V4_8;
}
void
trace_V4_8(FILE *outf, int **wbx, int **vx, int j, int d, int mid, 
	   struct tracekn_s *curr_tr, struct traceknstack_s *dolist, int traceback)
{  
  int i;

  i = j - d;

  if (mid < d-6){
    if (traceback) fprintf(outf," (v4_8) multiloop, trace vx  %d %d %d\n", 
			   i+2+mid, mid, vx[i+2+mid][mid]);
    if (traceback) fprintf(outf," (v4_8) multiloop, trace wbx %d %d %d\n", 
			   j-2, d-mid-6, wbx[j-2][d-mid-6]);
    
    PushTraceknstack(dolist, AttachTracekn(curr_tr, i+2, i+2+mid, i+2+(int)((mid)/2),
					   i+3+(int)((mid)/2), dpcP, dpcS));
    PushTraceknstack(dolist, AttachTracekn(curr_tr, i+4+mid, j-2, i+4+mid+(int)((d-mid-6)/2),
					   i+5+mid+(int)((d-mid-6)/2), dpcB, dpcS));  
  }
}

/*  FUNCTIONS: V5__/ 1 WBX 1 VX 
 */

/* (V5_1, V5_2, V5_3, V5_4)
 *   
 * contiguous coaxial pairs: (i,j) (i+2+mid,j-1)
 *                    ...add stack[j-1][i+2+mid][j][i] (+ 1 if no danglings)
 *   
 *  Return: int V5_X          
 */

/*  (V5_1) __/ no danglings / coaxial = stack + 1 */
int 
V5_1(int *s, int len, struct rnapar_2 *rnapar, int **icfg, int **wbx, int **vx, int j, int d, int mid)
{
  int V5_1;
  int    i;

  i = j - d;

  if (mid < d-3)
    V5_1 =  2*rnapar->P10 + rnapar->P5 
      + wsf*(icfg[idxPS(s[j-1],s[i+2+mid])][idxPS(s[j],s[i])] + INTSCALE)
      + wbx[i+1+mid][mid] + vx[j-1][d-mid-3];
  else
    V5_1 = -BIGINT;
  
  return V5_1;
}
void
trace_V5_1(FILE *outf, int **wbx, int **vx, int j, int d, int mid, 
	   struct tracekn_s *curr_tr, struct traceknstack_s *dolist, int traceback)
{  
  int i;

  i = j - d;

  if (mid < d-3){
    if (traceback) fprintf(outf," (v5_1) multiloop wbx %d %d %d\n", i+1+mid, mid, wbx[i+1+mid][mid]);
    if (traceback) fprintf(outf," (v5_1) multiloop vx  %d %d %d\n", j-1, d-mid-3, vx[j-1][d-mid-3]);
    
    PushTraceknstack(dolist, AttachTracekn(curr_tr, i+1, i+1+mid, i+1+(int)((mid)/2),
					   i+2+(int)((mid)/2), dpcB, dpcS));
    PushTraceknstack(dolist, AttachTracekn(curr_tr, i+2+mid, j-1, i+2+mid+(int)((d-mid-3)/2),
					   i+3+mid+(int)((d-mid-3)/2), dpcP, dpcS));
  }
}

/*  (V5_2)__/ i+1 dangles off i,j */
int 
V5_2(int *s, int len, struct rnapar_2 *rnapar, int **icfg, int **wbx, int **vx, int j, int d, int mid)
{
  int V5_2;
  int    i;

  i = j - d;
  
  if (mid > 1 && mid < d-3)
    V5_2 = 2*rnapar->P10 + rnapar->P5 + rnapar->P6 
      + wsf*icfg[idxR(s[i+1])][idxP(s[j],s[i])]
      + wsf*icfg[idxPS(s[j-1],s[i+2+mid])][idxPS(s[j],s[i])]
      + wbx[i+1+mid][mid-1] + vx[j-1][d-mid-3];
  else
    V5_2 = -BIGINT;
  
  return V5_2;
}
void
trace_V5_2(FILE *outf, int **wbx, int **vx, int j, int d, int mid, 
	   struct tracekn_s *curr_tr, struct traceknstack_s *dolist, int traceback)
{  
  int i;

  i = j - d;

  if (mid > 1 && mid < d-3){
    if (traceback) fprintf(outf," (v5_2) multiloop wbx %d %d %d\n", i+1+mid, mid-1, wbx[i+1+mid][mid-1]);
    if (traceback) fprintf(outf," (v5_2) multiloop vx  %d %d %d\n", j-1, d-mid-3, vx[j-1][d-mid-3]);
    
    PushTraceknstack(dolist, AttachTracekn(curr_tr, i+2, i+1+mid, i+2+(int)((mid-1)/2),
					   i+3+(int)((mid-1)/2), dpcB, dpcS));
    PushTraceknstack(dolist, AttachTracekn(curr_tr, i+2+mid, j-1, i+2+mid+(int)((d-mid-3)/2),
					   i+3+mid+(int)((d-mid-3)/2), dpcP, dpcS)); 
  }
}

/*  (V5_3)__/ i+mid+1 dangles off i+2+mid,j-1 */
int 
V5_3(int *s, int len, struct rnapar_2 *rnapar, int **icfg, int **wbx, int **vx, int j, int d, int mid)
{
  int V5_3;
  int    i;

  i = j - d;

  if (mid > 1 && mid < d-3)
    V5_3 = 2*rnapar->P10 + rnapar->P5 + rnapar->P6  
      + wsf*icfg[idxL(s[i+mid+1])][idxP(s[i+2+mid],s[j-1])]
      + wsf*icfg[idxPS(s[j-1],s[i+2+mid])][idxPS(s[j],s[i])] 
      + wbx[i+mid][mid-1] + vx[j-1][d-mid-3];
  else
    V5_3 = -BIGINT;

  return V5_3;
}
void
trace_V5_3(FILE *outf, int **wbx, int **vx, int j, int d, int mid, 
	   struct tracekn_s *curr_tr, struct traceknstack_s *dolist, int traceback)
{  
  int i;

  i = j - d;

  if (mid > 1 && mid < d-3){
    if (traceback) fprintf(outf," (v5_3) multiloop wbx %d %d %d\n", i+mid, mid-1, wbx[i+mid][mid-1]);
    if (traceback) fprintf(outf," (v5_3) multiloop vx  %d %d %d\n", j-1, d-mid-3, vx[j-1][d-mid-3]);
    
    PushTraceknstack(dolist, AttachTracekn(curr_tr, i+1, i+mid, i+1+(int)((mid-1)/2),
					   i+2+(int)((mid-1)/2), dpcB, dpcS));
    PushTraceknstack(dolist, AttachTracekn(curr_tr, i+2+mid, j-1, i+2+mid+(int)((d-mid-3)/2),
					   i+3+mid+(int)((d-mid-3)/2), dpcP, dpcS));  
  }
}

/*  (V5_4)__/  i+1 dangles off i,j /
            /  i+mid+1 dangles off i+2+mid,j-1 */
int 
V5_4(int *s, int len, struct rnapar_2 *rnapar, int **icfg, int **wbx, int **vx, int j, int d, int mid)
{
  int V5_4;
  int    i;

  i = j - d;

  if (mid > 2 && mid < d-3)
    V5_4 = 2*rnapar->P10 + rnapar->P5 + 2*rnapar->P6 
      + wsf*icfg[idxR(s[i+1])][idxP(s[j],s[i])]
      + wsf*icfg[idxL(s[i+mid+1])][idxP(s[i+2+mid],s[j-1])]
      + wsf*icfg[idxPS(s[j-1],s[i+2+mid])][idxPS(s[j],s[i])]
      + wbx[i+mid][mid-2] + vx[j-1][d-mid-3];
  else
    V5_4 = -BIGINT;
  
  return V5_4;
}
void
trace_V5_4(FILE *outf, int **wbx, int **vx, int j, int d, int mid, 
	   struct tracekn_s *curr_tr, struct traceknstack_s *dolist, int traceback)
{  
  int i;

  i = j - d;

  if (mid > 2 && mid < d-3){
    if (traceback) fprintf(outf," (v5_4) multiloop wbx %d %d %d\n", i+mid, mid-2, wbx[i+mid][mid-2]);
    if (traceback) fprintf(outf," (v5_4) multiloop vx  %d %d %d\n", j-1, d-mid-3, vx[j-1][d-mid-3]);
    
    PushTraceknstack(dolist, AttachTracekn(curr_tr, i+2, i+mid, i+2+(int)((mid-2)/2),
					   i+3+(int)((mid-2)/2), dpcB, dpcS));
    PushTraceknstack(dolist, AttachTracekn(curr_tr, i+2+mid, j-1, i+2+mid+(int)((d-mid-3)/2),
					   i+3+mid+(int)((d-mid-3)/2), dpcP, dpcS)); 
  }
}
 
/* (V5_5, V5_6, V5_7, V5_8)
 *   
 * non-contiguous coaxial pairs: (i,j) (i+2+mid,j-2)
 *                       ...add stack[j-2][i+2+mid][j][i]
 *   
 *  Return: int V5_X          
 */

/*  (V5_5) __/ no danglings */
int 
V5_5(int *s, int len, struct rnapar_2 *rnapar, int **icfg, int **wbx, int **vx, int j, int d, int mid)
{
  int V5_5;
  int    i;

  i = j - d;

  if (mid < d-4)
    V5_5 = 2*rnapar->P10 + rnapar->P5 + rnapar->P6
      + wsf*icfg[idxPS(s[j-2],s[i+2+mid])][idxPS(s[j],s[i])]
      + wbx[i+1+mid][mid] + vx[j-2][d-mid-4];
  else
    V5_5 = -BIGINT;
  
  return V5_5;
}
void
trace_V5_5(FILE *outf, int **wbx, int **vx, int j, int d, int mid, 
	   struct tracekn_s *curr_tr, struct traceknstack_s *dolist, int traceback)
{  
  int i;

  i = j - d;

  if (mid < d-4){
    if (traceback) fprintf(outf," (v5_5) multiloop wbx %d %d %d\n", i+1+mid, mid, wbx[i+1+mid][mid]);
    if (traceback) fprintf(outf," (v5_5) multiloop vx  %d %d %d\n", j-2, d-mid-4, vx[j-2][d-mid-4]);
    
    PushTraceknstack(dolist, AttachTracekn(curr_tr, i+1, i+1+mid, i+1+(int)((mid)/2),
					   i+2+(int)((mid)/2), dpcB, dpcS));
    PushTraceknstack(dolist, AttachTracekn(curr_tr, i+2+mid, j-2, i+2+mid+(int)((d-mid-4)/2),
					   i+3+mid+(int)((d-mid-4)/2), dpcP, dpcS));
  }
}

/*  (V5_6)__/ i+1 dangles off i,j */
int 
V5_6(int *s, int len, struct rnapar_2 *rnapar, int **icfg, int **wbx, int **vx, int j, int d, int mid)
{
  int V5_6;
  int    i;

  i = j - d;
  
  if (mid > 1 && mid < d-4)
    V5_6 = 2*rnapar->P10 + rnapar->P5 + 2*rnapar->P6 
      + wsf*icfg[idxR(s[i+1])][idxP(s[j],s[i])]
      + wsf*icfg[idxPS(s[j-2],s[i+2+mid])][idxPS(s[j],s[i])]
      + wbx[i+1+mid][mid-1] + vx[j-2][d-mid-4];
  else
    V5_6 = -BIGINT;
  
  return V5_6;
}
void
trace_V5_6(FILE *outf, int **wbx, int **vx, int j, int d, int mid, 
	   struct tracekn_s *curr_tr, struct traceknstack_s *dolist, int traceback)
{  
  int i;

  i = j - d;

  if (mid > 1 && mid < d-4){
    if (traceback) fprintf(outf," (v5_6)multiloop wbx %d %d %d\n", i+1+mid, mid-1, wbx[i+1+mid][mid-1]);
    if (traceback) fprintf(outf," (v5_6)multiloop vx  %d %d %d\n", j-2, d-mid-4, vx[j-2][d-mid-4]);
    
    PushTraceknstack(dolist, AttachTracekn(curr_tr, i+2, i+1+mid, i+2+(int)((mid-1)/2),
					   i+3+(int)((mid-1)/2), dpcB, dpcS));
    PushTraceknstack(dolist, AttachTracekn(curr_tr, i+2+mid, j-2, i+2+mid+(int)((d-mid-4)/2),
				   i+3+mid+(int)((d-mid-4)/2), dpcP, dpcS)); 
  }
}

/*  (V5_7)__/ i+mid+1 dangles off i+2+mid,j-2 */
int 
V5_7(int *s, int len, struct rnapar_2 *rnapar, int **icfg, int **wbx, int **vx, int j, int d, int mid)
{
  int V5_7;
  int    i;

  i = j - d;

  if (mid > 1 && mid < d-4)
    V5_7 = 2*rnapar->P10 + rnapar->P5 + 2*rnapar->P6  
      + wsf*icfg[idxL(s[i+mid+1])][idxP(s[i+2+mid],s[j-2])]
      + wsf*icfg[idxPS(s[j-2],s[i+2+mid])][idxPS(s[j],s[i])] 
      + wbx[i+mid][mid-1] + vx[j-2][d-mid-4];
  else
    V5_7 = -BIGINT;

  return V5_7;
}
void
trace_V5_7(FILE *outf, int **wbx, int **vx, int j, int d, int mid, 
	   struct tracekn_s *curr_tr, struct traceknstack_s *dolist, int traceback)
{  
  int i;

  i = j - d;

  if (mid > 1 && mid < d-4){
    if (traceback) fprintf(outf," (v5_7) multiloop wbx %d %d %d\n", i+mid, mid-1, wbx[i+mid][mid-1]);
    if (traceback) fprintf(outf," (v5_7) multiloop vx  %d %d %d\n", j-2, d-mid-4, vx[j-2][d-mid-4]);
    
    PushTraceknstack(dolist, AttachTracekn(curr_tr, i+1, i+mid, i+1+(int)((mid-1)/2),
					   i+2+(int)((mid-1)/2), dpcB, dpcS));
    PushTraceknstack(dolist, AttachTracekn(curr_tr, i+2+mid, j-2, i+2+mid+(int)((d-mid-4)/2),
					   i+3+mid+(int)((d-mid-4)/2), dpcP, dpcS));  
  }
}

/*  (V5_8)__/  i+1 dangles off i,j /
            /  i+mid+1 dangles off i+2+mid,j-2 */
int 
V5_8(int *s, int len, struct rnapar_2 *rnapar, int **icfg, int **wbx, int **vx, int j, int d, int mid)
{
  int V5_8;
  int    i;

  i = j - d;

  if (mid > 2 && mid < d-4)
    V5_8 = 2*rnapar->P10 + rnapar->P5 + 3*rnapar->P6 
      + wsf*icfg[idxR(s[i+1])][idxP(s[j],s[i])]
      + wsf*icfg[idxL(s[i+mid+1])][idxP(s[i+2+mid],s[j-2])]
      + wsf*icfg[idxPS(s[j-2],s[i+2+mid])][idxPS(s[j],s[i])]
      + wbx[i+mid][mid-2] + vx[j-2][d-mid-4];
  else
    V5_8 = -BIGINT;
  
  return V5_8;
}
void
trace_V5_8(FILE *outf, int **wbx, int **vx, int j, int d, int mid, 
	   struct tracekn_s *curr_tr, struct traceknstack_s *dolist, int traceback)
{  
  int i;

  i = j - d;

  if (mid > 2 && mid < d-4){
    if (traceback) fprintf(outf," (v5_8) multiloop wbx %d %d %d\n", i+mid, mid-2, wbx[i+mid][mid-2]);
    if (traceback) fprintf(outf," (v5_8) multiloop vx  %d %d %d\n", j-2, d-mid-4, vx[j-2][d-mid-4]);
    
    PushTraceknstack(dolist, AttachTracekn(curr_tr, i+2, i+mid, i+2+(int)((mid-2)/2),
					   i+3+(int)((mid-2)/2), dpcB, dpcS));
    PushTraceknstack(dolist, AttachTracekn(curr_tr, i+2+mid, j-2, i+2+mid+(int)((d-mid-4)/2),
					   i+3+mid+(int)((d-mid-4)/2), dpcP, dpcS)); 
  }
}

/*  FUNCTIONS: V6__/ 2 VX
 * 
 *  one connects i+mid1, i+mid1+mid   (i+mid1+mid, mid)
 *  one connects i+mid1+mid+1, j-mid2 (j-mid2, d-mid-mid1-mid2-1)

 */

/* V6_1-V6_4 (16 diagrams)
 *
 * contiguous coaxial  pairs (i+mid1, i+mid1+mid) (i+mid1+mid+1,j-mid2)
 * ...add stack[i+mid1+mid][i+mid1][i+mid1+mid+1][j-mid2]
 *   
 *  Return: int V6_X          
 */

/*  (V6_1) __/ no danglings / coaxial = stack + 1 */
int 
V6_1 (int *s, int len, struct rnapar_2 *rnapar, int **icfg, int **vx, int j, int d, int mid, int mid1, int mid2)
{
  int V6_1;
  int    i;

  i = j - d;

  if (mid < d-3 && mid1 == 1 && mid2 == 1)
    V6_1 = 3*rnapar->P10 + rnapar->P5
      + wsf*(icfg[idxPS(s[i+mid+mid1],s[i+mid1])]
	     [idxPS(s[i+mid+mid1+1],s[j-mid2])] + INTSCALE)
      + vx[i+mid+mid1][mid] + vx[j-mid2][d-mid-mid1-mid2-1];
  else
    V6_1 = -BIGINT;
  
  return V6_1;
}

/*  (V6_2a) __/  mid1 = 2 / i+1 dangles off i,j / 
 */
int 
V6_2a(int *s, int len, struct rnapar_2 *rnapar, int **icfg, int **vx, int j, int d, int mid, int mid1, int mid2)
{
  int V6_2a;
  int     i;

  i = j - d;

  if (mid < d-3 && mid1 == 2 && mid2 == 1)
    V6_2a = 3*rnapar->P10 + rnapar->P5 + rnapar->P6*(mid1-1) 
      + wsf*icfg[idxR(s[i+1])][idxP(s[j],s[i])]
      + wsf*icfg[idxPS(s[i+mid+mid1],s[i+mid1])]
                [idxPS(s[i+mid+mid1+1],s[j-mid2])]
      + vx[i+mid+mid1][mid] + vx[j-mid2][d-mid-mid1-mid2-1];
  else
    V6_2a = -BIGINT;

  return V6_2a;
}

/*  (V6_2b) __/  mid1 = 2 / i+1 dangles off i+mid1 i+mid1+mid */
int 
V6_2b(int *s, int len, struct rnapar_2 *rnapar, int **icfg, int **vx, int j, int d, int mid, int mid1, int mid2)
{
  int V6_2b;
  int     i;

  i = j - d;

  if (mid < d-3 && mid1 == 2 && mid2 == 1)
    V6_2b = 3*rnapar->P10 + rnapar->P5 + rnapar->P6*(mid1-1) 
      + wsf*icfg[idxL(s[i+1])][idxP(s[i+mid1],s[i+mid1+mid])]
      + wsf*icfg[idxPS(s[i+mid+mid1],s[i+mid1])]
                [idxPS(s[i+mid+mid1+1],s[j-mid2])]
      + vx[i+mid+mid1][mid] + vx[j-mid2][d-mid-mid1-mid2-1];
  else
    V6_2b = -BIGINT;
  
  return V6_2b;
}

/*  (V6_2c) __/   mid1 > 2 / i+1 dangles off i,j / 
 *	      / i+mid1-1 dangles off i+mid1 i+mid1+mid /
 */
int 
V6_2c(int *s, int len, struct rnapar_2 *rnapar, int **icfg, int **vx, int j, int d, int mid, int mid1, int mid2)
{
  int V6_2c;
  int     i;

  i = j - d;

  if (mid < d-3 && mid1 > 2 && mid2 == 1)
    V6_2c = 3*rnapar->P10 + rnapar->P5 + rnapar->P6*(mid1-1) 
      + wsf*icfg[idxL(s[i+1])][idxP(s[i+mid1],s[i+mid1+mid])]
      + wsf*icfg[idxL(s[i+mid1-1])][idxP(s[i+mid1],s[i+mid1+mid])]
      + wsf*icfg[idxPS(s[i+mid+mid1],s[i+mid1])]
                [idxPS(s[i+mid+mid1+1],s[j-mid2])]
      + vx[i+mid+mid1][mid] + vx[j-mid2][d-mid-mid1-mid2-1];
  else
    V6_2c = -BIGINT;

  return V6_2c;
}

/*  (V6_3a) __/ mid2 = 2 / j-1 dangles off i,j / 
 */
int 
V6_3a(int *s, int len, struct rnapar_2 *rnapar, int **icfg, int **vx, int j, int d, int mid, int mid1, int mid2)
{
  int V6_3a;
  int     i;

  i = j - d;

  if (mid < d-3 && mid1 == 1 && mid2 == 2)
    V6_3a = 3*rnapar->P10 + rnapar->P5 + rnapar->P6*(mid2-1) 
      + wsf*icfg[idxL(s[j-1])][idxP(s[j],s[i])]
      + wsf*icfg[idxPS(s[i+mid+mid1],s[i+mid1])]
                [idxPS(s[i+mid+mid1+1],s[j-mid2])]
      + vx[i+mid+mid1][mid] + vx[j-mid2][d-mid-mid1-mid2-1];
  else
    V6_3a = -BIGINT;
  
  return V6_3a;
}

/*  (V6_3b) __/  mid2 = 2 / j-1 dangles off i+mid+mid1+1,j-mid2 /
 */
int 
V6_3b(int *s, int len, struct rnapar_2 *rnapar, int **icfg, int **vx, int j, int d, int mid, int mid1, int mid2)
{
  int V6_3b;
  int     i;

  i = j - d;

  if (mid < d-3 && mid1 == 1 && mid2 == 2)
    V6_3b = 3*rnapar->P10 + rnapar->P5 + rnapar->P6*(mid2-1) 
      + wsf*icfg[idxR(s[j-1])][idxP(s[i+mid+mid1+1],s[j-mid2])]
      + wsf*icfg[idxPS(s[i+mid+mid1],s[i+mid1])]
                [idxPS(s[i+mid+mid1+1],s[j-mid2])]
      + vx[i+mid+mid1][mid] + vx[j-mid2][d-mid-mid1-mid2-1];
  else
    V6_3b = -BIGINT;
  
  return V6_3b;
}

/*  (V6_3c) __/ mid2 > 2 / j-1 dangles off i,j /
 *            / j-mid2+1 dangles off i+mid+mid1+1,j-mid2 /
 */
int 
V6_3c(int *s, int len, struct rnapar_2 *rnapar, int **icfg, int **vx, int j, int d, int mid, int mid1, int mid2)
{
  int V6_3c;
  int     i;

  i = j - d;

 if (mid < d-3 && mid1 == 1 && mid2 > 2)
   V6_3c = 3*rnapar->P10 + rnapar->P5 + rnapar->P6*(mid2-1) 
     + wsf*icfg[idxL(s[j-1])][idxP(s[j],s[i])]
     + wsf*icfg[idxR(s[j-mid2+1])][idxP(s[i+mid+mid1+1],s[j-mid2])]
     + wsf*icfg[idxPS(s[i+mid+mid1],s[i+mid1])]
               [idxPS(s[i+mid+mid1+1],s[j-mid2])]
     + vx[i+mid+mid1][mid] + vx[j-mid2][d-mid-mid1-mid2-1];
 else
    V6_3c = -BIGINT;
  
  return V6_3c;
}

/*  (V6_4a) __/ mid1 = 2 / i+1 dangles off i,j /
 *            / mid2 = 2 / j-1 dangles off i,j /
 */
int 
V6_4a(int *s, int len, struct rnapar_2 *rnapar, int **icfg, int **vx, int j, int d, int mid, int mid1, int mid2)
{
  int V6_4a;
  int     i;

  i = j - d;

  if (mid < d-3 && mid1 == 2 && mid2 == 2)
    V6_4a = 3*rnapar->P10 + rnapar->P5 + rnapar->P6*(mid1+mid2-2) 
      + wsf*icfg[idxR(s[i+1])][idxP(s[j],s[i])] 
      + wsf*icfg[idxL(s[j-1])][idxP(s[j],s[i])]
      + wsf*icfg[idxPS(s[i+mid+mid1],s[i+mid1])]
                [idxPS(s[i+mid+mid1+1],s[j-mid2])]
      + vx[i+mid+mid1][mid] + vx[j-mid2][d-mid-mid1-mid2-1];
 else
    V6_4a = -BIGINT;
  
  return V6_4a;
}

/*  (V6_4b) __/ mid1 = 2 / i+1 dangles off i,j /
 *            / mid2 = 2 / j-1 dangles off i+mid+mid1+1,j-mid2 /
 */
int 
V6_4b(int *s, int len, struct rnapar_2 *rnapar, int **icfg, int **vx, int j, int d, int mid, int mid1, int mid2)
{
  int V6_4b;
  int    i;

  i = j - d;

 if (mid < d-3 && mid1 == 2 && mid2 == 2)
   V6_4b = 3*rnapar->P10 + rnapar->P5 + rnapar->P6*(mid1+mid2-2) 
     + wsf*icfg[idxR(s[i+1])][idxP(s[j],s[i])]
     + wsf*icfg[idxR(s[j-1])][idxP(s[i+mid+mid1+1],s[j-mid2])]
     + wsf*icfg[idxPS(s[i+mid+mid1],s[i+mid1])]
               [idxPS(s[i+mid+mid1+1],s[j-mid2])]
     + vx[i+mid+mid1][mid] + vx[j-mid2][d-mid-mid1-mid2-1];
 else
    V6_4b = -BIGINT;

  return V6_4b;
}

/*  (V6_4c) __/ mid1 = 2 / i+1 dangles off i+mid1 i+mid1+mid /
 *	      / mid2 = 2 / j-1 dangles off i,j /
 */
int 
V6_4c(int *s, int len, struct rnapar_2 *rnapar, int **icfg, int **vx, int j, int d, int mid, int mid1, int mid2)
{
  int V6_4c;
  int     i;

  i = j - d;

  if (mid < d-3 && mid1 == 2 && mid2 == 2)
    V6_4c = 3*rnapar->P10 + rnapar->P5 + rnapar->P6*(mid1+mid2-2) 
      + wsf*icfg[idxL(s[j-1])][idxP(s[j],s[i])]
      + wsf*icfg[idxL(s[i+1])][idxP(s[i+mid1],s[i+mid1+mid])]
      + wsf*icfg[idxPS(s[i+mid+mid1],s[i+mid1])]
                [idxPS(s[i+mid+mid1+1],s[j-mid2])]
      + vx[i+mid+mid1][mid] + vx[j-mid2][d-mid-mid1-mid2-1];
 else
    V6_4c = -BIGINT;
  
  return V6_4c;
}

/*  (V6_4d) __/ mid1 = 2 / i+1 dangles off i+mid1 i+mid1+mid / 
 *            / mid2 = 2 / j-1 dangles off i+mid+mid1+1,j-mid2 /
 */
int 
V6_4d(int *s, int len, struct rnapar_2 *rnapar, int **icfg, int **vx, int j, int d, int mid, int mid1, int mid2)
{
  int V6_4d;
  int     i;

  i = j - d;

  if (mid < d-3 && mid1 == 2 && mid2 == 2)
    V6_4d = 3*rnapar->P10 + rnapar->P5 + rnapar->P6*(mid1+mid2-2) 
      + wsf*icfg[idxL(s[i+1])][idxP(s[i+mid1],s[i+mid1+mid])]
      + wsf*icfg[idxR(s[j-1])][idxP(s[i+mid+mid1+1],s[j-mid2])]
      + wsf*icfg[idxPS(s[i+mid+mid1],s[i+mid1])]
                [idxPS(s[i+mid+mid1+1],s[j-mid2])]
      + vx[i+mid+mid1][mid] + vx[j-mid2][d-mid-mid1-mid2-1];
 else
    V6_4d = -BIGINT;
  
  return V6_4d;
}

/*  (V6_4e) __/ mid1 > 2 / i+1 dangles off i,j /
 *            / i+mid1-1 dangles off i+mid1 i+mid1+mid / 
 */
int 
V6_4e(int *s, int len, struct rnapar_2 *rnapar, int **icfg, int **vx, int j, int d, int mid, int mid1, int mid2)
{
  int V6_4e;
  int     i;

  i = j - d;

  if (mid < d-3 && mid1 > 2 && mid2 == 2)
    V6_4e = 3*rnapar->P10 + rnapar->P5 + rnapar->P6*(mid1+mid2-2) 
      + wsf*icfg[idxR(s[i+1])][idxP(s[j],s[i])] 
      + wsf*icfg[idxL(s[j-1])][idxP(s[j],s[i])]
      + wsf*icfg[idxL(s[i+mid1-1])][idxP(s[i+mid1],s[i+mid1+mid])]
      + wsf*icfg[idxPS(s[i+mid+mid1],s[i+mid1])]
                [idxPS(s[i+mid+mid1+1],s[j-mid2])]
      + vx[i+mid+mid1][mid] + vx[j-mid2][d-mid-mid1-mid2-1];
 else
    V6_4e = -BIGINT;
  
  return V6_4e;
}

/*  (V6_4f) __/ mid1 > 2 / i+1 dangles off i,j /
 *            / i+mid1-1 dangles off i+mid1 i+mid1+mid / 
 *	      / mid2 = 2 / j-1 dangles off i+mid+mid1+1,j-mid2 / 
 */
int 
V6_4f(int *s, int len, struct rnapar_2 *rnapar, int **icfg, int **vx, int j, int d, int mid, int mid1, int mid2)
{
  int V6_4f;
  int     i;

  i = j - d;

  if (mid < d-3 && mid1 > 2 && mid2 == 2)
    V6_4f = 3*rnapar->P10 + rnapar->P5 + rnapar->P6*(mid1+mid2-2) 
      + wsf*icfg[idxR(s[i+1])][idxP(s[j],s[i])]
      + wsf*icfg[idxL(s[i+mid1-1])][idxP(s[i+mid1],s[i+mid1+mid])]
      + wsf*icfg[idxR(s[j-1])][idxP(s[i+mid+mid1+1],s[j-mid2])]
      + wsf*icfg[idxPS(s[i+mid+mid1],s[i+mid1])]
                [idxPS(s[i+mid+mid1+1],s[j-mid2])]
      + vx[i+mid+mid1][mid] + vx[j-mid2][d-mid-mid1-mid2-1];
 else
    V6_4f = -BIGINT;
  
  return V6_4f;
}

/*  (V6_4g) __/ mid1 = 2 / i+1 dangles off i,j /
 *            / mid2 > 2 / j-1 dangles off i,j /
 *            / j-mid2+1 dangles off i+mid+mid1+1,j-mid2 /
 */
int 
V6_4g(int *s, int len, struct rnapar_2 *rnapar, int **icfg, int **vx, int j, int d, int mid, int mid1, int mid2)
{
  int V6_4g;
  int     i;

  i = j - d;

  if (mid < d-3 && mid1 == 2 && mid2 > 2)
    V6_4g = 3*rnapar->P10 + rnapar->P5 + rnapar->P6*(mid1+mid2-2) 
      + wsf*icfg[idxR(s[i+1])][idxP(s[j],s[i])]
      + wsf*icfg[idxL(s[j-1])][idxP(s[j],s[i])]
      + wsf*icfg[idxR(s[j-mid2+1])][idxP(s[i+mid+mid1+1],s[j-mid2])]
      + wsf*icfg[idxPS(s[i+mid+mid1],s[i+mid1])]
                [idxPS(s[i+mid+mid1+1],s[j-mid2])]
      + vx[i+mid+mid1][mid] + vx[j-mid2][d-mid-mid1-mid2-1];
 else
    V6_4g = -BIGINT;
  
  return V6_4g;
}

/*  (V6_4h) __/ mid1 = 2 / i+1 dangles off i+mid1 i+mid1+mid /
 *            / mid2 > 2 / j-1 dangles off i,j /
 *            / j-mid2+1 dangles off i+mid+mid1+1,j-mid2 /
 */
int 
V6_4h(int *s, int len, struct rnapar_2 *rnapar, int **icfg, int **vx, int j, int d, int mid, int mid1, int mid2)
{
  int V6_4h;
  int     i;

  i = j - d;

  if (mid < d-3 && mid1 == 2 && mid2 > 2)
    V6_4h = 3*rnapar->P10 + rnapar->P5 + rnapar->P6*(mid1+mid2-2) 
      + wsf*icfg[idxL(s[j-1])][idxP(s[j],s[i])]
      + wsf*icfg[idxL(s[i+1])][idxP(s[i+mid1],s[i+mid1+mid])]
      + wsf*icfg[idxR(s[j-mid2+1])][idxP(s[i+mid+mid1+1],s[j-mid2])]
      + wsf*icfg[idxPS(s[i+mid+mid1],s[i+mid1])]
                [idxPS(s[i+mid+mid1+1],s[j-mid2])]
      + vx[i+mid+mid1][mid] + vx[j-mid2][d-mid-mid1-mid2-1];
 else
    V6_4h = -BIGINT;
  
  return V6_4h;
}

/*  (V6_4i) __/ mid1 > 2 / i+1 dangles off i,j /
 *            / i+mid1-1 dangles off i+mid1 i+mid1+mid / 
 *	      / mid2 > 2 / j-1 dangles off i,j /
 *            / j-mid2+1 dangles off i+mid+mid1+1,j-mid2 / 
 */
int 
V6_4i(int *s, int len, struct rnapar_2 *rnapar, int **icfg, int **vx, int j, int d, int mid, int mid1, int mid2)
{
  int V6_4i;
  int     i;

  i = j - d;

  if (mid < d-3 && mid1 > 2 && mid2 > 2)
    V6_4i = 3*rnapar->P10 + rnapar->P5 + rnapar->P6*(mid1+mid2-2) 
      + wsf*icfg[idxR(s[i+1])][idxP(s[j],s[i])]
      + wsf*icfg[idxL(s[j-1])][idxP(s[j],s[i])]
      + wsf*icfg[idxL(s[i+mid1-1])][idxP(s[i+mid1],s[i+mid1+mid])]
      + wsf*icfg[idxR(s[j-mid2+1])][idxP(s[i+mid+mid1+1],s[j-mid2])]
      + wsf*icfg[idxPS(s[i+mid+mid1],s[i+mid1])]
                [idxPS(s[i+mid+mid1+1],s[j-mid2])]
      + vx[i+mid+mid1][mid] + vx[j-mid2][d-mid-mid1-mid2-1];
 else
    V6_4i = -BIGINT;
  
  return V6_4i;
}
void
trace_V6_CC(FILE *outf, int **vx, int j, int d, int mid, int mid1, int mid2, 
	 struct tracekn_s *curr_tr, struct traceknstack_s *dolist, int traceback)
{  
  int i;

  i = j - d;

  if (mid < d-3){
    if (traceback) fprintf(outf," (v6_CC) multiloop vx %d %d %d\n", 
			   i+mid+mid1, mid, vx[i+mid+mid1][mid]);
    if (traceback) fprintf(outf," (v6_NCC) multiloop vx %d %d %d\n", j-mid2, d-mid-mid1-mid2-1, 
			   vx[j-mid2][d-mid-mid1-mid2-1]);
    
    PushTraceknstack(dolist, AttachTracekn(curr_tr, i+mid1, i+mid1+mid,
					   i+mid1+(int)((mid)/2),
					   i+mid1+(int)((mid)/2)+1, dpcP, dpcS));
    PushTraceknstack(dolist, AttachTracekn(curr_tr, i+mid1+mid+1, j-mid2,
					   i+mid1+mid+1+(int)((d-mid-mid1-mid2-1)/2),
					   i+mid1+mid+2+(int)((d-mid-mid1-mid2-1)/2), 
					   dpcP, dpcS));
  }
}

/* V6_5-V6_8 (16 diagrams)
 *
 * non-contiguous coaxial  pairs (i+mid1, i+mid1+mid) (i+mid1+mid+1,j-mid2)
 * ...add stack[i+mid1+mid][i+mid1][i+mid1+mid+1][j-mid2]
 *   
 *  Return: int V6_X          
 */

/*  (V6_5) __/ no danglings  */
int 
V6_5 (int *s, int len, struct rnapar_2 *rnapar, int **icfg, int **vx, int j, int d, int mid, int mid1, int mid2)
{
  int V6_5;
  int    i;

  i = j - d;

  if (mid < d-4 && mid1 == 1 && mid2 == 1)
    V6_5 = 3*rnapar->P10 + rnapar->P5 + rnapar->P6
      + wsf*icfg[idxPS(s[i+mid+mid1],s[i+mid1])]
	        [idxPS(s[i+mid+mid1+2],s[j-mid2])] 
      + vx[i+mid+mid1][mid] + vx[j-mid2][d-mid-mid1-mid2-2];
  else
    V6_5 = -BIGINT;
  
  return V6_5;
}

/*  (V6_6a) __/  mid1 = 2 / i+1 dangles off i,j / 
 */
int 
V6_6a(int *s, int len, struct rnapar_2 *rnapar, int **icfg, int **vx, int j, int d, int mid, int mid1, int mid2)
{
  int V6_6a;
  int     i;

  i = j - d;

  if (mid < d-4 && mid1 == 2 && mid2 == 1)
    V6_6a = 3*rnapar->P10 + rnapar->P5 + rnapar->P6*mid1 
      + wsf*icfg[idxR(s[i+1])][idxP(s[j],s[i])]
      + wsf*icfg[idxPS(s[i+mid+mid1],s[i+mid1])]
                [idxPS(s[i+mid+mid1+2],s[j-mid2])]
      + vx[i+mid+mid1][mid] + vx[j-mid2][d-mid-mid1-mid2-2];
  else
    V6_6a = -BIGINT;

  return V6_6a;
}

/*  (V6_6b) __/  mid1 = 2 / i+1 dangles off i+mid1 i+mid1+mid */
int 
V6_6b(int *s, int len, struct rnapar_2 *rnapar, int **icfg, int **vx, int j, int d, int mid, int mid1, int mid2)
{
  int V6_6b;
  int     i;

  i = j - d;

  if (mid < d-4 && mid1 == 2 && mid2 == 1)
    V6_6b = 3*rnapar->P10 + rnapar->P5 + rnapar->P6*mid1 
      + wsf*icfg[idxL(s[i+1])][idxP(s[i+mid1],s[i+mid1+mid])]
      + wsf*icfg[idxPS(s[i+mid+mid1],s[i+mid1])]
                [idxPS(s[i+mid+mid1+2],s[j-mid2])]
      + vx[i+mid+mid1][mid] + vx[j-mid2][d-mid-mid1-mid2-2];
  else
    V6_6b = -BIGINT;
  
  return V6_6b;
}

/*  (V6_6c) __/   mid1 > 2 / i+1 dangles off i,j / 
 *	      / i+mid1-1 dangles off i+mid1 i+mid1+mid /
 */
int 
V6_6c(int *s, int len, struct rnapar_2 *rnapar, int **icfg, int **vx, int j, int d, int mid, int mid1, int mid2)
{
  int V6_6c;
  int     i;

  i = j - d;

  if (mid < d-4 && mid1 > 2 && mid2 == 1)
    V6_6c = 3*rnapar->P10 + rnapar->P5 + rnapar->P6*mid1 
      + wsf*icfg[idxL(s[i+1])][idxP(s[i+mid1],s[i+mid1+mid])]
      + wsf*icfg[idxL(s[i+mid1-1])][idxP(s[i+mid1],s[i+mid1+mid])]
      + wsf*icfg[idxPS(s[i+mid+mid1],s[i+mid1])]
                [idxPS(s[i+mid+mid1+2],s[j-mid2])]
      + vx[i+mid+mid1][mid] + vx[j-mid2][d-mid-mid1-mid2-2];
  else
    V6_6c = -BIGINT;

  return V6_6c;
}

/*  (V6_7a) __/ mid2 = 2 / j-1 dangles off i,j / 
 */
int 
V6_7a(int *s, int len, struct rnapar_2 *rnapar, int **icfg, int **vx, int j, int d, int mid, int mid1, int mid2)
{
  int V6_7a;
  int     i;

  i = j - d;

  if (mid < d-4 && mid1 == 1 && mid2 == 2)
    V6_7a = 3*rnapar->P10 + rnapar->P5 + rnapar->P6*mid2
      + wsf*icfg[idxL(s[j-1])][idxP(s[j],s[i])]
      + wsf*icfg[idxPS(s[i+mid+mid1],s[i+mid1])]
                [idxPS(s[i+mid+mid1+2],s[j-mid2])]
      + vx[i+mid+mid1][mid] + vx[j-mid2][d-mid-mid1-mid2-2];
  else
    V6_7a = -BIGINT;
  
  return V6_7a;
}

/*  (V6_7b) __/  mid2 = 2 / j-1 dangles off i+mid+mid1+1,j-mid2 /
 */
int 
V6_7b(int *s, int len, struct rnapar_2 *rnapar, int **icfg, int **vx, int j, int d, int mid, int mid1, int mid2)
{
  int V6_7b;
  int     i;

  i = j - d;

  if (mid < d-4 && mid1 == 1 && mid2 == 2)
    V6_7b = 3*rnapar->P10 + rnapar->P5 + rnapar->P6*mid2 
      + wsf*icfg[idxR(s[j-1])][idxP(s[i+mid+mid1+1],s[j-mid2])]
      + wsf*icfg[idxPS(s[i+mid+mid1],s[i+mid1])]
                [idxPS(s[i+mid+mid1+2],s[j-mid2])]
      + vx[i+mid+mid1][mid] + vx[j-mid2][d-mid-mid1-mid2-2];
  else
    V6_7b = -BIGINT;
  
  return V6_7b;
}

/*  (V6_7c) __/ mid2 > 2 / j-1 dangles off i,j /
 *            / j-mid2+1 dangles off i+mid+mid1+1,j-mid2 /
 */
int 
V6_7c(int *s, int len, struct rnapar_2 *rnapar, int **icfg, int **vx, int j, int d, int mid, int mid1, int mid2)
{
  int V6_7c;
  int     i;

  i = j - d;

 if (mid < d-4 && mid1 == 1 && mid2 > 2)
   V6_7c = 3*rnapar->P10 + rnapar->P5 + rnapar->P6*mid2 
     + wsf*icfg[idxL(s[j-1])][idxP(s[j],s[i])]
     + wsf*icfg[idxR(s[j-mid2+1])][idxP(s[i+mid+mid1+1],s[j-mid2])]
     + wsf*icfg[idxPS(s[i+mid+mid1],s[i+mid1])]
               [idxPS(s[i+mid+mid1+2],s[j-mid2])]
     + vx[i+mid+mid1][mid] + vx[j-mid2][d-mid-mid1-mid2-2];
 else
    V6_7c = -BIGINT;
  
  return V6_7c;
}

/*  (V6_8a) __/ mid1 = 2 / i+1 dangles off i,j /
 *            / mid2 = 2 / j-1 dangles off i,j /
 */
int 
V6_8a(int *s, int len, struct rnapar_2 *rnapar, int **icfg, int **vx, int j, int d, int mid, int mid1, int mid2)
{
  int V6_8a;
  int     i;

  i = j - d;

  if (mid < d-4 && mid1 == 2 && mid2 == 2)
    V6_8a = 3*rnapar->P10 + rnapar->P5 + rnapar->P6*(mid1+mid2-1) 
      + wsf*icfg[idxR(s[i+1])][idxP(s[j],s[i])] 
      + wsf*icfg[idxL(s[j-1])][idxP(s[j],s[i])]
      + wsf*icfg[idxPS(s[i+mid+mid1],s[i+mid1])]
                [idxPS(s[i+mid+mid1+2],s[j-mid2])]
      + vx[i+mid+mid1][mid] + vx[j-mid2][d-mid-mid1-mid2-2];
 else
    V6_8a = -BIGINT;
  
  return V6_8a;
}

/*  (V6_8b) __/ mid1 = 2 / i+1 dangles off i,j /
 *            / mid2 = 2 / j-1 dangles off i+mid+mid1+1,j-mid2 /
 */
int 
V6_8b(int *s, int len, struct rnapar_2 *rnapar, int **icfg, int **vx, int j, int d, int mid, int mid1, int mid2)
{
  int V6_8b;
  int    i;

  i = j - d;

 if (mid < d-4 && mid1 == 2 && mid2 == 2)
   V6_8b = 3*rnapar->P10 + rnapar->P5 + rnapar->P6*(mid1+mid2-1) 
     + wsf*icfg[idxR(s[i+1])][idxP(s[j],s[i])]
     + wsf*icfg[idxR(s[j-1])][idxP(s[i+mid+mid1+1],s[j-mid2])]
     + wsf*icfg[idxPS(s[i+mid+mid1],s[i+mid1])]
               [idxPS(s[i+mid+mid1+2],s[j-mid2])]
     + vx[i+mid+mid1][mid] + vx[j-mid2][d-mid-mid1-mid2-2];
 else
    V6_8b = -BIGINT;

  return V6_8b;
}

/*  (V6_8c) __/ mid1 = 2 / i+1 dangles off i+mid1 i+mid1+mid /
 *	      / mid2 = 2 / j-1 dangles off i,j /
 */
int 
V6_8c(int *s, int len, struct rnapar_2 *rnapar, int **icfg, int **vx, int j, int d, int mid, int mid1, int mid2)
{
  int V6_8c;
  int     i;

  i = j - d;

  if (mid < d-4 && mid1 == 2 && mid2 == 2)
    V6_8c = 3*rnapar->P10 + rnapar->P5 + rnapar->P6*(mid1+mid2-1) 
      + wsf*icfg[idxL(s[j-1])][idxP(s[j],s[i])]
      + wsf*icfg[idxL(s[i+1])][idxP(s[i+mid1],s[i+mid1+mid])]
      + wsf*icfg[idxPS(s[i+mid+mid1],s[i+mid1])]
                [idxPS(s[i+mid+mid1+2],s[j-mid2])]
      + vx[i+mid+mid1][mid] + vx[j-mid2][d-mid-mid1-mid2-2];
 else
    V6_8c = -BIGINT;
  
  return V6_8c;
}

/*  (V6_8d) __/ mid1 = 2 / i+1 dangles off i+mid1 i+mid1+mid / 
 *            / mid2 = 2 / j-1 dangles off i+mid+mid1+1,j-mid2 /
 */
int 
V6_8d(int *s, int len, struct rnapar_2 *rnapar, int **icfg, int **vx, int j, int d, int mid, int mid1, int mid2)
{
  int V6_8d;
  int     i;

  i = j - d;

  if (mid < d-4 && mid1 == 2 && mid2 == 2)
    V6_8d = 3*rnapar->P10 + rnapar->P5 + rnapar->P6*(mid1+mid2-1) 
      + wsf*icfg[idxL(s[i+1])][idxP(s[i+mid1],s[i+mid1+mid])]
      + wsf*icfg[idxR(s[j-1])][idxP(s[i+mid+mid1+1],s[j-mid2])]
      + wsf*icfg[idxPS(s[i+mid+mid1],s[i+mid1])]
                [idxPS(s[i+mid+mid1+2],s[j-mid2])]
      + vx[i+mid+mid1][mid] + vx[j-mid2][d-mid-mid1-mid2-2];
 else
    V6_8d = -BIGINT;
  
  return V6_8d;
}

/*  (V6_8e) __/ mid1 > 2 / i+1 dangles off i,j /
 *            / i+mid1-1 dangles off i+mid1 i+mid1+mid / 
 */
int 
V6_8e(int *s, int len, struct rnapar_2 *rnapar, int **icfg, int **vx, int j, int d, int mid, int mid1, int mid2)
{
  int V6_8e;
  int     i;

  i = j - d;

  if (mid < d-4 && mid1 > 2 && mid2 == 2)
    V6_8e = 3*rnapar->P10 + rnapar->P5 + rnapar->P6*(mid1+mid2-1) 
      + wsf*icfg[idxR(s[i+1])][idxP(s[j],s[i])] 
      + wsf*icfg[idxL(s[j-1])][idxP(s[j],s[i])]
      + wsf*icfg[idxL(s[i+mid1-1])][idxP(s[i+mid1],s[i+mid1+mid])]
      + wsf*icfg[idxPS(s[i+mid+mid1],s[i+mid1])]
                [idxPS(s[i+mid+mid1+2],s[j-mid2])]
      + vx[i+mid+mid1][mid] + vx[j-mid2][d-mid-mid1-mid2-2];
 else
    V6_8e = -BIGINT;
  
  return V6_8e;
}

/*  (V6_8f) __/ mid1 > 2 / i+1 dangles off i,j /
 *            / i+mid1-1 dangles off i+mid1 i+mid1+mid / 
 *	      / mid2 = 2 / j-1 dangles off i+mid+mid1+1,j-mid2 / 
 */
int 
V6_8f(int *s, int len, struct rnapar_2 *rnapar, int **icfg, int **vx, int j, int d, int mid, int mid1, int mid2)
{
  int V6_8f;
  int     i;

  i = j - d;

  if (mid < d-4 && mid1 > 2 && mid2 == 2)
    V6_8f = 3*rnapar->P10 + rnapar->P5 + rnapar->P6*(mid1+mid2-1) 
      + wsf*icfg[idxR(s[i+1])][idxP(s[j],s[i])]
      + wsf*icfg[idxL(s[i+mid1-1])][idxP(s[i+mid1],s[i+mid1+mid])]
      + wsf*icfg[idxR(s[j-1])][idxP(s[i+mid+mid1+1],s[j-mid2])]
      + wsf*icfg[idxPS(s[i+mid+mid1],s[i+mid1])]
                [idxPS(s[i+mid+mid1+2],s[j-mid2])]
      + vx[i+mid+mid1][mid] + vx[j-mid2][d-mid-mid1-mid2-2];
 else
    V6_8f = -BIGINT;
  
  return V6_8f;
}

/*  (V6_8g) __/ mid1 = 2 / i+1 dangles off i,j /
 *            / mid2 > 2 / j-1 dangles off i,j /
 *            / j-mid2+1 dangles off i+mid+mid1+1,j-mid2 /
 */
int 
V6_8g(int *s, int len, struct rnapar_2 *rnapar, int **icfg, int **vx, int j, int d, int mid, int mid1, int mid2)
{
  int V6_8g;
  int     i;

  i = j - d;

  if (mid < d-4 && mid1 == 2 && mid2 > 2)
    V6_8g = 3*rnapar->P10 + rnapar->P5 + rnapar->P6*(mid1+mid2-1) 
      + wsf*icfg[idxR(s[i+1])][idxP(s[j],s[i])]
      + wsf*icfg[idxL(s[j-1])][idxP(s[j],s[i])]
      + wsf*icfg[idxR(s[j-mid2+1])][idxP(s[i+mid+mid1+1],s[j-mid2])]
      + wsf*icfg[idxPS(s[i+mid+mid1],s[i+mid1])]
                [idxPS(s[i+mid+mid1+2],s[j-mid2])]
      + vx[i+mid+mid1][mid] + vx[j-mid2][d-mid-mid1-mid2-2];
 else
    V6_8g = -BIGINT;
  
  return V6_8g;
}

/*  (V6_8h) __/ mid1 = 2 / i+1 dangles off i+mid1 i+mid1+mid /
 *            / mid2 > 2 / j-1 dangles off i,j /
 *            / j-mid2+1 dangles off i+mid+mid1+1,j-mid2 /
 */
int 
V6_8h(int *s, int len, struct rnapar_2 *rnapar, int **icfg, int **vx, int j, int d, int mid, int mid1, int mid2)
{
  int V6_8h;
  int     i;

  i = j - d;

  if (mid < d-4 && mid1 == 2 && mid2 > 2)
    V6_8h = 3*rnapar->P10 + rnapar->P5 + rnapar->P6*(mid1+mid2-1) 
      + wsf*icfg[idxL(s[j-1])][idxP(s[j],s[i])]
      + wsf*icfg[idxL(s[i+1])][idxP(s[i+mid1],s[i+mid1+mid])]
      + wsf*icfg[idxR(s[j-mid2+1])][idxP(s[i+mid+mid1+1],s[j-mid2])]
      + wsf*icfg[idxPS(s[i+mid+mid1],s[i+mid1])]
                [idxPS(s[i+mid+mid1+2],s[j-mid2])]
      + vx[i+mid+mid1][mid] + vx[j-mid2][d-mid-mid1-mid2-2];
 else
    V6_8h = -BIGINT;
  
  return V6_8h;
}

/*  (V6_8i) __/ mid1 > 2 / i+1 dangles off i,j /
 *            / i+mid1-1 dangles off i+mid1 i+mid1+mid / 
 *	      / mid2 > 2 / j-1 dangles off i,j /
 *            / j-mid2+1 dangles off i+mid+mid1+1,j-mid2 / 
 */
int 
V6_8i(int *s, int len, struct rnapar_2 *rnapar, int **icfg, int **vx, int j, int d, int mid, int mid1, int mid2)
{
  int V6_8i;
  int     i;

  i = j - d;

  if (mid < d-4 && mid1 > 2 && mid2 > 2)
    V6_8i = 3*rnapar->P10 + rnapar->P5 + rnapar->P6*(mid1+mid2-1) 
      + wsf*icfg[idxR(s[i+1])][idxP(s[j],s[i])]
      + wsf*icfg[idxL(s[j-1])][idxP(s[j],s[i])]
      + wsf*icfg[idxL(s[i+mid1-1])][idxP(s[i+mid1],s[i+mid1+mid])]
      + wsf*icfg[idxR(s[j-mid2+1])][idxP(s[i+mid+mid1+1],s[j-mid2])]
      + wsf*icfg[idxPS(s[i+mid+mid1],s[i+mid1])]
                [idxPS(s[i+mid+mid1+2],s[j-mid2])]
      + vx[i+mid+mid1][mid] + vx[j-mid2][d-mid-mid1-mid2-2];
 else
    V6_8i = -BIGINT;
  
  return V6_8i;
}
void
trace_V6_NCC(FILE *outf, int **vx, int j, int d, int mid, int mid1, int mid2, 
	 struct tracekn_s *curr_tr, struct traceknstack_s *dolist, int traceback)
{  
  int i;

  i = j - d;

  if (mid < d-4){
    if (traceback) fprintf(outf," (v6_NCC) multiloop vx %d %d %d\n", 
			   i+mid+mid1, mid, vx[i+mid+mid1][mid]);
    if (traceback) fprintf(outf," (v6_NCC) multiloop vx %d %d %d\n", j-mid2, d-mid-mid1-mid2-2, 
			   vx[j-mid2][d-mid-mid1-mid2-2]);
    
    PushTraceknstack(dolist, AttachTracekn(curr_tr, i+mid1, i+mid1+mid,
					   i+mid1+(int)((mid)/2),
					   i+mid1+(int)((mid)/2)+1, dpcP, dpcS));
    PushTraceknstack(dolist, AttachTracekn(curr_tr, i+mid1+mid+2, j-mid2,
					   i+mid1+mid+2+(int)((d-mid-mid1-mid2-2)/2),
					   i+mid1+mid+3+(int)((d-mid-mid1-mid2-2)/2), 
					   dpcP, dpcS));
  }
}

/*  FUNCTIONS: V7 (V7_1, V7_2, V7_3, V7_4)
 *   
 *  Purpose: calculate diagrams with 2 WHX 
 *   
 *  Return: int V7_X          
 */

/*  (V7_1) __/ no danglings  */
int 
V7_1 (int *s, int len, struct rnapar_2 *rnapar, int **icfg, int ****whx, int j, int d, int mid, int mid1, int mid2)
{
  int V7_1;
  int    i;

  i = j - d;

  V7_1 = rnapar->P10P + rnapar->P5P + rnapar->P13 + 3*rnapar->P6P 
    + whx[i+1+mid][mid][mid1][mid2] 
    + whx[j-1][d-mid1-4][mid-mid1-mid2-3][d-mid-5];
  
  return V7_1;
}
void
trace_V7_1(FILE *outf, int ****whx, int j, int d, int mid, int mid1, int mid2, 
	   struct tracekn_s *curr_tr, struct traceknstack_s *dolist, int traceback)
{  
  int i;

  i = j - d;

  if (traceback) fprintf(outf," (V7.1) multiloop whx %d %d %d %d %d\n", 
			 i+1+mid, mid, mid1, mid2,  
			 whx[i+1+mid][mid][mid1][mid2]);
  if (traceback) fprintf(outf," (V7.1) multiloop whx %d %d %d %d %d\n", 
			 j-1, d-mid1-4, mid-mid1-mid2-3, d-mid-5,
			 whx[j-1][d-mid1-4][mid-mid1-mid2-3][d-mid-5]);
  
  PushTraceknstack(dolist, AttachTracekn(curr_tr, i+1, i+1+mid,
					 i+1+mid1, i+1+mid-mid2, dpcS, dpcS));
  PushTraceknstack(dolist, AttachTracekn(curr_tr, i+3+mid1, j-1,
					 i+mid-mid2, i+mid+4, dpcS, dpcS));
}

/*  (V7_2) __/ i+1 dangles off i,j  */
int 
V7_2 (int *s, int len, struct rnapar_2 *rnapar, int **icfg, int ****whx, int j, int d, int mid, int mid1, int mid2)
{
  int V7_2;
  int    i;

  i = j - d;

  if (mid1 > 0)
    V7_2 = rnapar->P10P + rnapar->P5P + rnapar->P13 + 4*rnapar->P6P 
      + wkn*icfg[idxR(s[i+1])][idxP(s[j],s[i])]
      + whx[i+1+mid][mid-1][mid1-1][mid2] 
      + whx[j-1][d-mid1-4][mid-mid1-mid2-3][d-mid-5];
  else
    V7_2 = -BIGINT;
  
  return V7_2;
}
void
trace_V7_2(FILE *outf, int ****whx, int j, int d, int mid, int mid1, int mid2, 
	   struct tracekn_s *curr_tr, struct traceknstack_s *dolist, int traceback)
{
  int i;

  i = j - d;

  if (mid1 > 0){
    if (traceback) fprintf(outf," (V7.2) multiloop whx %d %d %d %d %d\n", 
			   i+1+mid, mid-1, mid1-1, mid2,  
			   whx[i+1+mid][mid-1][mid1-1][mid2]);
    if (traceback) fprintf(outf," (V7.2) multiloop whx %d %d %d %d %d\n", 
			   j-1, d-mid1-4, mid-mid1-mid2-3, d-mid-5,
			   whx[j-1][d-mid1-4][mid-mid1-mid2-3][d-mid-5]);
    
    PushTraceknstack(dolist, AttachTracekn(curr_tr, i+2, i+1+mid,
					   i+1+mid1, i+1+mid-mid2, dpcS, dpcS));
    PushTraceknstack(dolist, AttachTracekn(curr_tr, i+3+mid1, j-1,
					   i+mid-mid2, i+mid+4, dpcS, dpcS));
  }
}

/*  (V7_3) __/ j-1 dangles off i,j  */
int 
V7_3 (int *s, int len, struct rnapar_2 *rnapar, int **icfg, int ****whx, int j, int d, int mid, int mid1, int mid2)
{
  int V7_3;
  int    i;

  i = j - d;

  if (d-mid > 5)
    V7_3 = rnapar->P10P + rnapar->P5P + rnapar->P13 + 4*rnapar->P6P 
      + wkn*icfg[idxL(s[j-1])][idxP(s[j],s[i])]
      + whx[i+1+mid][mid][mid1][mid2] 
      + whx[j-2][d-mid1-5][mid-mid1-mid2-3][d-mid-6];
  else
    V7_3 = -BIGINT;
  
  return V7_3;
}
void
trace_V7_3(FILE *outf, int ****whx, int j, int d, int mid, int mid1, int mid2, 
	   struct tracekn_s *curr_tr, struct traceknstack_s *dolist, int traceback)
{
  int i;

  i = j - d;

  if (d-mid > 5){
    if (traceback) fprintf(outf," (V7.3) multiloop whx %d %d %d %d %d\n", 
			   i+1+mid, mid, mid1, mid2,  
			   whx[i+1+mid][mid][mid1][mid2]);
    if (traceback) fprintf(outf," (V7.3) multiloop whx %d %d %d %d %d\n", 
			   j-2, d-mid1-4, mid-mid1-mid2-2, d-mid-4,
			   whx[j-2][d-mid1-5][mid-mid1-mid2-3][d-mid-6]);
    
    PushTraceknstack(dolist, AttachTracekn(curr_tr, i+1, i+1+mid,
					   i+1+mid1, i+1+mid-mid2, dpcS, dpcS));
    PushTraceknstack(dolist, AttachTracekn(curr_tr, i+3+mid1, j-2,
					   i+mid-mid2, i+mid+4, dpcS, dpcS));
  }
}

/*  (V7_4) __/ i+1 dangles off i,j / j-1 dangles off i,j  */
int 
V7_4 (int *s, int len, struct rnapar_2 *rnapar, int **icfg, int ****whx, int j, int d, int mid, int mid1, int mid2)
{
  int V7_4;
  int    i;

  i = j - d;

  if (mid1 > 0 && d-mid > 5)
    V7_4 = rnapar->P10P + rnapar->P5P + rnapar->P13 + 5*rnapar->P6P 
      + wkn*icfg[idxR(s[i+1])][idxP(s[j],s[i])]
      + wkn*icfg[idxL(s[j-1])][idxP(s[j],s[i])]
      + whx[i+1+mid][mid-1][mid1-1][mid2] 
      + whx[j-2][d-mid1-5][mid-mid1-mid2-3][d-mid-6];
  else
    V7_4 = -BIGINT;
  
  return V7_4;
}
void
trace_V7_4(FILE *outf, int ****whx, int j, int d, int mid, int mid1, int mid2, 
	   struct tracekn_s *curr_tr, struct traceknstack_s *dolist, int traceback)
{
  int i;

  i = j - d;

  if (mid1 > 0 && d-mid > 5){
    if (traceback) fprintf(outf," (V7.4) multiloop whx %d %d %d %d %d\n", 
			   i+1+mid, mid-1, mid1-1, mid2,  
			   whx[i+1+mid][mid-1][mid1-1][mid2]);
    if (traceback) fprintf(outf," (V7.4) multiloop whx %d %d %d %d %d\n", 
			   j-2, d-mid1-5, mid-mid1-mid2-3, d-mid-6,
			   whx[j-2][d-mid1-5][mid-mid1-mid2-3][d-mid-6]);
    
    PushTraceknstack(dolist, AttachTracekn(curr_tr, i+2, i+1+mid,
					   i+1+mid1, i+1+mid-mid2, dpcS, dpcS));
    PushTraceknstack(dolist, AttachTracekn(curr_tr, i+3+mid1, j-2,
					   i+mid-mid2, i+mid+4, dpcS, dpcS));
  }
}

/*  FUNCTIONS: V8__/ 1 ZHX 1 WHX 
 */

/* (V8_1, V8_2)
 *   
 * contiguous coaxial pairs (i,j) (i+1,i+1+mid), 
 * ....add stack[i][j][i+1][i+1+mid1]
 * check for j-1 dangling off the pair (i,j)
 *   
 *  Return: int V8_n          
 */

/*  (V8_1) __/ i+mid+2 dangles off i+1,i+mid+1   */
int 
V8_1 (int *s, int len, struct rnapar_2 *rnapar, int **icfg, int ****whx, int ****zhx,
      int j, int d, int mid, int mid1, int mid2)
{
  int V8_1;
  int    i;

  i = j - d;

  V8_1 = 2*rnapar->P10P + rnapar->P5P + rnapar->P13 + 3*rnapar->P6P 
    + wkn*icfg[idxPS(s[i],s[j])][idxPS(s[i+1],s[i+mid+1])] 
    + zhx[i+1+mid][mid][mid1][mid2] 
    + whx[j-1][d-mid1-4][mid-mid1-mid2-3][d-mid-5];
  
  return V8_1;
}
void
trace_V8_1(FILE *outf, int ****whx, int ****zhx, int j, int d, int mid, int mid1, 
	   int mid2, struct tracekn_s *curr_tr, struct traceknstack_s *dolist, int traceback)
{
  int i;

  i = j - d;

  if (traceback) fprintf(outf," (V8.1) multiloop zhx %d %d %d %d %d\n", 
			 i+1+mid, mid, mid1, mid2,  
			 zhx[i+1+mid][mid][mid1][mid2]);
  if (traceback) fprintf(outf," (V8.1) multiloop whx %d %d %d %d %d\n", 
			 j-1, d-mid1-4, mid-mid1-mid2-3, d-mid-5,
			 whx[j-1][d-mid1-4][mid-mid1-mid2-3][d-mid-5]);
  
  PushTraceknstack(dolist, AttachTracekn(curr_tr, i+1, i+1+mid,
					 i+1+mid1, i+1+mid-mid2, dpcP, dpcS));
  PushTraceknstack(dolist, AttachTracekn(curr_tr, i+3+mid1, j-1,
					 i+mid-mid2, i+mid+4, dpcS, dpcS));
}
  
/*  (V8_2) __/ i+mid+2 dangles off i+1,i+mid+1 /  
             / j-1 dangles off i,j             */
int 
V8_2 (int *s, int len, struct rnapar_2 *rnapar, int **icfg, int ****whx, int ****zhx,
      int j, int d, int mid, int mid1, int mid2)
{
  int V8_2;
  int    i;

  i = j - d;

  if (d-mid > 5)
    V8_2 = 2*rnapar->P10P + rnapar->P5P + rnapar->P13 + 4*rnapar->P6P 
      + wkn*icfg[idxL(s[j-1])][idxP(s[j],s[i])]
      + wkn*icfg[idxPS(s[i],s[j])][idxPS(s[i+1],s[i+mid+1])]
      + zhx[i+1+mid][mid][mid1][mid2] 
      + whx[j-2][d-mid1-5][mid-mid1-mid2-3][d-mid-6];
  else
    V8_2 = -BIGINT;
  
  return V8_2;
}
void
trace_V8_2(FILE *outf, int ****whx, int ****zhx, int j, int d, int mid, int mid1, 
	   int mid2, struct tracekn_s *curr_tr, struct traceknstack_s *dolist, int traceback)
{
  int i;

  i = j - d;

  if (d-mid > 5){
    if (traceback) fprintf(outf," (V8.2) multiloop zhx %d %d %d %d %d\n", 
			   i+1+mid, mid, mid1, mid2,  
			   zhx[i+1+mid][mid][mid1][mid2]);
    if (traceback) fprintf(outf," (V8.2) multiloop whx %d %d %d %d %d\n", 
			   j-2, d-mid1-5, mid-mid1-mid2-3, d-mid-6,
			   whx[j-2][d-mid1-5][mid-mid1-mid2-3][d-mid-6]);
    
    PushTraceknstack(dolist, AttachTracekn(curr_tr, i+1, i+1+mid,
					   i+1+mid1, i+1+mid-mid2, dpcP, dpcS));
    PushTraceknstack(dolist, AttachTracekn(curr_tr, i+3+mid1, j-2,
					   i+mid-mid2, i+mid+4, dpcS, dpcS));
  }
}

/* (V8_3, V8_4)
 *   
 * non-contiguous coaxial pairs (i,j) (i+2,i+2+mid), 
 * ....add stack[i][j][i+2][i+2+mid1]
 * check for j-1 dangling off the pair (i,j)
 *   
 *  Return: int V8_n          
 */

/*  (V8_3) __/ i+mid+3 dangles off i+2,i+mid+2   */
int 
V8_3(int *s, int len, struct rnapar_2 *rnapar, int **icfg, int ****whx, int ****zhx,
     int j, int d, int mid, int mid1, int mid2)
{
  int V8_3;
  int    i;

  i = j - d;

  if (d-mid > 5)
    V8_3 = 2*rnapar->P10P + rnapar->P5P + rnapar->P13 + 4*rnapar->P6P 
      + wkn*icfg[idxPS(s[i],s[j])][idxPS(s[i+2],s[i+mid+2])] 
      + zhx[i+2+mid][mid][mid1][mid2] 
      + whx[j-1][d-mid1-5][mid-mid1-mid2-3][d-mid-6];
  else
    V8_3 = -BIGINT;
  
  return V8_3;
}
void
trace_V8_3(FILE *outf, int ****whx, int ****zhx, int j, int d, int mid, int mid1, 
	   int mid2, struct tracekn_s *curr_tr, struct traceknstack_s *dolist, int traceback)
{
  int i;

  i = j - d;

  if (d-mid > 5){
    if (traceback) fprintf(outf," (V8.3) multiloop zhx %d %d %d %d %d\n", 
			   i+2+mid, mid, mid1, mid2,  
			   zhx[i+2+mid][mid][mid1][mid2]);
    if (traceback) fprintf(outf," (V8.3) multiloop whx %d %d %d %d %d\n", 
			   j-1, d-mid1-5, mid-mid1-mid2-3, d-mid-6,
			   whx[j-1][d-mid1-5][mid-mid1-mid2-3][d-mid-6]);
    
    PushTraceknstack(dolist, AttachTracekn(curr_tr, i+2, i+2+mid,
					   i+2+mid1, i+2+mid-mid2, dpcP, dpcS));
    PushTraceknstack(dolist, AttachTracekn(curr_tr, i+4+mid1, j-1,
					   i+mid-mid2+1, i+mid+5, dpcS, dpcS));
  }
}

/*  (V8_4) __/ i+mid+3 dangles off i+2,i+mid+2 /  
             / j-1 dangles off i,j             */
int 
V8_4(int *s, int len, struct rnapar_2 *rnapar, int **icfg, int ****whx, int ****zhx,
     int j, int d, int mid, int mid1, int mid2)
{
  int V8_4;
  int    i;

  i = j - d;

  if (d-mid > 6)
    V8_4 = 2*rnapar->P10P + rnapar->P5P + rnapar->P13 + 5*rnapar->P6P 
      + wkn*icfg[idxL(s[j-1])][idxP(s[j],s[i])]
      + wkn*icfg[idxPS(s[i],s[j])][idxPS(s[i+2],s[i+mid+2])]
      + zhx[i+2+mid][mid][mid1][mid2] 
      + whx[j-2][d-mid1-6][mid-mid1-mid2-3][d-mid-7];
  else
    V8_4 = -BIGINT;
  
  return V8_4;
}
void
trace_V8_4(FILE *outf, int ****whx, int ****zhx, int j, int d, int mid, int mid1, 
	   int mid2, struct tracekn_s *curr_tr, struct traceknstack_s *dolist, int traceback)
{
  int i;

  i = j - d;

  if (d-mid > 6){
    if (traceback) fprintf(outf," (V8.4) multiloop zhx %d %d %d %d %d\n", 
			   i+2+mid, mid, mid1, mid2,  
			   zhx[i+2+mid][mid][mid1][mid2]);
    if (traceback) fprintf(outf," (V8.4) multiloop whx %d %d %d %d %d\n", 
			   j-2, d-mid1-6, mid-mid1-mid2-3, d-mid-7,
			   whx[j-2][d-mid1-6][mid-mid1-mid2-3][d-mid-7]);
    
    PushTraceknstack(dolist, AttachTracekn(curr_tr, i+2, i+2+mid,
					   i+2+mid1, i+2+mid-mid2, dpcP, dpcS));
    PushTraceknstack(dolist, AttachTracekn(curr_tr, i+4+mid1, j-2,
					   i+mid-mid2+1, i+mid+5, dpcS, dpcS));
  }
}


/*  FUNCTIONS: V9__/ 1 WHX 1 ZHX
 */

/* (V9_1, V9_2)
 *   
 * contiguous coaxial pairs: (i,j) (i+mid1+2, j-1)), 
 * ....add: stack[j-1][i+2+mid1][j][i]
 * check for i+1 dangling off the pair (i,j)
 *   
 *  Return: int V9_X          
 */

/*  (V9_1) __/ i+mid1+1 dangling off i+mid1+2,j-1  */
int 
V9_1(int *s, int len, struct rnapar_2 *rnapar, int **icfg, int ****whx, int ****zhx,
     int j, int d, int mid, int mid1, int mid2)
{
  int V9_1;
  int    i;

  i = j - d;

  V9_1 = 2*rnapar->P10P + rnapar->P5P + rnapar->P13 + 3*rnapar->P6P
    + wkn*icfg[idxPS(s[j-1],s[i+3+mid1])][idxPS(s[j],s[i])]
    + whx[i+1+mid][mid][mid1][mid2] 
    + zhx[j-1][d-mid1-4][mid-mid1-mid2-3][d-mid-5];
  
  return V9_1;
}
void 
trace_V9_1(FILE *outf, int ****whx, int ****zhx, 
	   int j, int d, int mid, int mid1, int mid2, 
	   struct tracekn_s *curr_tr, struct traceknstack_s *dolist, int traceback)
{
  int i;

  i = j - d;

  if (traceback) fprintf(outf," (V9.1) multiloop whx %d %d %d %d %d\n", 
			 i+1+mid, mid, mid1, mid2,  
			 whx[i+1+mid][mid][mid1][mid2]);
  if (traceback) fprintf(outf," (V9.1) multiloop zhx %d %d %d %d %d\n", 
			 j-1, d-mid1-4, mid-mid1-mid2-3, d-mid-5,
			 zhx[j-1][d-mid1-4][mid-mid1-mid2-3][d-mid-5]);
  
  PushTraceknstack(dolist, AttachTracekn(curr_tr, i+1, i+1+mid,
					 i+1+mid1, i+1+mid-mid2, dpcS, dpcS));
  PushTraceknstack(dolist, AttachTracekn(curr_tr, i+3+mid1, j-1,
					 i+mid-mid2, i+mid+4, dpcP, dpcS));
}

/*  (V9_2) __/ i+mid1+1 dangling off i+mid1+2,j-1 /
             / i+1 dangles off i,j                */
int 
V9_2(int *s, int len, struct rnapar_2 *rnapar, int **icfg, int ****whx, int ****zhx,
     int j, int d, int mid, int mid1, int mid2)
{
  int V9_2;
  int    i;

  i = j - d;

  if (mid1 > 0)
    V9_2 = 2*rnapar->P10P + rnapar->P5P + rnapar->P13 + 4*rnapar->P6P 
      + wkn*icfg[idxR(s[i+1])][idxP(s[j],s[i])]
      + wkn*icfg[idxPS(s[j-1],s[i+3+mid1])][idxPS(s[j],s[i])]
      + whx[i+1+mid][mid-1][mid1-1][mid2] 
      + zhx[j-1][d-mid1-4][mid-mid1-mid2-3][d-mid-5];
  else
    V9_2 = -BIGINT;
  
  return V9_2;
}
void 
trace_V9_2(FILE *outf, int ****whx, int ****zhx, 
	   int j, int d, int mid, int mid1, int mid2, 
	   struct tracekn_s *curr_tr, struct traceknstack_s *dolist, int traceback)
{
  int i;

  i = j - d;

  if (mid1 > 0){
    if (traceback) fprintf(outf," (V9.2) multiloop whx %d %d %d %d %d\n", 
			   i+1+mid, mid-1, mid1-1, mid2,  
			   whx[i+1+mid][mid-1][mid1-1][mid2]);
    if (traceback) fprintf(outf," (V9.2) multiloop zhx %d %d %d %d %d\n", 
			   j-1, d-mid1-4, mid-mid1-mid2-3, d-mid-5,
			   zhx[j-1][d-mid1-4][mid-mid1-mid2-3][d-mid-5]);
    
    PushTraceknstack(dolist, AttachTracekn(curr_tr, i+2, i+1+mid,
					   i+1+mid1, i+1+mid-mid2, dpcS, dpcS));
    PushTraceknstack(dolist, AttachTracekn(curr_tr, i+3+mid1, j-1,
					   i+mid-mid2, i+mid+4, dpcP, dpcS));
  }
}

/* (V9_3, V9_4)
 *   
 *  non-contiguous coaxial pairs: (i,j) (i+mid1+2, j-2)), 
 *   ....add: stack[j-2][i+2+mid1][j][i]
 *  check for i+1 dangling off the pair (i,j)
 *   
 *  Return: int V9_X          
 */

/*  (V9_3) __/ i+mid1+1 dangling off i+mid1+2,j-2  */
int 
V9_3(int *s, int len, struct rnapar_2 *rnapar, int **icfg, int ****whx, int ****zhx,
     int j, int d, int mid, int mid1, int mid2)
{
  int V9_3;
  int    i;

  i = j - d;

  if (d-mid > 5)
    V9_3 = 2*rnapar->P10P + rnapar->P5P + rnapar->P13 + 4*rnapar->P6P
      + wkn*icfg[idxPS(s[j-2],s[i+3+mid1])][idxPS(s[j],s[i])]
      + whx[i+1+mid][mid][mid1][mid2] 
      + zhx[j-2][d-mid1-5][mid-mid1-mid2-3][d-mid-6];
  else
    V9_3 = -BIGINT;
  
  return V9_3;
}
void 
trace_V9_3(FILE *outf, int ****whx, int ****zhx, 
	   int j, int d, int mid, int mid1, int mid2, 
	   struct tracekn_s *curr_tr, struct traceknstack_s *dolist, int traceback)
{
  int i;

  i = j - d;

  if (d-mid > 5){
    if (traceback) fprintf(outf," (V9.3) multiloop whx %d %d %d %d %d\n", 
			   i+1+mid, mid, mid1, mid2,  
			   whx[i+1+mid][mid][mid1][mid2]);
    if (traceback) fprintf(outf," (V9.3) multiloop zhx %d %d %d %d %d\n", 
			   j-2, d-mid1-5, mid-mid1-mid2-3, d-mid-6,
			   zhx[j-2][d-mid1-5][mid-mid1-mid2-3][d-mid-6]);
    
    PushTraceknstack(dolist, AttachTracekn(curr_tr, i+1, i+1+mid,
					   i+1+mid1, i+1+mid-mid2, dpcS, dpcS));
    PushTraceknstack(dolist, AttachTracekn(curr_tr, i+3+mid1, j-2,
					   i+mid-mid2, i+mid+4, dpcP, dpcS));
  }
}

/*  (V9_4) __/ i+mid1+1 dangling off i+mid1+2,j-2 /
             / i+1 dangles off i,j                */
int 
V9_4(int *s, int len, struct rnapar_2 *rnapar, int **icfg, int ****whx, int ****zhx,
     int j, int d, int mid, int mid1, int mid2)
{
  int V9_4;
  int    i;

  i = j - d;

  if (mid1 > 0 && d-mid > 5)
    V9_4 = 2*rnapar->P10P + rnapar->P5P + rnapar->P13 + 5*rnapar->P6P 
      + wkn*icfg[idxR(s[i+1])][idxP(s[j],s[i])]
      + wkn*icfg[idxPS(s[j-2],s[i+3+mid1])][idxPS(s[j],s[i])]
      + whx[i+1+mid][mid-1][mid1-1][mid2] 
      + zhx[j-2][d-mid1-5][mid-mid1-mid2-3][d-mid-6];
  else
    V9_4 = -BIGINT;
  
  return V9_4;
}
void 
trace_V9_4(FILE *outf, int ****whx, int ****zhx, 
	   int j, int d, int mid, int mid1, int mid2, 
	   struct tracekn_s *curr_tr, struct traceknstack_s *dolist, int traceback)
{
  int i;

  i = j - d;

  if (mid1 > 0 && d-mid > 5){
    if (traceback) fprintf(outf," (V9.4) multiloop whx %d %d %d %d %d\n", 
			   i+1+mid, mid-1, mid1-1, mid2,  
			   whx[i+1+mid][mid-1][mid1-1][mid2]);
    if (traceback) fprintf(outf," (V9.4) multiloop zhx %d %d %d %d %d\n", 
			   j-2, d-mid1-5, mid-mid1-mid2-3, d-mid-6,
			   zhx[j-2][d-mid1-5][mid-mid1-mid2-3][d-mid-6]);
    
    PushTraceknstack(dolist, AttachTracekn(curr_tr, i+2, i+1+mid,
					   i+1+mid1, i+1+mid-mid2, dpcS, dpcS));
    PushTraceknstack(dolist, AttachTracekn(curr_tr, i+3+mid1, j-2,
					   i+mid-mid2, i+mid+4, dpcP, dpcS));
  }
}

/*  FUNCTIONS: V10__/ 2 YHX 
 */

/*(V10_1-V10_4)
 *   
 * contiguous coaxial pairs (i+mid1+1, i+1+mid-mid2), (i+mid-mid2, i+mid+4)
 * ....add: stack[i+mid-mid2][i+2+mid][i+mid-mid2+1][i+mid1+1]
 * check for i+1 or j-1 dangling off the pair (i,j)
 *   
 *  Return: int V10_X          
 */

/*  (V10_1) __/ i+mid1+2 dangles off i+1+mid-mid2,i+1+mid1 /
              / i+mid+3 dangles off i+4+mid,i+mid-mid2     */
int 
V10_1(int *s, int len, struct rnapar_2 *rnapar, int **icfg, int ****yhx, int j, int d, int mid, int mid1, int mid2)
{
  int V10_1;
  int     i;

  i = j - d;

  V10_1 = 3*rnapar->P10P + rnapar->P5P + rnapar->P13 + 3*rnapar->P6P
    + wkn*icfg[idxR(s[i+mid1+2])][idxP(s[i+1+mid-mid2],s[i+1+mid1])]
    + wkn*icfg[idxL(s[i+mid+3])][idxP(s[i+mid+4],s[i+mid-mid2])]
    + wkn*icfg[idxPS(s[i+mid-mid2],s[i+mid+4])]
              [idxPS(s[i+1+mid-mid2],s[i+1+mid1])] 
    + yhx[i+1+mid][mid][mid1][mid2]
    + yhx[j-1][d-mid1-4][mid-mid1-mid2-3][d-mid-5];
  
  return V10_1;
}
void 
trace_V10_1(FILE *outf, int ****yhx, int j, int d, int mid, int mid1, int mid2, 
	    struct tracekn_s *curr_tr, struct traceknstack_s *dolist, int traceback)
{
  int i;

  i = j - d;

  if (traceback) fprintf(outf," (V10.1) multiloop yhx %d %d %d %d %d\n", 
			 i+1+mid, mid, mid1, mid2,  
			 yhx[i+1+mid][mid][mid1][mid2]);
  if (traceback) fprintf(outf," (V10.1) multiloop zhx %d %d %d %d %d\n", 
			 j-1, d-mid1-4, mid-mid1-mid2-3, d-mid-5,
			 yhx[j-1][d-mid1-4][mid-mid1-mid2-3][d-mid-5]);
  
  PushTraceknstack(dolist, AttachTracekn(curr_tr, i+1, i+1+mid,
					 i+1+mid1, i+1+mid-mid2, dpcS, dpcP));
  PushTraceknstack(dolist, AttachTracekn(curr_tr, i+3+mid1, j-1,
					 i+mid-mid2, i+mid+4, dpcS, dpcP));
}

/*  (V10_2) __/ i+mid1+2 dangles off i+1+mid-mid2,i+1+mid1  /
              / i+mid+3 dangles off i+4+mid,i+mid-mid2      /
              / i+1 dangles off i,j                         */
int 
V10_2(int *s, int len, struct rnapar_2 *rnapar, int **icfg, int ****yhx, int j, int d, int mid, int mid1, int mid2)
{
  int V10_2;
  int     i;

  i = j - d;

  if (mid > 2 && mid1 > 0)
    V10_2 = 3*rnapar->P10P + rnapar->P5P + rnapar->P13 + 4*rnapar->P6P 
      + wkn*icfg[idxR(s[i+1])][idxP(s[j],s[i])]
      + wkn*icfg[idxR(s[i+mid1+2])][idxP(s[i+1+mid-mid2],s[i+1+mid1])]
      + wkn*icfg[idxL(s[i+mid+3])][idxP(s[i+mid+4],s[i+mid-mid2])]
      + wkn*icfg[idxPS(s[i+mid-mid2],s[i+mid+4])]
	        [idxPS(s[i+1+mid-mid2],s[i+1+mid1])]
      + yhx[i+1+mid][mid-1][mid1-1][mid2] 
      + yhx[j-1][d-mid1-4][mid-mid1-mid2-3][d-mid-5];
  else
    V10_2 = -BIGINT;
  
  return V10_2;
}
void 
trace_V10_2(FILE *outf, int ****yhx, int j, int d, int mid, int mid1, int mid2, 
	    struct tracekn_s *curr_tr, struct traceknstack_s *dolist, int traceback)
{
  int i;

  i = j - d;

  if (mid > 2 && mid1 > 0){
    if (traceback) fprintf(outf," (V10.2) multiloop yhx %d %d %d %d %d\n", 
			   i+1+mid, mid-1, mid1-1, mid2,  
			   yhx[i+1+mid][mid-1][mid1-1][mid2]);
    if (traceback) fprintf(outf," (V10.2) multiloop zhx %d %d %d %d %d\n", 
			   j-1, d-mid1-4, mid-mid1-mid2-3, d-mid-5,
			   yhx[j-1][d-mid1-4][mid-mid1-mid2-3][d-mid-5]);
    
    PushTraceknstack(dolist, AttachTracekn(curr_tr, i+2, i+1+mid,
					   i+1+mid1, i+1+mid-mid2, dpcS, dpcP));
    PushTraceknstack(dolist, AttachTracekn(curr_tr, i+3+mid1, j-1,
					   i+mid-mid2, i+mid+4, dpcS, dpcP));
  }
}

/*  (V10_3) __/ i+mid1+2 dangles off i+1+mid-mid2,i+1+mid1  /
              / i+mid+3 dangles off i+4+mid,i+mid-mid2      /
              / j-1 dangles off i,j                         */
int 
V10_3(int *s, int len, struct rnapar_2 *rnapar, int **icfg, int ****yhx, int j, int d, int mid, int mid1, int mid2)
{
  int V10_3;
  int     i;

  i = j - d;

  if (mid > 1 && d-mid > 5)
    V10_3 = 3*rnapar->P10P + rnapar->P5P + rnapar->P13 + 4*rnapar->P6P 
      + wkn*icfg[idxL(s[j-1])][idxP(s[j],s[i])]
      + wkn*icfg[idxR(s[i+mid1+2])][idxP(s[i+1+mid-mid2],s[i+1+mid1])]
      + wkn*icfg[idxL(s[i+mid+3])][idxP(s[i+mid+4],s[i+mid-mid2])]
      + wkn*icfg[idxPS(s[i+mid-mid2],s[i+mid+4])]
	        [idxPS(s[i+1+mid-mid2],s[i+1+mid1])] 
      + yhx[i+1+mid][mid][mid1][mid2] 
      + yhx[j-2][d-mid1-5][mid-mid1-mid2-3][d-mid-6];
  else
    V10_3 = -BIGINT;
  
  return V10_3;
}
void 
trace_V10_3(FILE *outf, int ****yhx, int j, int d, int mid, int mid1, int mid2, 
	    struct tracekn_s *curr_tr, struct traceknstack_s *dolist, int traceback)

{
  int i;

  i = j - d;

  if (mid > 1 && d-mid > 5){
    if (traceback) fprintf(outf," (V10.3) multiloop yhx %d %d %d %d %d\n", 
			   i+1+mid, mid, mid1, mid2,  
			   yhx[i+1+mid][mid][mid1][mid2]);
    if (traceback) fprintf(outf," (V10.3) multiloop zhx %d %d %d %d %d\n", 
			   j-2, d-mid1-5, mid-mid1-mid2-3, d-mid-6,
			   yhx[j-2][d-mid1-5][mid-mid1-mid2-3][d-mid-6]);
    
    PushTraceknstack(dolist, AttachTracekn(curr_tr, i+1, i+1+mid,
					   i+1+mid1, i+1+mid-mid2, dpcS, dpcP));
    PushTraceknstack(dolist, AttachTracekn(curr_tr, i+3+mid1, j-2,
					   i+mid-mid2, i+mid+4, dpcS, dpcP));
  }
}

/*  (V10_4) __/ i+mid1+2 dangles off i+1+mid-mid2,i+1+mid1  /
              / i+mid+3 dangles off i+4+mid,i+mid-mid2      /
              / i+1 dangles off i,j                         /
              / j-1 dangles off i,j                         */
int 
V10_4(int *s, int len, struct rnapar_2 *rnapar, int **icfg, int ****yhx, int j, int d, int mid, int mid1, int mid2)
{
  int V10_4;
  int     i;

  i = j - d;

  if (mid > 2 && mid1 > 0 && d-mid > 5)
    V10_4 = 3*rnapar->P10P + rnapar->P5P + rnapar->P13 + 5*rnapar->P6P 
      + wkn*icfg[idxR(s[i+1])][idxP(s[j],s[i])]
      + wkn*icfg[idxL(s[j-1])][idxP(s[j],s[i])]
      + wkn*icfg[idxR(s[i+mid1+2])][idxP(s[i+1+mid-mid2],s[i+1+mid1])]
      + wkn*icfg[idxL(s[i+mid+3])][idxP(s[i+mid+4],s[i+mid-mid2])]
      + wkn*icfg[idxPS(s[i+mid-mid2],s[i+mid+4])]
                [idxPS(s[i+1+mid-mid2],s[i+1+mid1])]
      + yhx[i+1+mid][mid-1][mid1-1][mid2]
      + yhx[j-2][d-mid1-5][mid-mid1-mid2-3][d-mid-6];
  else
    V10_4 = -BIGINT;
  
  return V10_4;
}
void 
trace_V10_4(FILE *outf, int ****yhx, int j, int d, int mid, int mid1, int mid2, 
	    struct tracekn_s *curr_tr, struct traceknstack_s *dolist, int traceback)
{
  int i;

  i = j - d;

  if (mid > 2 && mid1 > 0 && d-mid > 5){
    if (traceback) fprintf(outf," (V10.4) multiloop yhx %d %d %d %d %d\n", 
			   i+1+mid, mid-1, mid1-1, mid2,  
			   yhx[i+1+mid][mid-1][mid1-1][mid2]);
    if (traceback) fprintf(outf," (V10.4) multiloop zhx %d %d %d %d %d\n", 
			   j-2, d-mid1-5, mid-mid1-mid2-3, d-mid-6,
			   yhx[j-2][d-mid1-5][mid-mid1-mid2-3][d-mid-6]);
    
    PushTraceknstack(dolist, AttachTracekn(curr_tr, i+2, i+1+mid,
					   i+1+mid1, i+1+mid-mid2, dpcS, dpcP));
    PushTraceknstack(dolist, AttachTracekn(curr_tr, i+4+mid1, j-2,
					   i+mid-mid2, i+mid+4, dpcS, dpcP));
  }
}

/*(V10_5-V10_8)
 *   
 * non-contiguous coaxial pairs (i+mid1+1, i+1+mid-mid2), (i+mid-mid2-1, i+mid+4)
 * ....add: stack[i+mid-mid2-1][i+2+mid][i+mid-mid2+1][i+mid1+1]
 * check for i+1 or j-1 dangling off the pair (i,j)
 *   
 *  Return: int V10_X          
 */

/*  (V10_5) __/ i+mid1+2 dangles off i+1+mid-mid2,i+1+mid1 /
              / i+mid+3 dangles off i+4+mid,i+mid-mid2-1     */
int 
V10_5(int *s, int len, struct rnapar_2 *rnapar, int **icfg, int ****yhx, int j, int d, int mid, int mid1, int mid2)
{
  int V10_5;
  int     i;

  i = j - d;
  if (mid-mid1-mid2 > 3)
    V10_5 = 3*rnapar->P10P + rnapar->P5P + rnapar->P13 + 4*rnapar->P6P
      + wkn*icfg[idxR(s[i+mid1+2])][idxP(s[i+1+mid-mid2],s[i+1+mid1])]
      + wkn*icfg[idxL(s[i+mid+3])][idxP(s[i+mid+4],s[i+mid-mid2-1])]
      + wkn*icfg[idxPS(s[i+mid-mid2-1],s[i+mid+4])]
                [idxPS(s[i+1+mid-mid2],s[i+1+mid1])] 
      + yhx[i+1+mid][mid][mid1][mid2]
      + yhx[j-1][d-mid1-4][mid-mid1-mid2-4][d-mid-5];
  else
    V10_5 = -BIGINT;
  
  return V10_5;
}
void 
trace_V10_5(FILE *outf, int ****yhx, int j, int d, int mid, int mid1, int mid2, 
	    struct tracekn_s *curr_tr, struct traceknstack_s *dolist, int traceback)
{
  int i;

  i = j - d;

  if (mid-mid1-mid2 > 3){
    if (traceback) fprintf(outf," (V10.5) multiloop yhx %d %d %d %d %d\n", 
			   i+1+mid, mid, mid1, mid2,  
			   yhx[i+1+mid][mid][mid1][mid2]);
    if (traceback) fprintf(outf," (V10.5) multiloop zhx %d %d %d %d %d\n", 
			   j-1, d-mid1-4, mid-mid1-mid2-4, d-mid-5,
			   yhx[j-1][d-mid1-4][mid-mid1-mid2-4][d-mid-5]);
    
    PushTraceknstack(dolist, AttachTracekn(curr_tr, i+1, i+1+mid,
					   i+1+mid1, i+1+mid-mid2, dpcS, dpcP));
    PushTraceknstack(dolist, AttachTracekn(curr_tr, i+3+mid1, j-1,
					   i+mid-mid2-1, i+mid+4, dpcS, dpcP));
  }
}

/*  (V10_6) __/ i+mid1+2 dangles off i+1+mid-mid2,i+1+mid1  /
              / i+mid+3 dangles off i+4+mid,i+mid-mid2-1      /
              / i+1 dangles off i,j                         */
int 
V10_6(int *s, int len, struct rnapar_2 *rnapar, int **icfg, int ****yhx, int j, int d, int mid, int mid1, int mid2)
{
  int V10_6;
  int     i;

  i = j - d;

  if (mid > 2 && mid1 > 0 && mid-mid1-mid2 > 3)
    V10_6 = 3*rnapar->P10P + rnapar->P5P + rnapar->P13 + 5*rnapar->P6P 
      + wkn*icfg[idxR(s[i+1])][idxP(s[j],s[i])]
      + wkn*icfg[idxR(s[i+mid1+2])][idxP(s[i+1+mid-mid2],s[i+1+mid1])]
      + wkn*icfg[idxL(s[i+mid+3])][idxP(s[i+mid+4],s[i+mid-mid2-1])]
      + wkn*icfg[idxPS(s[i+mid-mid2-1],s[i+mid+4])]
	        [idxPS(s[i+1+mid-mid2],s[i+1+mid1])]
      + yhx[i+1+mid][mid-1][mid1-1][mid2] 
      + yhx[j-1][d-mid1-4][mid-mid1-mid2-4][d-mid-5];
  else
    V10_6 = -BIGINT;
  
  return V10_6;
}
void 
trace_V10_6(FILE *outf, int ****yhx, int j, int d, int mid, int mid1, int mid2, 
	    struct tracekn_s *curr_tr, struct traceknstack_s *dolist, int traceback)
{
  int i;

  i = j - d;

  if (mid > 2 && mid1 > 0 && mid-mid1-mid2 > 3){
    if (traceback) fprintf(outf," (V10.6) multiloop yhx %d %d %d %d %d\n", 
			   i+1+mid, mid-1, mid1-1, mid2,  
			   yhx[i+1+mid][mid-1][mid1-1][mid2]);
    if (traceback) fprintf(outf," (V10.6) multiloop zhx %d %d %d %d %d\n", 
			   j-1, d-mid1-4, mid-mid1-mid2-4, d-mid-5,
			   yhx[j-1][d-mid1-4][mid-mid1-mid2-4][d-mid-5]);
    
    PushTraceknstack(dolist, AttachTracekn(curr_tr, i+2, i+1+mid,
					   i+1+mid1, i+1+mid-mid2, dpcS, dpcP));
    PushTraceknstack(dolist, AttachTracekn(curr_tr, i+3+mid1, j-1,
					   i+mid-mid2-1, i+mid+4, dpcS, dpcP));
  }
}

/*  (V10_7) __/ i+mid1+2 dangles off i+1+mid-mid2,i+1+mid1  /
              / i+mid+3 dangles off i+4+mid,i+mid-mid2      /
              / j-1 dangles off i,j                         */
int 
V10_7(int *s, int len, struct rnapar_2 *rnapar, int **icfg, int ****yhx, int j, int d, int mid, int mid1, int mid2)
{
  int V10_7;
  int     i;

  i = j - d;

  if (mid > 1 && d-mid > 5 && mid-mid1-mid2 > 3)
    V10_7 = 3*rnapar->P10P + rnapar->P5P + rnapar->P13 + 5*rnapar->P6P 
      + wkn*icfg[idxL(s[j-1])][idxP(s[j],s[i])]
      + wkn*icfg[idxR(s[i+mid1+2])][idxP(s[i+1+mid-mid2],s[i+1+mid1])]
      + wkn*icfg[idxL(s[i+mid+3])][idxP(s[i+mid+4],s[i+mid-mid2-1])]
      + wkn*icfg[idxPS(s[i+mid-mid2-1],s[i+mid+4])]
	        [idxPS(s[i+1+mid-mid2],s[i+1+mid1])] 
      + yhx[i+1+mid][mid][mid1][mid2] 
      + yhx[j-2][d-mid1-5][mid-mid1-mid2-4][d-mid-6];
  else
    V10_7 = -BIGINT;
  
  return V10_7;
}
void 
trace_V10_7(FILE *outf, int ****yhx, int j, int d, int mid, int mid1, int mid2, 
	    struct tracekn_s *curr_tr, struct traceknstack_s *dolist, int traceback)

{

  int i;

  i = j - d;

  if (mid > 1 && d-mid > 5 && mid-mid1-mid2 > 3){
    if (traceback) fprintf(outf," (V10.7) multiloop yhx %d %d %d %d %d\n", 
			   i+1+mid, mid, mid1, mid2,  
			   yhx[i+1+mid][mid][mid1][mid2]);
    if (traceback) fprintf(outf," (V10.7) multiloop zhx %d %d %d %d %d\n", 
			   j-2, d-mid1-5, mid-mid1-mid2-4, d-mid-6,
			   yhx[j-2][d-mid1-5][mid-mid1-mid2-4][d-mid-6]);
    
    PushTraceknstack(dolist, AttachTracekn(curr_tr, i+1, i+1+mid,
					   i+1+mid1, i+1+mid-mid2, dpcS, dpcP));
    PushTraceknstack(dolist, AttachTracekn(curr_tr, i+3+mid1, j-2,
					   i+mid-mid2-1, i+mid+4, dpcS, dpcP));
  }
}

/*  (V10_8) __/ i+mid1+2 dangles off i+1+mid-mid2,i+1+mid1  /
              / i+mid+3 dangles off i+4+mid,i+mid-mid2-1    /
              / i+1 dangles off i,j                         /
              / j-1 dangles off i,j                         */
int 
V10_8(int *s, int len, struct rnapar_2 *rnapar, int **icfg, int ****yhx, int j, int d, int mid, int mid1, int mid2)
{
  int V10_8;
  int     i;

  i = j - d;

  if (mid > 2 && mid1 > 0 && d-mid > 5 && mid-mid1-mid2 > 3)
    V10_8 = 3*rnapar->P10P + rnapar->P5P + rnapar->P13 + 6*rnapar->P6P 
      + wkn*icfg[idxR(s[i+1])][idxP(s[j],s[i])]
      + wkn*icfg[idxL(s[j-1])][idxP(s[j],s[i])]
      + wkn*icfg[idxR(s[i+mid1+2])][idxP(s[i+1+mid-mid2],s[i+1+mid1])]
      + wkn*icfg[idxL(s[i+mid+3])][idxP(s[i+mid+4],s[i+mid-mid2-1])]
      + wkn*icfg[idxPS(s[i+mid-mid2-1],s[i+mid+4])]
                [idxPS(s[i+1+mid-mid2],s[i+1+mid1])]
      + yhx[i+1+mid][mid-1][mid1-1][mid2]
      + yhx[j-2][d-mid1-5][mid-mid1-mid2-4][d-mid-6];
  else
    V10_8 = -BIGINT;
  
  return V10_8;
}
void 
trace_V10_8(FILE *outf, int ****yhx, int j, int d, int mid, int mid1, int mid2, 
	    struct tracekn_s *curr_tr, struct traceknstack_s *dolist, int traceback)
{
  int i;

  i = j - d;

  if (mid > 2 && mid1 > 0 && d-mid > 5 && mid-mid1-mid2 > 3){
    if (traceback) fprintf(outf," (V10.8) multiloop yhx %d %d %d %d %d\n", 
			   i+1+mid, mid-1, mid1-1, mid2,  
			   yhx[i+1+mid][mid-1][mid1-1][mid2]);
    if (traceback) fprintf(outf," (V10.8) multiloop zhx %d %d %d %d %d\n", 
			   j-2, d-mid1-5, mid-mid1-mid2-4, d-mid-6,
			   yhx[j-2][d-mid1-5][mid-mid1-mid2-4][d-mid-6]);
    
    PushTraceknstack(dolist, AttachTracekn(curr_tr, i+2, i+1+mid,
					   i+1+mid1, i+1+mid-mid2, dpcS, dpcP));
    PushTraceknstack(dolist, AttachTracekn(curr_tr, i+4+mid1, j-2,
					   i+mid-mid2-1, i+mid+4, dpcS, dpcP));
  }
}
