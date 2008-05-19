/* vhxgraphs.c 
 *
 * includes functions to calculate all the diagrams 
 * that fill and traceback the no-hole matrix:
 *      vhx[j][d][d1][d2] (len x len x len x len)
 *                        [(j-d,j) are base pared]
 *                        [(j-d+d1,j-d2) are base pared]
 *
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>

#include "cfg.h"
#include "proto.h"
#include "protovhx.h"
#include "squid.h"

/* (VH1)__/ IRREDUCIBLE SURFACES of  O(2). 
 * (stems, bulges, and internal loops) 
 */
int 
VH1(int *s, int len, int **icfg, int j, int d, int d1, int d2)
{ 
  int VH1;
  
  VH1 = wkn*F2(s, len, icfg, j, d, d1, d2);
  
  return VH1;
}
void
trace_VH1(FILE *outf, int ****vhx, int j, int d, int d1, int d2,
	  struct tracekn_s *curr_tr, struct traceknstack_s *dolist, int traceback)
{
  if (traceback) fprintf(outf," (vh1) IS2 (%d,%d) %d %d %d \n", d1, d2, j, d, vhx[j][d][d1][d2]); 
}

/* (VH2)__/ VHX + VHX */
int 
VH2(int *s, int len, int **icfg, int ****vhx, int j, int d, 
    int d1, int d2, int mid1, int mid2)
{ 
  int VH2;
  int i, k, l;

  i = j - d;
  k = i + d1;
  l = j - d2;
  
  VH2 = vhx[j][d][d1-mid1][d2-mid2]
    + vhx[l+mid2][d-d1-d2+mid1+mid2][mid1][mid2];
  
  return VH2;
}
void
trace_VH2(FILE *outf, int ****vhx, int j, int d, int d1, int d2, int mid1, int mid2,
	  struct tracekn_s *curr_tr, struct traceknstack_s *dolist, int traceback)
{
  int i, k, l;

  i = j - d;
  k = i + d1;
  l = j - d2;

  if (traceback) fprintf(outf," (vh2) trace vhx %d %d %d %d %d \n", 
			 j, d, d1-mid1, d2-mid2,
			 vhx[j][d][d1-mid1][d2-mid2]); 
  if (traceback) fprintf(outf," (vh2) trace vhx %d %d %d %d %d \n", 
			 l+mid2, d-d1-d2+mid1+mid2, mid1, mid2,
			 vhx[l+mid2][d-d1-d2+mid1+mid2][mid1][mid2]); 
  
  PushTraceknstack(dolist, AttachTracekn(curr_tr, i, j, k-mid1, l+mid2, dpcP, dpcP));
  PushTraceknstack(dolist, AttachTracekn(curr_tr, k-mid1, l+mid2, k, l, dpcP, dpcP));
}

/* (VH3)__/ REST of VHX (multiloops).  
 */
/* (VH3.1)__/ no danglings */
int 
VH3_1(int *s, int len, int **icfg, int ****whx, int j, int d, int d1, int d2)
{ 
  int VH3_1;
  int i, k, l;

  i = j - d;
  k = i + d1;
  l = j - d2;
  
  if (d1 > 1 && d2 > 1)
    VH3_1 = 2*P10P + P5P 
      + whx[j-1][d-2][d1-2][d2-2];
  else
    VH3_1 = -BIGINT;

  return VH3_1;
}
void
trace_VH3_1(FILE *outf, int ****whx, int j, int d, int d1, int d2,
	  struct tracekn_s *curr_tr, struct traceknstack_s *dolist, int traceback)
{
  int i, k, l;

  i = j - d;
  k = i + d1;
  l = j - d2;

  if (d1 > 1 && d2 > 1){
    if (traceback) fprintf(outf," (vh3.1) multiloop, trace whx  %d %d %d %d %d \n", 
			   j-1, d-2, d1-2, d2-2, whx[j-1][d-2][d1-2][d2-2]);
    
    PushTraceknstack(dolist, AttachTracekn(curr_tr, i+1, j-1, k-1, l+1, dpcS, dpcS));
  }
}

/*  (VH3.2)__/ i+1  dangles off j,i */
int 
VH3_2(int *s, int len, int **icfg, int ****whx, int j, int d, int d1, int d2)
{ 
  int VH3_2;
  int i, k, l;

  i = j - d;
  k = i + d1;
  l = j - d2;
  
  if(d1 > 2 && d2 > 1)
    VH3_2 = 2*P10P + P5P + P6P 
      + wkn*icfg[idxR(s[i+1])][idxP(s[j],s[i])]
      + whx[j-1][d-3][d1-3][d2-2];
  else
    VH3_2 = -BIGINT;
  
  return VH3_2;
}
void
trace_VH3_2(FILE *outf, int ****whx, int j, int d, int d1, int d2, 
	    struct tracekn_s *curr_tr, struct traceknstack_s *dolist, int traceback)
{
  int i, k, l;

  i = j - d;
  k = i + d1;
  l = j - d2;

  if(d1 > 2 && d2 > 1){
    if (traceback) fprintf(outf," (vh3.2) multiloop, trace whx  %d %d %d %d %d \n", 
			   j-1, d-3, d1-3, d2-2, whx[j-1][d-3][d1-3][d2-2]);
    
    PushTraceknstack(dolist, AttachTracekn(curr_tr, i+2, j-1, k-1, l+1, dpcS, dpcS));
  }
}

/*  (VH3.3)__/ j-1  dangles off i,j */
int 
VH3_3(int *s, int len, int **icfg, int ****whx, int j, int d, int d1, int d2)
{ 
  int VH3_3;
  int i, k, l;

  i = j - d;
  k = i + d1;
  l = j - d2;
  
  if(d1 > 1 && d2 > 2)
    VH3_3 = 2*P10P + P5P + P6P 
      + wkn*icfg[idxL(s[j-1])][idxP(s[j],s[i])]
      + whx[j-2][d-3][d1-2][d2-3];
  else
    VH3_3 = -BIGINT;
  
  return VH3_3;
}
void
trace_VH3_3(FILE *outf, int ****whx, int j, int d, int d1, int d2,
	    struct tracekn_s *curr_tr, struct traceknstack_s *dolist, int traceback)
{
  int i, k, l;

  i = j - d;
  k = i + d1;
  l = j - d2;

  if(d1 > 1 && d2 > 2){
    if (traceback) fprintf(outf," (vh3.3) multiloop, trace whx  %d %d %d %d %d \n", 
			   j-2, d-3, d1-2, d2-3, whx[j-2][d-3][d1-2][d2-3]);
    
    PushTraceknstack(dolist, AttachTracekn(curr_tr, i+1, j-2, k-1, l+1, dpcS, dpcS));
  }
}

/*  (VH3.4)__/ k-1  dangles off l,k */
int 
VH3_4(int *s, int len, int **icfg, int ****whx, int j, int d, int d1, int d2)
{ 
  int VH3_4;
  int i, k, l;

  i = j - d;
  k = i + d1;
  l = j - d2;
  
  if(d1 > 2 && d2 > 1)
    VH3_4 = 2*P10P + P5P + P6P 
      + wkn*icfg[idxL(s[k-1])][idxP(s[l],s[k])]
      + whx[j-1][d-2][d1-3][d2-2];
  else
    VH3_4 = -BIGINT;
  
  return VH3_4;
}
void
trace_VH3_4(FILE *outf, int ****whx, int j, int d, int d1, int d2, 
	    struct tracekn_s *curr_tr, struct traceknstack_s *dolist, int traceback)
{
  int i, k, l;

  i = j - d;
  k = i + d1;
  l = j - d2;

  if(d1 > 2 && d2 > 1){
    if (traceback) fprintf(outf," (vh3.4) multiloop, trace whx  %d %d %d %d %d \n", 
			   j-1, d-2, d1-3, d2-2, whx[j-1][d-2][d1-3][d2-2]);
    
    PushTraceknstack(dolist, AttachTracekn(curr_tr, i+1, j-1, k-2, l+1, dpcS, dpcS));
  }
}

/*  (VH3.5)__/ l+1  dangles off l,k */
int 
VH3_5(int *s, int len, int **icfg, int ****whx, int j, int d, int d1, int d2)
{ 
  int VH3_5;
  int i, k, l;

  i = j - d;
  k = i + d1;
  l = j - d2;
  
  if(d1 > 1 && d2 > 2)
    VH3_5 = 2*P10P + P5P + P6P 
      + wkn*icfg[idxR(s[l+1])][idxP(s[l],s[k])] 
      + whx[j-1][d-2][d1-2][d2-3];
  else
    VH3_5 = -BIGINT;
  
  return VH3_5;
}
void
trace_VH3_5(FILE *outf, int ****whx, int j, int d, int d1, int d2, 
	    struct tracekn_s *curr_tr, struct traceknstack_s *dolist, int traceback)
{
  int i, k, l;

  i = j - d;
  k = i + d1;
  l = j - d2;

  if(d1 > 1 && d2 > 2){
    if (traceback) fprintf(outf," (vh3.5) multiloop, trace whx  %d %d %d %d %d \n", 
			   j-1, d-2, d1-2, d2-3, whx[j-1][d-2][d1-2][d2-3]);
    
    PushTraceknstack(dolist, AttachTracekn(curr_tr, i+1, j-1, k-1, l+2, dpcS, dpcS));
  }
}

/*  (VH3.6)__/ i+1  dangles off j,i / j-1  dangles off j,i */
int 
VH3_6(int *s, int len, int **icfg, int ****whx, int j, int d, int d1, int d2)
{ 
  int VH3_6;
  int i, k, l;

  i = j - d;
  k = i + d1;
  l = j - d2;
  
  if(d1 > 3 && d2 > 1)
    VH3_6 = 2*P10P + P5P + 2*P6P 
      + wkn*icfg[idxR(s[i+1])][idxP(s[j],s[i])] 
      + wkn*icfg[idxL(s[k-1])][idxP(s[l],s[k])] 
      + whx[j-1][d-3][d1-4][d2-2];
  else
    VH3_6 = -BIGINT;
  
  return VH3_6;
}
void
trace_VH3_6(FILE *outf, int ****whx, int j, int d, int d1, int d2,
	    struct tracekn_s *curr_tr, struct traceknstack_s *dolist, int traceback)
{
  int i, k, l;

  i = j - d;
  k = i + d1;
  l = j - d2;

  if(d1 > 3 && d2 > 1){
    if (traceback) fprintf(outf," (vh3.6) multiloop, trace whx  %d %d %d %d %d \n", 
			   j-1, d-3, d1-4, d2-2, whx[j-1][d-3][d1-4][d2-2]);
    
    PushTraceknstack(dolist, AttachTracekn(curr_tr, i+2, j-1, k-2, l+1, dpcS, dpcS));
  }
}

/*  (VH3.7)__/ i+1  dangles off j,i / l+1  dangles off l,k */
int 
VH3_7(int *s, int len, int **icfg, int ****whx, int j, int d, int d1, int d2)
{ 
  int VH3_7;
  int i, k, l;

  i = j - d;
  k = i + d1;
  l = j - d2;
  
  if(d1 > 2 && d2 > 2)
    VH3_7 = 2*P10P + P5P + 2*P6P 
      + wkn*icfg[idxR(s[i+1])][idxP(s[j],s[i])]
      + wkn*icfg[idxR(s[l+1])][idxP(s[l],s[k])]
      + whx[j-1][d-3][d1-3][d2-3];
  else
    VH3_7 = -BIGINT;
  
  return VH3_7;
}
void
trace_VH3_7(FILE *outf, int ****whx, int j, int d, int d1, int d2, 
	    struct tracekn_s *curr_tr, struct traceknstack_s *dolist, int traceback)
{
  int i, k, l;

  i = j - d;
  k = i + d1;
  l = j - d2;

  if(d1 > 2 && d2 > 2){
    if (traceback) fprintf(outf," (vh3.7) multiloop, trace whx  %d %d %d %d %d \n", 
			   j-1, d-3, d1-3, d2-3, whx[j-1][d-3][d1-3][d2-3]);
    
    PushTraceknstack(dolist, AttachTracekn(curr_tr, i+2, j-1, k-1, l+2, dpcS, dpcS));
  }
}

/*  (VH3.8)__/ i+1  dangles off j,i / j-1  dangles off j,i */
int 
VH3_8(int *s, int len, int **icfg, int ****whx, int j, int d, int d1, int d2)
{ 
  int VH3_8;
  int i, k, l;

  i = j - d;
  k = i + d1;
  l = j - d2;
  
  if(d1 > 2 && d2 > 2)
    VH3_8 = 2*P10P + P5P + 2*P6P 
      + wkn*icfg[idxR(s[i+1])][idxP(s[j],s[i])]
      + wkn*icfg[idxL(s[j-1])][idxP(s[j],s[i])]
      + whx[j-2][d-4][d1-3][d2-3];
  else
    VH3_8 = -BIGINT;
  
  return VH3_8;
}
void
trace_VH3_8(FILE *outf, int ****whx, int j, int d, int d1, int d2, 
	    struct tracekn_s *curr_tr, struct traceknstack_s *dolist, int traceback)
{
  int i, k, l;

  i = j - d;
  k = i + d1;
  l = j - d2;

  if(d1 > 2 && d2 > 2){
    if (traceback) fprintf(outf," (vh3.8) multiloop, trace whx  %d %d %d %d %d \n", 
			   j-2, d-4, d1-3, d2-3, whx[j-2][d-4][d1-3][d2-3]);
    
    PushTraceknstack(dolist, AttachTracekn(curr_tr, i+2, j-2, k-1, l+1, dpcS, dpcS));
  }
}

/*  (VH3.9)__/ k-1  dangles off l,k / l+1  dangles off l,k */
int 
VH3_9(int *s, int len, int **icfg, int ****whx, int j, int d, int d1, int d2)
{ 
  int VH3_9;
  int i, k, l;

  i = j - d;
  k = i + d1;
  l = j - d2;
  
  if(d1 > 2 && d2 > 2)
    VH3_9 = 2*P10P + P5P + 2*P6P 
      + wkn*icfg[idxL(s[k-1])][idxP(s[l],s[k])] 
      + wkn*icfg[idxR(s[l+1])][idxP(s[l],s[k])]
      + whx[j-1][d-2][d1-3][d2-3];
  else
    VH3_9 = -BIGINT;
  
  return VH3_9;
}
void
trace_VH3_9(FILE *outf, int ****whx, int j, int d, int d1, int d2, 
	    struct tracekn_s *curr_tr, struct traceknstack_s *dolist, int traceback)
{
  int i, k, l;

  i = j - d;
  k = i + d1;
  l = j - d2;

  if(d1 > 2 && d2 > 2){
    if (traceback) fprintf(outf," (vh3.9) multiloop, trace whx  %d %d %d %d %d \n", 
			   j-1, d-2, d1-3, d2-3, whx[j-1][d-2][d1-3][d2-3]);
    
    PushTraceknstack(dolist, AttachTracekn(curr_tr, i+1, j-1, k-2, l+2, dpcS, dpcS));
  }
}

/*  (VH3.10)__/ k-1  dangles off l,k / j-1  dangles off j,i */
int 
VH3_10(int *s, int len, int **icfg, int ****whx, int j, int d, int d1, int d2)
{ 
  int VH3_10;
  int i, k, l;

  i = j - d;
  k = i + d1;
  l = j - d2;
  
  if(d1 > 2 && d2 > 2)
    VH3_10 = 2*P10P + P5P + 2*P6P 
      + wkn*icfg[idxL(s[k-1])][idxP(s[l],s[k])] 
      + wkn*icfg[idxL(s[j-1])][idxP(s[j],s[i])]
      + whx[j-2][d-3][d1-3][d2-3];
  else
    VH3_10 = -BIGINT;
  
  return VH3_10;
}
void
trace_VH3_10(FILE *outf, int ****whx, int j, int d, int d1, int d2, 
	     struct tracekn_s *curr_tr, struct traceknstack_s *dolist, int traceback)
{
  int i, k, l;

  i = j - d;
  k = i + d1;
  l = j - d2;

  if(d1 > 2 && d2 > 2){
    if (traceback) fprintf(outf," (vh3.10) multiloop, trace whx  %d %d %d %d %d \n", 
			   j-2, d-3, d1-3, d2-3, whx[j-2][d-3][d1-3][d2-3]);
    
    PushTraceknstack(dolist, AttachTracekn(curr_tr, i+1, j-2, k-2, l+1, dpcS, dpcS));
  }
}

/*  (VH3.11)__/ l+1  dangles off l,k / j-1  dangles off j,i */
int 
VH3_11(int *s, int len, int **icfg, int ****whx, int j, int d, int d1, int d2)
{ 
  int VH3_11;
  int i, k, l;

  i = j - d;
  k = i + d1;
  l = j - d2;
  
  if(d1 > 1 && d2 > 3)
    VH3_11 = 2*P10P + P5P + 2*P6P 
      + wkn*icfg[idxR(s[l+1])][idxP(s[l],s[k])]
      + wkn*icfg[idxL(s[j-1])][idxP(s[j],s[i])]
      + whx[j-2][d-3][d1-2][d2-4];
  else
    VH3_11 = -BIGINT;
  
  return VH3_11;
}
void
trace_VH3_11(FILE *outf, int ****whx, int j, int d, int d1, int d2, 
	     struct tracekn_s *curr_tr, struct traceknstack_s *dolist, int traceback)
{
  int i, k, l;

  i = j - d;
  k = i + d1;
  l = j - d2;

  if(d1 > 1 && d2 > 3){
    if (traceback) fprintf(outf," (vh3.11) multiloop, trace whx  %d %d %d %d %d \n", 
			   j-2, d-3, d1-2, d2-4, whx[j-2][d-3][d1-2][d2-4]);
    
    PushTraceknstack(dolist, AttachTracekn(curr_tr, i+1, j-2, k-1, l+2, dpcS, dpcS));
  }
}

/*  (VH3.12)__/ i+1  dangles off j,i / l+1  dangles off l,k / j-1  dangles off j,i */
int 
VH3_12(int *s, int len, int **icfg, int ****whx, int j, int d, int d1, int d2)
{ 
  int VH3_12;
  int i, k, l;

  i = j - d;
  k = i + d1;
  l = j - d2;
  
  if(d1 > 2 && d2 > 3)
    VH3_12 = 2*P10P + P5P + 3*P6P 
      + wkn*icfg[idxR(s[i+1])][idxP(s[j],s[i])] 
      + wkn*icfg[idxR(s[l+1])][idxP(s[l],s[k])] 
      + wkn*icfg[idxL(s[j-1])][idxP(s[j],s[i])]
      + whx[j-2][d-4][d1-3][d2-4];
  else
    VH3_12 = -BIGINT;
  
  return VH3_12;
}
void
trace_VH3_12(FILE *outf, int ****whx, int j, int d, int d1, int d2, 
	     struct tracekn_s *curr_tr, struct traceknstack_s *dolist, int traceback)
{
  int i, k, l;

  i = j - d;
  k = i + d1;
  l = j - d2;

  if(d1 > 2 && d2 > 3){
    if (traceback) fprintf(outf," (vh3.12) multiloop, trace whx  %d %d %d %d %d \n", 
			   j-2, d-4, d1-3, d2-3, whx[j-2][d-4][d1-3][d2-4]);
    
    PushTraceknstack(dolist, AttachTracekn(curr_tr, i+2, j-2, k-1, l+2, dpcS, dpcS));
  }
}

/*  (VH3.13)__/ i+1  dangles off j,i / k-1  dangles off l,k / l+1  dangles off l,k */
int 
VH3_13(int *s, int len, int **icfg, int ****whx, int j, int d, int d1, int d2)
{ 
  int VH3_13;
  int i, k, l;

  i = j - d;
  k = i + d1;
  l = j - d2;
  
  if(d1 > 3 && d2 > 2)
    VH3_13 = 2*P10P + P5P + 3*P6P 
      + wkn*icfg[idxR(s[i+1])][idxP(s[j],s[i])] 
      + wkn*icfg[idxL(s[k-1])][idxP(s[l],s[k])] 
      + wkn*icfg[idxR(s[l+1])][idxP(s[l],s[k])] 
      + whx[j-1][d-3][d1-4][d2-3];
  else
    VH3_13 = -BIGINT;
  
  return VH3_13;
}
void
trace_VH3_13(FILE *outf, int ****whx, int j, int d, int d1, int d2, 
	     struct tracekn_s *curr_tr, struct traceknstack_s *dolist, int traceback)
{
  int i, k, l;

  i = j - d;
  k = i + d1;
  l = j - d2;

  if(d1 > 3 && d2 > 2){
    if (traceback) fprintf(outf," (vh3.13) multiloop, trace whx  %d %d %d %d %d \n", 
			   j-1, d-3, d1-4, d2-3, whx[j-1][d-3][d1-4][d2-3]);
    
    PushTraceknstack(dolist, AttachTracekn(curr_tr, i+2, j-1, k-2, l+2, dpcS, dpcS));
  }
}

/*  (VH3.14)__/ i+1  dangles off j,i / k-1  dangles off l,k / j-1  dangles off j,i */
int 
VH3_14(int *s, int len, int **icfg, int ****whx, int j, int d, int d1, int d2)
{ 
  int VH3_14;
  int i, k, l;

  i = j - d;
  k = i + d1;
  l = j - d2;
  
  if(d1 > 3 && d2 > 2)
    VH3_14 = 2*P10P + P5P + 3*P6P 
      + wkn*icfg[idxR(s[i+1])][idxP(s[j],s[i])] 
      + wkn*icfg[idxL(s[k-1])][idxP(s[l],s[k])]
      + wkn*icfg[idxL(s[j-1])][idxP(s[j],s[i])]
      + whx[j-2][d-4][d1-4][d2-3];
  else
    VH3_14 = -BIGINT;
  
  return VH3_14;
}
void
trace_VH3_14(FILE *outf, int ****whx, int j, int d, int d1, int d2,
	     struct tracekn_s *curr_tr, struct traceknstack_s *dolist, int traceback)
{
  int i, k, l;

  i = j - d;
  k = i + d1;
  l = j - d2;

  if(d1 > 3 && d2 > 2){
    if (traceback) fprintf(outf," (vh3.14) multiloop, trace whx  %d %d %d %d %d \n", 
			   j-2, d-4, d1-4, d2-3, whx[j-2][d-4][d1-4][d2-3]);
    
    PushTraceknstack(dolist, AttachTracekn(curr_tr, i+2, j-2, k-2, l+1, dpcS, dpcS));
  }
}

/*  (VH3.15)__/ k-1  dangles off l,k / l+1  dangles off l,k / j-1  dangles off j,i */
int 
VH3_15(int *s, int len, int **icfg, int ****whx, int j, int d, int d1, int d2)
{ 
  int VH3_15;
  int i, k, l;

  i = j - d;
  k = i + d1;
  l = j - d2;
  
  if(d1 > 2 && d2 > 3)
    VH3_15 = 2*P10P + P5P + 3*P6P 
      + wkn*icfg[idxL(s[k-1])][idxP(s[l],s[k])]
      + wkn*icfg[idxR(s[l+1])][idxP(s[l],s[k])] 
      + wkn*icfg[idxL(s[j-1])][idxP(s[j],s[i])] 
      + whx[j-2][d-3][d1-3][d2-4];
  else
    VH3_15 = -BIGINT;
  
  return VH3_15;
}
void
trace_VH3_15(FILE *outf, int ****whx, int j, int d, int d1, int d2, 
	     struct tracekn_s *curr_tr, struct traceknstack_s *dolist, int traceback)
{
  int i, k, l;

  i = j - d;
  k = i + d1;
  l = j - d2;

  if(d1 > 2 && d2 > 3){
    if (traceback) fprintf(outf," (vh3.15) multiloop, trace whx  %d %d %d %d %d \n", 
			   j-2, d-3, d1-3, d2-4, whx[j-2][d-3][d1-3][d2-4]);
    
    PushTraceknstack(dolist, AttachTracekn(curr_tr, i+1, j-2, k-2, l+2, dpcS, dpcS));
  }
}

/*  (VH3.16)__/ i+1  dangles off j,i / k-1  dangles off l,k / 
              / l+1  dangles off l,k / j-1  dangles off j,i */
int 
VH3_16(int *s, int len, int **icfg, int ****whx, int j, int d, int d1, int d2)
{ 
  int VH3_16;
  int i, k, l;

  i = j - d;
  k = i + d1;
  l = j - d2;
  
  if(d1 > 3 && d2 > 3)
    VH3_16 = 2*P10P + P5P + 4*P6P 
      + wkn*icfg[idxR(s[i+1])][idxP(s[j],s[i])] 
      + wkn*icfg[idxL(s[k-1])][idxP(s[l],s[k])] 
      + wkn*icfg[idxR(s[l+1])][idxP(s[l],s[k])] 
      + wkn*icfg[idxL(s[j-1])][idxP(s[j],s[i])]
      + whx[j-2][d-4][d1-4][d2-4];
  else
    VH3_16 = -BIGINT;
  
  return VH3_16;
}
void
trace_VH3_16(FILE *outf, int ****whx, int j, int d, int d1, int d2, 
	     struct tracekn_s *curr_tr, struct traceknstack_s *dolist, int traceback)
{
  int i, k, l;

  i = j - d;
  k = i + d1;
  l = j - d2;

  if(d1 > 3 && d2 > 3){
    if (traceback) fprintf(outf,"(vh3.16)  multiloop, trace whx  %d %d %d %d %d \n", 
			   j-2, d-4, d1-4, d2-4, whx[j-2][d-4][d1-4][d2-4]);
    
    PushTraceknstack(dolist, AttachTracekn(curr_tr, i+2, j-2, k-2, l+2, dpcS, dpcS));
  }
}


