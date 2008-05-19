/* whxgraphs.c 
 *
 * includes functions to calculate all the diagrams 
 * that fill and traceback the no-hole matrix:
 *      whx[j][d][d1][d2] (len x len x len x len)
 *
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>

#include "cfg.h"
#include "proto.h"
#include "protowhx.h"
#include "squid.h"

/* (WH1)__/ STRUCTURES WITH ONE VHX (sixteen).
 */ 

/* (WH1.1)__/ no danglings  */
int 
WH1_1(int *s, int len, int **icfg, int ****vhx, int j, int d, int d1, int d2)
{ 
  int   WH1_1;
  int i, k, l;

  i = j - d;
  k = i + d1;
  l = j - d2;

  WH1_1 = 2*P10P + vhx[j][d][d1][d2];

  return WH1_1;
}
void
trace_WH1_1(FILE *outf, int ****vhx, int j, int d, int d1, int d2, 
	    struct tracekn_s *curr_tr, struct traceknstack_s *dolist, int traceback)
{
  int i, k, l;

  i = j - d;
  k = i + d1;
  l = j - d2;

  if (traceback) fprintf(outf," (WH1.1) trace vhx %d %d %d %d %d \n", 
			 j, d, d1, d2, vhx[j][d][d1][d2]); 
  
  curr_tr->type1 = dpcP;              
  curr_tr->type2 = dpcP;              
  PushTraceknstack(dolist, AttachTracekn(curr_tr, i, j, k, l, dpcP, dpcP)); 
}

/* TWO danglings (five) */
 
/* (WH1.2)__/ i and k dangling (j,d-1,d1-2,d2) */
int 
WH1_2(int *s, int len, int **icfg, int ****vhx, int j, int d, int d1, int d2)
{ 
  int   WH1_2;
  int i, k, l;

  i = j - d;
  k = i + d1;
  l = j - d2;

  if (d1 > 1)
    WH1_2 = 2*P10P + 2*P6P 
      + wkn*icfg[idxL(s[i])][idxP(s[i+1],s[j])] 
      + wkn*icfg[idxR(s[k])][idxP(s[k-1],s[l])] 
      + vhx[j][d-1][d1-2][d2];
  else
    WH1_2 = -BIGINT;

  return WH1_2;
}
void
trace_WH1_2(FILE *outf, int ****vhx, int j, int d, int d1, int d2, 
	    struct tracekn_s *curr_tr, struct traceknstack_s *dolist, int traceback)
{
  int i, k, l;

  i = j - d;
  k = i + d1;
  l = j - d2;

  if (traceback) fprintf(outf," (WH1.2) SLE, SRI, trace vhx %d %d %d %d %d \n", j, d-1, d1-2, d2,
			 vhx[j][d-1][d1-2][d2]); 
  
  curr_tr->type1 = dpcL;              
  curr_tr->type2 = dpcL;              
  PushTraceknstack(dolist, AttachTracekn(curr_tr, i+1, j, k-1, l, dpcP, dpcP)); 
}

/* (WH1.3)__/ l and j dangling (j-1,d-1,d1,d2-2) */
int 
WH1_3(int *s, int len, int **icfg, int ****vhx, int j, int d, int d1, int d2)
{ 
  int   WH1_3;
  int i, k, l;

  i = j - d;
  k = i + d1;
  l = j - d2;

  if (d2 > 1)
    WH1_3 = 2*P10P + 2*P6P 
      + wkn*icfg[idxL(s[l])][idxP(s[k],s[l+1])]
      + wkn*icfg[idxR(s[j])][idxP(s[i],s[j-1])]
      + vhx[j-1][d-1][d1][d2-2];
  else
    WH1_3 = -BIGINT;

  return WH1_3;
}
void
trace_WH1_3(FILE *outf, int ****vhx, int j, int d, int d1, int d2, 
	    struct tracekn_s *curr_tr, struct traceknstack_s *dolist, int traceback)
{
  int i, k, l;

  i = j - d;
  k = i + d1;
  l = j - d2;

  if (traceback) fprintf(outf," (WH1.3) SLI, SRE, trace vhx %d %d %d %d %d \n", j-1, d-1, d1, d2-2,
			 vhx[j-1][d-1][d1][d2-2]); 
  
  curr_tr->type1 = dpcL;              
  curr_tr->type2 = dpcL;              
  PushTraceknstack(dolist, AttachTracekn(curr_tr, i, j-1, k, l+1, dpcP, dpcP)); 
}

/* (WH1.4)__/ i and l dangling.  (j,d-1,d1-1,d2-1) */
 int 
WH1_4(int *s, int len, int **icfg, int ****vhx, int j, int d, int d1, int d2)
{ 
  int   WH1_4;
  int i, k, l;

  i = j - d;
  k = i + d1;
  l = j - d2;

  if (d1 > 0 && d2 > 0)
    WH1_4 = 2*P10P + 2*P6P 
      + wkn*icfg[idxL(s[i])][idxP(s[i+1],s[j])]
      + wkn*icfg[idxL(s[l])][idxP(s[k],s[l+1])]
      + vhx[j][d-1][d1-1][d2-1];
  else
    WH1_4 = -BIGINT;

  return WH1_4;
}
void
trace_WH1_4(FILE *outf, int ****vhx, int j, int d, int d1, int d2, 
	    struct tracekn_s *curr_tr, struct traceknstack_s *dolist, int traceback)
{
  int i, k, l;

  i = j - d;
  k = i + d1;
  l = j - d2;

  if (traceback) fprintf(outf," (WH1.4) SLE, SLI, trace vhx %d %d %d %d %d  \n", j, d-1, d1-1, d2-1,
			 vhx[j][d-1][d1-1][d2-1]); 
  
  curr_tr->type1 = dpcL;              
  curr_tr->type2 = dpcL;              
  PushTraceknstack(dolist, AttachTracekn(curr_tr, i+1, j, k, l+1, dpcP, dpcP)); 
}

/* (WH1.5)__/ k and j dangling (j-1,d-1,d1-1,d2-1) */
int 
WH1_5(int *s, int len, int **icfg, int ****vhx, int j, int d, int d1, int d2)
{ 
  int   WH1_5;
  int i, k, l;

  i = j - d;
  k = i + d1;
  l = j - d2;

  if (d1 > 0 && d2 > 0)
    WH1_5 = 2*P10P + 2*P6P 
      + wkn*icfg[idxR(s[k])][idxP(s[k-1],s[l])]
      + wkn*icfg[idxR(s[j])][idxP(s[i],s[j-1])]
      + vhx[j-1][d-1][d1-1][d2-1];
  else
    WH1_5 = -BIGINT;

  return WH1_5;
}
void
trace_WH1_5(FILE *outf, int ****vhx, int j, int d, int d1, int d2, 
	    struct tracekn_s *curr_tr, struct traceknstack_s *dolist, int traceback)
{
  int i, k, l;

  i = j - d;
  k = i + d1;
  l = j - d2;

  if (traceback) fprintf(outf," (WH1.5) SRI, SRE, trace vhx %d %d %d %d %d \n", j-1, d-1, d1-1, d2-1, 
			 vhx[j-1][d-1][d1-1][d2-1]); 
  
  curr_tr->type1 = dpcL;              
  curr_tr->type2 = dpcL;              
  PushTraceknstack(dolist, AttachTracekn(curr_tr, i, j-1, k-1, l, dpcP, dpcP)); 
}


/* THREE danglings (four) */

/* (WH1.6)__/ i, k, l dangling.   (j,d-1,d1-2,d2-1) */
int 
WH1_6(int *s, int len, int **icfg, int ****vhx, int j, int d, int d1, int d2)
{ 
  int   WH1_6;
  int i, k, l;

  i = j - d;
  k = i + d1;
  l = j - d2;

  if (d1 > 1 && d2 > 0)
    WH1_6 = 2*P10P + 3*P6P 
      + wkn*icfg[idxL(s[i])][idxP(s[i+1],s[j])]
      + wkn*icfg[idxR(s[k])][idxP(s[k-1],s[l+1])]
      + wkn*icfg[idxL(s[l])][idxP(s[k-1],s[l+1])] 
      + vhx[j][d-1][d1-2][d2-1];
  else
    WH1_6 = -BIGINT;

  return WH1_6;
}
void
trace_WH1_6(FILE *outf, int ****vhx, int j, int d, int d1, int d2, 
	    struct tracekn_s *curr_tr, struct traceknstack_s *dolist, int traceback)
{
  int i, k, l;

  i = j - d;
  k = i + d1;
  l = j - d2;

  if (traceback) fprintf(outf," (WH1.6) SLE, SRI, SLI, trace vhx %d %d %d %d %d \n", j, d-1, d1-2, d2-1, 
			 vhx[j][d-1][d1-2][d2-1]); 
  
  curr_tr->type1 = dpcL;              
  curr_tr->type2 = dpcL;              
  PushTraceknstack(dolist, AttachTracekn(curr_tr, i+1, j, k-1, l+1, dpcP, dpcP)); 
}

/* (WH1.7)__/ i, k, j dangling (j-1,d-2,d1-2,d2-1) */
int 
WH1_7(int *s, int len, int **icfg, int ****vhx, int j, int d, int d1, int d2)
{ 
  int   WH1_7;
  int i, k, l;

  i = j - d;
  k = i + d1;
  l = j - d2;

  if (d1 > 1 && d2 > 0)
    WH1_7 = 2*P10P + 3*P6P 
      + wkn*icfg[idxL(s[i])][idxP(s[i+1],s[j-1])]
      + wkn*icfg[idxR(s[j])][idxP(s[i+1],s[j-1])]
      + wkn*icfg[idxR(s[k])][idxP(s[k-1],s[l])] 
      + vhx[j-1][d-2][d1-2][d2-1];
  else
    WH1_7 = -BIGINT;

  return WH1_7;
}
void
trace_WH1_7(FILE *outf, int ****vhx, int j, int d, int d1, int d2, 
	    struct tracekn_s *curr_tr, struct traceknstack_s *dolist, int traceback)
{
  int i, k, l;

  i = j - d;
  k = i + d1;
  l = j - d2;

  if (traceback) fprintf(outf," (WH1.7) SLE, SRI, SRE, trace vhx %d %d %d %d %d \n", j-1, d-2, d1-2, d2-1,
			 vhx[j-1][d-2][d1-2][d2-1]); 
  
  curr_tr->type1 = dpcL;              
  curr_tr->type2 = dpcL;              
  PushTraceknstack(dolist, AttachTracekn(curr_tr, i+1, j-1, k-1, l, dpcP, dpcP)); 
}

/* (WH1.8)__/  i, l, j dangling (j-1,d-2,d1-1,d2-2) */
int 
WH1_8(int *s, int len, int **icfg, int ****vhx, int j, int d, int d1, int d2)
{ 
  int   WH1_8;
  int i, k, l;

  i = j - d;
  k = i + d1;
  l = j - d2;

  if (d1 > 0 && d2 > 1)
    WH1_8 = 2*P10P + 3*P6P 
      + wkn*icfg[idxL(s[i])][idxP(s[i+1],s[j-1])]
      + wkn*icfg[idxR(s[j])][idxP(s[i+1],s[j-1])] 
      + wkn*icfg[idxL(s[l])][idxP(s[k],s[l+1])]
      + vhx[j-1][d-2][d1-1][d2-2];
  else
    WH1_8 = -BIGINT;

  return WH1_8;
}
void
trace_WH1_8(FILE *outf, int ****vhx, int j, int d, int d1, int d2, 
	    struct tracekn_s *curr_tr, struct traceknstack_s *dolist, int traceback)
{
  int i, k, l;

  i = j - d;
  k = i + d1;
  l = j - d2;

  if (traceback) fprintf(outf," (WH1.8) SLE, SLI, SRE, trace vhx %d %d %d %d %d \n", j-1, d-2, d1-1, d2-2,
			 vhx[j-1][d-2][d1-1][d2-2]); 
  
  curr_tr->type1 = dpcL;              
  curr_tr->type2 = dpcL;              
  PushTraceknstack(dolist, AttachTracekn(curr_tr, i+1, j-1, k, l+1, dpcP, dpcP)); 
}

/* (WH1.9)__/  k, l, j dangling (j-1,d-1,d1-1,d2-2) */
int 
WH1_9(int *s, int len, int **icfg, int ****vhx, int j, int d, int d1, int d2)
{ 
  int   WH1_9;
  int i, k, l;

  i = j - d;
  k = i + d1;
  l = j - d2;

  if (d1 > 0 && d2 > 1)
    WH1_9 = 2*P10P + 3*P6P 
      + wkn*icfg[idxR(s[k])][idxP(s[k-1],s[l+1])]
      + wkn*icfg[idxL(s[l])][idxP(s[k-1],s[l+1])] 
      + wkn*icfg[idxR(s[j])][idxP(s[i],s[j-1])]
      + vhx[j-1][d-1][d1-1][d2-2];
  else
    WH1_9 = -BIGINT;

  return WH1_9;
}
void
trace_WH1_9(FILE *outf, int ****vhx, int j, int d, int d1, int d2, 
	    struct tracekn_s *curr_tr, struct traceknstack_s *dolist, int traceback)
{
  int i, k, l;

  i = j - d;
  k = i + d1;
  l = j - d2;

  if (traceback) fprintf(outf," (WH1.9) SRI, SLI, SRE, trace vhx %d %d %d %d %d \n", j-1, d-1, d1-1, d2-2,
			 vhx[j-1][d-1][d1-1][d2-2]); 
  
  curr_tr->type1 = dpcL;              
  curr_tr->type2 = dpcL;              
  PushTraceknstack(dolist, 
		   AttachTracekn(curr_tr, i, j-1, k-1, l+1, dpcP, dpcP)); 
}

/* (WH1.10)__/  FOUR  danglings (one diagram) */
int 
WH1_10(int *s, int len, int **icfg, int ****vhx, int j, int d, int d1, int d2)
{ 
  int   WH1_10;
  int i, k, l;

  i = j - d;
  k = i + d1;
  l = j - d2;

  if (d1 > 1 && d2 > 1)
    WH1_10 = 2*P10P + 4*P6P 
      + wkn*icfg[idxL(s[i])][idxP(s[i+1],s[j-1])]
      + wkn*icfg[idxR(s[j])][idxP(s[i+1],s[j-1])]
      + wkn*icfg[idxR(s[k])][idxP(s[k-1],s[l+1])]
      + wkn*icfg[idxL(s[l])][idxP(s[k-1],s[l+1])] 
      + vhx[j-1][d-2][d1-2][d2-2];
  else
    WH1_10 = -BIGINT;

  return WH1_10;
}
void
trace_WH1_10(FILE *outf, int ****vhx, int j, int d, int d1, int d2, 
	     struct tracekn_s *curr_tr, struct traceknstack_s *dolist, int traceback)
{
  int i, k, l;

  i = j - d;
  k = i + d1;
  l = j - d2;

  if (traceback) fprintf(outf," (WH1.10) SLE, SRI, SLI, SRE,  %d %d %d %d %d \n", j-1, d-2, d1-2, d2-2, 
			 vhx[j-1][d-2][d1-2][d2-2]);  
  
  curr_tr->type1 = dpcL;              
  curr_tr->type2 = dpcL;              
  PushTraceknstack(dolist, AttachTracekn(curr_tr, i+1, j-1, k-1, l+1, dpcP, dpcP)); 
}

/*(WH2)__/ STRUCTURES WITH ONE ZHX (four) */ 
  
/* (WH2.1)__/ Pair (i-j) */
int 
WH2_1(int *s, int len, int **icfg, int ****zhx, int j, int d, int d1, int d2)
{ 
  int   WH2_1;
  int i, k, l;

  i = j - d;
  k = i + d1;
  l = j - d2;

  WH2_1 = P10P + zhx[j][d][d1][d2];

  return WH2_1;
}
void
trace_WH2_1(FILE *outf, int ****zhx, int j, int d, int d1, int d2, 
	    struct tracekn_s *curr_tr, struct traceknstack_s *dolist, int traceback)
{
  int i, k, l;

  i = j - d;
  k = i + d1;
  l = j - d2;

  if (traceback) fprintf(outf," (WH2.1) trace zhx %d %d %d %d %d \n", j, d, d1, d2,
			 zhx[j][d][d1][d2]); 
  
  curr_tr->type1 = dpcP;              
  PushTraceknstack(dolist, AttachTracekn(curr_tr, i, j, k, l, dpcP, dpcS)); 
}

/* (WH2.2)__/ Pair (i+1-j-1)/ i and j dangle off i+1,j-1 */
int 
WH2_2(int *s, int len, int **icfg, int ****zhx, int j, int d, int d1, int d2)
{ 
  int   WH2_2;
  int i, k, l;

  i = j - d;
  k = i + d1;
  l = j - d2;

  if (d1 > 0 && d2 > 0)
    WH2_2 = P10P + 2*P6P 
      + wkn*icfg[idxL(s[i])][idxP(s[i+1],s[j-1])]
      + wkn*icfg[idxR(s[j])][idxP(s[i+1],s[j-1])] 
      + zhx[j-1][d-2][d1-1][d2-1];
  else
    WH2_2 = -BIGINT;

  return WH2_2;
}
void
trace_WH2_2(FILE *outf, int ****zhx, int j, int d, int d1, int d2, 
	    struct tracekn_s *curr_tr, struct traceknstack_s *dolist, int traceback)
{
  int i, k, l;

  i = j - d;
  k = i + d1;
  l = j - d2;

  if (traceback) fprintf(outf," (WH2.2) SLE, SRE, trace zhx %d %d %d %d %d \n", j-1, d-2, d1-1, d2-1,
			 zhx[j-1][d-2][d1-1][d2-1]); 
  
  curr_tr->type1 = dpcL;              
  PushTraceknstack(dolist, AttachTracekn(curr_tr, i+1, j-1, k, l, dpcP, dpcS)); 
}

/* (WH2.3)__/ Pair (i+1-j) / i dangles off i+1,j */
int 
WH2_3(int *s, int len, int **icfg, int ****zhx, int j, int d, int d1, int d2)
{ 
  int   WH2_3;
  int i, k, l;

  i = j - d;
  k = i + d1;
  l = j - d2;

  if (d1 > 0)
    WH2_3 = P10P + P6P 
      + wkn*icfg[idxL(s[i])][idxP(s[i+1],s[j])]
      + zhx[j][d-1][d1-1][d2];
  else
    WH2_3 = -BIGINT;

  return WH2_3;
}
void
trace_WH2_3(FILE *outf, int ****zhx, int j, int d, int d1, int d2, 
	    struct tracekn_s *curr_tr, struct traceknstack_s *dolist, int traceback)
{
  int i, k, l;

  i = j - d;
  k = i + d1;
  l = j - d2;

  if (traceback) fprintf(outf," (WH2.3) SLE, trace zhx  %d %d %d %d %d \n", j, d-1, d1-1, d2,
			 zhx[j][d-1][d1-1][d2]); 
  
  curr_tr->type1 = dpcL;              
  PushTraceknstack(dolist, AttachTracekn(curr_tr, i+1, j, k, l, dpcP, dpcS)); 
}

/* (WH2.4)__/ Pair (i-j-1) / j dangles off i,j-1 */
int 
WH2_4(int *s, int len, int **icfg, int ****zhx, int j, int d, int d1, int d2)
{ 
  int   WH2_4;
  int i, k, l;

  i = j - d;
  k = i + d1;
  l = j - d2;

  if (d2 > 0)
    WH2_4 = P10P + P6P 
      + wkn*icfg[idxR(s[j])][idxP(s[i],s[j-1])]
      + zhx[j-1][d-1][d1][d2-1];
  else
    WH2_4 = -BIGINT;

  return WH2_4;
}
void
trace_WH2_4(FILE *outf, int ****zhx, int j, int d, int d1, int d2, 
	    struct tracekn_s *curr_tr, struct traceknstack_s *dolist, int traceback)
{
  int i, k, l;

  i = j - d;
  k = i + d1;
  l = j - d2;

  if (traceback) fprintf(outf," (WH2.4) SLE, trace zhx  %d %d %d %d %d \n", j-1, d-1, d1, d2-1,
			 zhx[j-1][d-1][d1][d2-1]); 
  
  curr_tr->type1 = dpcL;              
  PushTraceknstack(dolist, AttachTracekn(curr_tr, i, j-1, k, l, dpcP, dpcS)); 
}

/* (WH3)__/ STRUCTURES WITH ONE YHX (four) */ 

/* (WH3.1)__/ no danglings / Pair (k-l) */
int 
WH3_1(int *s, int len, int **icfg, int ****yhx, int j, int d, int d1, int d2)
{ 
  int   WH3_1;
  int i, k, l;

  i = j - d;
  k = i + d1;
  l = j - d2;

  WH3_1 = P10P + yhx[j][d][d1][d2];

  return WH3_1;
}
void
trace_WH3_1(FILE *outf, int ****yhx, int j, int d, int d1, int d2, 
	    struct tracekn_s *curr_tr, struct traceknstack_s *dolist, int traceback)
{
  int i, k, l;

  i = j - d;
  k = i + d1;
  l = j - d2;

  if (traceback) fprintf(outf," (WH3.1) trace yhx %d %d %d %d %d \n", j, d, d1, d2,
			 yhx[j][d][d1][d2]); 
  
  curr_tr->type2 = dpcP;              
  PushTraceknstack(dolist, AttachTracekn(curr_tr, i, j, k, l, dpcS, dpcP)); 
}

/* (WH3.2)__/ Pair (k-1-l+1)/ k and l dangle off l+1,k-1 */
int 
WH3_2(int *s, int len, int **icfg, int ****yhx, int j, int d, int d1, int d2)
{ 
  int   WH3_2;
  int i, k, l;

  i = j - d;
  k = i + d1;
  l = j - d2;

  if (d1 > 0 && d2 > 0)
    WH3_2 = P10P + 2*P6P 
      + wkn*icfg[idxR(s[k])][idxP(s[k-1],s[l+1])]
      + wkn*icfg[idxL(s[l])][idxP(s[k-1],s[l+1])] 
      + yhx[j][d][d1-1][d2-1];
  else
    WH3_2 = -BIGINT;

  return WH3_2;
}
void
trace_WH3_2(FILE *outf, int ****yhx, int j, int d, int d1, int d2, 
	    struct tracekn_s *curr_tr, struct traceknstack_s *dolist, int traceback)
{
  int i, k, l;

  i = j - d;
  k = i + d1;
  l = j - d2;

  if (traceback) fprintf(outf," (W43.2) SRI, SLI, trace yhx %d %d %d %d %d \n", j, d, d1-1, d2-1,
			 yhx[j][d][d1-1][d2-1]); 
  
  curr_tr->type2 = dpcL;              
  PushTraceknstack(dolist, AttachTracekn(curr_tr, i, j, k-1, l+1, dpcS, dpcP)); 
}

/* (WH3.3)__/ Pair (k-l+1)/ l dangles off l+1,k */
int 
WH3_3(int *s, int len, int **icfg, int ****yhx, int j, int d, int d1, int d2)
{ 
  int   WH3_3;
  int i, k, l;

  i = j - d;
  k = i + d1;
  l = j - d2;

  if (d2 > 0)
    WH3_3 = P10P + P6P 
      + wkn*icfg[idxL(s[l])][idxP(s[k],s[l+1])] 
      + yhx[j][d][d1][d2-1];
  else
    WH3_3 = -BIGINT;

  return WH3_3;
}
void
trace_WH3_3(FILE *outf, int ****yhx, int j, int d, int d1, int d2, 
	  struct tracekn_s *curr_tr, struct traceknstack_s *dolist, int traceback)
{
  int i, k, l;

  i = j - d;
  k = i + d1;
  l = j - d2;

  if (traceback) fprintf(outf," (WH3.3) SLI, trace yhx  %d %d %d %d %d  \n", j, d, d1, d2-1,
			 yhx[j][d][d1][d2-1]); 
  
  curr_tr->type2 = dpcL;               
  PushTraceknstack(dolist, AttachTracekn(curr_tr, i, j, k, l+1, dpcS, dpcP)); 
}

/* (WH3.4)__/ Pair (k-1-l)/ k dangles off l,k-1 */
int 
WH3_4(int *s, int len, int **icfg, int ****yhx, int j, int d, int d1, int d2)
{ 
  int   WH3_4;
  int i, k, l;

  i = j - d;
  k = i + d1;
  l = j - d2;

  if (d1 > 0)
    WH3_4 = P10P + P6P 
      + wkn*icfg[idxR(s[k])][idxP(s[k-1],s[l])]
      + yhx[j][d][d1-1][d2];
  else
    WH3_4 = -BIGINT;

  return WH3_4;
}
void
trace_WH3_4(FILE *outf, int ****yhx, int j, int d, int d1, int d2, 
	  struct tracekn_s *curr_tr, struct traceknstack_s *dolist, int traceback)
{
  int i, k, l;

  i = j - d;
  k = i + d1;
  l = j - d2;

  if (traceback) fprintf(outf," (WH3.4) SRI, trace yhx %d %d %d %d %d \n", j, d, d1-1, d2,
			 yhx[j][d][d1-1][d2]); 
  
  curr_tr->type2 = dpcL;              
  PushTraceknstack(dolist, AttachTracekn(curr_tr, i, j, k-1, l, dpcS, dpcP));
}

/* (WH4)__/ STRUCTURES WITH ONE WHX (four) */ 

/* (WH4.1)__/ Singlet-left exterior / i dangles off i+1,j  */
int 
WH4_1(int *s, int len, int **icfg, int ****whx, int j, int d, int d1, int d2)
{ 
  int   WH4_1;
  int i, k, l;

  i = j - d;
  k = i + d1;
  l = j - d2;

  if (d1 > 0)
    WH4_1 = P6P + whx[j][d-1][d1-1][d2];
  else
    WH4_1 = -BIGINT;

  return WH4_1;
}
void
trace_WH4_1(FILE *outf, int ****whx, int j, int d, int d1, int d2, 
	    struct tracekn_s *curr_tr, struct traceknstack_s *dolist, int traceback)
{
  int i, k, l;

  i = j - d;
  k = i + d1;
  l = j - d2;

  if (traceback) fprintf(outf," (WH4.1) SLE, trace whx %d %d %d %d %d \n",  j, d-1, d1-1, d2,
			 whx[j][d-1][d1-1][d2]); 
  
  curr_tr->type1 = dpcL;              
  PushTraceknstack(dolist, AttachTracekn(curr_tr, i+1, j, k, l, dpcS, dpcS)); 
}

/* (WH4.2)__/ Singlet-right exterior / j dangles off i,j-1 */
int 
WH4_2(int *s, int len, int **icfg, int ****whx, int j, int d, int d1, int d2)
{ 
  int   WH4_2;
  int i, k, l;

  i = j - d;
  k = i + d1;
  l = j - d2;

  if (d2 > 0)
    WH4_2 = P6P + whx[j-1][d-1][d1][d2-1];
  else
    WH4_2 = -BIGINT;

  return WH4_2;
}
void
trace_WH4_2(FILE *outf, int ****whx, int j, int d, int d1, int d2, 
	    struct tracekn_s *curr_tr, struct traceknstack_s *dolist, int traceback)
{
  int i, k, l;

  i = j - d;
  k = i + d1;
  l = j - d2;

  if (traceback) fprintf(outf," (WH4.2) SRE, trace whx %d %d %d %d %d \n",  j-1, d-1, d1, d2-1,
			 whx[j-1][d-1][d1][d2-1]); 
  
  curr_tr->type1 = dpcL;              
  PushTraceknstack(dolist, AttachTracekn(curr_tr, i, j-1, k, l, dpcS, dpcS)); 
}

/* (WH4.3)__/ Singlet-left interior / l dangles off k,l+1 */ 
int 
WH4_3(int *s, int len, int **icfg, int ****whx, int j, int d, int d1, int d2)
{ 
  int   WH4_3;
  int i, k, l;

  i = j - d;
  k = i + d1;
  l = j - d2;

  if (d2 > 0)
   WH4_3 = P6P + whx[j][d][d1][d2-1];
  else
    WH4_3 = -BIGINT;

  return WH4_3;
}
void
trace_WH4_3(FILE *outf, int ****whx, int j, int d, int d1, int d2, 
	    struct tracekn_s *curr_tr, struct traceknstack_s *dolist, int traceback)
{
  int i, k, l;

  i = j - d;
  k = i + d1;
  l = j - d2;

  if (traceback) fprintf(outf," (WH4.3) SLI, trace whx %d %d %d %d %d \n", j, d, d1, d2-1,
			 whx[j][d][d1][d2-1]); 
  
  curr_tr->type2 = dpcL;               
  PushTraceknstack(dolist, AttachTracekn(curr_tr, i, j, k, l+1, dpcS, dpcS)); 
}

/* (WH4.4)__/ Singlet-right interior / k dangles off k-1,l */
int 
WH4_4(int *s, int len, int **icfg, int ****whx, int j, int d, int d1, int d2)
{ 
  int   WH4_4;
  int i, k, l;

  i = j - d;
  k = i + d1;
  l = j - d2;

  if (d1 > 0)
    WH4_4 = P6P + whx[j][d][d1-1][d2];
  else
    WH4_4 = -BIGINT;

  return WH4_4;
}
void
trace_WH4_4(FILE *outf, int ****whx, int j, int d, int d1, int d2, 
	    struct tracekn_s *curr_tr, struct traceknstack_s *dolist, int traceback)
{
  int i, k, l;

  i = j - d;
  k = i + d1;
  l = j - d2;

  if (traceback) fprintf(outf," (WH4.4) SRI, trace whx %d %d %d %d %d \n", j, d, d1-1, d2,
			 whx[j][d][d1-1][d2]); 
  
  curr_tr->type2 = dpcL;              
  PushTraceknstack(dolist, AttachTracekn(curr_tr, i, j, k-1, l, dpcS, dpcS));
}

/* (WH5)__/ STRUCTURES WITH TWO WBX (one). 
 * Connects i with k and l with j (k,d1) (j,d2) */
int 
WH5(int *s, int len, int **icfg, int **wbx, int j, int d, int d1, int d2)
{ 
  int   WH5;
  int i, k, l;

  i = j - d;
  k = i + d1;
  l = j - d2;

  WH5 = wbx[j][d2] + wbx[k][d1];

  return WH5;
}
void
trace_WH5(FILE *outf, int **wbx, int j, int d, int d1, int d2, 
	  struct tracekn_s *curr_tr, struct traceknstack_s *dolist, int traceback)
{
  int i, k, l;

  i = j - d;
  k = i + d1;
  l = j - d2;

  if (traceback) fprintf(outf," (WH5) trace wbx %d %d %d \n", j, d2, wbx[j][d2]); 
  if (traceback) fprintf(outf," (WH5) trace wbx %d %d %d \n", k, d1, wbx[k][d1]); 

  PushTraceknstack(dolist, AttachTracekn(curr_tr, l, j, l+(int)(d2/2), 
					 l+(int)(d2/2)+1, dpcB, dpcS)); 
  PushTraceknstack(dolist, AttachTracekn(curr_tr, i, k, i+(int)(d1/2), 
					 i+(int)(d1/2)+1, dpcB, dpcS));
}

/* STRUCTURES WITH ONE WBX AND ONE WHX: (WH6, WH8, WH10, WH12)
 * or         WITH ONE VX  AND ONE ZHX: (WH7, WH9)
 * or         WITH ONE VX  AND ONE YHX: (WH11, WH13)
 */

/* In the LEFT ARM (i,k).
 */

/* EXTERIOR LEFT */

/* (WH6)__/  1 WBX 1WHX
 * WBX connects i with i+mid (i+mid,mid)
 * WHX connects (i+mid+1,k) and (l,j) (j,d-mid-1,d1-mid-1,d2).
 */                                       
int 
WH6(int *s, int len, int **icfg, int **wbx, int ****whx, int j, int d, int d1, int d2, int mid)
{ 
  int   WH6;
  int i, k, l;

  i = j - d;
  k = i + d1;
  l = j - d2;
  
  WH6 = wbx[i+mid][mid] + whx[j][d-mid-1][d1-mid-1][d2];

  return WH6;
}
void
trace_WH6(FILE *outf, int **wbx, int ****whx, int j, int d, int d1, int d2, int mid,
	  struct tracekn_s *curr_tr, struct traceknstack_s *dolist, int traceback)
{
  int i, k, l;

  i = j - d;
  k = i + d1;
  l = j - d2;

  if (traceback) fprintf(outf," (WH6) EL, trace wbx %d %d %d \n", i+mid, mid, wbx[i+mid][mid]); 
  if (traceback) fprintf(outf," (WH6) EL, trace whx %d %d %d %d %d \n", j, d-mid-1, d1-mid-1, d2, 
			 whx[j][d-mid-1][d1-mid-1][d2]); 
  
  PushTraceknstack(dolist, AttachTracekn(curr_tr, i+mid+1, j, k, l, dpcS, dpcS));
  PushTraceknstack(dolist, AttachTracekn(curr_tr, i, i+mid, i+(int)(mid/2), 
					 i+(int)(mid/2)+1, dpcB, dpcS));
}

/* (WH7)__/  1 VX 1 ZHX
 * VX connects i with i+mid (i+mid,mid)
 * ZHX connects (i+mid+1,k) and (l,j) (j,d-mid-1,d1-mid-1,d2).
 * coaxial pairs (i, i+mid) (i+mid+1, j)
 * ....add stack[i+mid][i][i+mid+1][j]
 */   

/* (WH7.1)__/ no danglings / coaxial = stack + 1 */                                   
int 
WH7_1(int *s, int len, int **icfg, int **vx, int ****zhx, int j, int d, int d1, int d2, int mid)
{ 
  int   WH7_1;
  int i, k, l;

  i = j - d;
  k = i + d1;
  l = j - d2;
  
   WH7_1 = 2*P10P
     + wkn*(icfg[idxPS(s[i+mid],s[i])]
	    [idxPS(s[i+mid+1],s[j])] + INTSCALE)
     + vx[i+mid][mid] + zhx[j][d-mid-1][d1-mid-1][d2];
   
   return WH7_1;
}
void
trace_WH7_1(FILE *outf, int **vx, int ****zhx, int j, int d, int d1, int d2, int mid,
	    struct tracekn_s *curr_tr, struct traceknstack_s *dolist, int traceback)
{
  int i, k, l;

  i = j - d;
  k = i + d1;
  l = j - d2;

  if (traceback) fprintf(outf," (WH7.1) EL, trace vx  %d %d %d \n", i+mid, mid, vx[i+mid][mid]); 
  if (traceback) fprintf(outf," (WH7.1) EL, trace zhx %d %d %d %d %d \n", j, d-mid-1, d1-mid-1, d2, 
			 zhx[j][d-mid-1][d1-mid-1][d2]); 
  
  PushTraceknstack(dolist, AttachTracekn(curr_tr, i+mid+1, j, k, l, dpcP, dpcS));
  PushTraceknstack(dolist, AttachTracekn(curr_tr, i, i+mid, i+(int)(mid/2), 
					 i+(int)(mid/2)+1, dpcP, dpcS));
}

/* (WH7.2)__/ i dangles off i+1,i+mid */                                   
int 
WH7_2(int *s, int len, int **icfg, int **vx, int ****zhx, int j, int d, int d1, int d2, int mid)
{ 
  int   WH7_2;
  int i, k, l;

  i = j - d;
  k = i + d1;
  l = j - d2;
  
  if (mid > 0)
    WH7_2 = 2*P10P + P6P
      + wkn*icfg[idxL(s[i])][idxP(s[i+1],s[i+mid])]
      + wkn*icfg[idxPS(s[i+mid],s[i+1])][idxPS(s[i+mid+1],s[j])]
      + vx[i+mid][mid-1] + zhx[j][d-mid-1][d1-mid-1][d2];
  else
    WH7_2 = -BIGINT;

  return WH7_2;
}
void
trace_WH7_2(FILE *outf, int **vx, int ****zhx, int j, int d, int d1, int d2, int mid,
	    struct tracekn_s *curr_tr, struct traceknstack_s *dolist, int traceback)
{
  int i, k, l;

  i = j - d;
  k = i + d1;
  l = j - d2;

  if (traceback) fprintf(outf," (WH7.2) EL, dangle i, trace vx  %d %d %d \n", i+mid, mid-1, vx[i+mid][mid-1]); 
  if (traceback) fprintf(outf," (WH7.2) EL, trace zhx %d %d %d %d %d \n", j, d-mid-1, d1-mid-1, d2, 
			 zhx[j][d-mid-1][d1-mid-1][d2]); 
  
  PushTraceknstack(dolist, AttachTracekn(curr_tr, i+mid+1, j, k, l, dpcP, dpcS));
  PushTraceknstack(dolist, AttachTracekn(curr_tr, i+1, i+mid, i+1+(int)((mid-1)/2), 
					 i+(int)((mid-1)/2)+2, dpcP, dpcS));
}

/* (WH7.3)__/ j dangles off i+mid+1,j-1 */                                   
int 
WH7_3(int *s, int len, int **icfg, int **vx, int ****zhx, int j, int d, int d1, int d2, int mid)
{ 
  int   WH7_3;
  int i, k, l;

  i = j - d;
  k = i + d1;
  l = j - d2;
  
  if (d2 > 0)
    WH7_3 = 2*P10P + P6P
      + wkn*icfg[idxR(s[j])][idxP(s[i+mid+1],s[j-1])]
      + wkn*icfg[idxPS(s[i+mid],s[i])][idxPS(s[i+mid+1],s[j-1])]
      + vx[i+mid][mid] + zhx[j-1][d-mid-2][d1-mid-1][d2-1];
  else
    WH7_3 = -BIGINT;

  return WH7_3;
}
void
trace_WH7_3(FILE *outf, int **vx, int ****zhx, int j, int d, int d1, int d2, int mid,
	  struct tracekn_s *curr_tr, struct traceknstack_s *dolist, int traceback)
{
  int i, k, l;

  i = j - d;
  k = i + d1;
  l = j - d2;

  if (traceback) fprintf(outf," (WH7.3) EL, trace vx  %d %d %d \n", i+mid, mid, vx[i+mid][mid]); 
  if (traceback) fprintf(outf," (WH7.3) EL, dangle j, trace zhx %d %d %d %d %d \n", j-1, d-mid-2, d1-mid-1, d2-1, 
	  zhx[j-1][d-mid-2][d1-mid-1][d2-1]); 
  
  PushTraceknstack(dolist, AttachTracekn(curr_tr, i+mid+1, j-1, k, l, dpcP, dpcS));
  PushTraceknstack(dolist, AttachTracekn(curr_tr, i, i+mid, i+(int)(mid/2), 
					 i+(int)(mid/2)+1, dpcP, dpcS));
}

/* (WH7.4)__/ i dangles off i+1,i+mid / j dangles off i+mid+1,j-1*/                             int 
WH7_4(int *s, int len, int **icfg, int **vx, int ****zhx, int j, int d, int d1, int d2, int mid)
{ 
  int   WH7_4;
  int i, k, l;

  i = j - d;
  k = i + d1;
  l = j - d2;
  
  if (mid > 0 && d2 > 0)
    WH7_4 = 2*P10P + 2*P6P
      + wkn*icfg[idxL(s[i])][idxP(s[i+1],s[i+mid])]
      + wkn*icfg[idxR(s[j])][idxP(s[i+mid+1],s[j-1])]
      + wkn*icfg[idxPS(s[i+mid],s[i+1])][idxPS(s[i+mid+1],s[j-1])]
      + vx[i+mid][mid-1] + zhx[j-1][d-mid-2][d1-mid-1][d2-1];
  else
    WH7_4 = -BIGINT;

  return WH7_4;
}
void
trace_WH7_4(FILE *outf, int **vx, int ****zhx, int j, int d, int d1, int d2, int mid,
	  struct tracekn_s *curr_tr, struct traceknstack_s *dolist, int traceback)
{
  int i, k, l;

  i = j - d;
  k = i + d1;
  l = j - d2;

  if (traceback) fprintf(outf," (WH7.4) EL, dangle i, trace vx  %d %d %d \n", i+mid, mid-1, vx[i+mid][mid-1]); 
  if (traceback) fprintf(outf," (WH7.4) EL, dangle j, trace zhx %d %d %d %d %d \n", j-1, d-mid-2, d1-mid-1, d2-1, 
			 zhx[j-1][d-mid-2][d1-mid-1][d2-1]); 
  
  PushTraceknstack(dolist, AttachTracekn(curr_tr, i+mid+1, j-1, k, l, dpcP, dpcS));
  PushTraceknstack(dolist, AttachTracekn(curr_tr, i+1, i+mid, i+1+(int)((mid-1)/2), 
					 i+(int)((mid-1)/2)+2, dpcP, dpcS));
}

/* INTERIOR RIGHT */

/* (WH8)__/ 1 WBX 1 WHX
 * WBX connects k-mid with k (k,mid)
 * WHX connects (i,k-mid-1) and (l,j) (j,d,d1-mid-1,d2).
 */                                       
int 
WH8(int *s, int len, int **icfg, int **wbx, int ****whx, int j, int d, int d1, int d2, int mid)
{ 
  int   WH8;
  int i, k, l;

  i = j - d;
  k = i + d1;
  l = j - d2;
  
  WH8 = wbx[k][mid] + whx[j][d][d1-mid-1][d2];

  return WH8;
}
void
trace_WH8(FILE *outf, int **wbx, int ****whx, int j, int d, int d1, int d2, int mid,
	  struct tracekn_s *curr_tr, struct traceknstack_s *dolist, int traceback)
{
  int i, k, l;

  i = j - d;
  k = i + d1;
  l = j - d2;

  if (traceback) fprintf(outf," (WH8) IR, trace wbx  %d  %d %d \n", k, mid, wbx[k][mid]); 
  if (traceback) fprintf(outf," (WH8) IR, trace whx %d  %d  %d  %d %d \n", 
			 j, d, d1-mid-1, d2, whx[j][d][d1-mid-1][d2]); 
  
  PushTraceknstack(dolist, AttachTracekn(curr_tr, i, j, k-mid-1, l, dpcS, dpcS));
  PushTraceknstack(dolist, AttachTracekn(curr_tr, k-mid, k, k-mid+(int)(mid/2), 
					 k-mid+(int)(mid/2)+1, dpcB, dpcS));
}

/* (WH9)__/ 1 VX 1 YHX
 * VX  connects k-mid with k (k,mid)
 * YHX connects (i,k-mid-1) and (l,j) (j,d,d1-mid-1,d2).
 * coaxial pairs (k-mid-1,l) (k-mid,k)
 * ....add stack[k-mid-1][l][k-mid][k]
 */

/* (WH9.1)__/ no danglings / coaxial = stack + 1 */                                   

int 
WH9_1(int *s, int len, int **icfg, int **vx, int ****yhx, 
      int j, int d, int d1, int d2, int mid)
{ 
  int   WH9_1;
  int i, k, l;

  i = j - d;
  k = i + d1;
  l = j - d2;
  
  WH9_1 = 2*P10P
    + wkn*(icfg[idxPS(s[k-mid-1],s[l])]
	   [idxPS(s[k-mid],s[k])] + INTSCALE)
    + vx[k][mid] + yhx[j][d][d1-mid-1][d2];
  
  return WH9_1;
}
void
trace_WH9_1(FILE *outf, int **vx, int ****yhx, int j, int d, int d1, int d2, int mid,
	    struct tracekn_s *curr_tr, struct traceknstack_s *dolist, int traceback)
{
  int i, k, l;

  i = j - d;
  k = i + d1;
  l = j - d2;

  if (traceback) fprintf(outf," (WH9.1) IR, trace vx  %d  %d %d \n", k, mid, vx[k][mid]); 
  if (traceback) fprintf(outf," (WH9.1) IR, trace yhx %d  %d  %d  %d %d \n", 
			 j, d, d1-mid-1, d2, yhx[j][d][d1-mid-1][d2]); 
  
  PushTraceknstack(dolist, AttachTracekn(curr_tr, i, j, k-mid-1, l, dpcS, dpcP));
  PushTraceknstack(dolist, AttachTracekn(curr_tr, k-mid, k, k-mid+(int)(mid/2), 
					 k-mid+(int)(mid/2)+1, dpcP, dpcS));
}

/* (WH9.2)__/ k dangles off k-mid,k-1 */                                   
int 
WH9_2(int *s, int len, int **icfg, int **vx, int ****yhx, 
      int j, int d, int d1, int d2, int mid)
{ 
  int   WH9_2;
  int i, k, l;

  i = j - d;
  k = i + d1;
  l = j - d2;
  
  if (mid > 0)
    WH9_2 = 2*P10P + P6P
      + wkn*icfg[idxR(s[k])][idxP(s[k-mid],s[k-1])]
      + wkn*icfg[idxPS(s[k-mid-1],s[l])][idxPS(s[k-mid],s[k-1])]
      + vx[k-1][mid-1] + yhx[j][d][d1-mid-1][d2];
  else
    WH9_2 = -BIGINT;

  return WH9_2;
}
void
trace_WH9_2(FILE *outf, int **vx, int ****yhx, int j, int d, int d1, int d2, int mid,
	    struct tracekn_s *curr_tr, struct traceknstack_s *dolist, int traceback)
{
  int i, k, l;

  i = j - d;
  k = i + d1;
  l = j - d2;

  if (traceback) fprintf(outf," (WH9.2) IR, dangle k, trace vx %d  %d %d \n", k-1, mid-1, vx[k-1][mid-1]); 
  if (traceback) fprintf(outf," (WH9.2) IR, trace yhx %d  %d  %d  %d %d \n", 
			 j, d, d1-mid-1, d2, yhx[j][d][d1-mid-1][d2]); 
  
  PushTraceknstack(dolist,  AttachTracekn(curr_tr, i, j, k-mid-1, l, dpcS, dpcP));
  PushTraceknstack(dolist, AttachTracekn(curr_tr, k-mid, k-1, k-mid+(int)((mid-1)/2), 
					 k-mid+(int)((mid-1)/2)+ 1, dpcP, dpcS));
}

/* (WH9.3)__/ l dangles off l+1,k-mid-1 */                                   
int 
WH9_3(int *s, int len, int **icfg, int **vx, int ****yhx, 
      int j, int d, int d1, int d2, int mid)
{ 
  int   WH9_3;
  int i, k, l;

  i = j - d;
  k = i + d1;
  l = j - d2;
  
  if (d2 > 0)
    WH9_3 = 2*P10P + P6P
      + wkn*icfg[idxL(s[l])][idxP(s[l+1],s[k-mid-1])]
      + wkn*icfg[idxPS(s[k-mid-1],s[l+1])][idxPS(s[k-mid],s[k])]
      + vx[k][mid] + yhx[j][d][d1-mid-1][d2-1];
  else
    WH9_3 = -BIGINT;

  return WH9_3;
}
void
trace_WH9_3(FILE *outf, int **vx, int ****yhx, int j, int d, int d1, int d2, int mid,
	    struct tracekn_s *curr_tr, struct traceknstack_s *dolist, int traceback)
{
  int i, k, l;

  i = j - d;
  k = i + d1;
  l = j - d2;

  if (traceback) fprintf(outf," (WH9.3) IR, trace vx  %d  %d %d \n", k, mid, vx[k][mid]); 
  if (traceback) fprintf(outf," (WH9.3) IR, dangle l, trace yhx %d  %d  %d  %d %d \n", 
			 j, d, d1-mid-1, d2-1, yhx[j][d][d1-mid-1][d2-1]); 
  
  PushTraceknstack(dolist, AttachTracekn(curr_tr, i, j, k-mid-1, l+1, dpcS, dpcP));
  PushTraceknstack(dolist, AttachTracekn(curr_tr, k-mid, k, k-mid+(int)(mid/2), 
					 k-mid+(int)(mid/2)+1, dpcP, dpcS));
}

/* (WH9.4)__/ k dangles off k-mid,k-1 / l dangles offl+1,k-mid-1*/                  
int 
WH9_4(int *s, int len, int **icfg, int **vx, int ****yhx, 
      int j, int d, int d1, int d2, int mid)
{ 
  int   WH9_4;
  int i, k, l;

  i = j - d;
  k = i + d1;
  l = j - d2;
  
  if (mid > 0 && d2 > 0)
    WH9_4 = 2*P10P + 2*P6P
      + wkn*icfg[idxR(s[k])][idxP(s[k-mid],s[k-1])]
      + wkn*icfg[idxL(s[l])][idxP(s[l+1],s[k-mid-1])]
      + wkn*icfg[idxPS(s[k-mid-1],s[l+1])][idxPS(s[k-mid],s[k-1])]
      + vx[k-1][mid-1] + yhx[j][d][d1-mid-1][d2-1];
  else
    WH9_4 = -BIGINT;

  return WH9_4;
}
void
trace_WH9_4(FILE *outf, int **vx, int ****yhx, int j, int d, int d1, int d2, int mid,
	    struct tracekn_s *curr_tr, struct traceknstack_s *dolist, int traceback)
{
  int i, k, l;

  i = j - d;
  k = i + d1;
  l = j - d2;

  if (traceback) fprintf(outf," (WH9.4) IR, dangle k, trace vx %d  %d %d \n", k-1, mid-1, vx[k-1][mid-1]); 
  if (traceback) fprintf(outf," (WH9.4) IR, dangle l, trace yhx %d  %d  %d  %d %d \n", 
			 j, d, d1-mid-1, d2-1, yhx[j][d][d1-mid-1][d2-1]); 
  
  PushTraceknstack(dolist, AttachTracekn(curr_tr, i, j, k-mid-1, l+1, dpcS, dpcP));
  PushTraceknstack(dolist, AttachTracekn(curr_tr, k-mid, k-1, k-mid+(int)((mid-1)/2), 
					 k-mid+(int)((mid-1)/2)+ 1, dpcP, dpcS));
}

/* In the RIGHT ARM (l,j).
 */      
/* EXTERIOR RIGHT */

/* (WH10)__/ 1 WBX 1 WHX
 * WBX connects j-mid with j  (j,mid)
 * WHX connects (i,k) and (l,j-mid-1) 
 *               (j-mid-1, d-mid-1, d1, d2-mid-1).
 */                                       
int 
WH10(int *s, int len, int **icfg, int **wbx, int ****whx, int j, int d, int d1, int d2, int mid)
{ 
  int   WH10;
  int i, k, l;

  i = j - d;
  k = i + d1;
  l = j - d2;
  
  WH10 = wbx[j][mid] + whx[j-mid-1][d-mid-1][d1][d2-mid-1];

  return WH10;
}
void
trace_WH10(FILE *outf, int **wbx, int ****whx, int j, int d, int d1, int d2, int mid,
	   struct tracekn_s *curr_tr, struct traceknstack_s *dolist, int traceback)
{
  int i, k, l;

  i = j - d;
  k = i + d1;
  l = j - d2;

  if (traceback) fprintf(outf," (WH10) ER, trace wbx %d %d %d \n", j, mid, wbx[j][mid]); 
  if (traceback) fprintf(outf," (WH10) ER, trace whx %d %d %d %d %d \n", j-mid-1, d-mid-1, d1, d2-mid-1,
			 whx[j-mid-1][d-mid-1][d1][d2-mid-1]); 
  
  PushTraceknstack(dolist, AttachTracekn(curr_tr, j-mid, j, j-mid+(int)(mid/2), 
					 j-mid+(int)(mid/2)+1, dpcB, dpcS));
  PushTraceknstack(dolist, AttachTracekn(curr_tr, i, j-mid-1, k, l, dpcS, dpcS));
}

/* (WH11)__/ 1 VX 1 ZHX
 * VX  connects j-mid with j  (j,mid)
 * ZHX connects (i,k) and (l,j-mid-1) 
 *              (j-mid-1, d-mid-1, d1, d2-mid-1).
 * contiguos pairs (i,j-mid-1) (j-mid,j)
 * ....add stack[j-mid-1][i][j-mid][j]
 */                                       

/* (WH11.1)__/ no danglings / coaxial = stack + 1 */                                   
int 
WH11_1(int *s, int len, int **icfg, int **vx, int ****zhx, int j, int d, int d1, int d2, int mid)
{ 
  int  WH11_1;
  int i, k, l;

  i = j - d;
  k = i + d1;
  l = j - d2;
  
  WH11_1 = 2*P10P
    + wkn*(icfg[idxPS(s[j-mid-1],s[i])]
	   [idxPS(s[j-mid],s[j])] + INTSCALE)
    + vx[j][mid] + zhx[j-mid-1][d-mid-1][d1][d2-mid-1];
  
  return WH11_1;
}
void
trace_WH11_1(FILE *outf, int **vx, int ****zhx, int j, int d, int d1, int d2, int mid,
	     struct tracekn_s *curr_tr, struct traceknstack_s *dolist, int traceback)
{
  int i, k, l;

  i = j - d;
  k = i + d1;
  l = j - d2;

  if (traceback) fprintf(outf," (WH11) ER, trace vx  %d %d %d \n", j, mid, vx[j][mid]); 
  if (traceback) fprintf(outf," (WH11) ER, trace zhx %d %d %d %d %d \n", j-mid-1, d-mid-1, d1, d2-mid-1,
			 zhx[j-mid-1][d-mid-1][d1][d2-mid-1]); 
  
  PushTraceknstack(dolist, AttachTracekn(curr_tr, j-mid, j, j-mid+(int)(mid/2), 
					 j-mid+(int)(mid/2)+1, dpcP, dpcS));
  PushTraceknstack(dolist, AttachTracekn(curr_tr, i, j-mid-1, k, l, dpcP, dpcS));
}

/* (WH11.2)__/ i dangles off i+1,j-mid-1 */                                   
int 
WH11_2(int *s, int len, int **icfg, int **vx, int ****zhx, int j, int d, int d1, int d2, int mid)
{ 
  int  WH11_2;
  int i, k, l;

  i = j - d;
  k = i + d1;
  l = j - d2;
  
  if (d1 > 0)
    WH11_2 = 2*P10P + P6P
      + wkn*icfg[idxL(s[i])][idxP(s[i+1],s[j-mid-1])]
      + wkn*icfg[idxPS(s[j-mid-1],s[i+1])][idxPS(s[j-mid],s[j])]
      + vx[j][mid] + zhx[j-mid-1][d-mid-2][d1-1][d2-mid-1];
  else
    WH11_2 = -BIGINT;

  return WH11_2;
}
void
trace_WH11_2(FILE *outf, int **vx, int ****zhx, int j, int d, int d1, int d2, int mid,
	     struct tracekn_s *curr_tr, struct traceknstack_s *dolist, int traceback)
{
  int i, k, l;

  i = j - d;
  k = i + d1;
  l = j - d2;

  if (traceback) fprintf(outf," (WH11.2) ER, trace vx  %d %d %d \n", j, mid, vx[j][mid]); 
  if (traceback) fprintf(outf," (WH11.2) ER, dangle i, trace zhx %d %d %d %d %d \n", 
			 j-mid-1, d-mid-2, d1-1, d2-mid-1,
			 zhx[j-mid-1][d-mid-2][d1-1][d2-mid-1]); 
  
  PushTraceknstack(dolist, AttachTracekn(curr_tr, j-mid, j, j-mid+(int)(mid/2), 
					 j-mid+(int)(mid/2)+1, dpcP, dpcS));
  PushTraceknstack(dolist, AttachTracekn(curr_tr, i+1, j-mid-1, k, l, dpcP, dpcS));
}

/* (WH11.3)__/ j dangles off j-mid,j-1 */                                   
int 
WH11_3(int *s, int len, int **icfg, int **vx, int ****zhx, int j, int d, int d1, int d2, int mid)
{ 
  int   WH11_3;
  int i, k, l;

  i = j - d;
  k = i + d1;
  l = j - d2;
  
  if (mid > 0)
    WH11_3 = 2*P10P + P6P
      + wkn*icfg[idxR(s[j])][idxP(s[j-mid],s[j-1])]
      + wkn*icfg[idxPS(s[j-mid-1],s[i])][idxPS(s[j-mid],s[j-1])]
      + vx[j-1][mid-1] + zhx[j-mid-1][d-mid-1][d1][d2-mid-1];
  else
    WH11_3 = -BIGINT;

  return WH11_3;
}
void
trace_WH11_3(FILE *outf, int **vx, int ****zhx, int j, int d, int d1, int d2, int mid,
	     struct tracekn_s *curr_tr, struct traceknstack_s *dolist, int traceback)
{
  int i, k, l;

  i = j - d;
  k = i + d1;
  l = j - d2;

  if (traceback) fprintf(outf," (WH11.3) ER, dangle j, trace vx  %d %d %d \n", j-1, mid-1, vx[j-1][mid-1]); 
  if (traceback) fprintf(outf," (WH11.3) ER, trace zhx %d %d %d %d %d \n", j-mid-1, d-mid-1, d1, d2-mid-1,
			 zhx[j-mid-1][d-mid-1][d1][d2-mid-1]); 
  
  PushTraceknstack(dolist, AttachTracekn(curr_tr, j-mid, j-1, j-mid+(int)((mid-1)/2), 
					 j-mid+(int)((mid-1)/2)+1, dpcP, dpcS));
  PushTraceknstack(dolist, AttachTracekn(curr_tr, i, j-mid-1, k, l, dpcP, dpcS));
}

/* (WH11.4)__/ i dangles off i+1,j-mid-1 / j dangles off j-mid,j-1*/                           int 
WH11_4(int *s, int len, int **icfg, int **vx, int ****zhx, int j, int d, int d1, int d2, int mid)
{ 
  int   WH11_4;
  int i, k, l;

  i = j - d;
  k = i + d1;
  l = j - d2;
  
  if (d1 > 0 && mid > 0)
    WH11_4 = 2*P10P + 2*P6P
      + wkn*icfg[idxL(s[i])][idxP(s[i+1],s[j-mid-1])]
      + wkn*icfg[idxR(s[j])][idxP(s[j-mid],s[j-1])]
      + wkn*icfg[idxPS(s[j-mid-1],s[i+1])][idxPS(s[j-mid],s[j-1])]
      + vx[j-1][mid-1] + zhx[j-mid-1][d-mid-2][d1-1][d2-mid-1];
  else
    WH11_4 = -BIGINT;

  return WH11_4;
}
void
trace_WH11_4(FILE *outf, int **vx, int ****zhx, int j, int d, int d1, int d2, int mid,
	     struct tracekn_s *curr_tr, struct traceknstack_s *dolist, int traceback)
{
  int i, k, l;

  i = j - d;
  k = i + d1;
  l = j - d2;

  if (traceback) fprintf(outf," (WH11.4) ER, dangel j, trace vx  %d %d %d \n", j-1, mid-1, vx[j-1][mid-1]); 
  if (traceback) fprintf(outf," (WH11.4) ER, dangle i, trace zhx %d %d %d %d %d \n", 
			 j-mid-1, d-mid-2, d1-1, d2-mid-1,
			 zhx[j-mid-1][d-mid-2][d1-1][d2-mid-1]); 
  
  PushTraceknstack(dolist, AttachTracekn(curr_tr, j-mid, j-1, j-mid+(int)((mid-1)/2), 
					 j-mid+(int)((mid-1)/2)+1, dpcP, dpcS));
  PushTraceknstack(dolist, AttachTracekn(curr_tr, i+1, j-mid-1, k, l, dpcP, dpcS));
}

/* INTERIOR LEFT */

/* (WH12)__/ 1 WBX 1 WHX
 * WBX connects l with l+mid (l+mid,mid)
 * WHX connects (i,k) and (l+mid+1,j) (j, d, d1, d2-mid-1).
 */                                       
int 
WH12(int *s, int len, int **icfg, int **wbx, int ****whx, int j, int d, int d1, int d2, int mid)
{ 
  int    WH12;
  int i, k, l;

  i = j - d;
  k = i + d1;
  l = j - d2;
  
  WH12 = wbx[l+mid][mid] + whx[j][d][d1][d2-mid-1];

  return WH12;
}
void
trace_WH12(FILE *outf, int **wbx, int ****whx, int j, int d, int d1, int d2, int mid,
	  struct tracekn_s *curr_tr, struct traceknstack_s *dolist, int traceback)
{
  int i, k, l;

  i = j - d;
  k = i + d1;
  l = j - d2;

  if (traceback) fprintf(outf," (WH12) IL, trace wbx %d %d %d \n", l+mid, mid, wbx[l+mid][mid]); 
  if (traceback) fprintf(outf," (WH12) IL, trace whx %d %d %d %d %d \n", 
			 j, d, d1, d2-mid-1, whx[j][d][d1][d2-mid-1]); 
  
  PushTraceknstack(dolist, AttachTracekn(curr_tr, i, j, k, l+mid+1, dpcS, dpcS));
  PushTraceknstack(dolist, AttachTracekn(curr_tr, l, l+mid, l + (int)(mid/2), 
					 l+(int)(mid/2)+1, dpcB, dpcS));
}

/* (WH13)__/ 1 VX 1 YHX
 * VX  connects l with l+mid (l+mid,mid)
 * YHX connects (i,k) and (l+mid+1,j) (j, d, d1, d2-mid-1).
 * coaxial pairs (k,l+mid+1) (l,l+mid)
 * ....add stack[l+mid][l][l+mid+1][k]
 */ 

/* (WH13.1)__/ no danglings / coaxial = stack + 1 */                                   
int 
WH13_1(int *s, int len, int **icfg, int **vx, int ****yhx, int j, int d, int d1, int d2, int mid)
{ 
  int  WH13_1;
  int i, k, l;

  i = j - d;
  k = i + d1;
  l = j - d2;
  
    WH13_1 = 2*P10P
      + wkn*(icfg[idxPS(s[l+mid],s[l])]
	     [idxPS(s[l+mid+1],s[k])] + INTSCALE)
      + vx[l+mid][mid] + yhx[j][d][d1][d2-mid-1];
    
    return WH13_1;
}
void
trace_WH13_1(FILE *outf, int **vx, int ****yhx, int j, int d, int d1, int d2, int mid,
	     struct tracekn_s *curr_tr, struct traceknstack_s *dolist, int traceback)
{
  int i, k, l;

  i = j - d;
  k = i + d1;
  l = j - d2;

  if (traceback) fprintf(outf," (WH13.1) IL, trace vx  %d %d %d \n", l+mid, mid, vx[l+mid][mid]); 
  if (traceback) fprintf(outf," (WH13.1) IL, trace yhx %d %d %d %d %d \n", 
			 j, d, d1, d2-mid-1, yhx[j][d][d1][d2-mid-1]); 
  
  PushTraceknstack(dolist, AttachTracekn(curr_tr, i, j, k, l+mid+1, dpcS, dpcP));
  PushTraceknstack(dolist, AttachTracekn(curr_tr, l, l+mid, l + (int)(mid/2), 
					 l+(int)(mid/2)+1, dpcP, dpcS));
}

/* (WH13.2)__/ k dangles off l+mid+1,k-1 */                                   
int 
WH13_2(int *s, int len, int **icfg, int **vx, int ****yhx, int j, int d, int d1, int d2, int mid)
{ 
  int  WH13_2;
  int i, k, l;

  i = j - d;
  k = i + d1;
  l = j - d2;
  
  if (d1 > 0)
    WH13_2 = 2*P10P + P6P
      + wkn*icfg[idxL(s[k])][idxP(s[l+mid+1],s[k-1])]
      + wkn*icfg[idxPS(s[l+mid],s[l])][idxPS(s[l+mid+1],s[k-1])]
      + vx[l+mid][mid] + yhx[j][d][d1-1][d2-mid-1];
  else
    WH13_2 = -BIGINT;

  return WH13_2;
}
void
trace_WH13_2(FILE *outf, int **vx, int ****yhx, int j, int d, int d1, int d2, int mid,
	     struct tracekn_s *curr_tr, struct traceknstack_s *dolist, int traceback)
{
  int i, k, l;

  i = j - d;
  k = i + d1;
  l = j - d2;

  if (traceback) fprintf(outf," (WH13.2) IL, trace vx  %d %d %d \n", l+mid, mid, vx[l+mid][mid]); 
  if (traceback) fprintf(outf," (WH13.2) IL, dangle k, trace yhx %d %d %d %d %d \n", 
			 j, d, d1-1, d2-mid-1, yhx[j][d][d1-1][d2-mid-1]); 
  
  PushTraceknstack(dolist, AttachTracekn(curr_tr, i, j, k-1, l+mid+1, dpcS, dpcP));
  PushTraceknstack(dolist, AttachTracekn(curr_tr, l, l+mid, l + (int)(mid/2), 
					 l+(int)(mid/2)+1, dpcP, dpcS));
}

/* (WH13.3)__/ l dangles off l+1,l+mid */                                   
int 
WH13_3(int *s, int len, int **icfg, int **vx, int ****yhx, int j, int d, int d1, int d2, int mid)
{ 
  int  WH13_3;
  int i, k, l;

  i = j - d;
  k = i + d1;
  l = j - d2;
  
  if (mid > 0)
    WH13_3 = 2*P10P + P6P
      + wkn*icfg[idxR(s[l])][idxP(s[l+1],s[l+mid])]
      + wkn*icfg[idxPS(s[l+mid],s[l+1])][idxPS(s[l+mid+1],s[k])]
      + vx[l+mid][mid-1] + yhx[j][d][d1][d2-mid-1];
  else
    WH13_3 = -BIGINT;

  return WH13_3;
}
void
trace_WH13_3(FILE *outf, int **vx, int ****yhx, int j, int d, int d1, int d2, int mid,
	     struct tracekn_s *curr_tr, struct traceknstack_s *dolist, int traceback)
{
  int i, k, l;

  i = j - d;
  k = i + d1;
  l = j - d2;

  if (traceback) fprintf(outf," (WH13.3) IL, dangle l, trace vx  %d %d %d \n", 
			 l+mid, mid-1, vx[l+mid][mid-1]); 
  if (traceback) fprintf(outf," (WH13.3) IL, trace yhx %d %d %d %d %d \n", 
			 j, d, d1, d2-mid-1, yhx[j][d][d1][d2-mid-1]); 
  
  PushTraceknstack(dolist, AttachTracekn(curr_tr, i, j, k, l+mid+1, dpcS, dpcP));
  PushTraceknstack(dolist, AttachTracekn(curr_tr, l+1, l+mid, l+1+(int)((mid-1)/2), 
					 l+(int)((mid-1)/2)+2, dpcP, dpcS));
}

/* (WH13.4)__/ k dangles off l+mid+1,k-1 / l dangles off l+1,l+mid*/                           
int 
WH13_4(int *s, int len, int **icfg, int **vx, int ****yhx, int j, int d, int d1, int d2, int mid)
{ 
  int   WH13_4;
  int i, k, l;

  i = j - d;
  k = i + d1;
  l = j - d2;
  
  if (d1 > 0 && mid > 0)
    WH13_4 = 2*P10P + 2*P6P
      + wkn*icfg[idxL(s[k])][idxP(s[l+mid+1],s[k-1])]
      + wkn*icfg[idxR(s[l])][idxP(s[l+1],s[l+mid])]
      + wkn*icfg[idxPS(s[l+mid],s[l+1])][idxPS(s[l+mid+1],s[k-1])]
      + vx[l+mid][mid-1] + yhx[j][d][d1-1][d2-mid-1];
  else
    WH13_4 = -BIGINT;

  return WH13_4;
}
void
trace_WH13_4(FILE *outf, int **vx, int ****yhx, int j, int d, int d1, int d2, int mid,
	     struct tracekn_s *curr_tr, struct traceknstack_s *dolist, int traceback)
{
  int i, k, l;

  i = j - d;
  k = i + d1;
  l = j - d2;

  if (traceback) fprintf(outf," (WH13.4) IL, dangle l, trace vx  %d %d %d \n", 
			 l+mid, mid-1, vx[l+mid][mid-1]); 
  if (traceback) fprintf(outf," IL, dangle k, trace yhx %d %d %d %d %d \n", 
			 j, d, d1-1, d2-mid-1, yhx[j][d][d1-1][d2-mid-1]); 
  
  PushTraceknstack(dolist, AttachTracekn(curr_tr, i, j, k-1, l+mid+1, dpcS, dpcP));
  PushTraceknstack(dolist, AttachTracekn(curr_tr, l+1, l+mid, l+1+(int)((mid-1)/2), 
					 l+(int)((mid-1)/2)+2, dpcP, dpcS));
}

/* 
 * WITH TWO HOLE STRUCTURES. (WH14-WH24)
 */                    
/* (WH14)__/  STRUCTURE WITH ONE YHX AND ONE ZHX. (only one pair)
 *
 * yhx connects i to j and k-mid1 to l+mid2. 
 *              (j,d,d1-mid1,d2-mid2)
 *
 * zhx connect k-mid1 to l+mid2 and k to l. 
 *             (l+mid2,d-d1-d2+mid1+mid2,mid1,mid2)
 */ 

int 
WH14(int *s, int len, int **icfg, int ****zhx, int ****yhx, 
     int j, int d, int d1, int d2, int mid1, int mid2)
{ 
  int   WH14;
  int i, k, l;

  i = j - d;
  k = i + d1;
  l = j - d2;
  
    WH14 = yhx[j][d][d1-mid1][d2-mid2]
      + zhx[l+mid2][d-d1-d2+mid1+mid2][mid1][mid2];

  return WH14;
}
void
trace_WH14(FILE *outf, int ****zhx, int ****yhx, int j, int d, int d1, int d2, int mid1, int mid2,
	   struct tracekn_s *curr_tr, struct traceknstack_s *dolist, int traceback)
{
  int i, k, l;

  i = j - d;
  k = i + d1;
  l = j - d2;

  if (traceback) fprintf(outf," (WH14) trace YHX %d %d %d %d %d\n", 
			 j, d, d1-mid1, d2-mid2,
			 yhx[j][d][d1-mid1][d2-mid2]); 
  if (traceback) fprintf(outf," (WH14) trace ZHX %d %d %d %d %d\n", 
			 l+mid2, d-d1-d2+mid1+mid2, mid1, mid2,
			 zhx[l+mid2][d-d1-d2+mid1+mid2][mid1][mid2]); 
  
  PushTraceknstack(dolist, AttachTracekn(curr_tr, i, j, k-mid1, l+mid2, dpcS, dpcP));
  PushTraceknstack(dolist, AttachTracekn(curr_tr, k-mid1, l+mid2, k, l, dpcP, dpcS));
}

/* (WH15)__/  2 WHX (inclusive bifurcation)
 * One connects     (i,i+mid1)   to (j-mid2, j)  
 *                  (j,d,mid1,mid2)
 * Another connects (i+mid1+1,k) to (l,j-mid2-1) 
 *                   (j-mid2-1, d-mid1-mid2-2, d1-mid1-1, d2-mid2-1)
 */
int 
WH15(int *s, int len, int **icfg, int ****whx, 
     int j, int d, int d1, int d2, int mid1, int mid2)
{ 
  int   WH15;
  int i, k, l;
  
  i = j - d;
  k = i + d1;
  l = j - d2;

  if (d1-mid1 > 0 && d2-mid2 > 0)
    WH15 = P5P 
      + whx[j][d][mid1][mid2] 
      + whx[j-mid2-1][d-mid1-mid2-2][d1-mid1-1][d2-mid2-1];
  else
    WH15 = -BIGINT;
    
  return WH15;
}
void
trace_WH15(FILE *outf, int ****whx, int j, int d, int d1, int d2, int mid1, int mid2,
	   struct tracekn_s *curr_tr, struct traceknstack_s *dolist, int traceback)
{
  int i, k, l;

  i = j - d;
  k = i + d1;
  l = j - d2;

  if (traceback) fprintf(outf," (WH15) trace WHX %d %d %d %d %d\n", j, d, mid1, mid2, 
			 whx[j][d][mid1][mid2]); 
  if (traceback) fprintf(outf," (WH15) trace WHX %d %d %d %d %d\n", 
			 j-mid2-1, d-mid1-mid2-2, d1-mid1-1, d2-mid2-1,
			 whx[j-mid2-1][d-mid1-mid2-2][d1-mid1-1][d2-mid2-1]); 
  
  PushTraceknstack(dolist, AttachTracekn(curr_tr, i, j, i+mid1, j-mid2, dpcS, dpcS));
  PushTraceknstack(dolist, AttachTracekn(curr_tr, i+mid1+1, j-mid2-1, k, l, dpcS, dpcS));
}

/* (WH16)__/  2 WHX (croosed bifurcation)
 * One connects     (i,i+mid1)   to (l, l+mid2)  
 *                  (l+mid2,d-d2+mid2,mid1,mid2)
 * Another connects (i+mid1+2,k) to (l+mid2+3,j) 
 *                   (j,d-mid1-1,d1-mi1-1,d2-mid2-2)
 */
int 
WH16(int *s, int len, int **icfg, int ****whx, 
     int j, int d, int d1, int d2, int mid1, int mid2)
{ 
  int   WH16;
  int i, k, l;
  
  i = j - d;
  k = i + d1;
  l = j - d2;

  if (d1-mid1 > 1 && d2-mid2 > 2)
    WH16 = P12 
      + whx[l+mid2][d-d2+mid2][mid1][mid2] 
      + whx[j][d-mid1-2][d1-mid1-2][d2-mid2-3];
  else
    WH16 = -BIGINT;
    
  return WH16;
}
void
trace_WH16(FILE *outf, int ****whx, int j, int d, int d1, int d2, int mid1, int mid2,
	   struct tracekn_s *curr_tr, struct traceknstack_s *dolist, int traceback)
{
  int i, k, l;

  i = j - d;
  k = i + d1;
  l = j - d2;

  if (traceback) fprintf(outf," (WH16) trace WHX %d %d %d %d %d\n", 
			 l+mid2, d-d2+mid2, mid1, mid2,
			 whx[l+mid2][d-d2+mid2][mid1][mid2]); 
  if (traceback) fprintf(outf," (WH16) trace WHX %d %d %d %d %d\n", 
			 j, d-mid1-2, d1-mid1-2, d2-mid2-3,
			 whx[j][d-mid1-2][d1-mid1-2][d2-mid2-3]); 
  
  PushTraceknstack(dolist, AttachTracekn(curr_tr, i+mid1+2, j, k, l+mid2+3, dpcS, dpcS));
  PushTraceknstack(dolist, AttachTracekn(curr_tr, i, l+mid2, i+mid1, l, dpcS, dpcS));
}

/* Structures with pseudoknot on one of the arms (WH17-WH22)
 * One     connects (i,k)      to (l+mid1+1,j-mid2-1) 
 *                  (j-mid2-1,d-mid2-1,d1,d2-mid1-mid2-2)
 * Another connects (l,l+mid1) to (j-mid2,j)          
 *                  (j,d2,mid1,mid2) 
 */
/* (WH17)__/ 2 WHX */
int 
WH17(int *s, int len, int **icfg, int ****whx, 
     int j, int d, int d1, int d2, int mid1, int mid2)
{ 
  int    WH17;
  int i, k, l;
  
  i = j - d;
  k = i + d1;
  l = j - d2;
  
  WH17 = P12 
    + whx[j-mid2-1][d-mid2-1][d1][d2-mid1-mid2-2] 
    + whx[j][d2][mid1][mid2];
  
  return WH17;
}
void
trace_WH17(FILE *outf, int ****whx, int j, int d, int d1, int d2, int mid1, int mid2,
	   struct tracekn_s *curr_tr, struct traceknstack_s *dolist, int traceback)
{
  int i, k, l;

  i = j - d;
  k = i + d1;
  l = j - d2;

  if (traceback) fprintf(outf," (WH17) trace WHX %d %d %d %d %d \n", 
			 j-mid2-1, d-mid2-1, d1, d2-mid1-mid2-2, 
			 whx[j-mid2-1][d-mid2-1][d1][d2-mid1-mid2-2]); 
  if (traceback) fprintf(outf," (WH17) trace WHX %d %d %d %d %d \n", j, d2, mid1, mid2,
			 whx[j][d2][mid1][mid2]); 
  
  PushTraceknstack(dolist, AttachTracekn(curr_tr, 
					 l, j, l+mid1, j-mid2, dpcS, dpcS));
  PushTraceknstack(dolist, AttachTracekn(curr_tr, 
					 i, j-mid2-1, k, l+mid1+1, dpcS, dpcS));
}

/* (WH18)__/  1 ZHX 1 YHX
 * coaxial pairs (i,j-mid2-1) (l+mid1,j-mid2)
 * ....add stack[j-mid2-1][i][j-mid2][l+mid1]
 */

/* (WH18.1)__/ no danglings / coaxial = stack - 1 */                                   
int 
WH18_1(int *s, int len, int **icfg, int ****zhx, int ****yhx, 
       int j, int d, int d1, int d2, int mid1, int mid2)
{ 
  int   WH18_1;
  int i, k, l;
  
  i = j - d;
  k = i + d1;
  l = j - d2;
  
  WH18_1 = P12 + 2*P10P 
    + wkn*(icfg[idxPS(s[j-mid2-1],s[i])]
	   [idxPS(s[j-mid2],s[l+mid1])] + INTSCALE)
    + zhx[j-mid2-1][d-mid2-1][d1][d2-mid1-mid2-2] 
    + yhx[j][d2][mid1][mid2];
  
  return WH18_1;
}
void
trace_WH18_1(FILE *outf, int ****zhx, int ****yhx, 
	     int j, int d, int d1, int d2, int mid1, int mid2,
	     struct tracekn_s *curr_tr, struct traceknstack_s *dolist, int traceback)
{
  int i, k, l;

  i = j - d;
  k = i + d1;
  l = j - d2;

  if (traceback) fprintf(outf," (WH18.1) trace ZHX %d %d %d %d %d \n", 
			 j-mid2-1, d-mid2-1, d1, d2-mid1-mid2-2, 
			 zhx[j-mid2-1][d-mid2-1][d1][d2-mid1-mid2-2]); 
  if (traceback) fprintf(outf," (WH18.1) trace YHX %d %d %d %d %d \n", j, d2, mid1, mid2,
			 yhx[j][d2][mid1][mid2]); 
  
  PushTraceknstack(dolist, AttachTracekn(curr_tr, 
					 i, j-mid2-1, k, l+mid1+1, dpcP, dpcS));
  PushTraceknstack(dolist, AttachTracekn(curr_tr, 
					 l, j, l+mid1, j-mid2, dpcS, dpcP));
}

/* (WH18.2)__/ i dangles off i+1,j-mid2-1 */
int 
WH18_2(int *s, int len, int **icfg, int ****zhx, int ****yhx, 
       int j, int d, int d1, int d2, int mid1, int mid2)
{ 
  int   WH18_2;
  int i, k, l;
  
  i = j - d;
  k = i + d1;
  l = j - d2;
  
  if (d1 > 0)
    WH18_2 = P12 + 2*P10P + P6P
      + wkn*icfg[idxL(s[i])][idxP(s[i+1],s[j-mid2-1])]
      + wkn*icfg[idxPS(s[j-mid2-1],s[i+1])][idxPS(s[j-mid2],s[l+mid1])]
      + zhx[j-mid2-1][d-mid2-2][d1-1][d2-mid1-mid2-2] 
      + yhx[j][d2][mid1][mid2];
  else
    WH18_2 = -BIGINT;

  return WH18_2;
}
void
trace_WH18_2(FILE *outf, int ****zhx, int ****yhx, 
	     int j, int d, int d1, int d2, int mid1, int mid2,
	     struct tracekn_s *curr_tr, struct traceknstack_s *dolist, int traceback)
{
  int i, k, l;

  i = j - d;
  k = i + d1;
  l = j - d2;

  if (traceback) fprintf(outf," (WH18.2) trace ZHX %d %d %d %d %d \n", 
			 j-mid2-1, d-mid2-2, d1-1, d2-mid1-mid2-2, 
			 zhx[j-mid2-1][d-mid2-2][d1-1][d2-mid1-mid2-2]); 
  if (traceback) fprintf(outf," (WH18.2) trace YHX %d %d %d %d %d \n", j, d2, mid1, mid2,
			 yhx[j][d2][mid1][mid2]); 
  
  PushTraceknstack(dolist, AttachTracekn(curr_tr, 
					 i+1, j-mid2-1, k, l+mid1+1, dpcP, dpcS));
  PushTraceknstack(dolist, AttachTracekn(curr_tr, 
					 l, j, l+mid1, j-mid2, dpcS, dpcP));
}

/* (WH18.3)__/ l+mid1+1 dangles off j-mid2,l+mid1 */
int 
WH18_3(int *s, int len, int **icfg, int ****zhx, int ****yhx, 
       int j, int d, int d1, int d2, int mid1, int mid2)
{ 
  int  WH18_3;
  int i, k, l;
  
  i = j - d;
  k = i + d1;
  l = j - d2;
  
  if (d2-mid1-mid2 > 2)
    WH18_3 = P12 + 2*P10P + P6P
      + wkn*icfg[idxR(s[l+mid1+1])][idxP(s[j-mid2],s[l+mid1])]
      + wkn*icfg[idxPS(s[j-mid2-1],s[i])][idxPS(s[j-mid2],s[l+mid1])]
      + zhx[j-mid2-1][d-mid2-1][d1][d2-mid1-mid2-3] 
      + yhx[j][d2][mid1][mid2];
  else
    WH18_3 = -BIGINT;

  return WH18_3;
}
void
trace_WH18_3(FILE *outf, int ****zhx, int ****yhx, 
	     int j, int d, int d1, int d2, int mid1, int mid2,
	     struct tracekn_s *curr_tr, struct traceknstack_s *dolist, int traceback)
{
  int i, k, l;

  i = j - d;
  k = i + d1;
  l = j - d2;

  if (traceback) fprintf(outf," (WH18.3) trace ZHX %d %d %d %d %d \n", 
			 j-mid2-1, d-mid2-1, d1, d2-mid1-mid2-3, 
			 zhx[j-mid2-1][d-mid2-1][d1][d2-mid1-mid2-3]); 
  if (traceback) fprintf(outf," (WH18.3) trace YHX %d %d %d %d %d \n", 
			 j, d2, mid1, mid2,
			 yhx[j][d2][mid1][mid2]); 
  
  PushTraceknstack(dolist, AttachTracekn(curr_tr, 
					 i, j-mid2-1, k, l+mid1+2, dpcP, dpcS));
  PushTraceknstack(dolist, AttachTracekn(curr_tr, 
					 l, j, l+mid1, j-mid2, dpcS, dpcP));
}

/* (WH18.4)__/ i dangles off i+1,j-mid2-1 / l+mid1+1 dangles off j-mid2,l+mid1 */
int 
WH18_4(int *s, int len, int **icfg, int ****zhx, int ****yhx, 
       int j, int d, int d1, int d2, int mid1, int mid2)
{ 
  int  WH18_4;
  int i, k, l;
  
  i = j - d;
  k = i + d1;
  l = j - d2;
  
  if (d1 > 0 && d2-mid1-mid2 > 2)
    WH18_4 = P12 + 2*P10P + 2*P6P
      + wkn*icfg[idxL(s[i])][idxP(s[i+1],s[j-mid2-1])]
      + wkn*icfg[idxR(s[l+mid1+1])][idxP(s[j-mid2],s[l+mid1])]
      + wkn*icfg[idxPS(s[j-mid2-1],s[i+1])][idxPS(s[j-mid2],s[l+mid1])]
      + zhx[j-mid2-1][d-mid2-2][d1-1][d2-mid1-mid2-3] 
      + yhx[j][d2][mid1][mid2];
  else
    WH18_4 = -BIGINT;
  
  return WH18_4;
}
void
trace_WH18_4(FILE *outf, int ****zhx, int ****yhx, 
	     int j, int d, int d1, int d2, int mid1, int mid2,
	     struct tracekn_s *curr_tr, struct traceknstack_s *dolist, int traceback)
{
  int i, k, l;

  i = j - d;
  k = i + d1;
  l = j - d2;

  if (traceback) fprintf(outf," (WH18.4) trace ZHX %d %d %d %d %d \n", 
			 j-mid2-1, d-mid2-2, d1-1, d2-mid1-mid2-3, 
			 zhx[j-mid2-1][d-mid2-2][d1-1][d2-mid1-mid2-3]); 
  if (traceback) fprintf(outf," (WH18.4) trace YHX %d %d %d %d %d \n", 
			 j, d2, mid1, mid2,
			 yhx[j][d2][mid1][mid2]); 
  
  PushTraceknstack(dolist, AttachTracekn(curr_tr, 
					 i+1, j-mid2-1, k, l+mid1+2, dpcP, dpcS));
  PushTraceknstack(dolist, AttachTracekn(curr_tr, 
					 l, j, l+mid1, j-mid2, dpcS, dpcP));
}

/* (WH19)__/ 2 YHX 
 * coaxial pairs (k,l+mid1+1) (l+mid1,j-mid2)
 * ....add stack[l+mid1][j-mid2][l+mid1+1][k]
 */

/* (WH19.1)__/ no danglings / coaxial = stack - 1 */                                   
int 
WH19_1(int *s, int len, int **icfg, int ****yhx, 
       int j, int d, int d1, int d2, int mid1, int mid2)
{ 
  int   WH19_1;
  int i, k, l;
  
  i = j - d;
  k = i + d1;
  l = j - d2;
  
  WH19_1 = P12 + 2*P10P 
    + wkn*(icfg[idxPS(s[l+mid1],s[j-mid2])]
	   [idxPS(s[l+mid1+1],s[k])] + INTSCALE)
    + yhx[j-mid2-1][d-mid2-1][d1][d2-mid1-mid2-2] 
    + yhx[j][d2][mid1][mid2];
  
  return WH19_1;
}
void
trace_WH19_1(FILE *outf, int ****yhx, int j, int d, int d1, int d2, int mid1, int mid2,
	     struct tracekn_s *curr_tr, struct traceknstack_s *dolist, int traceback)
{
  int i, k, l;

  i = j - d;
  k = i + d1;
  l = j - d2;

  if (traceback) fprintf(outf," (WH19.1) trace YHX %d %d %d %d %d \n", 
			 j-mid2-1, d-mid2-1, d1, d2-mid1-mid2-2, 
			 yhx[j-mid2-1][d-mid2-1][d1][d2-mid1-mid2-2]); 
  if (traceback) fprintf(outf," (WH19.1) trace YHX %d %d %d %d %d \n", 
			 j, d2, mid1, mid2,
			 yhx[j][d2][mid1][mid2]); 
  
  PushTraceknstack(dolist, AttachTracekn(curr_tr, 
					 i, j-mid2-1, k, l+mid1+1, dpcS, dpcP));
  PushTraceknstack(dolist, AttachTracekn(curr_tr, 
					 l, j, l+mid1, j-mid2, dpcS, dpcP));
}

/* (WH19.2)__/ k dangles off l+mid1+1,k-1 */
int 
WH19_2(int *s, int len, int **icfg, int ****yhx, 
       int j, int d, int d1, int d2, int mid1, int mid2)
{ 
  int  WH19_2;
  int i, k, l;
  
  i = j - d;
  k = i + d1;
  l = j - d2;
  
  if (d1 > 0)
    WH19_2 = P12 + 2*P10P + P6P
      + wkn*icfg[idxR(s[k])][idxP(s[l+mid1+1],s[k-1])]
      + wkn*icfg[idxPS(s[l+mid1],s[j-mid2])][idxPS(s[l+mid1+1],s[k-1])]
      + yhx[j-mid2-1][d-mid2-1][d1-1][d2-mid1-mid2-2] 
      + yhx[j][d2][mid1][mid2];
  else
    WH19_2 = -BIGINT;

  return WH19_2;
}
void
trace_WH19_2(FILE *outf, int ****yhx, int j, int d, int d1, int d2, int mid1, int mid2,
	     struct tracekn_s *curr_tr, struct traceknstack_s *dolist, int traceback)
{
  int i, k, l;

  i = j - d;
  k = i + d1;
  l = j - d2;

  if (traceback) fprintf(outf," (WH19.2) trace YHX %d %d %d %d %d \n", 
			 j-mid2-1, d-mid2-1, d1-1, d2-mid1-mid2-2, 
			 yhx[j-mid2-1][d-mid2-1][d1-1][d2-mid1-mid2-2]); 
  if (traceback) fprintf(outf," (WH19.2) trace YHX %d %d %d %d %d \n", 
			 j, d2, mid1, mid2,
			 yhx[j][d2][mid1][mid2]); 
  
  PushTraceknstack(dolist, AttachTracekn(curr_tr, 
					 i, j-mid2-1, k-1, l+mid1+1, dpcS, dpcP));
  PushTraceknstack(dolist, AttachTracekn(curr_tr, 
					 l, j, l+mid1, j-mid2, dpcS, dpcP));
}

/* (WH19.3)__/ j-mid2-1 dangles off j-mid2,l+mid1 */
int 
WH19_3(int *s, int len, int **icfg, int ****yhx, 
       int j, int d, int d1, int d2, int mid1, int mid2)
{ 
  int  WH19_3;
  int i, k, l;
  
  i = j - d;
  k = i + d1;
  l = j - d2;
  
  if (d2-mid1-mid2 > 2)
    WH19_3 = P12 + 2*P10P + P6P
      + wkn*icfg[idxL(s[j-mid2-1])][idxP(s[j-mid2],s[l+mid1])]
      + wkn*icfg[idxPS(s[l+mid1],s[j-mid2])][idxPS(s[l+mid1+1],s[k])]
      + yhx[j-mid2-2][d-mid2-2][d1][d2-mid1-mid2-3] 
      + yhx[j][d2][mid1][mid2];
  else
    WH19_3 = -BIGINT;

  return WH19_3;
}
void
trace_WH19_3(FILE *outf, int ****yhx, int j, int d, int d1, int d2, int mid1, int mid2,
	     struct tracekn_s *curr_tr, struct traceknstack_s *dolist, int traceback)
{
  int i, k, l;

  i = j - d;
  k = i + d1;
  l = j - d2;

  if (traceback) fprintf(outf," (WH19.3) trace YHX %d %d %d %d %d \n", 
			 j-mid2-2, d-mid2-2, d1, d2-mid1-mid2-3, 
			 yhx[j-mid2-2][d-mid2-2][d1][d2-mid1-mid2-3]); 
  if (traceback) fprintf(outf," (WH19.3) trace YHX %d %d %d %d %d \n", 
			 j, d2, mid1, mid2,
			 yhx[j][d2][mid1][mid2]); 
  
  PushTraceknstack(dolist, AttachTracekn(curr_tr, 
					 i, j-mid2-2, k, l+mid1+1, dpcS, dpcP));
  PushTraceknstack(dolist, AttachTracekn(curr_tr, 
					 l, j, l+mid1, j-mid2, dpcS, dpcP));
}

/* (WH19.4)__/ k dangles off l+mid1+1,k-1 / j-mid2-1 dangles off j-mid2,l+mid */
int 
WH19_4(int *s, int len, int **icfg, int ****yhx, 
       int j, int d, int d1, int d2, int mid1, int mid2)
{ 
  int  WH19_4;
  int i, k, l;
  
  i = j - d;
  k = i + d1;
  l = j - d2;
  
  if (d1 > 0 && d2-mid1-mid2 > 2)
    WH19_4 = P12 + 2*P10P + 2*P6P
      + wkn*icfg[idxR(s[k])][idxP(s[l+mid1+1],s[k-1])]
      + wkn*icfg[idxL(s[j-mid2-1])][idxP(s[j-mid2],s[l+mid1])]
      + wkn*icfg[idxPS(s[l+mid1],s[j-mid2])][idxPS(s[l+mid1+1],s[k-1])]
      + yhx[j-mid2-2][d-mid2-2][d1-1][d2-mid1-mid2-3] 
      + yhx[j][d2][mid1][mid2];
  else
    WH19_4 = -BIGINT;
  
  return WH19_4;
}
void
trace_WH19_4(FILE *outf, int ****yhx, int j, int d, int d1, int d2, int mid1, int mid2,
	     struct tracekn_s *curr_tr, struct traceknstack_s *dolist, int traceback)
{
  int i, k, l;

  i = j - d;
  k = i + d1;
  l = j - d2;

  if (traceback) fprintf(outf," (WH19.4) trace YHX %d %d %d %d %d \n", 
			 j-mid2-2, d-mid2-2, d1-1, d2-mid1-mid2-3, 
			 yhx[j-mid2-2][d-mid2-2][d1-1][d2-mid1-mid2-3]); 
  if (traceback) fprintf(outf," (WH19.4) trace YHX %d %d %d %d %d \n", 
			 j, d2, mid1, mid2,
			 yhx[j][d2][mid1][mid2]); 
  
  PushTraceknstack(dolist, AttachTracekn(curr_tr, 
					 i, j-mid2-2, k-1, l+mid1+1, dpcS, dpcP));
  PushTraceknstack(dolist, AttachTracekn(curr_tr, 
					 l, j, l+mid1, j-mid2, dpcS, dpcP));
}

/* One     connects (i+mid1+1,k-mid2-1) to (l,j) 
 *                  (j,d-mid1-1,d1-mid1-mid2-2,d2)
 * Another connects (i,i+mid1) to (k-mid2,k)          
 *                  (k,d1,mid1,mid2) 
 */
/* (WH20)__/ 2 WHX */
int 
WH20(int *s, int len, int **icfg, int ****whx, 
     int j, int d, int d1, int d2, int mid1, int mid2)
{ 
  int   WH20;
  int i, k, l;
  
  i = j - d;
  k = i + d1;
  l = j - d2;
  
    WH20 = P12 
      + whx[j][d-mid1-1][d1-mid1-mid2-2][d2] 
      + whx[k][d1][mid1][mid2];

  return WH20;
}
void
trace_WH20(FILE *outf, int ****whx, int j, int d, int d1, int d2, int mid1, int mid2,
	   struct tracekn_s *curr_tr, struct traceknstack_s *dolist, int traceback)
{
  int i, k, l;

  i = j - d;
  k = i + d1;
  l = j - d2;

  if (traceback) fprintf(outf," (WH20) trace WHX %d %d %d %d %d \n", 
			 k, d1, mid1, mid2,
			 whx[k][d1][mid1][mid2]); 
  if (traceback) fprintf(outf," (WH20) trace WHX %d %d %d %d %d \n", 
			 j, d-mid1-1, d1-mid1-mid2-2, d2, 
			 whx[j][d-mid1-1][d1-mid1-mid2-2][d2]); 
  
  PushTraceknstack(dolist, AttachTracekn(curr_tr, 
					 i+mid1+1, j, k-mid2-1, l, dpcS, dpcS));
  PushTraceknstack(dolist, AttachTracekn(curr_tr, 
					 i, k, i+mid1, k-mid2, dpcS, dpcS));
}

/* (WH21)__/ 1 ZHX 1 YHX
 * coaxial pairs (i+mid1+1,j) (i+mid1,k-mid2)
 * ....add stack[i+mid1][k-mid2][i+mid1+1][j]
 */

/* (WH21.1)__/ no danglings / coaxial = stack - 1 */                                   
int 
WH21_1(int *s, int len, int **icfg, int ****zhx, int ****yhx, 
       int j, int d, int d1, int d2, int mid1, int mid2)
{ 
  int   WH21_1;
  int i, k, l;
  
  i = j - d;
  k = i + d1;
  l = j - d2;
  
  WH21_1 = P12 + 2*P10P 
    + wkn*(icfg[idxPS(s[i+mid1],s[k-mid2])]
	   [idxPS(s[i+mid1+1],s[j])] + INTSCALE)
    + zhx[j][d-mid1-1][d1-mid1-mid2-2][d2] 
    + yhx[k][d1][mid1][mid2];
  
  return WH21_1;
}
void
trace_WH21_1(FILE *outf, int ****zhx, int ****yhx, 
	     int j, int d, int d1, int d2, int mid1, int mid2,
	     struct tracekn_s *curr_tr, struct traceknstack_s *dolist, int traceback)
{
  int i, k, l;

  i = j - d;
  k = i + d1;
  l = j - d2;

  if (traceback) fprintf(outf," (WH21.1) trace YHX %d %d %d %d %d \n", 
			 k, d1, mid1, mid2,
			 yhx[k][d1][mid1][mid2]); 
  if (traceback) fprintf(outf," (WH21.1) trace ZHX %d %d %d %d %d \n", 
			 j, d-mid1-1, d1-mid1-mid2-2, d2, 
			 zhx[j][d-mid1-1][d1-mid1-mid2-2][d2]); 
  
  PushTraceknstack(dolist, AttachTracekn(curr_tr, 
					 i, k, i+mid1, k-mid2, dpcS, dpcP));
  PushTraceknstack(dolist, AttachTracekn(curr_tr, 
					 i+mid1+1, j, k-mid2-1, l, dpcP, dpcS));
}

/* (WH21.2)__/ j dangles off i+mid1+1,j-1 */
int 
WH21_2(int *s, int len, int **icfg, int ****zhx, int ****yhx, 
       int j, int d, int d1, int d2, int mid1, int mid2)
{ 
  int   WH21_2;
  int i, k, l;
  
  i = j - d;
  k = i + d1;
  l = j - d2;
  
  if (d2 > 0)
    WH21_2 = P12 + 2*P10P + P6P
      + wkn*icfg[idxR(s[j])][idxP(s[i+mid1+1],s[j-1])]
      + wkn*icfg[idxPS(s[i+mid1],s[k-mid2])][idxPS(s[i+mid1+1],s[j-1])]
      + zhx[j-1][d-mid1-2][d1-mid1-mid2-2][d2-1] 
      + yhx[k][d1][mid1][mid2];
  else
    WH21_2 = -BIGINT;

  return WH21_2;
}
void
trace_WH21_2(FILE *outf, int ****zhx, int ****yhx, 
	     int j, int d, int d1, int d2, int mid1, int mid2,
	     struct tracekn_s *curr_tr, struct traceknstack_s *dolist, int traceback)
{
  int i, k, l;

  i = j - d;
  k = i + d1;
  l = j - d2;

  if (traceback) fprintf(outf," (WH21.2) trace YHX %d %d %d %d %d \n", 
			 k, d1, mid1, mid2,
			 yhx[k][d1][mid1][mid2]); 
  if (traceback) fprintf(outf," (WH21.2) trace ZHX %d %d %d %d %d \n", 
			 j-1, d-mid1-2, d1-mid1-mid2-2, d2-1, 
			 zhx[j-1][d-mid1-2][d1-mid1-mid2-2][d2-1]); 
  
  PushTraceknstack(dolist, AttachTracekn(curr_tr, 
					 i, k, i+mid1, k-mid2, dpcS, dpcP));
  PushTraceknstack(dolist, AttachTracekn(curr_tr, 
					 i+mid1+1, j-1, k-mid2-1, l, dpcP, dpcS));
}

/* (WH21.3)__/ k-mid2-1 dangles off k-mid2,i+mid1 */
int 
WH21_3(int *s, int len, int **icfg, int ****zhx, int ****yhx, 
       int j, int d, int d1, int d2, int mid1, int mid2)
{ 
  int   WH21_3;
  int i, k, l;
  
  i = j - d;
  k = i + d1;
  l = j - d2;
  
  if (d1-mid1-mid2 > 2)
    WH21_3 = P12 + 2*P10P + P6P
      + wkn*icfg[idxL(s[k-mid2-1])][idxP(s[ k-mid2],s[i+mid1])]
      + wkn*icfg[idxPS(s[i+mid1],s[k-mid2])][idxPS(s[i+mid1+1],s[j])] 
      + zhx[j][d-mid1-1][d1-mid1-mid2-3][d2] 
      + yhx[k][d1][mid1][mid2];
  else
    WH21_3 = -BIGINT;

  return WH21_3;
}
void
trace_WH21_3(FILE *outf, int ****zhx, int ****yhx, 
	     int j, int d, int d1, int d2, int mid1, int mid2,
	     struct tracekn_s *curr_tr, struct traceknstack_s *dolist, int traceback)
{
  int i, k, l;

  i = j - d;
  k = i + d1;
  l = j - d2;

  if (traceback) fprintf(outf," (WH21.3) trace YHX %d %d %d %d %d \n", 
			 k, d1, mid1, mid2,
			 yhx[k][d1][mid1][mid2]); 
  if (traceback) fprintf(outf," (WH21.3) trace ZHX %d %d %d %d %d \n", 
			 j, d-mid1-1, d1-mid1-mid2-3, d2, 
			 zhx[j][d-mid1-1][d1-mid1-mid2-3][d2]); 
  
  PushTraceknstack(dolist, AttachTracekn(curr_tr, 
					 i, k, i+mid1, k-mid2, dpcS, dpcP));
  PushTraceknstack(dolist, AttachTracekn(curr_tr, 
					 i+mid1+1, j, k-mid2-2, l, dpcP, dpcS));
}

/* (WH21.4)__/ j dangles off i+mid1+1,j-1 
             / k-mid2-1 dangles off k-mid2,i+mid1 */
int 
WH21_4(int *s, int len, int **icfg, int ****zhx, int ****yhx, 
       int j, int d, int d1, int d2, int mid1, int mid2)
{ 
  int   WH21_4;
  int i, k, l;
  
  i = j - d;
  k = i + d1;
  l = j - d2;
  
  if (d2 > 0 && d1-mid1-mid2 > 2)
    WH21_4 = P12 + 2*P10P + 2*P6P
      + wkn*icfg[idxR(s[j])][idxP(s[i+mid1+1],s[j-1])]
      + wkn*icfg[idxL(s[k-mid2-1])][idxP(s[ k-mid2],s[i+mid1])]
      + wkn*icfg[idxPS(s[i+mid1],s[k-mid2])][idxPS(s[i+mid1+1],s[j-1])]
      + zhx[j-1][d-mid1-2][d1-mid1-mid2-3][d2-1] 
      + yhx[k][d1][mid1][mid2];
  else
    WH21_4 = -BIGINT;
  
  return WH21_4;
}
void
trace_WH21_4(FILE *outf, int ****zhx, int ****yhx, 
	     int j, int d, int d1, int d2, int mid1, int mid2,
	     struct tracekn_s *curr_tr, struct traceknstack_s *dolist, int traceback)
{
  int i, k, l;

  i = j - d;
  k = i + d1;
  l = j - d2;

  if (traceback) fprintf(outf," (WH21.4) trace YHX %d %d %d %d %d \n", 
			 k, d1, mid1, mid2,
			 yhx[k][d1][mid1][mid2]); 
  if (traceback) fprintf(outf," (WH21.4) trace ZHX %d %d %d %d %d \n", 
			 j-1, d-mid1-2, d1-mid1-mid2-3, d2-1, 
	  zhx[j-1][d-mid1-2][d1-mid1-mid2-3][d2-1]); 
  
  PushTraceknstack(dolist, AttachTracekn(curr_tr, 
					 i, k, i+mid1, k-mid2, dpcS, dpcP));
  PushTraceknstack(dolist, AttachTracekn(curr_tr, 
					 i+mid1+1, j-1, k-mid2-2, l, dpcP, dpcS));
}

/* (WH22)__/ 2 YHX 
 * coaxial pairs (k-mid2-1,l) (i+mid1,k-mid2)
 * ....add stack[k-mid2-1][l][k-mid2][i+mid1]
 */

/* (WH22.1)__/ no danglings / 2 coaxial = stack - 1 */                                   
int 
WH22_1(int *s, int len, int **icfg, int ****yhx, 
       int j, int d, int d1, int d2, int mid1, int mid2)
{ 
  int   WH22_1;
  int i, k, l;
  
  i = j - d;
  k = i + d1;
  l = j - d2;
  
  WH22_1 = P12 + 2*P10P 
    + wkn*(icfg[idxPS(s[k-mid2-1],s[l])]
	   [idxPS(s[k-mid2],s[i+mid1])] + INTSCALE)
    + yhx[j][d-mid1-1][d1-mid1-mid2-2][d2] 
    + yhx[k][d1][mid1][mid2];
  
  return WH22_1;
}
void
trace_WH22_1(FILE *outf, int ****yhx, int j, int d, int d1, int d2, int mid1, int mid2,
	     struct tracekn_s *curr_tr, struct traceknstack_s *dolist, int traceback)
{
  int i, k, l;

  i = j - d;
  k = i + d1;
  l = j - d2;

  if (traceback) fprintf(outf," (WH22.1) trace YHX %d %d %d %d %d \n", 
			 k, d1, mid1, mid2,
			 yhx[k][d1][mid1][mid2]); 
  if (traceback) fprintf(outf," (WH22.1) trace YHX %d %d %d %d %d \n", 
			 j, d-mid1-1, d1-mid1-mid2-2, d2, 
			 yhx[j][d-mid1-1][d1-mid1-mid2-2][d2]); 
  
  PushTraceknstack(dolist, AttachTracekn(curr_tr, 
					 i, k, i+mid1, k-mid2, dpcS, dpcP));
  PushTraceknstack(dolist, AttachTracekn(curr_tr, 
					 i+mid1+1, j, k-mid2-1, l, dpcS, dpcP));
}

/* (WH22.2)__/ l dangles off l+1,k-mid2-1 */
int 
WH22_2(int *s, int len, int **icfg, int ****yhx, 
       int j, int d, int d1, int d2, int mid1, int mid2)
{ 
  int   WH22_2;
  int i, k, l;
  
  i = j - d;
  k = i + d1;
  l = j - d2;
  
  if (d2 > 0)
    WH22_2 = P12 + 2*P10P + P6P
      + wkn*icfg[idxL(s[l])][idxP(s[l+1],s[k-mid2-1])]
      + wkn*icfg[idxPS(s[k-mid2-1],s[l+1])][idxPS(s[k-mid2],s[i+mid1])]
      + yhx[j][d-mid1-1][d1-mid1-mid2-2][d2-1] 
      + yhx[k][d1][mid1][mid2];
  else
    WH22_2 = -BIGINT;

  return WH22_2;
}
void
trace_WH22_2(FILE *outf, int ****yhx, int j, int d, int d1, int d2, int mid1, int mid2,
	     struct tracekn_s *curr_tr, struct traceknstack_s *dolist, int traceback)
{
  int i, k, l;

  i = j - d;
  k = i + d1;
  l = j - d2;

  if (traceback) fprintf(outf," (WH22.2) trace YHX %d %d %d %d %d \n", 
			 k, d1, mid1, mid2,
			 yhx[k][d1][mid1][mid2]); 
  if (traceback) fprintf(outf," (WH22.2) trace YHX %d %d %d %d %d \n", 
			 j, d-mid1-1, d1-mid1-mid2-2, d2-1, 
			 yhx[j][d-mid1-1][d1-mid1-mid2-2][d2-1]); 
  
  PushTraceknstack(dolist, AttachTracekn(curr_tr, 
					 i, k, i+mid1, k-mid2, dpcS, dpcP));
  PushTraceknstack(dolist, AttachTracekn(curr_tr, 
					 i+mid1+1, j, k-mid2-1, l, dpcS, dpcP));
}

/* (WH22.3)__/ i+mid1+1 dangles off k-mid2,i+mid1 */
int 
WH22_3(int *s, int len, int **icfg, int ****yhx, 
       int j, int d, int d1, int d2, int mid1, int mid2)
{ 
  int   WH22_3;
  int i, k, l;
  
  i = j - d;
  k = i + d1;
  l = j - d2;
  
  if (d1-mid1-mid2 > 2)
    WH22_3 = P12 + 2*P10P + P6P 
      + wkn*icfg[idxR(s[i+mid1+1])][idxP(s[k-mid2],s[i+mid1])]
      + wkn*icfg[idxPS(s[k-mid2-1],s[l])][idxPS(s[k-mid2],s[i+mid1])]
      + yhx[j][d-mid1-2][d1-mid1-mid2-3][d2] 
      + yhx[k][d1][mid1][mid2];
  else
    WH22_3 = -BIGINT;

  return WH22_3;
}
void
trace_WH22_3(FILE *outf, int ****yhx, int j, int d, int d1, int d2, int mid1, int mid2,
	     struct tracekn_s *curr_tr, struct traceknstack_s *dolist, int traceback)
{
  int i, k, l;

  i = j - d;
  k = i + d1;
  l = j - d2;

  if (traceback) fprintf(outf," (WH22.3) trace YHX %d %d %d %d %d \n", 
			 k, d1, mid1, mid2,
			 yhx[k][d1][mid1][mid2]); 
  if (traceback) fprintf(outf," (WH22.3) trace YHX %d %d %d %d %d \n", 
			 j, d-mid1-2, d1-mid1-mid2-3, d2, 
			 yhx[j][d-mid1-2][d1-mid1-mid2-3][d2]); 
  
  PushTraceknstack(dolist, AttachTracekn(curr_tr, 
					 i, k, i+mid1, k-mid2, dpcS, dpcP));
  PushTraceknstack(dolist, AttachTracekn(curr_tr, 
					 i+mid1+2, j, k-mid2-1, l, dpcS, dpcP));
}

/* (WH22.4)__/ l dangles off l+1,k-mid2-1 / i+mid1+1 dangles off k-mid2,i+mid1 */
int 
WH22_4(int *s, int len, int **icfg, int ****yhx, 
       int j, int d, int d1, int d2, int mid1, int mid2)
{ 
  int   WH22_4;
  int i, k, l;
  
  i = j - d;
  k = i + d1;
  l = j - d2;
  
  if (d2 > 0 && d1-mid1-mid2 > 2)
    WH22_4 = P12 + 2*P10P + 2*P6P
      + wkn*icfg[idxL(s[l])][idxP(s[l+1],s[k-mid2-1])]
      + wkn*icfg[idxR(s[i+mid1+1])][idxP(s[k-mid2],s[i+mid1])]
      + wkn*icfg[idxPS(s[k-mid2-1],s[l+1])][idxPS(s[k-mid2],s[i+mid1])]
      + yhx[j][d-mid1-2][d1-mid1-mid2-3][d2-1] 
      + yhx[k][d1][mid1][mid2];
  else
    WH22_4 = -BIGINT;
  
  return WH22_4;
}
void
trace_WH22_4(FILE *outf, int ****yhx, int j, int d, int d1, int d2, int mid1, int mid2,
	     struct tracekn_s *curr_tr, struct traceknstack_s *dolist, int traceback)
{
  int i, k, l;

  i = j - d;
  k = i + d1;
  l = j - d2;

  if (traceback) fprintf(outf," (WH22.4) trace YHX %d %d %d %d %d \n", 
			 k, d1, mid1, mid2,
			 yhx[k][d1][mid1][mid2]); 
  if (traceback) fprintf(outf," (WH22.4) trace YHX %d %d %d %d %d \n", 
			 j, d-mid1-2, d1-mid1-mid2-3, d2-1, 
			 yhx[j][d-mid1-2][d1-mid1-mid2-3][d2-1]); 
  
  PushTraceknstack(dolist, AttachTracekn(curr_tr, 
					 i, k, i+mid1, k-mid2, dpcS, dpcP));
  PushTraceknstack(dolist, AttachTracekn(curr_tr, 
					 i+mid1+2, j, k-mid2-1, l, dpcS, dpcP));
}




