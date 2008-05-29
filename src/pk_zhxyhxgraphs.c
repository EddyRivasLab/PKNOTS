/* zhxyhxgraphs.c 
 *
 * includes functions to calculate all the diagrams 
 * that fill and traceback the no-hole matrices:
 *      zhx[j][d][d1][d2] (len x len x len x len)
 *                        [(j-d,j) are base pared]
 *      yhx[j][d][d1][d2] (len x len x len x len)
 *                        [(j-d+d1,j-d2) are base pared]
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
#include "pk_zhxyhxgraphs.h"
#include "pk_trace.h"
#include "pk_util.h"


/* diagrams for ZHX (ZH1-ZH12)
 */

/* (ZH1)__/ pair k,l */
int 
ZH1(ESL_DSQ *s, int len, struct rnapar_2 *rnapar, int **icfg, int ****vhx, int j, int d, int d1, int d2)
{ 
  int     ZH1;
  int i, k, l;

  i = j - d;
  k = i + d1;
  l = j - d2;
  
  ZH1 = rnapar->P10P + vhx[j][d][d1][d2];

  return ZH1;
}
void
trace_ZH1(FILE *outf, int ****vhx, int j, int d, int d1, int d2,
	  struct tracekn_s *curr_tr, struct traceknstack_s *dolist, int traceback)
{
  int i, k, l;

  i = j - d;
  k = i + d1;
  l = j - d2;

  if (traceback) fprintf(outf," (ZH1) trace vhx %d %d %d %d %d \n", 
			 j, d, d1, d2, vhx[j][d][d1][d2]); 
  
  PushTraceknstack(dolist, AttachTracekn(curr_tr, i, j, k, l, dpcP, dpcP)); 
}

/* (ZH2)__/  pair k+1,l-1 */
int 
ZH2(ESL_DSQ *s, int len, struct rnapar_2 *rnapar, int **icfg, int ****vhx, int j, int d, int d1, int d2)
{ 
  int     ZH2;
  int i, k, l;

  i = j - d;
  k = i + d1;
  l = j - d2;
  
  if (d1 > 0 && d2 > 0)
    ZH2 = rnapar->P10P + 2*rnapar->P6P 
      + wkn*icfg[idxR(s[k])][idxP(s[k-1],s[l+1])]
      + wkn*icfg[idxL(s[l])][idxP(s[k-1],s[l+1])]
      + vhx[j][d][d1-1][d2-1];
  else
    ZH2 = -BIGINT;

  return ZH2;
}
void
trace_ZH2(FILE *outf, int ****vhx, int j, int d, int d1, int d2,
	  struct tracekn_s *curr_tr, struct traceknstack_s *dolist, int traceback)
{
  int i, k, l;

  i = j - d;
  k = i + d1;
  l = j - d2;

  if (traceback) fprintf(outf," (ZH2) dangle k,l, trace vhx %d %d %d %d %d \n", 
			 j, d, d1-1, d2-1, vhx[j][d][d1-1][d2-1]); 
  
  curr_tr->type2 = dpcL;           
  PushTraceknstack(dolist, AttachTracekn(curr_tr, i, j, k-1, l+1, dpcP, dpcP)); 
}


/* (ZH3)__/ k dangles off k-1,l base-paired  */
int 
ZH3(ESL_DSQ *s, int len, struct rnapar_2 *rnapar, int **icfg, int ****vhx, int j, int d, int d1, int d2)
{ 
  int     ZH3;
  int i, k, l;

  i = j - d;
  k = i + d1;
  l = j - d2;
  
  if (d1 > 0)
    ZH3 = rnapar->P10P + rnapar->P6P 
      + wkn*icfg[idxR(s[k])][idxP(s[k-1],s[l])]
      + vhx[j][d][d1-1][d2];
  else
    ZH3 = -BIGINT;

  return ZH3;
}
void
trace_ZH3(FILE *outf, int ****vhx, int j, int d, int d1, int d2,
	  struct tracekn_s *curr_tr, struct traceknstack_s *dolist, int traceback)
{
  int i, k, l;

  i = j - d;
  k = i + d1;
  l = j - d2;

  if (traceback) fprintf(outf," (ZH3) dangle k, trace vhx %d %d %d %d %d \n", 
			 j, d, d1-1, d2, vhx[j][d][d1-1][d2]); 
  
  curr_tr->type2 = dpcL;           
  PushTraceknstack(dolist, AttachTracekn(curr_tr, i, j, k-1, l, dpcP, dpcP)); 
}

/* (ZH4)__/ l dangles off k,l+1 base-paired 
 */
int
ZH4(ESL_DSQ *s, int len, struct rnapar_2 *rnapar, int **icfg, int ****vhx, int j, int d, int d1, int d2)
{ 
  int     ZH4;
  int i, k, l;

  i = j - d;
  k = i + d1;
  l = j - d2;
  
  if (d2 > 0)
    ZH4 = rnapar->P10P + rnapar->P6P 
      + wkn*icfg[idxL(s[l])][idxP(s[k],s[l+1])]
      + vhx[j][d][d1][d2-1];
  else
    ZH4 = -BIGINT;

  return ZH4;
}
void
trace_ZH4(FILE *outf, int ****vhx, int j, int d, int d1, int d2,
	  struct tracekn_s *curr_tr, struct traceknstack_s *dolist, int traceback)
{
  int i, k, l;

  i = j - d;
  k = i + d1;
  l = j - d2;

  if (traceback) fprintf(outf," (ZH4) dangle l, trace vhx %d %d %d %d %d \n", 
			 j, d, d1, d2-1, vhx[j][d][d1][d2-1]); 
  
  curr_tr->type2 = dpcL;           
  PushTraceknstack(dolist, AttachTracekn(curr_tr, i, j, k, l+1, dpcP, dpcP)); 
}

/* (ZH5)__/  k ss off k-1,l  */
int 
ZH5(ESL_DSQ *s, int len, struct rnapar_2 *rnapar, int **icfg, int ****zhx, int j, int d, int d1, int d2)
{ 
  int     ZH5;
  int i, k, l;

  i = j - d;
  k = i + d1;
  l = j - d2;
  
  if (d1 > 0)
    ZH5 = rnapar->P6P + zhx[j][d][d1-1][d2];
  else
    ZH5 = -BIGINT;

  return ZH5;
}
void
trace_ZH5(FILE *outf, int ****zhx, int j, int d, int d1, int d2,
	  struct tracekn_s *curr_tr, struct traceknstack_s *dolist, int traceback)
{
  int i, k, l;

  i = j - d;
  k = i + d1;
  l = j - d2;

  if (traceback) fprintf(outf," (ZH5) ss dangle k, trace zhx %d %d %d %d %d \n", 
			 j, d, d1-1, d2, zhx[j][d][d1-1][d2]); 
  
  curr_tr->type2 = dpcL;           
  PushTraceknstack(dolist, AttachTracekn(curr_tr, i, j, k-1, l, dpcP, dpcS)); 
}

/* (ZH6)__/ l dangles ss off k,l+1  */
int 
ZH6(ESL_DSQ *s, int len, struct rnapar_2 *rnapar, int **icfg, int ****zhx, int j, int d, int d1, int d2)
{ 
  int     ZH6;
  int i, k, l;

  i = j - d;
  k = i + d1;
  l = j - d2;
  
  if (d2 > 0)
    ZH6 = rnapar->P6P + zhx[j][d][d1][d2-1];
  else
    ZH6 = -BIGINT;

  return ZH6;
}
void
trace_ZH6(FILE *outf, int ****zhx, int j, int d, int d1, int d2,
	  struct tracekn_s *curr_tr, struct traceknstack_s *dolist, int traceback)
{
  int i, k, l;

  i = j - d;
  k = i + d1;
  l = j - d2;

  if (traceback) fprintf(outf," (ZH6) ss dangle l, trace zhx %d %d %d %d %d \n", 
			 j, d, d1, d2-1, zhx[j][d][d1][d2-1]); 
  
  curr_tr->type2 = dpcL;           
  PushTraceknstack(dolist, AttachTracekn(curr_tr, i, j, k, l+1, dpcP, dpcS)); 
}

/* 
 * STRUCTURES WITH ONE WBX AND ONE ZHX. (ZH7, ZH9)
 * or         WITH ONE VX  AND ONE VHX. (ZH8, ZH10)
 */

/* In the LEFT ARM (i,k).
 */
/* Interior RIGHT */

/* (ZH7)__/ 1 WBX 1 ZHX
 * WBX connects k-mid with k (k,mid)
 * ZHX connects (i,k-mid-1) and (l,j) (j,d,d1-mid-1,d2).
 */                                       
int 
ZH7(ESL_DSQ *s, int len, struct rnapar_2 *rnapar, int **icfg, int **wbx,int ****zhx, int j, int d, int d1, int d2,int mid)
{ 
  int     ZH7;
  int i, k, l;

  i = j - d;
  k = i + d1;
  l = j - d2;
  
  ZH7 = wbx[k][mid] + zhx[j][d][d1-mid-1][d2];

  return ZH7;
}
void
trace_ZH7(FILE *outf, int **wbx, int ****zhx, int j, int d, int d1, int d2, int mid,
	  struct tracekn_s *curr_tr, struct traceknstack_s *dolist, int traceback)
{
  int i, k, l;

  i = j - d;
  k = i + d1;
  l = j - d2;

  if (traceback) fprintf(outf," (ZH7) IR, trace  wbx %d  %d %d \n", k, mid, wbx[k][mid]); 
  if (traceback) fprintf(outf," (ZH7) IR, trace  zhx %d  %d  %d  %d %d \n", 
			 j, d, d1-mid-1, d2, zhx[j][d][d1-mid-1][d2]); 
  
  PushTraceknstack(dolist, AttachTracekn(curr_tr, i, j, k-mid-1, l, dpcP, dpcS));
  PushTraceknstack(dolist, AttachTracekn(curr_tr, k-mid, k, k-mid+(int)(mid/2), 
					 k-mid+(int)(mid/2)+1, dpcB, dpcS));
}

/* (ZH8)__/ 1 VX 1 VHX 
 * VX  connects k-mid with k (k,mid)
 * VHX connects (i,k-mid-1) and (l,j) (j,d,d1-mid-1,d2).
 * coaxial pairs (k-mid-1,l) (k-mid,k)
 * ....add stack[k-mid-1][l][k-mid][k]
 */

/* (ZH8.1)__/ no danglings / coaxial = stacking + 1 */
int 
ZH8_1(ESL_DSQ *s, int len, struct rnapar_2 *rnapar, int **icfg, int **vx,int ****vhx, int j, int d, int d1, int d2,int mid)
{ 
  int   ZH8_1;
  int i, k, l;

  i = j - d;
  k = i + d1;
  l = j - d2;
  
  ZH8_1 = 2*rnapar->P10P 
    + wkn*(icfg[idxPS(s[k-mid-1],s[l])]
	   [idxPS(s[k-mid],s[k])] + INTSCALE)
    + vx[k][mid] + vhx[j][d][d1-mid-1][d2];
  
  return ZH8_1;
}
void
trace_ZH8_1(FILE *outf, int **vx, int ****vhx, int j, int d, int d1, int d2, int mid,
	    struct tracekn_s *curr_tr, struct traceknstack_s *dolist, int traceback)
{
  int i, k, l;

  i = j - d;
  k = i + d1;
  l = j - d2;

  if (traceback) fprintf(outf," (ZH8.1) IR, trace vx  %d  %d %d \n", k, mid, vx[k][mid]); 
  if (traceback) fprintf(outf," (ZH8.1) IR, trace vhx %d  %d  %d  %d %d \n", 
			 j, d, d1-mid-1, d2, vhx[j][d][d1-mid-1][d2]); 
  
  PushTraceknstack(dolist, AttachTracekn(curr_tr, i, j, k-mid-1, l, dpcP, dpcP));
  PushTraceknstack(dolist, AttachTracekn(curr_tr, k-mid, k, k-mid+(int)(mid/2), 
					 k-mid+(int)(mid/2)+1, dpcP, dpcS));
}

/* (ZH8.2)__/ k dangles off k-mid,k-1 */
int 
ZH8_2(ESL_DSQ *s, int len, struct rnapar_2 *rnapar, int **icfg, int **vx,int ****vhx, int j, int d, int d1, int d2,int mid)
{ 
  int   ZH8_2;
  int i, k, l;

  i = j - d;
  k = i + d1;
  l = j - d2;
  
  if (mid > 0)
    ZH8_2 = 2*rnapar->P10P + rnapar->P6P 
      + wkn*icfg[idxR(s[k])][idxP(s[k-mid],s[k-1])]
      + wkn*icfg[idxPS(s[k-mid-1],s[l])][idxPS(s[k-mid],s[k-1])]
      + vx[k-1][mid-1] + vhx[j][d][d1-mid-1][d2];
  else
    ZH8_2 = -BIGINT;

  return ZH8_2;
}
void
trace_ZH8_2(FILE *outf, int **vx, int ****vhx, int j, int d, int d1, int d2, int mid,
	    struct tracekn_s *curr_tr, struct traceknstack_s *dolist, int traceback)
{
  int i, k, l;

  i = j - d;
  k = i + d1;
  l = j - d2;

  if (traceback) fprintf(outf," (ZH8.2) IR, dangle k, trace vx %d  %d %d \n", 
			 k-1, mid-1, vx[k-1][mid-1]); 
  if (traceback) fprintf(outf," (ZH8.2) IR, trace  vhx %d  %d  %d  %d %d \n", 
			 j, d, d1-mid-1, d2, vhx[j][d][d1-mid-1][d2]); 
  
  PushTraceknstack(dolist, AttachTracekn(curr_tr, i, j, k-mid-1, l, dpcP, dpcP));
  PushTraceknstack(dolist, AttachTracekn(curr_tr, k-mid, k-1, k-mid+(int)((mid-1)/2), 
					 k-mid+(int)((mid-1)/2)+1, dpcP, dpcS));
}

/* (ZH8.3)__/ l dangles off l+1,k-mid-1 */
int 
ZH8_3(ESL_DSQ *s, int len, struct rnapar_2 *rnapar, int **icfg, int **vx,int ****vhx, int j, int d, int d1, int d2,int mid)
{ 
  int   ZH8_3;
  int i, k, l;

  i = j - d;
  k = i + d1;
  l = j - d2;
  
  if (d2 > 0)
    ZH8_3 = 2*rnapar->P10P + rnapar->P6P 
      + wkn*icfg[idxL(s[l])][idxP(s[l+1],s[k-mid-1])]
      + wkn*icfg[idxPS(s[k-mid-1],s[l+1])][idxPS(s[k-mid],s[k])]
      + vx[k][mid] + vhx[j][d][d1-mid-1][d2-1];
  else
    ZH8_3 = -BIGINT;

  return ZH8_3;
}
void
trace_ZH8_3(FILE *outf, int **vx, int ****vhx, int j, int d, int d1, int d2, int mid,
	    struct tracekn_s *curr_tr, struct traceknstack_s *dolist, int traceback)
{
  int i, k, l;

  i = j - d;
  k = i + d1;
  l = j - d2;

  if (traceback) fprintf(outf," (ZH8.3) IR, trace vx %d  %d %d \n", k, mid, vx[k][mid]); 
  if (traceback) fprintf(outf," (ZH8.3) IR, dangle l, trace vhx %d  %d  %d  %d %d \n", 
			 j, d, d1-mid-1, d2-1, vhx[j][d][d1-mid-1][d2-1]); 
  
  PushTraceknstack(dolist, AttachTracekn(curr_tr, i, j, k-mid-1, l+1, dpcP, dpcP));
  PushTraceknstack(dolist, AttachTracekn(curr_tr, k-mid, k, k-mid+(int)(mid/2), 
					 k-mid+(int)(mid/2)+1, dpcP, dpcS));
}

/* (ZH8.4)__/ k dangles off k-mid,k-1 / l dangles off l+1,k-mid-1 */
int 
ZH8_4(ESL_DSQ *s, int len, struct rnapar_2 *rnapar, int **icfg, int **vx,int ****vhx, int j, int d, int d1, int d2,int mid)
{ 
  int   ZH8_4;
  int i, k, l;

  i = j - d;
  k = i + d1;
  l = j - d2;
  
  if (mid > 0 && d2 > 0)
    ZH8_4 = 2*rnapar->P10P + 2*rnapar->P6P 
      + wkn*icfg[idxR(s[k])][idxP(s[k-mid],s[k-1])]
      + wkn*icfg[idxL(s[l])][idxP(s[l+1],s[k-mid-1])]
      + wkn*icfg[idxPS(s[k-mid-1],s[l+1])][idxPS(s[k-mid],s[k-1])]
      + vx[k-1][mid-1] + vhx[j][d][d1-mid-1][d2-1];
  else
    ZH8_4 = -BIGINT;

  return ZH8_4;
}
void
trace_ZH8_4(FILE *outf, int **vx, int ****vhx, int j, int d, int d1, int d2, int mid,
	    struct tracekn_s *curr_tr, struct traceknstack_s *dolist, int traceback)
{
  int i, k, l;

  i = j - d;
  k = i + d1;
  l = j - d2;

  if (traceback) fprintf(outf," (ZH8.4) IR, dangle k, trace vx  %d  %d %d \n", 
			 k-1, mid-1, vx[k-1][mid-1]); 
  if (traceback) fprintf(outf," (ZH8.4) IR, dangle l, trace vhx %d  %d  %d  %d %d \n", 
			 j, d, d1-mid-1, d2-1, vhx[j][d][d1-mid-1][d2-1]); 
  
  PushTraceknstack(dolist, 
		   AttachTracekn(curr_tr, i, j, k-mid-1, l+1, dpcP, dpcP));
  PushTraceknstack(dolist, 
		   AttachTracekn(curr_tr, 
				 k-mid, k-1, k-mid+(int)((mid-1)/2), 
				 k-mid+(int)((mid-1)/2)+1, dpcP, dpcS));
}

/* In the RIGHT ARM (l,j).
 */      

/* Interior LEFT */

/* (ZH9)__/ 1 WBX 1 ZHX
 * WBX connects l with l+mid (l+mid,mid)
 * ZHX connects (i,k) and (l+mid+1,j) (j, d, d1, d2-mid-1).
 */  
int 
ZH9(ESL_DSQ *s, int len, struct rnapar_2 *rnapar, int **icfg, int **wbx,int ****zhx, int j, int d, int d1, int d2, int mid)
{ 
  int     ZH9;
  int i, k, l;

  i = j - d;
  k = i + d1;
  l = j - d2;
  
  ZH9 = wbx[l+mid][mid] + zhx[j][d][d1][d2-mid-1];

  return ZH9;
}
void
trace_ZH9(FILE *outf, int **wbx, int ****zhx, int j, int d, int d1, int d2, int mid,
	  struct tracekn_s *curr_tr, struct traceknstack_s *dolist, int traceback)
{
  int i, k, l;

  i = j - d;
  k = i + d1;
  l = j - d2;

  if (traceback) fprintf(outf," (ZH9) IL wbx        %d %d %d \n", 
			 l+mid, mid, wbx[l+mid][mid]); 
  if (traceback) fprintf(outf," (ZH9) IL zhx %d %d %d %d %d \n", 
			 j, d, d1, d2-mid-1, zhx[j][d][d1][d2-mid-1]); 
  
  PushTraceknstack(dolist, AttachTracekn(curr_tr, i, j, k, l+mid+1, dpcP, dpcS));
  PushTraceknstack(dolist, AttachTracekn(curr_tr, l, l+mid, l + (int)(mid/2), 
					 l + (int)(mid/2) + 1, dpcB, dpcS));
}

/* (ZH10)__/ 1 VX 1 VHX 
 * VX  connects l with l+mid (l+mid,mid)
 * VHX connects (i,k) and (l+mid+1,j) (j, d, d1, d2-mid-1).
 * coaxial pairs (k,l+mid+1) (l,l+mid)
 * ....add stack[l+mid][l][l+mid+1][k]
 */ 

/* (ZH10.1)__/ no danglings / coaxial = stacking + 1 */
int 
ZH10_1(ESL_DSQ *s, int len, struct rnapar_2 *rnapar, int **icfg, int **vx,int ****vhx, int j, int d, int d1, int d2,int mid)
{ 
  int  ZH10_1;
  int i, k, l;

  i = j - d;
  k = i + d1;
  l = j - d2;
  
    ZH10_1 = 2*rnapar->P10P 
      + wkn*(icfg[idxPS(s[l+mid],s[l])]
	     [idxPS(s[l+mid+1],s[k])] + INTSCALE)
      + vx[l+mid][mid] + vhx[j][d][d1][d2-mid-1];
    
    return ZH10_1;
}
void
trace_ZH10_1(FILE *outf, int **vx, int ****vhx, int j, int d, int d1, int d2, int mid,
	     struct tracekn_s *curr_tr, struct traceknstack_s *dolist, int traceback)
{
  int i, k, l;

  i = j - d;
  k = i + d1;
  l = j - d2;

  if (traceback) fprintf(outf," (ZH10.1) IL vx        %d %d %d \n", 
			 l+mid, mid, vx[l+mid][mid]); 
  if (traceback) fprintf(outf," (ZH10.1) IL vhx %d %d %d %d %d \n", 
			 j, d, d1, d2-mid-1, vhx[j][d][d1][d2-mid-1]); 
  
  PushTraceknstack(dolist, AttachTracekn(curr_tr, i, j, k, l+mid+1, dpcP, dpcP));
  PushTraceknstack(dolist, AttachTracekn(curr_tr, l, l+mid, l+(int)(mid/2), 
					 l+(int)(mid/2)+1, dpcP, dpcS));
}

/* (ZH10.2)__/ l dangles off  l+1,l+mid */
int 
ZH10_2(ESL_DSQ *s, int len, struct rnapar_2 *rnapar, int **icfg, int **vx,int ****vhx, int j, int d, int d1, int d2,int mid)
{ 
  int   ZH10_2;
  int i, k, l;

  i = j - d;
  k = i + d1;
  l = j - d2;
  
  if (mid > 0)
    ZH10_2 = 2*rnapar->P10P + rnapar->P6P 
      + wkn*icfg[idxL(s[l])][idxP(s[l+1],s[l+mid])]
      + wkn*icfg[idxPS(s[l+mid],s[l+1])][idxPS(s[l+mid+1],s[k])]
      + vx[l+mid][mid-1] + vhx[j][d][d1][d2-mid-1];
  else
    ZH10_2 = -BIGINT;

  return ZH10_2;
}
void
trace_ZH10_2(FILE *outf, int **vx, int ****vhx, int j, int d, int d1, int d2, int mid,
	     struct tracekn_s *curr_tr, struct traceknstack_s *dolist, int traceback)
{
  int i, k, l;

  i = j - d;
  k = i + d1;
  l = j - d2;

  if (traceback) fprintf(outf," (ZH10.2) IL, dangle l, trace vx %d %d %d \n", 
			 l+mid, mid-1, vx[l+mid][mid-1]); 
  if (traceback) fprintf(outf," (ZH10.2) IL, trace vhx %d %d %d %d %d \n", 
			 j, d, d1, d2-mid-1, vhx[j][d][d1][d2-mid-1]); 
  
  PushTraceknstack(dolist, AttachTracekn(curr_tr, i, j, k, l+mid+1, dpcP, dpcP));
  PushTraceknstack(dolist, AttachTracekn(curr_tr, l+1, l+mid, l+1+(int)((mid-1)/2), 
					 l+(int)((mid-1)/2)+2, dpcP, dpcS));
}

/* (ZH10.3)__/ k dangles off  l+mid+1,k-1 */
int 
ZH10_3(ESL_DSQ *s, int len, struct rnapar_2 *rnapar, int **icfg, int **vx,int ****vhx, int j, int d, int d1, int d2,int mid)
{ 
  int   ZH10_3;
  int i, k, l;

  i = j - d;
  k = i + d1;
  l = j - d2;
  
  if (d1 > 0)
    ZH10_3 = 2*rnapar->P10P + rnapar->P6P 
      + wkn*icfg[idxR(s[k])][idxP(s[l+mid+1],s[k-1])]
      + wkn*icfg[idxPS(s[l+mid],s[l])][idxPS(s[l+mid+1],s[k-1])]
      + vx[l+mid][mid] + vhx[j][d][d1-1][d2-mid-1];
  else
    ZH10_3 = -BIGINT;

  return ZH10_3;
}
void
trace_ZH10_3(FILE *outf, int **vx, int ****vhx, int j, int d, int d1, int d2, int mid,
	     struct tracekn_s *curr_tr, struct traceknstack_s *dolist, int traceback)
{
  int i, k, l;

  i = j - d;
  k = i + d1;
  l = j - d2;

  if (traceback) fprintf(outf," (ZH10.3) IL, trace vx %d %d %d \n", 
			 l+mid, mid, vx[l+mid][mid]); 
  if (traceback) fprintf(outf," (ZH10.3) IL, dangle k, trace vhx %d %d %d %d %d \n", 
			 j, d, d1-1, d2-mid-1, vhx[j][d][d1-1][d2-mid-1]); 
  
  PushTraceknstack(dolist, AttachTracekn(curr_tr, i, j, k-1, l+mid+1, dpcP, dpcP));
  PushTraceknstack(dolist, AttachTracekn(curr_tr, l, l+mid, l+(int)(mid/2), 
					 l+(int)(mid/2)+1, dpcP, dpcS));
}

/* (ZH10.4)__/ l dangles off  l+1,l+mid / k dangles off  l+mid+1,k-1*/
int 
ZH10_4(ESL_DSQ *s, int len, struct rnapar_2 *rnapar, int **icfg, int **vx,int ****vhx, int j, int d, int d1, int d2,int mid)
{ 
  int  ZH10_4;
  int i, k, l;

  i = j - d;
  k = i + d1;
  l = j - d2;
  
  if (mid > 0 && d1 > 0)
    ZH10_4 = 2*rnapar->P10P + 2*rnapar->P6P 
      + wkn*icfg[idxR(s[k])][idxP(s[l+mid+1],s[k-1])]
      + wkn*icfg[idxL(s[l])][idxP(s[l+1],s[l+mid])]
      + wkn*icfg[idxPS(s[l+mid],s[l+1])][idxPS(s[l+mid+1],s[k-1])]
      + vx[l+mid][mid-1] + vhx[j][d][d1-1][d2-mid-1];
  else
    ZH10_4 = -BIGINT;

  return ZH10_4;
}
void
trace_ZH10_4(FILE *outf, int **vx, int ****vhx, int j, int d, int d1, int d2, int mid,
	     struct tracekn_s *curr_tr, struct traceknstack_s *dolist, int traceback)
{
  int i, k, l;

  i = j - d;
  k = i + d1;
  l = j - d2;

  if (traceback) fprintf(outf," (ZH10.4) IL, dangle l, trace vx %d %d %d \n", 
			 l+mid, mid-1, vx[l+mid][mid-1]); 
  if (traceback) fprintf(outf," (ZH10.4) IL, dangle k, trace vhx %d %d %d %d %d \n", 
			 j, d, d1-1, d2-mid-1, vhx[j][d][d1-1][d2-mid-1]); 
  
  PushTraceknstack(dolist, AttachTracekn(curr_tr, i, j, k-1, l+mid+1, dpcP, dpcP));
  PushTraceknstack(dolist, AttachTracekn(curr_tr, l+1, l+mid, l+1+(int)((mid-1)/2), 
					 l+(int)((mid-1)/2)+2, dpcP, dpcS));
}

/* (ZH11)__/ 1 VHX 1 ZHX
 *
 * vhx connects i to j and mid1 to mid2. (j,d,mid1,mid2)
 *
 * zhx connects i+mid1 to j-mid2 and k to l. 
 *              (j-mid2,d-mid1-mid2,d1-mid1, d2-mid2)
 */
int 
ZH11(ESL_DSQ *s, int len, struct rnapar_2 *rnapar, int **icfg, int ****zhx,int ****vhx, 
     int j, int d, int d1, int d2, int mid1, int mid2)
{ 
  int    ZH11;
  int i, k, l;

  i = j - d;
  k = i + d1;
  l = j - d2;
  
  ZH11 = vhx[j][d][mid1][mid2]
    + zhx[j-mid2][d-mid1-mid2][d1-mid1][d2-mid2];

  return ZH11;
}
void
trace_ZH11(FILE *outf, int ****zhx, int ****vhx, int j, int d, int d1, int d2, int mid1, 
	   int mid2, struct tracekn_s *curr_tr, struct traceknstack_s *dolist, int traceback)
{
  int i, k, l;

  i = j - d;
  k = i + d1;
  l = j - d2;

  if (traceback) fprintf(outf," (ZH11) trace vhx  %d %d %d %d %d \n", 
			 j, d, mid1, mid2, 
			 vhx[j][d][mid1][mid2]); 
  if (traceback) fprintf(outf," (ZH11) trace zhx  %d %d %d %d %d \n", 
			 j-mid2, d-mid1-mid2-2, d1-mid1, d2-mid2, 
			 zhx[j-mid2][d-mid1-mid2][d1-mid1][d2-mid2]); 
  
  PushTraceknstack(dolist, AttachTracekn(curr_tr, i, j, i+mid1, j-mid2, 
					 dpcP, dpcP));
  PushTraceknstack(dolist, AttachTracekn(curr_tr, i+mid1, j-mid2, k, l, 
					 dpcP, dpcS));
}

/* (ZH12)__/ REST of ZHX (multiloops).  
 */
/* (ZH12.1)__/ no danglings */
int 
ZH12_1(ESL_DSQ *s, int len, struct rnapar_2 *rnapar, int **icfg, int ****whx, int j, int d, int d1, int d2)
{ 
  int ZH12_1;
  int i, k, l;

  i = j - d;
  k = i + d1;
  l = j - d2;
  
  if (d1 > 0 && d2 > 0)
    ZH12_1 = rnapar->P10P + rnapar->P5P 
      + whx[j-1][d-2][d1-1][d2-1];
  else
    ZH12_1 = -BIGINT;

  return ZH12_1;
}
void
trace_ZH12_1(FILE *outf, int ****whx, int j, int d, int d1, int d2,
	  struct tracekn_s *curr_tr, struct traceknstack_s *dolist, int traceback)
{
  int i, k, l;

  i = j - d;
  k = i + d1;
  l = j - d2;

  if (d1 > 0 && d2 > 0){
    if (traceback) fprintf(outf," (ZH12.1) multiloop, trace whx  %d %d %d %d %d \n", 
			   j-1, d-2, d1-1, d2-1, whx[j-1][d-2][d1-1][d2-1]);
    
    PushTraceknstack(dolist, AttachTracekn(curr_tr, i+1, j-1, k, l, dpcS, dpcS));
  }
}

/*  (ZH12.2)__/ i+1  dangles off j,i */
int 
ZH12_2(ESL_DSQ *s, int len, struct rnapar_2 *rnapar, int **icfg, int ****whx, int j, int d, int d1, int d2)
{ 
  int ZH12_2;
  int i, k, l;

  i = j - d;
  k = i + d1;
  l = j - d2;
  
  if(d1 > 1 && d2 > 0)
    ZH12_2 = rnapar->P10P + rnapar->P5P + rnapar->P6P 
      + wkn*icfg[idxR(s[i+1])][idxP(s[j],s[i])]
      + whx[j-1][d-3][d1-2][d2-1];
  else
    ZH12_2 = -BIGINT;
  
  return ZH12_2;
}
void
trace_ZH12_2(FILE *outf, int ****whx, int j, int d, int d1, int d2, 
	    struct tracekn_s *curr_tr, struct traceknstack_s *dolist, int traceback)
{
  int i, k, l;

  i = j - d;
  k = i + d1;
  l = j - d2;

  if(d1 > 1 && d2 > 0){
    if (traceback) fprintf(outf," (ZH12.2) multiloop, trace whx  %d %d %d %d %d \n", 
			   j-1, d-3, d1-2, d2-1, whx[j-1][d-3][d1-2][d2-1]);
    
    PushTraceknstack(dolist, AttachTracekn(curr_tr, i+2, j-1, k, l, dpcS, dpcS));
  }
}

/*  (ZH12.3)__/ j-1  dangles off i,j */
int 
ZH12_3(ESL_DSQ *s, int len, struct rnapar_2 *rnapar, int **icfg, int ****whx, int j, int d, int d1, int d2)
{ 
  int ZH12_3;
  int i, k, l;

  i = j - d;
  k = i + d1;
  l = j - d2;
  
  if(d1 > 0 && d2 > 1)
    ZH12_3 = rnapar->P10P + rnapar->P5P + rnapar->P6P 
      + wkn*icfg[idxL(s[j-1])][idxP(s[j],s[i])]
      + whx[j-2][d-3][d1-1][d2-2];
  else
    ZH12_3 = -BIGINT;
  
  return ZH12_3;
}
void
trace_ZH12_3(FILE *outf, int ****whx, int j, int d, int d1, int d2,
	    struct tracekn_s *curr_tr, struct traceknstack_s *dolist, int traceback)
{
  int i, k, l;

  i = j - d;
  k = i + d1;
  l = j - d2;

  if(d1 > 0 && d2 > 1){
    if (traceback) fprintf(outf," (ZH12.3) multiloop, trace whx  %d %d %d %d %d \n", 
			   j-2, d-3, d1-1, d2-2, whx[j-2][d-3][d1-1][d2-2]);
    
    PushTraceknstack(dolist, AttachTracekn(curr_tr, i+1, j-2, k, l, dpcS, dpcS));
  }
}

/*  (ZH12.4)__/ i+1 and j-1 dangle off i,j*/
int 
ZH12_4(ESL_DSQ *s, int len, struct rnapar_2 *rnapar, int **icfg, int ****whx, int j, int d, int d1, int d2)
{ 
  int ZH12_4;
  int i, k, l;

  i = j - d;
  k = i + d1;
  l = j - d2;
  
  if(d1 > 1 && d2 > 1)
    ZH12_4 = rnapar->P10P + rnapar->P5P + 2*rnapar->P6P 
      + wkn*icfg[idxR(s[i+1])][idxP(s[j],s[i])]
      + wkn*icfg[idxL(s[j-1])][idxP(s[j],s[i])]
      + whx[j-2][d-4][d1-2][d2-2];
  else
    ZH12_4 = -BIGINT;
  
  return ZH12_4;
}
void
trace_ZH12_4(FILE *outf, int ****whx, int j, int d, int d1, int d2, 
	    struct tracekn_s *curr_tr, struct traceknstack_s *dolist, int traceback)
{
  int i, k, l;

  i = j - d;
  k = i + d1;
  l = j - d2;

  if(d1 > 1 && d2 > 1){
    if (traceback) fprintf(outf," (ZH12.4) multiloop, trace whx  %d %d %d %d %d \n", 
			   j-2, d-4, d1-2, d2-2, whx[j-2][d-4][d1-2][d2-2]);
    
    PushTraceknstack(dolist, AttachTracekn(curr_tr, i+2, j-2, k, l, dpcS, dpcS));
  }
}


/* diagrams for YHX (YH1-YH12)
 */

/* (YH1)__/ pair k,l */
int 
YH1(ESL_DSQ *s, int len, struct rnapar_2 *rnapar, int **icfg, int ****vhx, int j, int d, int d1, int d2)
{ 
  int     YH1;
  int i, k, l;

  i = j - d;
  k = i + d1;
  l = j - d2;
  
    YH1 = rnapar->P10P + vhx[j][d][d1][d2];
  
  return YH1;
}
void
trace_YH1(FILE *outf, int ****vhx, int j, int d, int d1, int d2,
	  struct tracekn_s *curr_tr, struct traceknstack_s *dolist, int traceback)
{
  int i, k, l;

  i = j - d;
  k = i + d1;
  l = j - d2;

  if (traceback) fprintf(outf," (YH1) trace vhx %d %d %d %d %d \n", 
			 j, d, d1, d2, vhx[j][d][d1][d2]); 
  
  PushTraceknstack(dolist, AttachTracekn(curr_tr, i, j, k, l, dpcP, dpcP)); 
}

/* (YH2)__/ pair i+1,j-1  */
int 
YH2(ESL_DSQ *s, int len, struct rnapar_2 *rnapar, int **icfg, int ****vhx, int j, int d, int d1, int d2)
{ 
  int     YH2;
  int i, k, l;

  i = j - d;
  k = i + d1;
  l = j - d2;
  
  if (d1 > 0 && d2 > 0)
    YH2 = rnapar->P10P + 2*rnapar->P6P 
      + wkn*icfg[idxL(s[i])][idxP(s[i+1],s[j-1])]
      + wkn*icfg[idxR(s[j])][idxP(s[i+1],s[j-1])]
      + vhx[j-1][d-2][d1-1][d2-1];
  else
    YH2 = -BIGINT;

  return YH2;
}
void
trace_YH2(FILE *outf, int ****vhx, int j, int d, int d1, int d2,
	  struct tracekn_s *curr_tr, struct traceknstack_s *dolist, int traceback)
{
  int i, k, l;

  i = j - d;
  k = i + d1;
  l = j - d2;

  if (traceback) fprintf(outf," (YH2) SLE, SRE, trace vhx %d %d %d %d %d \n", 
			 j-1, d-2, d1-1, d2-1, vhx[j-1][d-2][d1-1][d2-1]); 
  
  curr_tr->type1 = dpcL;           
  PushTraceknstack(dolist, AttachTracekn(curr_tr, i+1, j-1, k, l, dpcP, dpcP)); 
}

/* (YH3)__/ i dangles off i+1,j base-paired */
int 
YH3(ESL_DSQ *s, int len, struct rnapar_2 *rnapar, int **icfg, int ****vhx, int j, int d, int d1, int d2)
{ 
  int     YH3;
  int i, k, l;

  i = j - d;
  k = i + d1;
  l = j - d2;
  
  if (d1 > 0)
    YH3 = rnapar->P10P + rnapar->P6P 
      + wkn*icfg[idxL(s[i])][idxP(s[i+1],s[j])] 
      + vhx[j][d-1][d1-1][d2];
  else
    YH3 = -BIGINT;

  return YH3;
}
void
trace_YH3(FILE *outf, int ****vhx, int j, int d, int d1, int d2,
	  struct tracekn_s *curr_tr, struct traceknstack_s *dolist, int traceback)
{
  int i, k, l;

  i = j - d;
  k = i + d1;
  l = j - d2;

  if (traceback) fprintf(outf," (YH3) SLE, trace vhx %d %d %d %d %d \n", 
			 j, d-1, d1-1, d2, vhx[j][d-1][d1-1][d2]); 
  
  curr_tr->type1 = dpcL;           
  PushTraceknstack(dolist, AttachTracekn(curr_tr, i+1, j, k, l, dpcP, dpcP)); 
}

/* (YH4)__/ j dangles off i,j-1. i,j-1  base-paired  */
int 
YH4(ESL_DSQ *s, int len, struct rnapar_2 *rnapar, int **icfg, int ****vhx, int j, int d, int d1, int d2)
{ 
  int     YH4;
  int i, k, l;

  i = j - d;
  k = i + d1;
  l = j - d2;
  
  if (d2 > 0)
    YH4 = rnapar->P10P + rnapar->P6P 
      + wkn*icfg[idxR(s[j])][idxP(s[i],s[j-1])]
      + vhx[j-1][d-1][d1][d2-1];
  else
    YH4 = -BIGINT;

  return YH4;
}
void
trace_YH4(FILE *outf, int ****vhx, int j, int d, int d1, int d2,
	  struct tracekn_s *curr_tr, struct traceknstack_s *dolist, int traceback)
{
  int i, k, l;

  i = j - d;
  k = i + d1;
  l = j - d2;

  if (traceback) fprintf(outf," (YH4) SRE, trace vhx %d %d %d %d %d \n", 
			 j-1, d-1, d1, d2-1, vhx[j-1][d-1][d1][d2-1]); 
  
  curr_tr->type1 = dpcL;           
  PushTraceknstack(dolist, AttachTracekn(curr_tr, i, j-1, k, l, dpcP, dpcP)); 
}

/* (YH5)__/ i dangles ss off i+1,j */
int 
YH5(ESL_DSQ *s, int len, struct rnapar_2 *rnapar, int **icfg, int ****yhx, int j, int d, int d1, int d2)
{ 
  int     YH5;
  int i, k, l;

  i = j - d;
  k = i + d1;
  l = j - d2;
  
  if (d1 > 0)
    YH5 = rnapar->P6P + yhx[j][d-1][d1-1][d2];
  else
    YH5 = -BIGINT;

  return YH5;
}
void
trace_YH5(FILE *outf, int ****yhx, int j, int d, int d1, int d2,
	  struct tracekn_s *curr_tr, struct traceknstack_s *dolist, int traceback)
{
  int i, k, l;

  i = j - d;
  k = i + d1;
  l = j - d2;

  if (traceback) fprintf(outf," (YH5) SLE, trace yhx %d %d %d %d %d \n", 
			 j, d-1, d1-1, d2, yhx[j][d-1][d1-1][d2]); 
  
  curr_tr->type1 = dpcL;           
  PushTraceknstack(dolist, AttachTracekn(curr_tr, i+1, j, k, l, dpcS, dpcP)); 
}

/* (YH6)__/ j dangles off i,j-1  */
int 
YH6(ESL_DSQ *s, int len, struct rnapar_2 *rnapar, int **icfg, int ****yhx, int j, int d, int d1, int d2)
{ 
  int     YH6;
  int i, k, l;

  i = j - d;
  k = i + d1;
  l = j - d2;
  
  if (d2 > 0)
    YH6 = rnapar->P6P + yhx[j-1][d-1][d1][d2-1];
  else
    YH6 = -BIGINT;

  return YH6;
}
void
trace_YH6(FILE *outf, int ****yhx, int j, int d, int d1, int d2,
	  struct tracekn_s *curr_tr, struct traceknstack_s *dolist, int traceback)
{
  int i, k, l;

  i = j - d;
  k = i + d1;
  l = j - d2;

  if (traceback) fprintf(outf," (YH6) SRE, trace yhx %d %d %d %d %d \n", 
			 j-1, d-1, d1, d2-1, yhx[j-1][d-1][d1][d2-1]); 
  
  curr_tr->type1 = dpcL;           
  PushTraceknstack(dolist, AttachTracekn(curr_tr, i, j-1, k, l, dpcS, dpcP)); 
}

/* 
 * STRUCTURES WITH ONE WBX AND ONE YHX. (YH7, YH9)
 * or         WITH ONE VX  AND ONE VHX. (YH8, YH10)
 */

/* In the LEFT ARM (i,k).
 */

/* EXTERIOR LEFT */

/* (YH7)__/ 1 WBX 1 YHX
 * WBX connects i with i+mid (i+mid,mid)
 * YHX connects (i+mid+1,k) and (l,j) (j,d-mid-1,d1-mid-1,d2).
 */                                       
int 
YH7(ESL_DSQ *s, int len, struct rnapar_2 *rnapar, int **icfg, int **wbx,int ****yhx, int j, int d, int d1, int d2,int mid)
{ 
  int     YH7;
  int i, k, l;

  i = j - d;
  k = i + d1;
  l = j - d2;
  
    YH7 = wbx[i+mid][mid] + yhx[j][d-mid-1][d1-mid-1][d2];

  return YH7;
}
void
trace_YH7(FILE *outf, int **wbx, int ****yhx, int j, int d, int d1, int d2, int mid,
	  struct tracekn_s *curr_tr, struct traceknstack_s *dolist, int traceback)
{
  int i, k, l;

  i = j - d;
  k = i + d1;
  l = j - d2;

  if (traceback) fprintf(outf," (YH7) EL wbx        %d %d %d \n", 
			 i+mid, mid, wbx[i+mid][mid]); 
  if (traceback) fprintf(outf," (YH7) EL yhx %d %d %d %d %d \n", j, d-mid-1, d1-mid-1, d2, 
			 yhx[j][d-mid-1][d1-mid-1][d2]); 
  
  PushTraceknstack(dolist, AttachTracekn(curr_tr, i+mid+1, j, k, l, dpcS, dpcP));
  PushTraceknstack(dolist, AttachTracekn(curr_tr, i, i+mid, i+(int)(mid/2), 
					 i+(int)(mid/2)+1, dpcB, dpcS));
}

/* (YH8)__/ 1 VX 1 ZHX 
 * VX connects i with i+mid (i+mid,mid)
 * VHX connects (i+mid+1,k) and (l,j) (j,d-mid-1,d1-mid-1,d2).
 * contiguos pairs (i, i+mid) (i+mid+1, j)
 * ....add stack[i+mid][i][i+mid+1][j]
 */                                       

/* (YH8.1)__/ no danglings / coaxial = stacking - 1 */
int 
YH8_1(ESL_DSQ *s, int len, struct rnapar_2 *rnapar, int **icfg, int **vx,int ****vhx, int j, int d, int d1, int d2,int mid)
{ 
  int   YH8_1;
  int i, k, l;

  i = j - d;
  k = i + d1;
  l = j - d2;
  
  YH8_1 = 2*rnapar->P10P 
    + wkn*(icfg[idxPS(s[i+mid],s[i])][idxPS(s[i+mid+1],s[j])] + INTSCALE)
    + vx[i+mid][mid] + vhx[j][d-mid-1][d1-mid-1][d2];
  
  return YH8_1;
}
void
trace_YH8_1(FILE *outf, int **vx, int ****vhx, int j, int d, int d1, int d2, int mid,
	    struct tracekn_s *curr_tr, struct traceknstack_s *dolist, int traceback)
{
  int i, k, l;

  i = j - d;
  k = i + d1;
  l = j - d2;

  if (traceback) fprintf(outf," (YH8.1) EL, trace vx  %d %d %d \n", 
			 i+mid, mid, vx[i+mid][mid]); 
  if (traceback) fprintf(outf," (YH8.1) EL, trace vhx %d %d %d %d %d \n", j, d-mid-1, d1-mid-1, d2, 
			 vhx[j][d-mid-1][d1-mid-1][d2]); 
  
  PushTraceknstack(dolist, AttachTracekn(curr_tr, i+mid+1, j, k, l, dpcP, dpcP));
  PushTraceknstack(dolist, AttachTracekn(curr_tr, i, i+mid, i+(int)(mid/2), 
					 i+(int)(mid/2)+1, dpcP, dpcS));
}

/* (YH8.2)__/ i dangles off i+1,i+mid */
int 
YH8_2(ESL_DSQ *s, int len, struct rnapar_2 *rnapar, int **icfg, int **vx,int ****vhx, int j, int d, int d1, int d2,int mid)
{ 
  int   YH8_2;
  int i, k, l;

  i = j - d;
  k = i + d1;
  l = j - d2;
  
  if (mid > 0)
    YH8_2 = 2*rnapar->P10P + rnapar->P6P
      + wkn*icfg[idxL(s[i])][idxP(s[i+1],s[i+mid])]
      + wkn*icfg[idxPS(s[i+mid],s[i+1])][idxPS(s[i+mid+1],s[j])]
      + vx[i+mid][mid-1] + vhx[j][d-mid-1][d1-mid-1][d2];
  else
    YH8_2 = -BIGINT;

  return YH8_2;
}
void
trace_YH8_2(FILE *outf, int **vx, int ****vhx, int j, int d, int d1, int d2, int mid,
	    struct tracekn_s *curr_tr, struct traceknstack_s *dolist, int traceback)
{
  int i, k, l;

  i = j - d;
  k = i + d1;
  l = j - d2;

  if (traceback) fprintf(outf," (YH8.2) EL, dangle i, trace vx  %d %d %d \n", 
			 i+mid, mid-1, vx[i+mid][mid-1]); 
  if (traceback) fprintf(outf," (YH8.2) EL, trace vhx %d %d %d %d %d \n", 
			 j, d-mid-1, d1-mid-1, d2, 
			 vhx[j][d-mid-1][d1-mid-1][d2]); 
  
  PushTraceknstack(dolist, AttachTracekn(curr_tr, i+mid+1, j, k, l, dpcP, dpcP));
  PushTraceknstack(dolist, AttachTracekn(curr_tr, i+1, i+mid, i+(int)((mid-1)/2), 
					 i+(int)((mid-1)/2)+1, dpcP, dpcS));
}

/* (YH8.3)__/ j dangles off i+mid+1,j-1 */
int 
YH8_3(ESL_DSQ *s, int len, struct rnapar_2 *rnapar, int **icfg, int **vx,int ****vhx, int j, int d, int d1, int d2,int mid)
{ 
  int   YH8_3;
  int i, k, l;

  i = j - d;
  k = i + d1;
  l = j - d2;
  
  if (d2 > 0)
    YH8_3 = 2*rnapar->P10P + rnapar->P6P
      + wkn*icfg[idxR(s[j])][idxP(s[i+mid+1],s[j-1])]
      + wkn*icfg[idxPS(s[i+mid],s[i])][idxPS(s[i+mid+1],s[j-1])]
      + vx[i+mid][mid] + vhx[j-1][d-mid-2][d1-mid-1][d2-1];
  else
    YH8_3 = -BIGINT;

  return YH8_3;
}
void
trace_YH8_3(FILE *outf, int **vx, int ****vhx, int j, int d, int d1, int d2, int mid,
	    struct tracekn_s *curr_tr, struct traceknstack_s *dolist, int traceback)
{
  int i, k, l;

  i = j - d;
  k = i + d1;
  l = j - d2;

  if (traceback) fprintf(outf," (YH8.3) EL, trace vx  %d %d %d \n", 
			 i+mid, mid, vx[i+mid][mid]); 
  if (traceback) fprintf(outf," (YH8.3) EL, dangle j, trace vhx %d %d %d %d %d \n",
			 j-1, d-mid-2, d1-mid-1, d2-1, 
			 vhx[j-1][d-mid-2][d1-mid-1][d2-1]); 
  
  PushTraceknstack(dolist, AttachTracekn(curr_tr, i+mid+1, j-1, k, l, dpcP, dpcP));
  PushTraceknstack(dolist, AttachTracekn(curr_tr, i, i+mid, i+(int)(mid/2), 
					 i+(int)(mid/2)+1, dpcP, dpcS));
}

/* (YH8.4)__/ i dangles off i+1,i+mid / j dangles off i+mid+1,j-1*/
int 
YH8_4(ESL_DSQ *s, int len, struct rnapar_2 *rnapar, int **icfg, int **vx,int ****vhx, int j, int d, int d1, int d2,int mid)
{ 
  int   YH8_4;
  int i, k, l;

  i = j - d;
  k = i + d1;
  l = j - d2;
  
  if (mid > 0 && d2 > 0)
    YH8_4 = 2*rnapar->P10P + 2*rnapar->P6P
      + wkn*icfg[idxL(s[i])][idxP(s[i+1],s[i+mid])]
      + wkn*icfg[idxR(s[j])][idxP(s[i+mid+1],s[j-1])]
      + wkn*icfg[idxPS(s[i+mid],s[i+1])][idxPS(s[i+mid+1],s[j-1])]
      + vx[i+mid][mid-1] + vhx[j-1][d-mid-2][d1-mid-1][d2-1];
  else
    YH8_4 = -BIGINT;

  return YH8_4;
}
void
trace_YH8_4(FILE *outf, int **vx, int ****vhx, int j, int d, int d1, int d2, int mid,
	    struct tracekn_s *curr_tr, struct traceknstack_s *dolist, int traceback)
{
  int i, k, l;

  i = j - d;
  k = i + d1;
  l = j - d2;

  if (traceback) fprintf(outf," (YH8.4) EL, dangle i, trace vx  %d %d %d \n", 
			 i+mid, mid-1, vx[i+mid][mid-1]); 
  if (traceback) fprintf(outf," (YH8.4) EL, dangle j, trace vhx %d %d %d %d %d \n", 
			 j-1, d-mid-2, d1-mid-1, d2-1, 
			 vhx[j][d-mid-1][d1-mid-1][d2]); 
  
  PushTraceknstack(dolist, AttachTracekn(curr_tr, i+mid+1, j-1, k, l, dpcP, dpcP));
  PushTraceknstack(dolist, AttachTracekn(curr_tr, i+1, i+mid, i+(int)((mid-1)/2), 
					 i+(int)((mid-1)/2)+1, dpcP, dpcS));
}

/* In the RIGHT ARM (l,j).
 */      
/* Exterior RIGHT */

/* (YH9)__/ 1 WBX 1 YHX
 * WBX connects j-mid with j  (j,mid)
 * YHX connects (i,k) and (l,j-mid-1) 
 *               (j-mid-1, d-mid-1, d1, d2-mid-1).
 */                                       
int 
YH9(ESL_DSQ *s, int len, struct rnapar_2 *rnapar, int **icfg, int **wbx,int ****yhx, int j, int d, int d1, int d2, int mid)
{ 
  int     YH9;
  int i, k, l;

  i = j - d;
  k = i + d1;
  l = j - d2;
  
  YH9 = wbx[j][mid] + yhx[j-mid-1][d-mid-1][d1][d2-mid-1];

  return YH9;
}
void
trace_YH9(FILE *outf, int **wbx, int ****yhx, int j, int d, int d1, int d2, int mid,
	  struct tracekn_s *curr_tr, struct traceknstack_s *dolist, int traceback)
{
  int i, k, l;

  i = j - d;
  k = i + d1;
  l = j - d2;

  if (traceback) fprintf(outf," (YH9) ER, trace wbx %d %d %d \n", 
			 j, mid, wbx[j][mid]); 
  if (traceback) fprintf(outf," (YH9) ER, trace yhx %d %d %d %d %d \n", 
			 j-mid-1, d-mid-1, d1, d2-mid-1,
			 yhx[j-mid-1][d-mid-1][d1][d2-mid-1]); 
  
  PushTraceknstack(dolist, AttachTracekn(curr_tr, j-mid, j, j-mid+(int)(mid/2), 
					 j-mid+(int)(mid/2)+1, dpcB, dpcS));
  PushTraceknstack(dolist, AttachTracekn(curr_tr, i, j-mid-1, k, l, dpcS, dpcP));
}

/* (YH10)__/ 1 VX 1 VHX 
 * VX  connects j-mid with j  (j,mid)
 * VHX connects (i,k) and (l,j-mid-1) 
 *              (j-mid-1, d-mid-1, d1, d2-mid-1).
 * coaxial pairs (i,j-mid-1) (j-mid,j)
 * ....add stack[j-mid-1][i][j-mid][j]
 */                                       
/* (YH10.1)__/ no danglings / coaxial = stacking + 1 */
int 
YH10_1(ESL_DSQ *s, int len, struct rnapar_2 *rnapar, int **icfg, int **vx,int ****vhx, int j, int d, int d1, int d2,int mid)
{ 
  int  YH10_1;
  int i, k, l;

  i = j - d;
  k = i + d1;
  l = j - d2;
  
  YH10_1 = 2*rnapar->P10P 
    + wkn*(icfg[idxPS(s[j-mid-1],s[i])]
	   [idxPS(s[j-mid],s[j])] + INTSCALE)
    + vx[j][mid] + vhx[j-mid-1][d-mid-1][d1][d2-mid-1];
  
  return YH10_1;
}
void
trace_YH10_1(FILE *outf, int **vx, int ****vhx, int j, int d, int d1, int d2, int mid,
	     struct tracekn_s *curr_tr, struct traceknstack_s *dolist, int traceback)
{
  int i, k, l;

  i = j - d;
  k = i + d1;
  l = j - d2;
  
  if (traceback) fprintf(outf," (YH10.1) ER, trace vx  %d %d %d \n", 
			 j, mid, vx[j][mid]); 
  if (traceback) fprintf(outf," (YH10.1) ER, trace vhx %d %d %d %d %d \n", 
			 j-mid-1, d-mid-1, d1, d2-mid-1,
			 vhx[j-mid-1][d-mid-1][d1][d2-mid-1]); 
  
  PushTraceknstack(dolist, AttachTracekn(curr_tr, j-mid, j, j-mid+(int)(mid/2), 
					 j-mid+(int)(mid/2)+1, dpcP, dpcS));
  PushTraceknstack(dolist, AttachTracekn(curr_tr, i, j-mid-1, k, l, dpcP, dpcP));
}

/* (YH10.2)__/ i dangles off i+1,j-mid-1 */
int 
YH10_2(ESL_DSQ *s, int len, struct rnapar_2 *rnapar, int **icfg, int **vx,int ****vhx, int j, int d, int d1, int d2,int mid)
{ 
  int   YH10_2;
  int i, k, l;

  i = j - d;
  k = i + d1;
  l = j - d2;
  
  if (d1 > 0)
    YH10_2 = 2*rnapar->P10P + rnapar->P6P
      + wkn*icfg[idxL(s[i])][idxP(s[i+1],s[j-mid-1])]
      + wkn*icfg[idxPS(s[j-mid-1],s[i+1])][idxPS(s[j-mid],s[j])]
      + vx[j][mid] + vhx[j-mid-1][d-mid-2][d1-1][d2-mid-1];
  else
    YH10_2 = -BIGINT;

  return YH10_2;
}
void
trace_YH10_2(FILE *outf, int **vx, int ****vhx, int j, int d, int d1, int d2, int mid,
	     struct tracekn_s *curr_tr, struct traceknstack_s *dolist, int traceback)
{
  int i, k, l;

  i = j - d;
  k = i + d1;
  l = j - d2;

  if (traceback) fprintf(outf," (YH10.2) ER, trace vx  %d %d %d \n", 
			 j, mid, vx[j][mid]); 
  if (traceback) fprintf(outf," (YH10.2) ER, dangle i, trace vhx %d %d %d %d %d \n", 
			 j-mid-1, d-mid-2, d1-1, d2-mid-1,
			 vhx[j-mid-1][d-mid-2][d1-1][d2-mid-1]); 
  
  PushTraceknstack(dolist, AttachTracekn(curr_tr, j-mid, j, j-mid+(int)(mid/2), 
					 j-mid+(int)(mid/2)+1, dpcP, dpcS));
  PushTraceknstack(dolist, AttachTracekn(curr_tr, i+1, j-mid-1, k, l, dpcP, dpcP));
}

/* (YH10.3)__/ j dangles off j-mid, j-1 */
int 
YH10_3(ESL_DSQ *s, int len, struct rnapar_2 *rnapar, int **icfg, int **vx,int ****vhx, int j, int d, int d1, int d2,int mid)
{ 
  int   YH10_3;
  int i, k, l;

  i = j - d;
  k = i + d1;
  l = j - d2;
  
  if (mid > 0)
    YH10_3 = 2*rnapar->P10P + rnapar->P6P
      + wkn*icfg[idxR(s[j])][idxP(s[j-mid],s[j-1])]
      + wkn*icfg[idxPS(s[j-mid-1],s[i])][idxPS(s[j-mid],s[j-1])]
      + vx[j-1][mid-1] + vhx[j-mid-1][d-mid-1][d1][d2-mid-1];
  else
    YH10_3 = -BIGINT;

  return YH10_3;
}
void
trace_YH10_3(FILE *outf, int **vx, int ****vhx, int j, int d, int d1, int d2, int mid,
	     struct tracekn_s *curr_tr, struct traceknstack_s *dolist, int traceback)
{
  int i, k, l;

  i = j - d;
  k = i + d1;
  l = j - d2;

  if (traceback) fprintf(outf," (YH10.3) ER, dangle j, trace vx  %d %d %d \n", 
			 j-1, mid-1, vx[j-1][mid-1]); 
  if (traceback) fprintf(outf," (YH10.3) ER, trace vhx %d %d %d %d %d \n", 
			 j-mid-1, d-mid-1, d1, d2-mid-1,
			 vhx[j-mid-1][d-mid-1][d1][d2-mid-1]); 
  
  PushTraceknstack(dolist, AttachTracekn(curr_tr, j-mid, j-1, j-mid+(int)((mid-1)/2), 
					 j-mid+(int)((mid-1)/2)+1, dpcP, dpcS));
  PushTraceknstack(dolist, AttachTracekn(curr_tr, i, j-mid-1, k, l, dpcP, dpcP));
}

/* (YH10.4)__/ i dangles off i+1,j-mid-1 / j dangles off j-mid, j-1*/
int 
YH10_4(ESL_DSQ *s, int len, struct rnapar_2 *rnapar, int **icfg, int **vx,int ****vhx, int j, int d, int d1, int d2,int mid)
{ 
  int  YH10_4;
  int i, k, l;

  i = j - d;
  k = i + d1;
  l = j - d2;
  
  if (d1 > 0 && mid > 0)
    YH10_4 = 2*rnapar->P10P + 2*rnapar->P6P
      + wkn*icfg[idxL(s[i])][idxP(s[i+1],s[j-mid-1])]
      + wkn*icfg[idxR(s[j])][idxP(s[j-mid],s[j-1])]
      + wkn*icfg[idxPS(s[j-mid-1],s[i+1])][idxPS(s[j-mid],s[j-1])]
      + vx[j-1][mid-1] + vhx[j-mid-1][d-mid-2][d1-1][d2-mid-1];
  else
    YH10_4 = -BIGINT;

  return YH10_4;
}
void
trace_YH10_4(FILE *outf, int **vx, int ****vhx, int j, int d, int d1, int d2, int mid,
	     struct tracekn_s *curr_tr, struct traceknstack_s *dolist, int traceback)
{
  int i, k, l;

  i = j - d;
  k = i + d1;
  l = j - d2;

  if (traceback) fprintf(outf," (YH10.4) ER, dangle j, trace vx  %d %d %d \n", 
			 j-1, mid-1, vx[j-1][mid-1]); 
  if (traceback) fprintf(outf," (YH10.4) ER, dangle i, trace vhx %d %d %d %d %d \n", 
			 j-mid-1, d-mid-2, d1-1, d2-mid-1,
			 vhx[j-mid-1][d-mid-2][d1-1][d2-mid-1]); 
  
  PushTraceknstack(dolist, AttachTracekn(curr_tr, j-mid, j-1, j-mid+(int)((mid-1)/2), 
					 j-mid+(int)((mid-1)/2)+1, dpcP, dpcS));
  PushTraceknstack(dolist, AttachTracekn(curr_tr, i+1, j-mid-1, k, l, dpcP, dpcP));
}

/* (YH11)__/ 1 VHX 1 YHX
 *
 * vhx connects k-mid1 to l+mid2 and k to l. 
 *              (l+mid2,d-d1-d2+mid1+mid2,mid1,mid2)
 *
 * yhx connects i to j and k-mid1 to l+mid2. 
 *              (j,d,d1-mid1, d2-mid2)
 */
int 
YH11(ESL_DSQ *s, int len, struct rnapar_2 *rnapar, int **icfg, int ****yhx,int ****vhx, 
     int j, int d, int d1, int d2, int mid1, int mid2)
{ 
  int    YH11;
  int i, k, l;

  i = j - d;
  k = i + d1;
  l = j - d2;
  
  YH11 = vhx[l+mid2][d-d1-d2+mid1+mid2][mid1][mid2]
    + yhx[j][d][d1-mid1][d2-mid2];
  
  return YH11;
}
void
trace_YH11(FILE *outf, int ****yhx, int ****vhx, int j, int d, int d1, int d2, int mid1, 
	   int mid2, struct tracekn_s *curr_tr, struct traceknstack_s *dolist, int traceback)
{
  int i, k, l;

  i = j - d;
  k = i + d1;
  l = j - d2;

  if (traceback) fprintf(outf," (YH11) trace yhx  %d %d %d %d %d \n", 
			 j, d, d1-mid1, d2-mid2, 
			 yhx[j][d][d1-mid1][d2-mid2]); 
  if (traceback) fprintf(outf," (YH11) trace vhx  %d %d %d %d %d \n", 
			 l+mid2, d-d1-d2+mid1+mid2, mid1, mid2, 
			 vhx[l+mid2][d-d1-d2+mid1+mid2][mid1][mid2]); 
  
  PushTraceknstack(dolist, AttachTracekn(curr_tr, i, j, k-mid1, l+mid2, 
					 dpcS, dpcP));
  PushTraceknstack(dolist, AttachTracekn(curr_tr, k-mid1, l+mid2, k, l, 
					 dpcP, dpcP));
}

/* (YH12)__/ REST of YHX (multiloops).  
 */
/* (YH12.1)__/ no danglings */
int 
YH12_1(ESL_DSQ *s, int len, struct rnapar_2 *rnapar, int **icfg, int ****whx, int j, int d, int d1, int d2)
{ 
  int YH12_1;
  int i, k, l;

  i = j - d;
  k = i + d1;
  l = j - d2;
  
  if (d1 > 0 && d2 > 0)
    YH12_1 = rnapar->P10P + rnapar->P5P 
      + whx[j][d][d1-1][d2-1];
  else
    YH12_1 = -BIGINT;

  return YH12_1;
}
void
trace_YH12_1(FILE *outf, int ****whx, int j, int d, int d1, int d2,
	  struct tracekn_s *curr_tr, struct traceknstack_s *dolist, int traceback)
{
  int i, k, l;

  i = j - d;
  k = i + d1;
  l = j - d2;

  if (d1 > 0 && d2 > 0){
    if (traceback) fprintf(outf," (YH12.1) multiloop, trace whx  %d %d %d %d %d \n", 
			   j, d, d1-1, d2-1, whx[j][d][d1-1][d2-1]);
    
    PushTraceknstack(dolist, AttachTracekn(curr_tr, i, j, k-1, l+1, dpcS, dpcS));
  }
}

/*  (YH12.2)__/ k-1  dangles off l,k */
int 
YH12_2(ESL_DSQ *s, int len, struct rnapar_2 *rnapar, int **icfg, int ****whx, int j, int d, int d1, int d2)
{ 
  int YH12_2;
  int i, k, l;

  i = j - d;
  k = i + d1;
  l = j - d2;
  
  if(d1 > 1 && d2 > 0)
    YH12_2 = rnapar->P10P + rnapar->P5P + rnapar->P6P 
      + wkn*icfg[idxL(s[k-1])][idxP(s[l],s[k])]
      + whx[j][d][d1-2][d2-1];
  else
    YH12_2 = -BIGINT;
  
  return YH12_2;
}
void
trace_YH12_2(FILE *outf, int ****whx, int j, int d, int d1, int d2, 
	    struct tracekn_s *curr_tr, struct traceknstack_s *dolist, int traceback)
{
  int i, k, l;

  i = j - d;
  k = i + d1;
  l = j - d2;

  if(d1 > 1 && d2 > 0){
    if (traceback) fprintf(outf," (YH12.2) multiloop, trace whx  %d %d %d %d %d \n", 
			   j, d, d1-2, d2-1, whx[j][d][d1-2][d2-1]);
    
    PushTraceknstack(dolist, AttachTracekn(curr_tr, i, j, k-2, l+1, dpcS, dpcS));
  }
}

/*  (YH12.3)__/ l+1  dangles off l,k */
int 
YH12_3(ESL_DSQ *s, int len, struct rnapar_2 *rnapar, int **icfg, int ****whx, int j, int d, int d1, int d2)
{ 
  int YH12_3;
  int i, k, l;

  i = j - d;
  k = i + d1;
  l = j - d2;
  
  if(d1 > 0 && d2 > 1)
    YH12_3 = rnapar->P10P + rnapar->P5P + rnapar->P6P 
      + wkn*icfg[idxR(s[l+1])][idxP(s[l],s[k])] 
      + whx[j][d][d1-1][d2-2];
  else
    YH12_3 = -BIGINT;
  
  return YH12_3;
}
void
trace_YH12_3(FILE *outf, int ****whx, int j, int d, int d1, int d2,
	    struct tracekn_s *curr_tr, struct traceknstack_s *dolist, int traceback)
{
  int i, k, l;

  i = j - d;
  k = i + d1;
  l = j - d2;

  if(d1 > 0 && d2 > 1){
    if (traceback) fprintf(outf," (YH12.3) multiloop, trace whx  %d %d %d %d %d \n", 
			   j, d, d1-1, d2-2, whx[j][d][d1-1][d2-2]);
    
    PushTraceknstack(dolist, AttachTracekn(curr_tr, i, j, k-1, l+2, dpcS, dpcS));
  }
}

/*  (YH12.4)__/ l+1 and k-1 dangle off l,k */
int 
YH12_4(ESL_DSQ *s, int len, struct rnapar_2 *rnapar, int **icfg, int ****whx, int j, int d, int d1, int d2)
{ 
  int YH12_4;
  int i, k, l;

  i = j - d;
  k = i + d1;
  l = j - d2;
  
  if(d1 > 1 && d2 > 1)
    YH12_4 = rnapar->P10P + rnapar->P5P + 2*rnapar->P6P 
      + wkn*icfg[idxL(s[k-1])][idxP(s[l],s[k])]
      + wkn*icfg[idxR(s[l+1])][idxP(s[l],s[k])] 
      + whx[j][d][d1-2][d2-2];
  else
    YH12_4 = -BIGINT;
  
  return YH12_4;
}
void
trace_YH12_4(FILE *outf, int ****whx, int j, int d, int d1, int d2, 
	    struct tracekn_s *curr_tr, struct traceknstack_s *dolist, int traceback)
{
  int i, k, l;

  i = j - d;
  k = i + d1;
  l = j - d2;

  if(d1 > 1 && d2 > 1){
    if (traceback) fprintf(outf," (YH12.4) multiloop, trace whx  %d %d %d %d %d \n", 
			   j, d, d1-2, d2-2, whx[j][d][d1-2][d2-2]);
    
    PushTraceknstack(dolist, AttachTracekn(curr_tr, i, j, k-2, l+2, dpcS, dpcS));
  }
}






