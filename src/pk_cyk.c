/* viterbi.c (mdified)
 * 
 * includes knots at second order in the zuker approximation
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>

#include <easel.h>
#include <esl_sq.h>
#include <esl_sqio.h>

#include "pknots.h"

#include "pk_cyk.h"
#include "pk_fillmtx.h"
#include "pk_tracemtx.h"
#include "pk_rnaparam.h"
#include "pk_vxgraphs.h"
#include "pk_wbxgraphs.h"
#include "pk_wxgraphs.h"
#include "pk_util.h"


/* Function: StructurePredictkn_2IS()
 * 
 * Purpose:  Given a sequence and an RNA structure SCFG, predict
 *           the secondary structure (best SCFG parse) of the 
 *           sequence. Uses a single-best-path (Viterbi) dynamic
 *           programming algorithm.
 *           
 * Args:     seq         - sequence to fold. 0..len-1;
 *                         allowable characters are A,C,G,U upper case 
 *           len         - length of sequence
 *           icfg        - grammar for RNA structure, integer log form
 *           verbose     - print debug info if TRUE.
 *           ret_trace   - RETURN: traceback of best alignment 
 *           ret_score   - RETURN: score of best alignment
 *         
 * Return:   (void)
 */          
int
StructurePredictkn_2IS(FILE *outf, ESL_ALPHABET *abc, ESL_SQ *sq, int len, struct rnapar_2 *rnapar, int **icfg,  
		       int verbose, int traceback, struct tracekn_s **ret_tr, 
		       int *ret_score, int allow_coaxials, int allow_pseudoknots, int approx)
{
  struct tracekn_s *tr;          /* traceback of optimal structure              */
  int     ****whx;               /* One-hole matrix, 4D (len x len x len x len) */
  int     ****vhx;               /* One-hole matrix, 4D (len x len x len x len) */
  int     ****zhx;               /* One-hole matrix, 4D (len x len x len x len) */
  int     ****yhx;               /* One-hole matrix, 4D (len x len x len x len) */
  int     **wx;                  /*  no-hole matrix, 2D (len x len)             */
  int     **wbx;                 /*  no-hole matrix, 2D (len x len)             */
  int     **vx;                  /*  no-hole matrix, 2D (len x len)             */
  int     *vp;                   /*  no-hole matrix, for internal loops         */
  ESL_DSQ *dsq;
  int     i;
  int     status;

  ESL_ALLOC(sq->dsq, sizeof(char) * (len+2));
  if (esl_abc_Digitize(abc, sq->seq, sq->dsq) != eslOK)
    pk_fatal("failed to digitize sequence");
     
  /* the dsq goes from 1..L, pknots arrays traditionally go from 0...L-1 */
  ESL_ALLOC(dsq, sizeof(char) * (len));
  for (i = 1; i <= len; i ++) 
    dsq[i-1] = sq->dsq[i];

  Alloc_Mtx(len, &wx, &wbx, &vx, &vp);
  Pattern_Mtx(len, wx, wbx, vx, vp); /* for debugging */

 if (allow_pseudoknots) {
    Alloc_Mgp(len, &whx, &vhx, &zhx, &yhx);
    Pattern_Mgp(len, whx, vhx, zhx, yhx); /* for debugging */
 
    FillMtx(dsq, len, rnapar, icfg,  wx, wbx, vx, vp, whx, vhx, zhx, yhx, allow_coaxials, approx);
  }
  else
    FillMtx_nested(dsq, len, rnapar, icfg,  wx, wbx, vx, vp, allow_coaxials);
  
  if (allow_pseudoknots) 
    TraceMtx(outf, dsq, len, rnapar, icfg, wx, wbx, vx,  whx, vhx, zhx, yhx,
	     len-1, len-1, (int)((len-1)/2), len-(int)((len-1)/2)-2,  
	     approx, &tr, traceback);
 else
    TraceMtx_nested(outf, dsq, len, rnapar, icfg, wx, wbx, vx, len-1, len-1, &tr, traceback);
    
  *ret_score = (float)(wx[len-1][len-1]);
  *ret_tr    = tr;
  
  Free_Mtx(len, wx, wbx, vx, vp);
  if (allow_pseudoknots)
    Free_Mgp(len, whx, vhx, zhx, yhx);

  free(dsq);
  return eslOK;

 ERROR:
  free(dsq);
  return status;
}
                       

/* Function: Alloc_Mtx()
 * 
 * Purpose:  Malloc space for the score matrices wx, wbx, vx.
 *           wx, wbx, vx are indexed as j, d.      
 *           j  = 0,...,len-1
 *           d  = 0,...,j
 *           
 *           i = j - d
 *           
 * Args:     len     - length of sequence
 *           ret_wx  - RETURN: score matrix wx
 *           ret_wbx - RETURN: score matrix wbx
 *           ret_vx  - RETURN: score matrix vx
 * 
 * Return:   Ptr to allocated scoring matrix, or
 *           dies and exits.
 */
void
Alloc_Mtx(int len,  int ***ret_wx, int ***ret_wbx, int ***ret_vx, int **ret_vp)
{
  int **wx, **wbx, **vx, *vp;
  int j;
                                                     
  /* no-hole matrices, wx (j,d), wbx (j,d), vx (j,d).
   *  fastest varying index is j.
   */
  /* This way of alloc'ing a 2D array keeps the CFG all in one
   * contiguous chunk of RAM and might keep us in cache.
   */
  
  if ((wx  = (int **) malloc (sizeof(int *) * len)) == NULL)
    pk_fatal("malloc failed in wx");
  if ((wx[0]  = (int *) malloc (sizeof(int) * Dim2len(len))) == NULL)
    pk_fatal("malloc failed in wx[0]");
  
  if ((wbx = (int **) malloc (sizeof(int *) * len)) == NULL)
    pk_fatal("malloc failed in wbx");
  if( (wbx[0] = (int *) malloc (sizeof(int) * Dim2len(len))) == NULL)
    pk_fatal("malloc failed in wbx[0]");
  
  if ((vx  = (int **) malloc (sizeof(int *) * len)) == NULL)
    pk_fatal("malloc failed in vx");
  if ((vx[0]  = (int *) malloc (sizeof(int) * Dim2len(len))) == NULL)
    pk_fatal("malloc failed in vx[0]");
  
  if ((vp = (int *) malloc (sizeof(int) * len)) == NULL)
    pk_fatal("malloc failed in vp");

  for (j = 1; j < len; j++)
    {
      vx[j]  =  vx[0]  + Dim2len(j);  
      wx[j]  =  wx[0]  + Dim2len(j);  
      wbx[j] =  wbx[0] + Dim2len(j);  
    }
  
  *ret_wx  = wx;
  *ret_wbx = wbx;
  *ret_vx  = vx;
  *ret_vp  = vp;
}


/* Function: Pattern_Mtx()
 * 
 * Purpose:  For debugging: write a pattern into the matrix
 *           so I can tell what's been written to.
 */
void
Pattern_Mtx(int len, int **wx, int **wbx, int **vx, int *vp)
{
  int j, d;
  
  for (j = 0; j < len; j++) {

    vp[j] = -123456789;
 
    for (d = 0; d <= j ; d++)
      {
	wbx[j][d] = -123456789; 
	wx[j][d] = -123456789; 
	vx[j][d] = -123456789; 
      }   
  }
}

void
Free_Mtx(int len, int **wx, int **wbx, int **vx, int *vp)
{
 /* Free matrices: wx, wbx and vx.
  */
 free( wx[0]);
 free(wbx[0]);
 free( vx[0]);
 
 free( wx);
 free(wbx);
 free( vx);
 free( vp);
}


/* Function: Alloc_Mgp()
 * 
 * Purpose:  Malloc space for the score matrices whx, vhx, zhx, yz.
 *           whx, vhx, zhx, yhx are indexed as j, d, d1, d2.
 *           j  = 0,...,len-1
 *           d  = 0,...,j
 *           d1 = 0,...,d
 *           d2 = 0,   ,d-d1
 *           
 *           i = j - d
 *           i = k - d1
 *           l = j - d2
 *           
 * Args:     len     - length of sequence
 *           ret_whx - RETURN: one-hole score matrix whx
 *           ret_vhx - RETURN: one-hole score matrix vhx
 *           ret_zhx - RETURN: one-hole score matrix zhx
 *           ret_yhx - RETURN: one-hole score matrix yhx
 * 
 * Return:   Ptr to allocated scoring matrix, or
 *           dies and exits.
 */
void
Alloc_Mgp(int len, int *****ret_whx, int *****ret_vhx, 
		int *****ret_zhx, int *****ret_yhx)
{
  int ****whx, ****vhx, ****zhx, ****yhx;
  int j, d, d1;
  
  /* gap matrices, whx(j,d,d1,d2), vhx(j,d,d1,d2), zhx(j,d,d1,d2), yhx(j,d,d1,d2).
   * fastest varying index is j.
   */
 if ((whx = (int ****) malloc (sizeof(int ***) * len)) == NULL)
    pk_fatal("malloc failed in whx");
  if ((whx[0] = (int ***) malloc (sizeof(int **) * Dim2len(len))) == NULL)
    pk_fatal("malloc failed in whx[0]");
  if ((whx[0][0] = (int **) malloc (sizeof(int *) * Dim3(len))) == NULL)
    pk_fatal("malloc failed in whx[0][0]");
  if ((whx[0][0][0] = (int *) malloc (sizeof(int) * Dim4(len))) == NULL)
    pk_fatal("malloc failed in whx[0][0][0]");
  
  if ((vhx = (int ****) malloc (sizeof(int ***) * len)) == NULL)
    pk_fatal("malloc failed in vhx");
  if ((vhx[0] = (int ***) malloc (sizeof(int **) * Dim2len(len))) == NULL)
    pk_fatal("malloc failed in vhx[0]");
  if ((vhx[0][0] = (int **) malloc (sizeof(int *) * Dim3(len))) == NULL)
    pk_fatal("malloc failed in vhx[0][0]");
  if ((vhx[0][0][0] = (int *) malloc (sizeof(int) * Dim4(len))) == NULL)
    pk_fatal("malloc failed in vhx[0][0][0]");
  
  if ((zhx = (int ****) malloc (sizeof(int ***) * len)) == NULL)
    pk_fatal("malloc failed in zhx");
  if ((zhx[0] = (int ***) malloc (sizeof(int **) * Dim2len(len))) == NULL)
    pk_fatal("malloc failed in zhx[0]");
  if ((zhx[0][0] = (int **) malloc (sizeof(int *) * Dim3(len))) == NULL)
    pk_fatal("malloc failed in zhx[0][0]");
  if ((zhx[0][0][0] = (int *) malloc (sizeof(int) * Dim4(len))) == NULL)
    pk_fatal("malloc failed in zhx[0][0][0]");
  
  if ((yhx = (int ****) malloc (sizeof(int ***) * len)) == NULL)
    pk_fatal("malloc failed in yhx");
  if ((yhx[0] = (int ***) malloc (sizeof(int **) * Dim2len(len))) == NULL)
    pk_fatal("malloc failed in yhx[0]");
  if ((yhx[0][0] = (int **) malloc (sizeof(int *) * Dim3(len))) == NULL)
    pk_fatal("malloc failed in yhx[0][0]");
  if ((yhx[0][0][0] = (int *) malloc (sizeof(int) * Dim4(len))) == NULL)
    pk_fatal("malloc failed in yhx[0][0][0]");
  
  for (j = 0; j < len; j++)
    {
      whx[j]  =  whx[0] + Dim2len(j);                          
      vhx[j]  =  vhx[0] + Dim2len(j);                          
      zhx[j]  =  zhx[0] + Dim2len(j);                          
      yhx[j]  =  yhx[0] + Dim2len(j);   
      
      for (d = 0; d <= j; d++)
        {
	  whx[j][d]  =  whx[0][0] + Dim3(j) + Dim2j(d-1);                          
	  vhx[j][d]  =  vhx[0][0] + Dim3(j) + Dim2j(d-1);                          
	  zhx[j][d]  =  zhx[0][0] + Dim3(j) + Dim2j(d-1);                          
	  yhx[j][d]  =  yhx[0][0] + Dim3(j) + Dim2j(d-1);   
	  
	  for (d1 = 0; d1 <= d; d1++)
	    {   
	      whx[j][d][d1]  =  whx[0][0][0] + Dim4(j) + Dim3j(d-1) + Dim2d(d1-1,d);              
	      vhx[j][d][d1]  =  vhx[0][0][0] + Dim4(j) + Dim3j(d-1) + Dim2d(d1-1,d);              
	      zhx[j][d][d1]  =  zhx[0][0][0] + Dim4(j) + Dim3j(d-1) + Dim2d(d1-1,d);            
	      yhx[j][d][d1]  =  yhx[0][0][0] + Dim4(j) + Dim3j(d-1) + Dim2d(d1-1,d);  
	    }
	}
    }

  *ret_whx = whx;
  *ret_vhx = vhx;
  *ret_zhx = zhx;
  *ret_yhx = yhx;
}


/* Function: Pattern_Mgp()
 * 
 * Purpose:  For debugging: write a pattern into the matrix
 *           so I can tell what's been written to.
 */
void
Pattern_Mgp(int len, int ****whx, int ****vhx, int ****zhx, int ****yhx)
{
  int j, d, d1, d2;
 
  for (j = 0; j < len; j++)
    for (d = 0; d <= j ; d++)
      for (d1 = 0; d1 <= d; d1++) 
	for (d2 = 0; d2 <= (d-d1); d2++)
	  {
 	    vhx[j][d][d1][d2] = -123456789; 
	    whx[j][d][d1][d2] = -123456789; 
	    zhx[j][d][d1][d2] = -123456789; 
	    yhx[j][d][d1][d2] = -123456789; 
	  }
}   

/* Function: Free_Mgp()
 * 
 * Purpose:  Free the space allocated to the gap scoring matrices.
 *           Precisely mirrors the allocations above in allocate_mgp().
 * 
 * Return:   (void)
 */
void
Free_Mgp(int len, int ****whx, int ****vhx, int ****zhx, int ****yhx)
{
 /* Free matrices: whx, vhx, zhx, yhx.
  */
 free(whx[0][0][0]);
 free(vhx[0][0][0]);
 free(zhx[0][0][0]);
 free(yhx[0][0][0]);
 
 free(whx[0][0]);
 free(vhx[0][0]);
 free(zhx[0][0]);
 free(yhx[0][0]);
 
 free(whx[0]);
 free(vhx[0]);
 free(zhx[0]);
 free(yhx[0]);
 
 free(whx);
 free(vhx);
 free(zhx);
 free(yhx);
}


/* Function: Print_2DMtx()
 * 
 * Purpose:  For debugging: print a matrix out.
 *
 * Args:     fp   - open FILE ptr to write to (stdout, perhaps)
 *           mtx   - DP matrix, len x len 
 *           len  - width/height of wx
 *
 * Return:   (void)                  
 */
void
Print_2DMtx(FILE *fp, int **mtx, int len, struct rnapar_2 *rnapar)
{
  int j, d;

  fprintf(fp, "### Here's the 2D matrix:\n");
  for (j = 0; j < len; j++)
    {
      fprintf(fp, "%d   \n", j);
      for (d = 0; d <= j; d++)
             if (mtx[j][d] == -123456789) /* pattern */
               fprintf(fp, "%9s ", "..");
        else if (mtx[j][d] <= (d+1) * rnapar->P6)  /* prohibited */
               fprintf(fp, "%9s ", "**");
        else
               fprintf(fp, "%9d ", mtx[j][d]);
      fputs("\n", fp);
     }
}
/* Function: Print_4DMtx()
 * 
 * Purpose:  For debugging: print a matrix out.
 *
 * Args:     fp   - open FILE ptr to write to (stdout, perhaps)
 *           whx   - DP matrix, len x len x len x len
 *           len  - width/height of whx
 *
 * Return:   (void)                  
 */
void
Print_4DMtx(FILE *fp, int ****mtx, int len, struct rnapar_2 *rnapar)
{
  int j, d, d1, d2;

  fprintf(fp, "### Here's the 4D matrix:\n");
  for (j = 0; j < len; j++)
    for (d = 0; d <= j; d++)
      for (d1 = 0; d1 <= d; d1++)
        {
          fprintf(fp, "%d %d %d  \n", j,d,d1);
          for (d2 = 0; d2 <= (d - d1); d2++)
                 if (mtx[j][d][d1][d2] == -123456789) /* pattern */
	           fprintf(fp, "%9s ", "..");
	    else if (d1 + d2 == d && 
		     mtx[j][d][d1][d2] <= (d + 1) * rnapar->P6) /* prohibited */
	           fprintf(fp, "%9s ", "**");
	    else if (d1 + d2 == d - 1 && 
		     mtx[j][d][d1][d2] <= (d1 + d2 + 1) * rnapar->P6) /* prohibited */
	           fprintf(fp, "%9s ", "**");
	    else if (mtx[j][d][d1][d2] <= (d1+d2+2) * rnapar->P6P) /* prohibited */
    	           fprintf(fp, "%9s ", "**");
	    else
	           fprintf(fp, "%9d ", mtx[j][d][d1][d2]);
          fputs("\n", fp);
	}
}	
