/* fillwindow.c 
 *
 * fill the no-hole matrices: VX, WX, WBX.         Dimension  (len x len)
 * and the hole matrices:     VHX, ZHX, YHX, WHX.  Dimension  (len x len x len x len)
 *
 */
                        
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>

#include "version.h"
#include "cfg.h"
#include "proto.h"
#include "squid.h" 

#ifdef MEMDEBUG
#include "dbmalloc.h"
#endif


/* Function: FillMtx_nestedScan()
 * 
 * Purpose:  The dynamic programming algorithm for folding
 *           an RNA using a stochastic context-free grammar.
 *           No pseudoknots
 *           
 * Arguments: s     - sequence (0..len-1) converted to integers
 *                    representing array indices of symbols
 *            len   - length of iseq
 *            cfg   - context-free grammar state transitions, integer log form
 *             vp   - DP matrix, already alloc'ed 
 *             vx   - DP matrix, already alloc'ed 
 *            wbx   - DP matrix, already alloc'ed 
 *             wx   - DP matrix, already alloc'ed 
 *            
 * Return:    (void)           
 *            vx, wx, wbx are filled (in that order).
 */       
void
FillMtx_nestedScan(FILE *outf, int *s, int seqlen, int len, int pos_end, int win, int off, 
		   int **icfg, int **wbx, int **wx, int **vx, int *vp, float **ret_score)
{    
  float   *score;     /* best scores backward for a given position       */
  int    jabs, j;     /* row index (sequence end position)               */
  int          d;     /* column indices ("distance from diagonal")       */
  int       jmod;
  int       dmax;  

  Alloc_Scores(len, &score);
  Pattern_Scores(len, score);
  
  for (j = 0; j < len; j++)
    {
      jabs = j + (pos_end - len + 1);

      IdxWindow(win, off, j, &jmod, &dmax);

      for (d = 0; d <= dmax; d++)        
	{
	  /*  filling of VX(j)(d) without pseudoknots
	   *                      
	   */
	  FillVPScan(s, len, win, icfg, vx, vp, jabs, jmod, d);
	  FillVX_nestedScan(s, len, win, icfg, wbx, vx, vp, jabs, jmod, d);
	  
	  /* filling  of  WX(j)(d) == No hole between (i,j). 
	   * filling  of WBX(j)(d) == No hole between (i,j) for multiloop loops. 
	   */
	  FillWBX_nestedScan(s, len, win, icfg, wbx, vx, jabs, jmod, d);

	  FillWX_nestedScan(s, len, win, icfg, wx, vx, jabs, jmod, d);
  
	}   /* while filling d */
      
      score[j] = (float)wx[jmod][dmax];
      
    }   /* while filling j */

  *ret_score = score;
}


 

/* Function: IdxWindow()
 * 
 * Purpose:  For a given position j, calculate the index j mod (win, off),
 *           and the maximum length of the interval dmax.
 *
 * Args:     j     - 
 *           jmod  - 
 *           dmax  - 
 *
 * Return:   (void) returns jmod and dmax               
 */
void
IdxWindow(int win, int off, int j, int *ret_jmod, int *ret_dmax)
{

  int jw;
  int joff;
  int jmod, dmax;

  jw = j % win;
  
  if (win >= off) 
    {
      jmod = jw;
      dmax = (j < win)? j : (win-1) - (off-1) + jw % off;
    }
  else if (win < off && off > 0)
    {
      joff = j % off;

      if (joff >= win) 
	{
	  jmod = 0; 
	  dmax = -1;
	}
      else 
	{
	  jmod = joff;
	  dmax = (j < win)? j : (joff < win)? joff : -1;
	}
    }
  else Die("you have to slide at least one nt at the time!");

  if (jmod >= win) Die("Bad assignment of jmod. IdxWindow()");
  if (dmax >= win) Die("Bad assignment of dmax. IdxWindow()");
  
  *ret_jmod = jmod;
  *ret_dmax = dmax;
}



/* Function: Alloc_Scores()
 * 
 * Purpose:  Malloc space for the score[j] matrix.
 *           j  = 0,...,len-1
 *           
 * Args:     len        - length of sequence
 *           ret_score  - RETURN: score matrix 
 * 
 * Return:   Ptr to allocated scoring matrix, or
 *           dies and exits.
 */
void
Alloc_Scores(int win, float **ret_score)
{
  float *score;
                                                     
  if ((score = (float *) malloc (sizeof(float) * win)) == NULL)
    Die("malloc failed in score");

  *ret_score  = score;
}

/* Function: pattern_scores()
 * 
 * Purpose:  For debugging: write a pattern into the matrix
 *           so I can tell what's been written to.
 */
void
Pattern_Scores(int win, float *score)
{
  int j;
  
  for (j = 0; j < win; j++)
    score[j] = -123456789.;
}   

