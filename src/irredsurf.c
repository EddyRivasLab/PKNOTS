/* irredsurf.c (mdified)
 * 
 * calculates the IS1 (hairpin loops) and IS2(bulges, internal loops, stacks(
 */
                           
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>

#include "cfg.h"
#include "proto.h"
#include "squid.h"
 
/* Function: F1()
 * 
 * Purpose:   calculates the IRREDUCIBLE SURFACES of  O(1). (hairpin loops) 
 *
 * Arguments: s    - sequence (0..len-1) converted to integers
 *                    representing array indices of symbols
 *            len  - length of iseq
 *            icfg - context-free grammar state transitions, integer log form
 *            j,d
 *         
 * Return:   F1
 */       
int
F1(int *s, int len, struct rnapar_2 *rnapar, int **icfg, int j, int d)
{
  int F1;
  int  i;     /* sequence start position */

  i = j - d;
  
  if (icfg[idxS][idxP(s[i],s[j])] == 1*INTSCALE)
    { 
      /*no stacking energies */
      if (d == 4)     
	F1 = rnapar->P4 
	  + wsf*icfg[idxS][idxX(d-1)];
      
      /* extrapolate for large loops */
      else if (d > 31)  
	F1 = rnapar->P4 
	  + wsf*icfg[idxPL(s[i],s[j])][idxPL(s[i+1],s[j-1])] 
	  + wsf*icfg[idxS][idxX(30)] 
	  + wsf*(int)(PRELOG*log((d-1)/30.));
      
      /* check for ultrastable tetraloops */
      else if (d == 5)  
	F1 = rnapar->P4 
	  + wsf*icfg[idxPL(s[i],s[j])][idxPL(s[i+1],s[j-1])] 
	  + wsf*icfg[idxS][idxX(d-1)]
	  + wsf*icfg[idxS][idxTL(s[i]+s[i+1]*4+s[i+2]*4*4+s[i+3]*4*4*4+s[i+4]*4*4*4*4+s[i+5]*4*4*4*4*4)];
      else
	F1 = rnapar->P4 
	  + wsf*icfg[idxPL(s[i],s[j])][idxPL(s[i+1],s[j-1])] 
	  + wsf*icfg[idxS][idxX(d-1)];
    }

  else F1 = -BIGINT;

  return F1;
  
}

/* Function: F2()
 * 
 * Purpose:  calculates the IRREDUCIBLE SURFACES of  O(2)
 *                          (internal loops, bulges, and stems)
 *           
 * Arguments: s    - sequence (0..len-1) converted to integers
 *                    representing array indices of symbols
 *            len  - length of iseq
 *            icfg - context-free grammar state transitions, integer log form
 *            
 * Return:   F2
 */       
int
F2(int *s, int len, struct rnapar_2 *rnapar, int **icfg, int j, int d, int d1, int d2)
{
  int F2;
  int i;			/* sequence start position                         */
  int k, l;                     /* sequence intermediate positions                 */
  int dmin, arg;                /* minimum of d1, d2                               */   
  int dpen;                     /* penalty for asymmetric internal loops           */
  int penalty;                  /* final penalty for asymmetric internal loops     */

  i = j - d;
  k = i + d1;
  l = j - d2;

  dmin = (d1 < d2) ? d1 : d2;
  arg  = (dmin < 4) ? dmin-1 : 3;
  
  dpen    = abs(d1-d2)*wsf*icfg[idxS][idxWA(arg)];
  penalty = (dpen > MAXPEN) ? dpen : MAXPEN;
  
  if (icfg[idxS][idxP(s[i],s[j])] == 1*INTSCALE
      && icfg[idxS][idxP(s[k],s[l])] == 1*INTSCALE)
    { 
      
      /* IRREDUCIBLE SURFACES of  O(2). F2[j][d][d1][d2].
       * (stems, bulges, and internal loops) 
       */
      
	                        /* STEMS */
      if (d1 == 1 && d2 == 1)
	F2 = rnapar->P1 
	  + wsf*icfg[idxPS(s[i],s[j])][idxPS(s[k],s[l])];
      
	                        /* BULGES */
      else if (d1 == 1 && d2 == 2)  /* add stacking term */
	F2 = rnapar->P2  
	  + wsf*icfg[idxPS(s[i],s[j])][idxPS(s[k],s[l])] 
	  + wsf*icfg[idxS][idxV(1)];
      
      else if (d1 == 2 && d2 == 1)  /* add stacking term */
	F2 = rnapar->P2  
	  + wsf*icfg[idxPS(s[i],s[j])][idxPS(s[k],s[l])] 
	  + wsf*icfg[idxS][idxV(1)];
      
      else if (d1 > 31 && d2 == 1)  /* extrapolate for large loops */
	F2 = rnapar->P2  
	  + wsf*icfg[idxS][idxV(30)] 
	  + wsf*(int)(PRELOG*log((d1-1)/30.));
      
      else if (d1 == 1 && d2 > 31)  /* extrapolate for large loops */
	F2 = rnapar->P2 
	  + wsf*icfg[idxS][idxV(30)] 
	  + wsf*(int)(PRELOG*log((d2-1)/30.));
      
      else if (d1 == 1 && d2 > 2 && d2 < 32)
	F2 = rnapar->P2  
	  + wsf*icfg[idxS][idxV(d2-1)];
      
      else if (d2 == 1 && d1 > 2 && d1 < 32)
	F2 = rnapar->P2  
	  + wsf*icfg[idxS][idxV(d1-1)];
      
	                        /* INTERNAL LOOPS */
      else if (d1 == 2 && d2 == 2 && /* Jame's rule (special case) */
	       icfg[idxS][idxP(s[i+1],s[j-1])] != INTSCALE)
	F2 = rnapar->P3 		
	  + wsf*icfg[idxPS(s[i],s[j])][idxPS(s[i+1],s[j-1])] 
	  + wsf*icfg[idxPS(s[l],s[k])][idxPS(s[l+1],s[k-1])];
      
      else if (d1 == 2 && d2 == 2 && 
	       icfg[idxS][idxP(s[i+1],s[j-1])] == INTSCALE)
	{
	  if (icfg[idxPS(s[k-1],s[l+1])][idxPS(s[k],s[l])] ==  /* this is the case so far */
	      icfg[idxPS(s[l],s[k])][idxPS(s[l+1],s[k-1])])    /* with current parameters */
	    
	    F2 = -BIGINT;  /* better read it as two stems */
	  
	  else
	    F2 = rnapar->P3 		
	      + wsf*icfg[idxPI(s[i],s[j])][idxPI(s[i+1],s[j-1])] 
	      + wsf*icfg[idxPI(s[l],s[k])][idxPI(s[l+1],s[k-1])];
	}
      
       else if (d1 + d2 > 32)
	F2 = rnapar->P3 
	  + wsf*icfg[idxPI(s[i],s[j])][idxPI(s[i+1],s[j-1])] 
	  + wsf*icfg[idxPI(s[l],s[k])][idxPI(s[l+1],s[k-1])] 
	  + wsf*icfg[idxS][idxW(30)] 
	  /*+ penalty*/
	  + wsf*(int)(PRELOG*log((d1+d2-2)/30.));
      
      else 
	F2 = rnapar->P3
	  + wsf*icfg[idxPI(s[i],s[j])][idxPI(s[i+1],s[j-1])] 
	  + wsf*icfg[idxPI(s[l],s[k])][idxPI(s[l+1],s[k-1])] 
	  + wsf*icfg[idxS][idxW(d1+d2-2)]
	  /*+ penalty*/;
    }
  else F2 = -BIGINT;

  return F2;

}



