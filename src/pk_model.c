/* model.c
 * SRE, Fri May 27 11:25:06 1994
 * 
 * Allocation, free'ing, initialization of the SCFG
 */

#include <stdio.h>
#include <stdlib.h>
#include <limits.h>

#include "cfg.h"
#include "proto.h"
#include "squid.h"


void
AllocBaseFreq(float **ret_basefreq)
{
  float *basefreq;
  int          nt;

  basefreq = (float *) MallocOrDie(sizeof(float) * 4);


  /* Initialize all probabilities to zero
   */
  for (nt = 0; nt < 4; nt++)
    basefreq[nt] = 0.;

  *ret_basefreq = basefreq;
}

/* Function: BaseFreq() 
 * 
 * Purpose:  to calculate the base freq of a sequence in characters
 *           
 * Args:     s            - sequence composition (A, C, G, U)
 *           j            - end position                 
 *           d            - length of sequence             
 *           ret_basefreq - vector of frequencies              
 *                  
 * Return:  (void) basefreq are filled up.                
 */
void
BaseFreq(char *s, int j, int d, float **ret_basefreq)
{
  int    mid;
  int    nt;
  float *basefreq;
  float  sum = 0.;

  AllocBaseFreq(&basefreq);

  for (mid = 0; mid <= d; mid++)
    if (s[j-mid] == 'A')      basefreq[0] += 1.;
    else if (s[j-mid] == 'C') basefreq[1] += 1.;
    else if (s[j-mid] == 'G') basefreq[2] += 1.;
    else if (s[j-mid] == 'U') basefreq[3] += 1.;
    else Die("Unrecognized character (%c, pos = %d) in sequence\n", s[j-mid], j-mid);

  for (nt = 0; nt < 4; nt++) {
    basefreq[nt] /= (float)(d+1);    /* normalize            */
    sum += basefreq[nt];           /* add for sumcheck     */
  }

  if (sum > 1.001 || sum < 0.999) 
    Die ("freqs of nts do not add up to one (sum = %f)", sum);

  *ret_basefreq = basefreq;
}

float **
AllocSCFG(void)
{
  float **cfg;
  int     i, j;

  /* This way of alloc'ing a 2D array keeps the CFG all in one
   * contiguous chunk of RAM and might keep us in cache.
   */
  if ((cfg    = (float **) malloc (sizeof(float *) * NSTATES)) == NULL ||
      (cfg[0] = (float *)  malloc (sizeof(float)   * NSTATES * NSTATES)) == NULL)
    Die("malloc failed");
  for (i = 1; i < NSTATES; i++)
    cfg[i] = cfg[0] + i * NSTATES;

  /* Initialize all transitions to zero
   */
  for (i = 0; i < NSTATES; i++)
    for (j = 0; j < NSTATES; j++)
      cfg[i][j] = 0.0;

  return cfg;
}

int **
AllocIntSCFG(void)
{
  int  **icfg;
  int    i;

  /* As above, we take care to keep the iSCFG in a contiguous
   * chunk of RAM, hoping to stay in cache.
   */
  if ((icfg    = (int **) malloc (sizeof(int *) * NSTATES)) == NULL ||
      (icfg[0] = (int *)  malloc (sizeof(int)   * NSTATES * NSTATES)) == NULL)
    Die("malloc failed");
  for (i = 1; i < NSTATES; i++)
    icfg[i] = icfg[0] + i * NSTATES;
  return icfg;
}

float **
DupSCFG(float **cfg)
{
  float **new;
  int i,j;

  new = AllocSCFG();
  for (i = 0; i < NSTATES; i++)
    for (j = 0; j < NSTATES; j++)
      new[i][j] = cfg[i][j];
  return new;
}
    
void
FreeSCFG(float **cfg)
{
  free(cfg[0]);
  free(cfg);
}


void 
FreeIntSCFG(int **icfg)
{
  free(icfg[0]);
  free(icfg);
}


/* Function: LogifySCFG()
 * 
 * Purpose:  Take an SCFG and convert it to integer log form
 *           for an alignment.
 *           
 * Args:     cfg  - the grammar, floating point form
 *                  
 * Return:   integer log form of the grammar
 *           Alloc'ed here; caller must free.                 
 */
int **
LogifySCFG(float **cfg)
{
  int **icfg;
  int   i,j;

  icfg = AllocIntSCFG();
  for (i = 0; i < NSTATES; i++)
    for (j = 0; j < NSTATES; j++)
      if (Connects[Ntype[i]][Ntype[j]])
	icfg[i][j] = IntizeScale(LOG2(cfg[i][j]));
      else
	icfg[i][j] = INT_MIN;

  return icfg;
}



/* Function: LogoddsifySCFG()
 * 
 * Purpose:  Take an SCFG and convert it to integer log odds form
 *           for an alignment.
 *           
 * Args:     cfg    - the grammar, floating point form
 *                  
 * Return:   integer log form of the grammar
 *           Alloc'ed here; caller must free.                 
 */
int **
LogoddsifySCFG(float **cfg)
{
  int **icfg;
  int   i,j;

  icfg = AllocIntSCFG();
  for (i = 0; i < NSTATES; i++)
    for (j = 0; j < NSTATES; j++)
      {
	if (Connects[Ntype[i]][Ntype[j]])
	  {
	    icfg[i][j] = IntizeScale(LOG2(cfg[i][j]));

	    /* Transitions into pairs; add ~LOG2(16) (4 bits)
	     * Transitions into singlets; add ~LOG2(4) (2 bits)
	     */
	    if (Ntype[j] == dpcP) 
	      icfg[i][j] += IntizeScale(4.0);
	    else if (Ntype[j] != dpcE && Ntype[j] != dpcS && Ntype[j] != dpcB)
	      icfg[i][j] += IntizeScale(2.0);
	  }
	else
	  icfg[i][j] = INT_MIN;
      }
  return icfg;
}

void
RandomSCFG(float **cfg)
{
  int i, j;

  for (i = 0; i < NSTATES; i++)
    for (j = 0; j < NSTATES; j++)
      cfg[i][j] = (Connects[Ntype[i]][Ntype[j]]) ? sre_random() : 0.0;
  NormalizeSCFG(cfg);
}

void
NormalizeSCFG(float **cfg)
{
  float sum;
  int   i, j;

  for (i = 0; i < NSTATES; i++)
    {
      for (sum = 0.0, j = 0; j < NSTATES; j++)
	sum += cfg[i][j];
      for (j = 0; j < NSTATES; j++)
	cfg[i][j] /= sum;
    }
}


/* Function: NussinovIntSCFG()
 * 
 * Purpose:  Construct an SCFG for simulating Nussinov-style
 *           RNA folding.
 *           
 *           The returned model is in integer form, but it
 *           is not based on probabilities. A score of
 *           1 is assigned to base pairs. Everything else
 *           scores 0. An alignment to this model is
 *           a maximization of base pairing.
 *           
 * Return:   an integer SCFG, which must be free'd by caller.
 */         
int **
NussinovIntSCFG()
{
  int   **icfg;
  int     i,j;
  int     symi, symj;

  icfg = AllocIntSCFG();

  /* Set everything to a score of zero (or prohibit)
   */
  for (i = 0; i < NSTATES; i++)
    for (j = 0; j < NSTATES; j++)
      icfg[i][j] = (Connects[Ntype[i]][Ntype[j]]) ? 0 : INT_MIN;
  
  /* Set transitions into Watson-Crick pairs to 1,
   * all transitions to other pairs to -50.
   */
  for (i = 0; i < NSTATES; i++)
    if (Connects[Ntype[i]][dpcP])
      {
	for (symi = 0; symi < 4; symi ++)
	  for (symj = 0; symj < 4; symj ++)
	    if (symi+symj == 3)	/* trick to determine base pairs (0,3), (3,0), (1,2), (2,1) */
	      icfg[i][idxP(symi,symj)] = 1 * INTSCALE;   /* Watson-Crick */
	    else
	      icfg[i][idxP(symi,symj)] = -1* BIGINT; /* PROHIBIT */
      }
  return icfg;
}


/* Function: ProbifySCFG()
 * 
 * Purpose:  Convert a model from counts to probabilities
 * 
 * Args:     cfg - model, in counts form. Returns in probability form.
 *                 
 * Return:   (void)
 */
void
ProbifySCFG(float **cfg)
{
  int   fs, ts;			/* from state, to state                      */
  float norm;			/* sum used for normalization in denominator */
  
  for (fs = 0; fs < NSTATES; fs++) 
    {
      norm = 0.0;
      for (ts = 0; ts < NSTATES; ts++)
	if (Connects[Ntype[fs]][Ntype[ts]])
	  norm += cfg[fs][ts] + 1.0; /* plus-one Laplace prior */
	else			/* paranoia never hurts. */
	  if (cfg[fs][ts] != 0.0)
	    Die("Somebody screwed up. %s -> %s transition is nonzero.\n",
		stNAME[fs], stNAME[ts]);
      
      for (ts = 0; ts < NSTATES; ts++)
	if (Connects[Ntype[fs]][Ntype[ts]])
	  cfg[fs][ts] = (cfg[fs][ts] + 1.0) / norm;
    }
}


/* Function: Stype()
 * 
 * Purpose:  Maps node/symi/symj to state type. Like the Ntype macro,
 *           but too complicated to be a macro.
 *           
 * Args:     node -- e.g. dpcP, etc. Index of a node type.
 *           symi -- 0..3 integer for ACGU. 
 *           symj -- 0..3 integer for ACGU.
 *           size -- if !=0 it is a bulge, hairpin or intloop
 *           tlp  -- number of ultrastable tetraloops
 *           
 * Return:   0..NSTATES-1 state index.
 *           If something goes wrong, returns -1.
 */
int
Stype(int node, int symi, int symj, int size, int asy, int tlp)
{          
  int st;

  switch (node) {
  case dpcP:  st = (symi<4 && symj<4) ? idxP(symi,symj) : -1; break;
  case dpcL:  st = (symi<4) ? idxL(symi)  : -1;               break;
  case dpcR:  st = (symj<4) ? idxR(symj)  : -1;               break;
  case dpcB:  st = idxB;                                      break;
  case dpcS:  st = idxS;                                      break;
  case dpcPL: st = (symi<4 && symj<4) ? idxPL(symi,symj) : -1; break;
  case dpcPI: st = (symi<4 && symj<4) ? idxPI(symi,symj) : -1; break;
  case dpcX:  st = (size<31) ? idxX(size) : -1;               break;
  case dpcV:  st = (size<31) ? idxV(size) : -1;               break;
  case dpcW:  st = (size<31) ? idxW(size) : -1;               break;
  case dpcWA: st = (asy<31)  ? idxWA(asy) : -1;               break;
  case dpcTL: st = (tlp<256)  ? idxTL(tlp) : -1;               break;
  case dpcE:  st = idxE;                                      break;
  default: Die("can't find state type for node, symi, symj, size, asy, tlp = %d,%d,%d,%d,%d,%d\n",
	       node, symi, symj, size, asy, tlp);
  }
  return st;
}
 

























