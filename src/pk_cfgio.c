/* cfgio.c
 * I/O of models to/from disk files
 * SRE, Sat May 28 14:42:02 1994
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "pknots.h"
#include "pk_cfgio.h"
#include "pk_model.h"

#ifdef MEMDEBUG
#include "dbmalloc.h"
#endif

/* The magic number is "yrn1" + 0x80808080 */
static int v10magic     = 0xf9f2eeb1;


/* Function: SaveSCFG()
 * 
 * Purpose:  Write an SCFG to disk in binary format.
 *           Return 1 on success, 0 on failure
 */
int
SaveSCFG(FILE *ofp, float **cfg)
{
  int fs;

  if (fwrite(&v10magic, sizeof(int), 1, ofp) < 1)
    return 0;
  for (fs = 0; fs < NSTATES; fs++)
    if (fwrite(cfg[fs], sizeof(float), NSTATES, ofp) < NSTATES)
      return 0;
  return 1;
}


/* Function: ReadSCFG()
 * 
 * Purpose:  Read an SCFG from disk in binary format.
 *           
 *           
 * Return:   Return 1 on success, 0 on failure.
 *           ret_cfg is freed by called with FreeSCFG().
 */
int 
ReadSCFG(FILE *ofp, float ***ret_cfg)
{
  int fs;
  float **cfg;
  int magic;

  cfg = AllocSCFG();
  if (fread(&magic, sizeof(int), 1, ofp) < 1)
    { printf("Failed to read magic number from SCFG save file"); return 0; }
  if (magic == v10magic)
    {
      for (fs = 0; fs < NSTATES; fs++)
	if (fread(cfg[fs], sizeof(float), NSTATES, ofp) < NSTATES)
	  { printf("fread failed on SCFG save file"); return 0; }
    }
  else
    { printf("bad magic on that SCFG save file"); return 0; }
  *ret_cfg = cfg;
  return 1;
}


/* Function: WriteRdbSCFG()
 * 
 * Purpose:  Write a SCFG to disk in RDB database format.
 *           Though the file has such long lines that it is
 *           essentially unreadable, the strategy is deliberate:
 *           as SCFGs change, this code remains completely general,
 *           and necessary tables can be generated by RDB commands.
 */
void
WriteRdbSCFG(FILE *ofp, float **cfg)
{
  int i, j;

  /* Comment section of RDB database file
   */
  fprintf(ofp, "# SCFG output\n");
  fprintf(ofp, "#\n");

  /* First line contains column names.
   */
  fprintf(ofp, "Node\tState\t");
  for (j = 0; j < NSTATES; j++)
    fprintf(ofp, "%s\t", stNAME[j]);
  fprintf(ofp, "\n");

  /* Second line contains formatting
   */
  fprintf(ofp, "9\t9\t");
  for (j = 0; j < NSTATES; j++)
    fprintf(ofp, "9n\t");
  fprintf(ofp, "\n");

  /* And the rest is data.
   */
  for (i = 0; i < NSTATES; i++)
    {
      fprintf(ofp, "%s\t%s\t", dpcNAME[Ntype[i]], stNAME[i]);
      for (j = 0; j < NSTATES; j++)
	fprintf(ofp, "%9.3f\t", cfg[i][j]);
      fprintf(ofp, "\n");
    }
}

/* Function: WriteRdbISCFG()
 * 
 * Purpose:  Write a integer-form SCFG to disk in RDB database format.
 */
void
WriteRdbISCFG(FILE *ofp, int **icfg)
{
  int i, j;

  /* Comment section of RDB database file
   */
  fprintf(ofp, "# SCFG output\n");
  fprintf(ofp, "#\n");

  /* First line contains column names.
   */
  fprintf(ofp, "Node\tState\t");
  for (j = 0; j < NSTATES; j++)
    fprintf(ofp, "%s\t", stNAME[j]);
  fprintf(ofp, "\n");

  /* Second line contains formatting
   */
  fprintf(ofp, "9\t9\t");
  for (j = 0; j < NSTATES; j++)
    fprintf(ofp, "9n\t");
  fprintf(ofp, "\n");

  /* And the rest is data.
   */
  for (i = 0; i < NSTATES; i++)
    {
      fprintf(ofp, "%s\t%s\t", dpcNAME[Ntype[i]], stNAME[i]);
      for (j = 0; j < NSTATES; j++)
	fprintf(ofp, "%9.3f\t", (float) icfg[i][j] / (float) INTSCALE);
      fprintf(ofp, "\n");
    }
}

/* Function: WriteRdbSummary()
 * 
 * Purpose:  Write a more compact view of a counts-based SCFG in
 *           RDB database format. The summary is generated
 *           by summing over node types. 
 */
void
WriteRdbSummary(FILE *ofp, float **cfg)
{
  int i, j;
  float sum[NNODES][NNODES];

  /* Calculate the by-node summary
   */
  for (i = 0; i < NNODES; i++)
    for (j = 0; j < NNODES; j++)
      sum[i][j] = 0.0;
  for (i = 0; i < NSTATES; i++)
    for (j = 0; j < NSTATES; j++)
      sum[Ntype[i]][Ntype[j]] += cfg[i][j];

  /* Comment section of RDB database file
   */
  fprintf(ofp, "# SCFG counts data, summarized by node type\n");
  fprintf(ofp, "#\n");

  /* First line contains column names.
   */
  fprintf(ofp, "Node\t");
  for (j = 0; j < NNODES; j++)
    fprintf(ofp, "%s\t", dpcNAME[j]);
  fprintf(ofp, "\n");

  /* Second line contains formatting
   */
  fprintf(ofp, "9\t");
  for (j = 0; j < NNODES; j++)
    fprintf(ofp, "6n\t");
  fprintf(ofp, "\n");

  /* And the rest is data.
   */
  for (i = 0; i < NNODES; i++)
    {
      fprintf(ofp, "%s\t", dpcNAME[i]);
      for (j = 0; j < NNODES; j++)
	fprintf(ofp, "%6.0f\t", sum[i][j]);
      fputs("\n", ofp);
    }
}
