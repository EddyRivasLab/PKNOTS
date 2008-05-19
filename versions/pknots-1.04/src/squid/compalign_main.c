/* SQUID - A C function library for biological sequence analysis
 * Copyright (C) 1992-1996 Sean R. Eddy	
 *
 *    This source code is distributed under terms of the
 *    GNU General Public License. See the files COPYING 
 *    and GNULICENSE for further details.
 *
 */

/* main for compalign
 * 
 * Compalign -- a program to compare two sequence alignments
 * SRE, Tue Nov  3 07:38:03 1992
 * incorporated into SQUID, Thu Jan 26 16:52:41 1995
 * 
 * Usage: compalign <trusted-alignment> <test-alignment>
 * 
 * Calculate the fractional "identity" between the trusted alignment
 * and the test alignment. The two files must contain exactly the same
 * sequences, in exactly the same order.
 * 
 * The identity of the multiple sequence alignments is defined as
 * the averaged identity over all N(N-1)/2 pairwise alignments. 
 * 
 * The fractional identity of two sets of pairwise alignments
 * is in turn defined as follows (for aligned known sequences k1 and k2,
 * and aligned test sequences t1 and t2):
 * 
 *           matched columns / total columns, 
 *       
 *       where total columns = the total number of columns in
 *        which there is a valid (nongap) symbol in k1 or k2;
 *        
 *       matched columns = the number of columns in which one of the
 *         following is true:
 *         
 *          k1 and k2 both have valid symbols at a given column; t1 and t2
 *             have the same symbols aligned in a column of the t1/t2
 *             alignment;
 *             
 *          k1 has a symbol aligned to a gap in k2; that symbol in t1
 *             is also aligned to a gap;
 *             
 *          k2 has a symbol aligned to a gap in k1; that symbol in t2
 *             is also aligned to a gap.
 * 
 * Because scores for all possible pairs are calculated, the
 * algorithm is of order (N^2)L for N sequences of length L;
 * large sequence sets will take a while.
 * 
 * Sean Eddy, Tue Nov  3 07:46:59 1992
 * 
 */

#include <stdio.h>
#include <string.h>
#include "squid.h"

struct opt_s OPTIONS[] = {
  { "-c", TRUE, sqdARG_NONE },     
  { "-h", TRUE, sqdARG_NONE },     
};
#define NOPTIONS (sizeof(OPTIONS) / sizeof(struct opt_s))

char usage[]  = "Usage: compalign [-options] <trusted.ali> <test.ali>\n\
  Compare two multiple alignments of the same sequences.\n\
  Available options are:\n\
   -c       : only compare under marked #=CS consensus structure\n\
   -h       : print short help and usage info\n";


int
main(int argc, char **argv)
{
  char  *kfile;                 /* name of file of trusted (known) alignment */
  char  *tfile;                 /* name of file of test alignment            */
  int    format;		/* alignment file format (kSelex or kMSF)    */
  char **kseqs;                 /* trusted alignment                         */
  char **kraw;			/* dealigned trusted seqs                    */
  AINFO  kinfo;			/* alignment info for trusted alignment      */
  char **tseqs;                 /* test alignment                            */
  char **traw;			/* dealigned test sequences                  */
  AINFO  tinfo;			/* alignment info for test alignment         */
  int    idx;			/* counter for sequences                     */
  int    apos;			/* position in alignment                     */
  float  score;			/* RESULT: score for the comparison          */

  int    cs_only;		/* TRUE to compare under #=CS annotation only */
  int   *ref;			

  char *optname;
  char *optarg;
  int   optind;

  /***********************************************
   * Parse command line
   ***********************************************/

  cs_only = FALSE;

  while (Getopt(argc, argv, OPTIONS, NOPTIONS, usage,
                &optind, &optname, &optarg))
    {
      if      (strcmp(optname, "-c") == 0) { cs_only = TRUE; }
      else if (strcmp(optname, "-h") == 0) {
        printf("compalign %s, %s\n%s\n", squid_version, squid_date, usage);
        exit(EXIT_SUCCESS);
      }
    }

  if (argc - optind != 2)
    Die("Incorrect number of command line arguments.\n%s\n", usage); 

  kfile = argv[optind++];
  tfile = argv[optind];
		 
  /***********************************************
   * Read in the alignments
   ***********************************************/

  if (! SeqfileFormat(kfile, &format, NULL))
    switch (squid_errno) {
    case SQERR_NOFILE: 
      Die("Trusted alignment %s could not be opened for reading", kfile);
      break;
    /*FALLTHRU*/
    case SQERR_FORMAT: 
    default:           Die("Failed to determine format of trusted alignment %s", kfile);
    }
  if (! ReadAlignment(kfile, format, &kseqs, &kinfo))
    Die("Failed to parse trusted alignment file %s", kfile);

  if (! SeqfileFormat(tfile, &format, NULL))
    switch (squid_errno) {
    case SQERR_NOFILE: Die("Test alignment %s could not be opened for reading", tfile);
    case SQERR_FORMAT: 
    default:           Die("Failed to determine format of test alignment %s", tfile);
    }
  if (! ReadAlignment(tfile, format, &tseqs, &tinfo))
    Die("Failed to parse test alignment file %s", tfile);

				/* test that they're the same! */
  if (kinfo.nseq != tinfo.nseq)
    Die("files %s and %s do not contain same number of seqs!\n", kfile, tfile);

  for (idx = 0; idx < kinfo.nseq; idx++)
    {
      s2upper(kseqs[idx]);
      s2upper(tseqs[idx]);
    }

  for (idx = 0; idx < kinfo.nseq; idx++)
    if (strcmp(kinfo.sqinfo[idx].name, tinfo.sqinfo[idx].name) != 0)
      Die("seqs in %s and %s don't seem to be in the same order\n  (%s != %s)",
	  kfile, tfile, kinfo.sqinfo[idx].name, tinfo.sqinfo[idx].name);

  DealignAseqs(kseqs, kinfo.nseq, &kraw);
  DealignAseqs(tseqs, tinfo.nseq, &traw);
  for (idx = 0; idx < kinfo.nseq; idx++)
    if (strcmp(kraw[idx], traw[idx]) != 0)
      Die("raw seqs in %s and %s are not the same (died at %s, number %d)\n",
	  kfile, tfile, kinfo.sqinfo[idx].name, idx);
  Free2DArray(kraw, kinfo.nseq);
  Free2DArray(traw, tinfo.nseq);

  if (cs_only)
    {
      if (! (kinfo.flags & AINFO_CS))
	Die("Trusted alignment %s has no consensus structure annotation\n  -- can't use -c!\n",
	    kfile);
      ref = (int *) MallocOrDie (sizeof(int) * kinfo.alen);
      for (apos = 0; apos < kinfo.alen; apos++)
	ref[apos] = (isgap(kinfo.cs[apos])) ? FALSE : TRUE;
    }	

  /***********************************************
   * Compare the alignments, print results
   ***********************************************/

  printf("compalign %s, %s\n", squid_version, squid_date);
  printf("\n");
  printf("Trusted alignment: \t%s\n", kfile);
  printf("Test alignment:    \t%s\n", tfile);
  printf("Total sequences:   \t%d\n", kinfo.nseq);
  printf("\n");

  if (cs_only)
    score = CompareRefMultAlignments(ref, kseqs, tseqs, kinfo.nseq);
  else
    score = CompareMultAlignments(kseqs, tseqs, kinfo.nseq);

  printf("Alignment identity = %.4f\n", score);
  
  /***********************************************
   * Garbage collection and exit
   ***********************************************/

  if (cs_only) free(ref);
  FreeAlignment(kseqs, &kinfo);
  FreeAlignment(tseqs, &tinfo);
  return 0;
}


