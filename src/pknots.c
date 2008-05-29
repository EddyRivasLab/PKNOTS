/* knotf_main.c
 *
 * Fold RNAs including knots.
 */
                                                                       
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>

#include <easel.h>
#include <esl_alphabet.h>
#include <esl_sqio.h>
#include <esl_wuss.h>

#include "version.h"
#include "cfg.h"
#include "proto.h"
#include "protovx.h"
#include "protowbx.h"
#include "protowx.h"

#ifdef MEMDEBUG
#include "dbmalloc.h"
#endif
              
static void ct_output(FILE *ofp, char *seq, SQINFO sqinfo, int *ss, int j, int d);
static void param_output(FILE *outf, struct rnapar_2 zkn_param, 
			 int format, int shuffleseq, int allow_pseudoknots, 
			 int approx, int score, float pairs, float cykpairs);

static ESL_OPTIONS options[] = {
 /* name                type             default  env_var range   toggles req   incompat help                                      docgroup */
  { "-c",               eslARG_NONE,     FALSE,   NULL,   NULL,   NULL,   NULL, NULL,    "add L^5 coaxials (V6)",                  0 },
  { "-g",               eslARG_STRING,   NULL,    NULL,   NULL,   NULL,   NULL, NULL,    "save ct-format files",                   0 },
  { "-h",               eslARG_NONE,     FALSE,   NULL,   NULL,   NULL,   NULL, NULL,    "show help and usage",                    0 },
  { "-k",               eslARG_NONE,     FALSE,   NULL,   NULL,   NULL,   NULL, NULL,    "allow pseudoknots",                      0 },
  { "-s",               eslARG_NONE,     FALSE,   NULL,   NULL,   NULL,   NULL, NULL,    "shuffle sequences",                      0 },
  { "-t",               eslARG_NONE,     FALSE,   NULL,   NULL,   NULL,   NULL, NULL,    "print traceback",                        0 },
  { "-v",               eslARG_NONE,     FALSE,   NULL,   NULL,   NULL,   NULL, NULL,    "be verbose?",                            0 },
  {  0, 0, 0, 0, 0, 0, 0, 0, 0, 0 },
};

static char usage[] = "\
Usage: ./pknots [-options] <infile> <outfile>\n\
 where <infile> is a fasta or Stockholm file with sequences to fold\n\
 results in Stockholm format are saved to <outfile>.\n\
 \n\
 Available options are:\n\
   -g      : ct output\n\
   -h      : print short help and usage info\n\
   -c      : add V6 the time consuming (N^5) coaxials\n\
   -k      : allow pseudoknots\n\
   -s      : shuffle the sequences\n\
   -t      : print traceback\n\
   -v      : verbose debugging output\n\
"
;

int
main(int argc, char **argv)
{ 
  struct rnapar_2   zkn_param;          /* rna parameters secon order + knots              */
  struct tracekn_s *tr; 	        /* traceback of a predicted RNA structure          */
  char             *seqfile;            /* input sequence file                             */
  char             *outfile;            /* where to send the output                        */
  ESL_GETOPTS      *go;
  ESL_SQFILE       *sqfp;      	        /* open sequence file                              */
  ESL_SQ           *sq;
  ESL_ALPHABET     *abc;
  int             **icfg;               /* integer log form grammar for alignment          */
  int               L;
  CYKVAL            sc;                 /* score of predicted structure                    */
  float             pairs;              /* number of base paired as input                  */
  float             cykpairs;           /* number of base paired of predicted structure    */

  int               approx;             /* TRUE  == external pseudoknot approximation      *
				         *          exclude diagrams V7-V10 and WB9-WB10   *
				         * FALSE == full pseudoknot model                  *
				           include diagrams V7-V10 and WB9-WB10   */
  int               allow_coaxials;     /* TRUE add V6 coaxial diagrams                    */
  int               allow_pseudoknots;	/* TRUE  == include pseudoknots                    *
				         * FALSE == no pseudoknots                         */
  int               ctoutput;           /* TRUE print ctoutput                             */
  int               shuffleseq;         /* TRUE to shuffle the sequences                   */
  int               traceback;          /* TRUE print traceback                            */
  int               verbose;            /* TRUE to be extremely verbose to debug           */
  int               i;
  int               status;
  
  /* Process command line
   */
  if ((go = esl_getopts_Create(options)) == NULL)      esl_fatal("Bad option structure\n");
  if (esl_opt_ProcessCmdline(go, argc, argv) != eslOK) esl_fatal("Failed to parse command line: %s\n", go->errbuf);
  if (esl_opt_VerifyConfig(go)               != eslOK) esl_fatal("Failed to parse command line: %s\n", go->errbuf);

  if (esl_opt_GetBoolean(go, "-h") == TRUE) {
    puts(usage); 
    puts("\n  where options are:");
    esl_opt_DisplayHelp(stdout, go, 0, 2, 80); /* 0=all docgroups; 2=indentation; 80=width */
    return 0;
  }
  if (esl_opt_ArgNumber(go) != 3) esl_fatal("Incorrect number of command line arguments.\n%s\n", usage);
 
  allow_coaxials    = esl_opt_GetBoolean(go, "-c");
  ctoutput          = esl_opt_GetBoolean(go, "-g");
  allow_pseudoknots = esl_opt_GetBoolean(go, "-k");
  shuffleseq        = esl_opt_GetBoolean(go, "-s");
  traceback         = esl_opt_GetBoolean(go, "-t");
  verbose           = esl_opt_GetBoolean(go, "-v");

  seqfile = esl_opt_GetArg(go, 1);
  outfile = esl_opt_GetArg(go, 2);
  esl_getopts_Destroy(go);

  /*********************************************** 
   * Create the starting model 
   ***********************************************/
  Parameters2_Zkn(&zkn_param);
  icfg = ParamIntSCFG(&zkn_param); 

  /*********************************************** 
   * Open sequence file and output file
   ***********************************************/

  if (! SeqfileFormat(seqfile, &format, NULL))
    Die("Failed to determine format of sequence file %s\n", seqfile);
  
  if ((sqfp = SeqfileOpen(seqfile, format, NULL)) == NULL)
    Die("Failed to open sequence file %s", seqfile);

  ofp = stdout;
  if (outfile != NULL && (ofp = fopen(outfile, "w")) == NULL)
    Die("Failed to open output file %s", outfile);
 
  /* Open the output file and the input seqfile.
   */
  if ((fp = fopen(outfile, "w")) == NULL) 
    pk_fatal("failed to open %s for output", outfile);
  if (esl_sqfile_Open(seqfile, infmt, NULL, &sqfp) != eslOK)
    pk_fatal("failed to open %s", seqfile);
  sq = esl_sq_Create();


 
  /*********************************************** 
   * Print banner
   ***********************************************/

  puts(banner);  
  printf("        PKNOTS %s (%s)", RELEASE, RELEASEDATE);
  printf(" using squid %s (%s)\n", squid_version, squid_date);
  printf("---------------------------------------------------\n");
  printf("Folding sequences from:  %s \n", seqfile);
  printf("---------------------------------------------------\n");
  puts("");

  /* OK, start folding.
   */
  while ((status = esl_sqio_Read(sqfp, sq)) == eslOK) 
    {
      L   = sq->n;

      if (verbose) {
	printf("SEQ:\n");
	for (i = 0; i < L; i++) printf("%c ", seq[i]);
	printf("\n");
      }
     
      ESL_ALLOC(sq->dsq, sizeof(char) * (L+2));

      if (esl_abc_Digitize(abc, sq->seq, sq->dsq) != eslOK)
	pk_fatal("failed to digitize sequence");


    }

  while (ReadSeq(sqfp, format, &seq, &sqinfo))
    {
      s2upper(seq);
      StripDegeneracy(seq);

      if (verbose) {
	printf("SEQ:\n");
	for (i = 0; i < sqinfo.len; i++) printf("%c ", seq[i]);
	printf("\n");
      }

       if ((format == kSquid || format == kSelex) && sqinfo.flags & SQINFO_SS) 
	 PrintCtSeq(ofp, &sqinfo, seq, 0, sqinfo.len, sqinfo.ss);
      
     if (shuffleseq) ShuffleSequence(seq, sqinfo.len, sqinfo.len-1, verbose);
      
     StructurePredictkn_2IS(ofp, seq, sqinfo.len, &zkn_param, icfg, verbose, traceback,
			     &tr, &score, allow_coaxials, allow_pseudoknots, approx);
      Tracekn(tr, seq, sqinfo.len, FALSE, &cc);
      Traceintkn(tr, seq, sqinfo.len, FALSE, &ss);
      sqinfo.flags |= SQINFO_SS;
      
      /* 
       * print output
       */
      WriteSeqkn(ofp, seq, &sqinfo, ss, format, &pairs, &cykpairs);
      
      /* If this was a kSquid or Selex file with a "given" secondary structure,
       * compare the two structures
       */
      if ((format == kSquid || format == kSelex) && sqinfo.flags & SQINFO_SS) { 
	CompareRNAStructures(ofp, 0, sqinfo.len, sqinfo.ss, ss);
      }

      param_output(ofp, zkn_param, format, shuffleseq, allow_pseudoknots, approx, score, pairs, cykpairs);

      if (ctoutput) 	  
	ct_output(ofp, seq, sqinfo, ss, sqinfo.len-1, sqinfo.len-1);

      FreeTracekn(tr);
      FreeSequence(seq, &sqinfo);
      free(cc);
      free(ss);
    }
  
  /*********************************************** 
   * Cleanup and exit
   ***********************************************/

  FreeIntSCFG(icfg);
  esl_sq_Destroy(sq);
  esl_sqfile_Close(sqfp);
  esl_alphabet_Destroy(abc);
  fclose(fp);
  exit (0);
}




/* Function: param_output()
 * 
 * Purpose:  Prints the info about the parameters used in the search.
 *           
 * Args:     outf     - where to send the output
 *           seq      - sequence to fold. 
 *           sqinfo   - info structures for seq                 
 *           ss       - secondary structure in integer form.
 *                      ss[i] = j  if position i paired to position j
 *                      ss[i] = -1 if position i is unpaired
 *           j        - final position of fragment to report. 
 *                      j\in [0,seqlen-1], where seqlen is the total length of the sequence.
 *           d        - length of the fragment.
 *         
 * Return:   (void)
 */          
static void
param_output(FILE *ofp, struct rnapar_2 zkn_param, int format, int shuffleseq, int allow_pseudoknots, int approx, 
	     int score, float pairs, float cykpairs)
{  
      fprintf(ofp,"----------------------------------------------------------------------\n");
      
      fprintf(ofp,"   Log odds score:    %8d \n",  score);
      fprintf(ofp,"energy (kcal/mol):  %10.2f \n", -(float)(score)/(INTSCALE*wsf));
      fprintf(ofp,"number of pairs found:      %4.2f \n", cykpairs);
      if (format == kSquid)
	fprintf(ofp,"number of pairs annotated:  %4.2f \n", pairs);
      fprintf(ofp,"----------------------------------------------------------------------\n");
      
      if (shuffleseq) { 
	fprintf(ofp,"----------------------------------------------------------------------\n");
	fprintf(ofp,"Shuffled sequence\n");
	fprintf(ofp,"----------------------------------------------------------------------\n");
      }
      
      if (!allow_pseudoknots) {
	fprintf(ofp,"Allow pseudoknots? NO\n");
	fprintf(ofp,"----------------------------------------------------------------------\n");
	fprintf(ofp,"Parameters (energy  units, kcal/mol)\n");
	fprintf(ofp," P1   parameter:    %10.2f \n", -(float)zkn_param.P1/(INTSCALE*wsf));
	fprintf(ofp," P2   parameter:    %10.2f \n", -(float)zkn_param.P2/(INTSCALE*wsf));
	fprintf(ofp," P3   parameter:    %10.2f \n", -(float)zkn_param.P3/(INTSCALE*wsf));
	fprintf(ofp," P4   parameter:    %10.2f \n", -(float)zkn_param.P4/(INTSCALE*wsf));
	fprintf(ofp," P5   parameter:    %10.2f \n", -(float)zkn_param.P5/(INTSCALE*wsf));
	fprintf(ofp," P6   parameter:    %10.2f \n", -(float)zkn_param.P6/(INTSCALE*wsf));
	fprintf(ofp," P10  parameter:    %10.2f \n", -(float)zkn_param.P10/(INTSCALE*wsf));
	fprintf(ofp,"----------------------------------------------------------------------\n");
      }
      else {
	fprintf(ofp,"Allow pseudoknots? YES\n");
	fprintf(ofp,"----------------------------------------------------------------------\n");
	fprintf(ofp,"Parameters (energy units, kcal/mol)\n");
	fprintf(ofp," wkn  parameter:    %10.2f (pseudoknots weight)\n", (float)wkn/wsf);
	fprintf(ofp," P1   parameter:    %10.2f \n", -(float)zkn_param.P1/(INTSCALE*wsf));
	fprintf(ofp," P2   parameter:    %10.2f \n", -(float)zkn_param.P2/(INTSCALE*wsf));
	fprintf(ofp," P3   parameter:    %10.2f \n", -(float)zkn_param.P3/(INTSCALE*wsf));
	fprintf(ofp," P4   parameter:    %10.2f \n", -(float)zkn_param.P4/(INTSCALE*wsf));
	fprintf(ofp," P5   parameter:    %10.2f \n", -(float)zkn_param.P5/(INTSCALE*wsf));
	fprintf(ofp," P5P  parameter:    %10.2f \n", -(float)zkn_param.P5P/(INTSCALE*wkn));
	fprintf(ofp," P6   parameter:    %10.2f \n", -(float)zkn_param.P6/(INTSCALE*wsf));
	fprintf(ofp," P6P  parameter:    %10.2f \n", -(float)zkn_param.P6P/(INTSCALE*wkn));
	fprintf(ofp," P10  parameter:    %10.2f \n", -(float)zkn_param.P10/(INTSCALE*wsf));
	fprintf(ofp," P10P parameter:    %10.2f \n", -(float)zkn_param.P10P/(INTSCALE*wkn));
	fprintf(ofp," P11  parameter:    %10.2f \n", -(float)zkn_param.P11/(INTSCALE*wsf));
	fprintf(ofp," P12  parameter:    %10.2f \n", -(float)zkn_param.P12/(INTSCALE*wsf));
	
	if (approx){
	  fprintf(ofp," P13  parameter:    infinity\n");
	  fprintf(ofp," (external pseudoknots approximation)\n");
	}
	else {
	  fprintf(ofp," P13  parameter:    %10.2f \n", -(float)zkn_param.P13/(INTSCALE*wsf));
	  fprintf(ofp," full pseudoknot model.\n");
	fprintf(ofp,"----------------------------------------------------------------------\n");
	}
      }
}      

/* Function: ct_output()
 * 
 * Purpose:  Converts the output of pknots program to the "connect" format.
 *           It allows to convert any fragment of the sequence.
 *           
 * Args:     outf     - where to send the output
 *           seq      - sequence to fold. 
 *           sqinfo   - info structures for seq                 
 *           ss       - secondary structure in integer form.
 *                      ss[i] = j  if position i paired to position j
 *                      ss[i] = -1 if position i is unpaired
 *           j        - final position of fragment to report. 
 *                      j\in [0,seqlen-1], where seqlen is the total length of the sequence.
 *           d        - length of the fragment.
 *         
 * Return:   (void)
 */          
static void
ct_output(FILE *ofp, char *seq, SQINFO sqinfo, int *ss, int j, int d)
{  
  int i;          /* initial position of fragment. 0,..,d-1 */
  int iabs;       /* absolute value of initial position     */

  fprintf(ofp,"\n ct_output \n");
  fprintf(ofp,"----------------------------------------------------------------------\n");

  for (i = 0; i <= d; i++) {
    iabs = i + j - d;
    fprintf(ofp, "%5d %c   %5d %4d %4d %4d\n",
            i+1, seq[iabs], i, i+2, (ss[iabs] != -1)? ss[iabs]+1-j+d:0, iabs+1);
  }
  
  fprintf(ofp, "\n");
  
}      

