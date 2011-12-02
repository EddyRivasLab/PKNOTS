/* pknots.c
 *
 * Fold RNAs including knots.
 */
                                                                       
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>

#include <easel.h>
#include <esl_alphabet.h>
#include <esl_getopts.h>
#include <esl_random.h>
#include <esl_randomseq.h>
#include <esl_sq.h>
#include <esl_sqio.h>
#include <esl_wuss.h>

#include "pknots.h"
#include "pk_cyk.h"
#include "pk_model.h"
#include "pk_rnaparam.h"
#include "pk_rnaoutput.h"
#include "pk_trace.h"
#include "pk_util.h"
#include "pk_version.h"


#ifdef MEMDEBUG
#include "dbmalloc.h"
#endif
              
static char banner[] = 
"PKNOTS: optimal minimum-energy RNA folding with pseudoknots and coaxial energies";

static ESL_OPTIONS options[] = {
 /* name                type             default  env_var range   toggles req   incompat help                                      docgroup */
  { "-a",               eslARG_NONE,     FALSE,   NULL,   NULL,   NULL,   NULL, NULL,    "pseudoknot approx, exclude V7-V10 and WB9-WB1", 0 },
  { "-c",               eslARG_NONE,     FALSE,   NULL,   NULL,   NULL,   NULL, NULL,    "add L^5 coaxials (V6)",                         0 },
  { "-g",               eslARG_NONE,     NULL,    NULL,   NULL,   NULL,   NULL, NULL,    "save as ct-format files",                       0 },
  { "-h",               eslARG_NONE,     FALSE,   NULL,   NULL,   NULL,   NULL, NULL,    "show help and usage",                           0 },
  { "-k",               eslARG_NONE,     FALSE,   NULL,   NULL,   NULL,   NULL, NULL,    "allow pseudoknots",                             0 },
  { "-s",               eslARG_NONE,     FALSE,   NULL,   NULL,   NULL,   NULL, NULL,    "shuffle sequences",                             0 },
  { "-t",               eslARG_NONE,     FALSE,   NULL,   NULL,   NULL,   NULL, NULL,    "print traceback",                               0 },
  { "-v",               eslARG_NONE,     FALSE,   NULL,   NULL,   NULL,   NULL, NULL,    "be verbose?",                                   0 },
  { "--infmt",          eslARG_STRING,   NULL,    NULL,   NULL,  NULL,   NULL, NULL,     "specify format",                                0 },
  {  0, 0, 0, 0, 0, 0, 0, 0, 0, 0 },
};

static char usage[] = "\
Usage: ./pknots [-options] <infile> <outfile>\n\
 where <infile> is a fasta or Stockholm file with sequences to fold\n\
 results in Stockholm format are saved to <outfile>.\n\
 \n\
"
;

int
main(int argc, char **argv)
{ 
  ESL_RANDOMNESS   *r;
  struct rnapar_2  *zkn_param = NULL;   /* rna parameters secon order + knots              */
  struct tracekn_s *tr = NULL; 	        /* traceback of a predicted RNA structure          */
  char             *seqfile;            /* input sequence file                             */
  char             *outfile;            /* where to send the output                        */
  FILE             *ofp;
  ESL_GETOPTS      *go;
  ESL_SQFILE       *sqfp;      	        /* open sequence file                              */
  ESL_SQ           *sq = NULL;
  ESL_ALPHABET     *abc = NULL;
  int               format = eslSQFILE_UNKNOWN;
  int             **icfg;               /* integer log form grammar for alignment          */
  long              L;
  CYKVAL            sc;                 /* score of predicted structure                    */

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
  char             *informat;
  int               infmt;
  int               verbose;            /* TRUE to be extremely verbose to debug           */
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
  if (esl_opt_ArgNumber(go) != 2) esl_fatal("Incorrect number of command line arguments.\n%s\n", usage);
 
  approx            = esl_opt_GetBoolean(go, "-a");
  allow_coaxials    = esl_opt_GetBoolean(go, "-c");
  allow_pseudoknots = esl_opt_GetBoolean(go, "-k");
  ctoutput          = esl_opt_GetBoolean(go, "-g");
  informat          = esl_opt_GetString(go, "--infmt");
  shuffleseq        = esl_opt_GetBoolean(go, "-s");
  traceback         = esl_opt_GetBoolean(go, "-t");
  verbose           = esl_opt_GetBoolean(go, "-v");

  seqfile = esl_opt_GetArg(go, 1);
  outfile = esl_opt_GetArg(go, 2);

  infmt = eslSQFILE_FASTA;
  if (informat != NULL)
    {
      infmt = esl_sqio_EncodeFormat(informat);
      if (infmt == eslSQFILE_UNKNOWN) 
	pk_fatal("Unrecognized file format %s\n", informat);
    }


  r = esl_randomness_CreateTimeseeded();

  /* Open the output file and the input seqfile.
   */
  abc = esl_alphabet_Create(eslRNA);
  if ((ofp = fopen(outfile, "w")) == NULL) 
    pk_fatal("failed to open %s for output", outfile);
  if (esl_sqfile_OpenDigital(abc, seqfile, infmt, NULL, &sqfp) != eslOK)
    pk_fatal("failed to open %s", seqfile);
  sq = esl_sq_CreateDigital(abc);

  /* Create the starting model 
   */
  Parameters2_Zkn(&zkn_param);
  icfg = ParamIntSCFG(zkn_param); 

  /* Print banner 
   */
  puts(banner);  
  printf("        PKNOTS %s (%s)", RELEASE, RELEASEDATE);
  printf(" using easel\n");
  printf("---------------------------------------------------\n");
  printf("Folding sequences from:  %s \n", seqfile);
  printf("---------------------------------------------------\n");
  puts("");
  
  /* OK, start folding.
   */
  while ((status = esl_sqio_Read(sqfp, sq)) == eslOK) 
    {
      L = (long)sq->n;
      
      if (verbose) {
	printf("[Original sequence:] len %ld\n", L);
	esl_sqio_Write(stdout, sq, eslSQFILE_FASTA, FALSE);
      }
      
      if (shuffleseq) esl_rsq_CShuffle(r, sq->seq, sq->seq); /* shuffle in place */
      
      if (StructurePredictkn_2IS(stdout, sq, L, zkn_param, icfg, verbose, traceback,
				 &tr, &sc, allow_coaxials, allow_pseudoknots, approx) != eslOK)
	pk_fatal("could not fold the sequence");
      
      /* convert traceback to a ss */
      if (Tracekn(tr, sq, FALSE) != eslOK) pk_fatal("could not convert to ss");
      
      /* print output */
      WriteSeqkn(ofp, abc, sq, ctoutput, zkn_param, format, shuffleseq, allow_pseudoknots, approx, sc);
      
      FreeTracekn(tr);
      esl_sq_Reuse(sq);
    }
  if      (status == eslEFORMAT) esl_fatal("Parse failed (sequence file %s)\n%s\n", sqfp->filename, sqfp->get_error(sqfp));     
  else if (status != eslEOF)     esl_fatal("Unexpected error %d reading sequence file %s", status, sqfp->filename);
 
  /* Cleanup and exit
   */
  esl_getopts_Destroy(go);
  free(zkn_param);
  FreeIntSCFG(icfg);
  esl_randomness_Destroy(r);
  esl_sq_Destroy(sq);
  esl_sqfile_Close(sqfp);
  esl_alphabet_Destroy(abc);
  fclose(ofp);
  exit (0);
}

 

