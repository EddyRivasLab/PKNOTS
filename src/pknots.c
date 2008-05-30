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
              
static void ct_output(FILE *ofp, char *seq, int *ct, int j, int d);
static void param_output(FILE *outf, struct rnapar_2 zkn_param, 
			 int format, int shuffleseq, int allow_pseudoknots, 
			 int approx, CYKVAL sc, float cykpairs);

static char banner[] = 
"PKNOTS: optimal minimum-energy RNA folding with pseudoknots and coaxial energies";

static ESL_OPTIONS options[] = {
 /* name                type             default  env_var range   toggles req   incompat help                                      docgroup */
  { "-c",               eslARG_NONE,     FALSE,   NULL,   NULL,   NULL,   NULL, NULL,    "add L^5 coaxials (V6)",                  0 },
  { "-g",               eslARG_NONE,     NULL,    NULL,   NULL,   NULL,   NULL, NULL,    "save as ct-format files",                0 },
  { "-h",               eslARG_NONE,     FALSE,   NULL,   NULL,   NULL,   NULL, NULL,    "show help and usage",                    0 },
  { "-k",               eslARG_NONE,     FALSE,   NULL,   NULL,   NULL,   NULL, NULL,    "allow pseudoknots",                      0 },
  { "-s",               eslARG_NONE,     FALSE,   NULL,   NULL,   NULL,   NULL, NULL,    "shuffle sequences",                      0 },
  { "-t",               eslARG_NONE,     FALSE,   NULL,   NULL,   NULL,   NULL, NULL,    "print traceback",                        0 },
  { "-v",               eslARG_NONE,     FALSE,   NULL,   NULL,   NULL,   NULL, NULL,    "be verbose?",                            0 },
  { "--infmt",          eslARG_STRING,   NULL,    NULL,   NULL,  NULL,   NULL, NULL,     "specify format",                         0 },
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
  struct rnapar_2   zkn_param;          /* rna parameters secon order + knots              */
  struct tracekn_s *tr; 	        /* traceback of a predicted RNA structure          */
  char             *seqfile;            /* input sequence file                             */
  char             *outfile;            /* where to send the output                        */
  FILE             *ofp;
  ESL_GETOPTS      *go;
  ESL_SQFILE       *sqfp;      	        /* open sequence file                              */
  ESL_SQ           *sq;
  ESL_ALPHABET     *abc;
  int               format = eslSQFILE_UNKNOWN;
  int              *ct = NULL;
  int             **icfg;               /* integer log form grammar for alignment          */
  int               L;
  CYKVAL            sc;                 /* score of predicted structure                    */
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
      infmt = esl_sqio_FormatCode(informat);
      if (infmt == eslSQFILE_UNKNOWN) 
	pk_fatal("Unrecognized file format %s\n", informat);
    }
  esl_getopts_Destroy(go);

  r = esl_randomness_CreateTimeseeded();

  /*********************************************** 
   * Create the starting model 
   ***********************************************/
  abc = esl_alphabet_Create(eslRNA);
  Parameters2_Zkn(&zkn_param);
  icfg = ParamIntSCFG(&zkn_param); 

  /* Open the output file and the input seqfile.
   */
  if ((ofp = fopen(outfile, "w")) == NULL) 
    pk_fatal("failed to open %s for output", outfile);
  if (esl_sqfile_Open(seqfile, infmt, NULL, &sqfp) != eslOK)
    pk_fatal("failed to open %s", seqfile);
  sq = esl_sq_Create();

  /*********************************************** 
   * Print banner
   ***********************************************/
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
     L = sq->n;
      
      if (verbose) {
	printf("[Original sequence:]\n");
	esl_sqio_Write(stdout, sq, eslSQFILE_FASTA);
;
      }

      if (shuffleseq) esl_rsq_CShuffle(r, sq->seq, sq->seq); /* shuffle in place */

      if (StructurePredictkn_2IS(ofp, abc, sq, L, &zkn_param, icfg, verbose, traceback,
				 &tr, &sc, allow_coaxials, allow_pseudoknots, approx) != eslOK)
	pk_fatal("could not fold the sequence");
     
     /* convert traceback to a ss */
     Tracekn(tr, sq->seq, L, FALSE, &(sq->ss));
          
     /* the CT array*/
     ESL_ALLOC(ct, sizeof(int) * (L+1));
     if (esl_wuss2ct(sq->ss, L, ct) != eslOK)
       pk_fatal("could not generate ctfile");

     /* write ss to a stockholm file or ct output*/
     if (ctoutput) ct_output(ofp, sq->seq, ct, L-1, L-1);
     else esl_sqio_Write(ofp, sq, eslMSAFILE_STOCKHOLM);

     /* print old style output to stdout*/
     WriteSeqkn(stdout, sq, ct, &cykpairs);
     param_output(stdout, zkn_param, format, shuffleseq, allow_pseudoknots, approx, sc, cykpairs);

     FreeTracekn(tr);
     free(sq->dsq); sq->dsq = NULL;
     free(sq->ss);  sq->ss  = NULL; 
      free(ct);
      esl_sq_Reuse(sq);
    }
  
  /*********************************************** 
   * Cleanup and exit
   ***********************************************/

  FreeIntSCFG(icfg);
  esl_randomness_Destroy(r);
  esl_sq_Destroy(sq);
  esl_sqfile_Close(sqfp);
  esl_alphabet_Destroy(abc);
  fclose(ofp);
  exit (0);

 ERROR:
  return status;
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
ct_output(FILE *ofp, char *seq, int *ct, int j, int d)
{  
  int i;          /* initial position of fragment. 0,..,d-1 */
  int iabs;       /* absolute value of initial position     */

  fprintf(ofp,"\n ct_output \n");
  fprintf(ofp,"----------------------------------------------------------------------\n");

  for (i = 0; i <= d; i++) {
    iabs = i + j - d;
    fprintf(ofp, "%5d %c   %5d %4d %4d %4d\n",
            i+1, seq[iabs], i, i+2, ct[iabs], iabs+1);
  }
  
  fprintf(ofp, "\n");
  
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
	     CYKVAL sc, float cykpairs)
{  
      fprintf(ofp,"----------------------------------------------------------------------\n");
      
      fprintf(ofp,"   Log odds score:    %8d \n",  sc);
      fprintf(ofp,"energy (kcal/mol):  %10.2f \n", -(float)(sc)/(INTSCALE*wsf));
      fprintf(ofp,"number of pairs found:      %4.2f \n", cykpairs);
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


