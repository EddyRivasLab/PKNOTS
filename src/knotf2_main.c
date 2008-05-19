/* knotf_main.c
 *
 * Fold RNAs including knots.
 */
                                                                       
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>

#include "version.h"
#include "cfg.h"
#include "proto.h"
#include "protovx.h"
#include "protowbx.h"
#include "protowx.h"
#include "squid.h" 

#ifdef MEMDEBUG
#include "dbmalloc.h"
#endif
              
static struct opt_s OPTIONS[] = {
  { "-g",        TRUE,  sqdARG_NONE  },
  { "-h",        TRUE,  sqdARG_NONE  },
  { "-c",        TRUE,  sqdARG_NONE  },
  { "-k",        TRUE,  sqdARG_NONE  },
  { "-o",        TRUE,  sqdARG_STRING},
  { "-s",        TRUE,  sqdARG_NONE  },
  { "-t",        TRUE,  sqdARG_NONE  },
  { "-v",        TRUE,  sqdARG_NONE  },
};
                
#define NOPTIONS (sizeof(OPTIONS) / sizeof(struct opt_s))

static char usage[]  = "\
Usage: pknots [-options] <seqfile in>\n\
where options are:\n\
   -g            : ct output\n\
   -h            : print short help and usage info\n\
   -c            : add V6 the time consuming (N^5) coaxials\n\
   -k            : allow pseudoknots\n\
   -o <outfile>  : direct structure-annotated sequence to <outfile>\n\
   -s            : shuffle the sequences\n\
   -t            : print traceback\n\
   -v            : verbose debugging output\n\
";

static char banner[] = 
"PKNOTS: optimal minimum-energy RNA folding with pseudoknots and coaxial energies";

static void ct_output(FILE *ofp, char *seq, SQINFO sqinfo, int *ss, int j, int d);
static void param_output(FILE *outf, struct rnapar_2 zkn_param, 
			 int format, int shuffleseq, int allow_pseudoknots, 
			 int approx, int score, float pairs, float cykpairs);

int
main(int argc, char **argv)
{ 
  struct rnapar_2 zkn_param;    /* rna parameters secon order + knots              */
  struct tracekn_s      *tr;	/* traceback of a predicted RNA structure          */
  char                 *seq;	/* sequence to fold                                */
  SQINFO             sqinfo;    /* info structures for seq                         */
  char                  *cc;    /* secondary structure (should go to SQINFO)       */
  int                   *ss;    /* secondary structure (should go to SQINFO)       */
  char             *seqfile;    /* input sequence file                             */
  SQFILE              *sqfp;	/* open sequence file                              */
  int                format;    /* format of sequence file                         */
  int                **icfg;    /* integer log form grammar for alignment          */
  int                 score;    /* score of predicted structure                    */
  float               pairs;    /* number of base paired as input                  */
  float            cykpairs;    /* number of base paired of predicted structure    */

  int               verbose;    /* TRUE to be extremely verbose to debug           */
  char             *outfile;    /* where to send the output                        */
  FILE                 *ofp;	/* open output file                                */

  char             *optname;
  char              *optarg; 
  int                optind;	
  int             traceback;    /* TRUE print traceback                            */
  int              ctoutput;    /* TRUE print ctoutput                             */
  int        allow_coaxials;    /* TRUE add V6 coaxial diagrams                    */
  int            shuffleseq;    /* TRUE to shuffle the sequences                   */
  int     allow_pseudoknots;	/* TRUE  == include pseudoknots                    *
				 * FALSE == no pseudoknots                         */
  int                approx;    /* TRUE  == external pseudoknot approximation      *
				 *          exclude diagrams V7-V10 and WB9-WB10   *
				 * FALSE == full pseudoknot model                  *
				            include diagrams V7-V10 and WB9-WB10   */
  int                     i;

#ifdef MEMDEBUG
  unsigned long histid1, histid2, orig_size, current_size;
  orig_size = malloc_size(&histid1);
#endif

  /*********************************************** 
   * Parse command line
   ***********************************************/
                                 
  verbose           = FALSE;
  outfile           = NULL;
  traceback         = FALSE;    /* TRUE  ==  print traceback         */
  ctoutput          = FALSE;    /* TRUE  ==  print ctoutput          */
  shuffleseq        = FALSE;    /* TRUE  ==  shuffle the sequence    */
  allow_coaxials    = FALSE;    /* TRUE  ==  include V6     diagrams */
  allow_pseudoknots = FALSE;    /* TRUE  ==  pseudoknots             */
  approx            = FALSE;    /* FALSE ==  includes V7-V10 and WB9-WB10 diagrams */

  while (Getopt(argc, argv, OPTIONS, NOPTIONS, usage,
		&optind, &optname, &optarg))
    {
           if (strcmp(optname, "-c") == 0) allow_coaxials    = TRUE;
      else if (strcmp(optname, "-g") == 0) ctoutput          = TRUE;
      else if (strcmp(optname, "-k") == 0) allow_pseudoknots = TRUE;
      else if (strcmp(optname, "-o") == 0) outfile           = optarg;
      else if (strcmp(optname, "-s") == 0) shuffleseq        = TRUE;
      else if (strcmp(optname, "-t") == 0) traceback         = TRUE;
      else if (strcmp(optname, "-v") == 0) verbose           = TRUE;
      
      else if (strcmp(optname, "-h") == 0) 
	{
	  puts(banner); 
	  printf("         PKNOT %s (%s)", RELEASE, RELEASEDATE);
	  printf(" using squid %s (%s)\n", squid_version, squid_date);
	  puts(usage);
	  exit(0);
	}
    }

  if (argc - optind != 1)
    Die("Incorrect number of command line arguments.\n%s\n", usage);
  seqfile  = argv[optind]; 

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
 
  /*********************************************** 
   * Create the starting model 
   ***********************************************/
  Parameters2_Zkn(&zkn_param);
  icfg = ParamIntSCFG(&zkn_param); 
 
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

  if (outfile != NULL) fclose(ofp);
  FreeIntSCFG(icfg);
  SeqfileClose(sqfp);
  return EXIT_SUCCESS;
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

