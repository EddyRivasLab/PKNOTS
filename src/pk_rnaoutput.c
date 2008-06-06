/* rnaoutput.c
 * derived from COVE's konings.c
 * SRE, Sun Aug 28 10:39:18 1994
 * 
 * Representation of secondary structure and secondary structural 
 * alignments using a variant of Danielle Konings' string notation.
 * 
 * See: Konings and Hogeweg, J. Mol. Biol. 207:597-614 1989
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>

#include <easel.h>
#include <esl_sqio.h>
#include <esl_wuss.h>

#include "pknots.h"
#include "pk_irredsurf.h"
#include "pk_rnaparam.h"
#include "pk_rnaoutput.h"
#include "pk_vxgraphs.h"
#include "pk_wbxgraphs.h"
#include "pk_wxgraphs.h"
#include "pk_trace.h"
#include "pk_tracemtx.h"
#include "pk_util.h"

#ifdef MEMDEBUG
#include "dbmalloc.h"
#endif

static void ct_output(FILE *ofp, char *seq, int *ct, int j, int d);
static void param_output(FILE *outf, struct rnapar_2 zkn_param, 
			 int format, int shuffleseq, int allow_pseudoknots, 
			 int approx, CYKVAL sc, float cykpairs);

/* Function: Tracekn()
 * 
 * Purpose:  From a traceback tree of seq, produce a
 *           secondary structure string. ">" and "<" are
 *           used for pairwise emissions; "." for single-stranded stuff.
 *           Note that structure is defined by pairwise emissions,
 *           not by Watson-Crick-isms and stacking rules.
 *
 *           Constructed for the one-hole algorithm.
 *           
 * Args:     tr          - traceback structure
 *           seq         - sequence, 0..rlen-1
 *           rlen        - length of seq and returned ss string
 *           watsoncrick - TRUE to annotate canonical pairs only
 *           ret_ss      - RETURN: alloc'ed secondary structure string
 *
 * Return:   
 */
void
Tracekn(struct tracekn_s *tr, ESL_DSQ *dsq, int rlen, int watsoncrick, char *ss)  
{ 
  struct traceknstack_s *dolist;
  struct tracekn_s      *curr;
  int                    j;

  /* initialize */
  ss[0] = '\0';
  for (j = 1; j <= rlen; j ++)
    ss[j] = '.';
  
  dolist = InitTraceknstack();
  PushTraceknstack(dolist, tr->nxtl);

  while ((curr = PopTraceknstack(dolist)) != NULL)
    {
       if (curr->type1 == dpcP || curr->type1 == dpcPS || curr->type1 == dpcPI || curr->type1 == dpcPL)
	{
	  if (! watsoncrick  ||
	      IsRNAComplementDigital(dsq[curr->emiti+1], dsq[curr->emitj+1], TRUE))
	    {
	      ss[curr->emiti+1] = '<';
	      ss[curr->emitj+1] = '>';
	    }
	}

      if (curr->type2 == dpcP || curr->type2 == dpcPS || curr->type2 == dpcPI || curr->type2 == dpcPL)
	{
	  if (! watsoncrick  ||
	      IsRNAComplementDigital(dsq[curr->emitk+1], dsq[curr->emitl+1], TRUE))
	    {
	      ss[curr->emitk+1] = '<';
	      ss[curr->emitl+1] = '>';
	    }
	}

      if (curr->nxtr) PushTraceknstack(dolist, curr->nxtr);
      if (curr->nxtl) PushTraceknstack(dolist, curr->nxtl);
    }

  if (esl_wuss_full(ss, ss) != eslOK) 
    pk_fatal("Tracekn(): could not convert structure to wuss_full format");
  
   FreeTraceknstack(dolist);
}

/* Function: WriteSeqkn()
 * 
 * Purpose:  writes in outf file the Tracekn
 *           
 *           Constructed for the one-hole algorithm.
 *
 *           "." are used for single-stranded stuff.
 *           Note that structure is defined by pairwise emissions,
 *           not by Watson-Crick-isms and stacking rules.
 *           
 */
int
WriteSeqkn(FILE *outf, ESL_ALPHABET *abc, ESL_SQ *sq, int *ct, int ctoutput, 
	   struct rnapar_2 zkn_param, int format, int shuffleseq, 
	   int allow_pseudoknots, int approx, CYKVAL sc)
{
  int   numline = 0;
  int   lines = 0, spacer = 4, width = 20, tab = 0;
  int   i, j, l, l1, ibase, m;
  char  endstr[10]; 
  char  s[100];			/* buffer for sequence  */
  int   pos[100];		/* buffer for structure */
  int   lss[100];		/* buffer for secondary structure */
  int   seqlen;   
  int   cykpairs = 0;

  if (esl_sq_Textize(sq)!= eslOK)  pk_fatal("coudnot textize %s", sq->name); /* sq is now text mode */

  seqlen = sq->n;

  strcpy( endstr,"");
  l1 = 0;

  printf("NAM  %s\n", sq->name);

  printf("SEQ\n");

  numline = 1;                /* number seq lines w/ coords  */
  strcpy(endstr, "\n");

  for (i=0, l=0, ibase = 1, lines = 0; i < seqlen; ) {

    if (l1 < 0) 
      l1 = 0;
    else if (l1 == 0) {
      if (numline) 
	printf("%8d ",ibase);
      for (j=0; j<tab; j++) 
	fputc(' ',stdout);
    }
      
    if (spacer != 0 && l%spacer == 1) {
      s[l] = ' '; 
      lss[l] = 1234; 
      l++;
    }
    
    if (spacer != 0 && l%spacer == 2) {
      s[l] = ' '; 
      lss[l] = 1234; 
      l++;
    }
      
    if (spacer != 0 && l%spacer == 3) {
      s[l] = ' '; 
      lss[l] = 1234;
      l++;
    }

    pos[l] = i+1;
    s[l]   = *(sq->seq+i);
      
    if (sq->ss[i+1] != 0) {
      lss[l]  = sq->ss[i+1];
      cykpairs += 1;
    }
    else 
      lss[l] = 56789;
    
    l++; i++;
    l1++;                 /* don't count spaces for width*/
    if (l1 == width || i == seqlen) {
      s[l]  = '\0';
      lss[l] = 888888;
      
      printf("%s\n", s);
      
      if (numline) 
	printf("         ");
      
      for (j=0; j<tab; j++) 
	fputc(' ',stdout);
      
      for (m=0; m<l; m++)
	  if (s[m] != ' '  &&  pos[m] <= 9) 
	    printf("%d   ", *(pos+m));
	  else if (s[m] != ' '  && pos[m] > 9 && pos[m] <= 99) 
	    printf("%d  ", *(pos+m));
	  else if (s[m] != ' ') 
	    printf("%d ", *(pos+m));
      
      printf("\n");
      
      if (numline) 
	printf("         ");
      
      for (j = 0; j < tab; j++) 
	fputc(' ',stdout);
      
      for (m = 0; m < l; m++)
	if (s[m] != ' '  && lss[m] <= 9 && lss[m] != 56789)
	  printf("%d   ", *(lss+m));
	else if (s[m] != ' '  && lss[m] > 9 && lss[m] <= 99 && lss[m] != 56789)
	  printf("%d  ", *(lss+m));
	else if (s[m] != ' ' && lss[m] != 56789)
	  printf("%d ", *(lss+m));
	else if (s[m] != ' ') 
	  printf(".   ");
      
      printf("%s\n",endstr);
 
      l = 0; l1 = 0;
      lines++;
      ibase = i+1;
    }
  }

  cykpairs = cykpairs/2;

  param_output(stdout, zkn_param, format, shuffleseq, allow_pseudoknots, approx, sc, cykpairs);
  
  /* write ss to a stockholm file or ctfile */
  if (ctoutput) ct_output(outf, sq->seq, ct, seqlen-1, seqlen-1);
  else          esl_sqio_Write(outf, sq, eslMSAFILE_STOCKHOLM);

  if (esl_sq_Digitize(abc, sq)!= eslOK)  pk_fatal("coudnot digitize %s", sq->name); /* sq is digital again */

  return lines;
} 



/* Function: IsRNAComplement()
 * 
 * Purpose:  Returns TRUE if sym1, sym2 are Watson-Crick complementary.
 *           If allow_gu is TRUE, GU pairs also return TRUE.
 */
int
IsRNAComplement(char sym1, char sym2, int allow_gu)
{
  sym1 = toupper(sym1);
  sym2 = toupper(sym2);

  if ((sym1 == 'A' && sym2 == 'U') ||
      (sym1 == 'C' && sym2 == 'G') ||
      (sym1 == 'G' && sym2 == 'C') ||
      (sym1 == 'U' && sym2 == 'A') ||
      (allow_gu && sym1 == 'G' && sym2 == 'U') ||
      (allow_gu && sym1 == 'U' && sym2 == 'G'))
    return TRUE;
  else
    return FALSE;
}

int
IsRNAComplementDigital(int sym1, int sym2, int allow_gu)
{
  if ((sym1 == 0 && sym2 == 3) ||
      (sym1 == 1 && sym2 == 2) ||
      (sym1 == 2 && sym2 == 1) ||
      (sym1 == 3 && sym2 == 1) ||
      (allow_gu && sym1 == 2 && sym2 == 3) ||
      (allow_gu && sym1 == 3 && sym2 == 2))
    return TRUE;
  else
    return FALSE;
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


