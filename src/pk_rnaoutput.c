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

static void ct_output(FILE *ofp, char *seq, char *seqname, int *ct, int j, int d);
static void param_output(FILE *outf, struct rnapar_2 *zkn_param, 
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
int
Tracekn(struct tracekn_s *tr, ESL_SQ *sq, int watsoncrick)  
{ 
  struct traceknstack_s *dolist;
  struct tracekn_s      *curr;
  int                    j;
  int                    status;

  ESL_ALLOC(sq->ss, sizeof(char) * (sq->salloc));

  /* initialize ss */
  sq->ss[0] = '\0';
  for (j = 1; j <= sq->n; j ++)
    sq->ss[j] = '.';
  sq->ss[sq->n+1] = '\0';
  
  dolist = InitTraceknstack();
  PushTraceknstack(dolist, tr->nxtl);

  while ((curr = PopTraceknstack(dolist)) != NULL)
    {
      if (curr->emiti < 0 ||curr->emiti >= sq->n)  
	pk_fatal("Tracekn(): bad traceback");
      if (curr->emitj < 0 ||curr->emitj >= sq->n)  
	pk_fatal("Tracekn(): bad traceback");
      if (curr->emitk < 0 ||curr->emitk >= sq->n)  
	pk_fatal("Tracekn(): bad traceback");
      if (curr->emitl < 0 ||curr->emitl >= sq->n)  
	pk_fatal("Tracekn(): bad traceback");
      
       if (curr->type1 == dpcP || curr->type1 == dpcPS || curr->type1 == dpcPI || curr->type1 == dpcPL)
	{
	  if (! watsoncrick  ||
	      IsRNAComplementDigital(sq->dsq[curr->emiti+1], sq->dsq[curr->emitj+1], TRUE))
	    {
	      sq->ss[curr->emiti+1] = '<';
	      sq->ss[curr->emitj+1] = '>';
	    }
	}

      if (curr->type2 == dpcP || curr->type2 == dpcPS || curr->type2 == dpcPI || curr->type2 == dpcPL)
	{
	  if (! watsoncrick  ||
	      IsRNAComplementDigital(sq->dsq[curr->emitk+1], sq->dsq[curr->emitl+1], TRUE))
	    {
	      sq->ss[curr->emitk+1] = '<';
	      sq->ss[curr->emitl+1] = '>';
	    }
	}

      if (curr->nxtr) PushTraceknstack(dolist, curr->nxtr);
      if (curr->nxtl) PushTraceknstack(dolist, curr->nxtl);
    }

  if (esl_wuss_full(sq->ss+1, sq->ss+1) != eslOK) 
    pk_fatal("Tracekn(): could not convert structure to wuss_full format");
  
   FreeTraceknstack(dolist);
   return eslOK;

 ERROR:
  return status;
}

/* Function: Traceintkn()
 * 
 * Purpose:  Convert a secondary structure string to an array of integers
 *           representing what position each position is base-paired 
 *           to (0..len-1), or -1 if none. 
 *           
 *           Constructed for the one-hole algorithm.
 *
 *           "." are used for single-stranded stuff.
 *           Note that structure is defined by pairwise emissions,
 *           not by Watson-Crick-isms and stacking rules.
 *           
 * Args:     tr          - traceback structure
 *           seq         - sequence, 0..rlen-1
 *           rlen        - length of seq and returned ss string
 *           watsoncrick - TRUE to annotate canonical pairs only
 *           ret_ss      - RETURN: alloc'ed secondary structure string
 *
 * Return:   ret_ss contains a string 0..rlen-1 containing the
 *           secondary structure. Must be free'd by caller.
 */
int
Traceintkn(struct tracekn_s *tr, ESL_SQ *sq, int watsoncrick, 
	   int **ret_ct)  
{ 
  struct traceknstack_s *dolist = NULL;
  struct tracekn_s      *curr = NULL;
  int                   *ct = NULL;
  int                    i;
  int                    status;

  ESL_ALLOC(ct, sizeof(int) * (sq->n+1));

  for (i = 0; i <= sq->n; i++)
    ct[i] = 0;

  dolist = InitTraceknstack();
  PushTraceknstack(dolist, tr->nxtl);

  while ((curr = PopTraceknstack(dolist)) != NULL)
    {
      if (curr->type1 == dpcP || curr->type1 == dpcPS || curr->type1 == dpcPI || curr->type1 == dpcPL)
	{
	  if (! watsoncrick  ||
	      IsRNAComplement(sq->dsq[curr->emiti+1], sq->dsq[curr->emitj+1], TRUE))
	    {
	      ct[curr->emiti+1] = curr->emitj+1;
	      ct[curr->emitj+1] = curr->emiti+1;
	    }
	}


      if (curr->type2 == dpcP || curr->type2 == dpcPS || curr->type2 == dpcPI || curr->type2 == dpcPL )
	{
	  if (! watsoncrick  ||
	      IsRNAComplement(sq->dsq[curr->emitk+1], sq->dsq[curr->emitl+1], TRUE))
	    {
	      ct[curr->emitk+1] = curr->emitl+1;
	      ct[curr->emitl+1] = curr->emitk+1;
	    }
	}

      if (curr->nxtr) PushTraceknstack(dolist, curr->nxtr);
      if (curr->nxtl) PushTraceknstack(dolist, curr->nxtl);
    }
  FreeTraceknstack(dolist);
  *ret_ct = ct;
  return eslOK;
  
 ERROR:
  if (dolist != NULL) FreeTraceknstack(dolist);
  if (ct != NULL) free(ct);
  return status;
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
WriteSeqkn(FILE *outf, ESL_ALPHABET *abc, ESL_SQ *sq, struct tracekn_s *tr, int ctoutput, 
	   struct rnapar_2 *zkn_param, int format, int shuffleseq, 
	   int allow_pseudoknots, int approx, CYKVAL sc)
{
  int    *ct = NULL;
  int     numline = 0;
  int     lines = 0, spacer = 4, width = 20, tab = 0;
  int     i, j, l, l1, ibase, m;
  char    endstr[10]; 
  char    s[100];			/* buffer for sequence  */
  int     pos[100];		/* buffer for structure */
  int     lss[100];		/* buffer for secondary structure */
  int     cykpairs = 0;
  int     L;
  int     status;

  /* the CT array (starts at 1) */
  L = strlen(sq->ss+1);
  if ((status = Traceintkn(tr, sq, FALSE, &ct)) != eslOK) goto ERROR;
  
  /* textize sequence for output */
  if (esl_sq_Textize(sq) != eslOK)  pk_fatal("couldnot textize %s", sq->name); /* sq is now text mode */

   /* write ss to output */ 
  if (ctoutput) ct_output(outf, sq->seq, sq->name, ct, sq->n-1, sq->n-1);
  else {
    
    strcpy( endstr,"");
    l1 = 0;
    
    fprintf(outf, "NAM  %s\n", sq->name);
    
    fprintf(outf, "SEQ\n");
    
    numline = 1;                /* number seq lines w/ coords  */
    strcpy(endstr, "\n");
    
    for (i=0, l=0, ibase = 1, lines = 0; i < sq->n; ) {
      
      if (l1 < 0) 
	l1 = 0;
      else if (l1 == 0) {
	if (numline) 
	  fprintf(outf, "%8d ",ibase);
	for (j=0; j<tab; j++) 
	  fputc(' ',outf);
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
      
      if (ct[i+1] != 0) {
	lss[l]  = ct[i+1];
	cykpairs += 1;
      }
      else 
	lss[l] = 56789;
      
      l++; i++;
      l1++;                 /* don't count spaces for width*/
      if (l1 == width || i == sq->n) {
	s[l]  = '\0';
	lss[l] = 888888;
	
	fprintf(outf, "%s\n", s);
	
	if (numline) 
	  fprintf(outf, "         ");
	
	for (j=0; j<tab; j++) 
	  fputc(' ',outf);
	
	for (m=0; m<l; m++)
	  if (s[m] != ' '  &&  pos[m] <= 9) 
	    fprintf(outf, "%d   ", *(pos+m));
	  else if (s[m] != ' '  && pos[m] > 9 && pos[m] <= 99) 
	    fprintf(outf, "%d  ", *(pos+m));
	  else if (s[m] != ' ') 
	    fprintf(outf, "%d ", *(pos+m));
	
	fprintf(outf, "\n");
	
	if (numline) 
	  fprintf(outf, "         ");
	
	for (j = 0; j < tab; j++) 
	  fputc(' ',outf);
	
	for (m = 0; m < l; m++)
	  if (s[m] != ' '  && lss[m] <= 9 && lss[m] != 56789)
	    fprintf(outf, "%d   ", *(lss+m));
	  else if (s[m] != ' '  && lss[m] > 9 && lss[m] <= 99 && lss[m] != 56789)
	    fprintf(outf, "%d  ", *(lss+m));
	  else if (s[m] != ' ' && lss[m] != 56789)
	    fprintf(outf, "%d ", *(lss+m));
	  else if (s[m] != ' ') 
	    fprintf(outf, ".   ");
	
	fprintf(outf, "%s\n",endstr);
	
	l = 0; l1 = 0;
	lines++;
	ibase = i+1;
      }
    }

    param_output(outf, zkn_param, format, shuffleseq, allow_pseudoknots, approx, sc, cykpairs);
  }


  
  cykpairs = cykpairs/2;

  /* digitize back */
  if (esl_sq_Digitize(abc, sq)!= eslOK)  pk_fatal("coudnot digitize %s", sq->name); /* sq is digital again */

  free(ct);
  return lines;

 ERROR:
  return status;
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
ct_output(FILE *ofp, char *seq, char *seqname, int *ct, int j, int d)
{  
  int i;          /* initial position of fragment. 0,..,d-1 */
  int iabs;       /* absolute value of initial position     */

  fprintf(ofp,"\n ct_output seq: %s\n", seqname);
  fprintf(ofp,"----------------------------------------------------------------------\n");

  for (i = 0; i <= d; i++) {
    iabs = i + j - d;
    fprintf(ofp, "%5d %c   %5d %4d %4d %4d\n",
            i+1, seq[iabs+1], i, i+2, ct[iabs+1], iabs+1);
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
param_output(FILE *ofp, struct rnapar_2 *zkn_param, int format, int shuffleseq, int allow_pseudoknots, int approx, 
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
	fprintf(ofp," P1   parameter:    %10.2f \n", -(float)zkn_param->P1/(INTSCALE*wsf));
	fprintf(ofp," P2   parameter:    %10.2f \n", -(float)zkn_param->P2/(INTSCALE*wsf));
	fprintf(ofp," P3   parameter:    %10.2f \n", -(float)zkn_param->P3/(INTSCALE*wsf));
	fprintf(ofp," P4   parameter:    %10.2f \n", -(float)zkn_param->P4/(INTSCALE*wsf));
	fprintf(ofp," P5   parameter:    %10.2f \n", -(float)zkn_param->P5/(INTSCALE*wsf));
	fprintf(ofp," P6   parameter:    %10.2f \n", -(float)zkn_param->P6/(INTSCALE*wsf));
	fprintf(ofp," P10  parameter:    %10.2f \n", -(float)zkn_param->P10/(INTSCALE*wsf));
	fprintf(ofp,"----------------------------------------------------------------------\n");
      }
      else {
	fprintf(ofp,"Allow pseudoknots? YES\n");
	fprintf(ofp,"----------------------------------------------------------------------\n");
	fprintf(ofp,"Parameters (energy units, kcal/mol)\n");
	fprintf(ofp," wkn  parameter:    %10.2f (pseudoknots weight)\n", (float)wkn/wsf);
	fprintf(ofp," P1   parameter:    %10.2f \n", -(float)zkn_param->P1/(INTSCALE*wsf));
	fprintf(ofp," P2   parameter:    %10.2f \n", -(float)zkn_param->P2/(INTSCALE*wsf));
	fprintf(ofp," P3   parameter:    %10.2f \n", -(float)zkn_param->P3/(INTSCALE*wsf));
	fprintf(ofp," P4   parameter:    %10.2f \n", -(float)zkn_param->P4/(INTSCALE*wsf));
	fprintf(ofp," P5   parameter:    %10.2f \n", -(float)zkn_param->P5/(INTSCALE*wsf));
	fprintf(ofp," P5P  parameter:    %10.2f \n", -(float)zkn_param->P5P/(INTSCALE*wkn));
	fprintf(ofp," P6   parameter:    %10.2f \n", -(float)zkn_param->P6/(INTSCALE*wsf));
	fprintf(ofp," P6P  parameter:    %10.2f \n", -(float)zkn_param->P6P/(INTSCALE*wkn));
	fprintf(ofp," P10  parameter:    %10.2f \n", -(float)zkn_param->P10/(INTSCALE*wsf));
	fprintf(ofp," P10P parameter:    %10.2f \n", -(float)zkn_param->P10P/(INTSCALE*wkn));
	fprintf(ofp," P11  parameter:    %10.2f \n", -(float)zkn_param->P11/(INTSCALE*wsf));
	fprintf(ofp," P12  parameter:    %10.2f \n", -(float)zkn_param->P12/(INTSCALE*wsf));
	
	if (approx){
	  fprintf(ofp," P13  parameter:    infinity\n");
	  fprintf(ofp," (external pseudoknots approximation)\n");
	}
	else {
	  fprintf(ofp," P13  parameter:    %10.2f \n", -(float)zkn_param->P13/(INTSCALE*wsf));
	  fprintf(ofp," full pseudoknot model.\n");
	fprintf(ofp,"----------------------------------------------------------------------\n");
	}
      }
}      


