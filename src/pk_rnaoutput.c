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
 * Return:   ret_ss contains a string 0..rlen-1 containing the
 *           secondary structure. Must be free'd by caller.
 */
void
Tracekn(struct tracekn_s *tr, char *seq, int rlen, int watsoncrick, 
        char **ret_ss)  
{ 
  struct traceknstack_s *dolist;
  struct tracekn_s      *curr;
  char                  *ss;

  if ((ss = (char *) malloc (sizeof(char) * rlen+1)) == NULL)
    pk_fatal("malloc failed");
  memset(ss, '.', rlen);
  ss[rlen] = '\0';

  dolist = InitTraceknstack();
  PushTraceknstack(dolist, tr->nxtl);

  while ((curr = PopTraceknstack(dolist)) != NULL)
    {
      if (curr->type1 == dpcP || curr->type1 == dpcPS || curr->type1 == dpcPI || curr->type1 == dpcPL)
	{
	  if (! watsoncrick  ||
	      IsRNAComplement(seq[curr->emiti], seq[curr->emitj], TRUE))
	    {
	      ss[curr->emiti] = '>';
	      ss[curr->emitj] = '<';
	    }
	}

      if (curr->type2 == dpcP || curr->type2 == dpcPS || curr->type2 == dpcPI || curr->type2 == dpcPL)
	{
	  if (! watsoncrick  ||
	      IsRNAComplement(seq[curr->emitk], seq[curr->emitl], TRUE))
	    {
	      ss[curr->emitk] = '>';
	      ss[curr->emitl] = '<';
	    }
	}

      if (curr->nxtr) PushTraceknstack(dolist, curr->nxtr);
      if (curr->nxtl) PushTraceknstack(dolist, curr->nxtl);
    }

  if (esl_kh2wuss(ss, ss) != eslOK) 
    pk_fatal("Tracekn(): could not convert structure to wuss format");
  if (esl_wuss_full(ss, ss) != eslOK) 
    pk_fatal("Tracekn(): could not convert structure to wuss_full format");
 
  FreeTraceknstack(dolist);
  *ret_ss = ss;
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
void
Traceintkn(struct tracekn_s *tr, char *seq, int rlen, int watsoncrick, 
	   int **ret_ss)  
{ 
  struct traceknstack_s *dolist;
  struct tracekn_s      *curr;
  int                   *ss;
  int                    i;

  if ((ss = (int *) malloc (sizeof(int) * rlen+1)) == NULL)
    pk_fatal("malloc failed");

  for (i = 0; i < rlen; i++)
    ss[i] = -1;

  dolist = InitTraceknstack();
  PushTraceknstack(dolist, tr->nxtl);

  while ((curr = PopTraceknstack(dolist)) != NULL)
    {
      if (curr->type1 == dpcP || curr->type1 == dpcPS || curr->type1 == dpcPI || curr->type1 == dpcPL)
	{
	  if (! watsoncrick  ||
	      IsRNAComplement(seq[curr->emiti], seq[curr->emitj], TRUE))
	    {
	      ss[curr->emiti] = curr->emitj;
	      ss[curr->emitj] = curr->emiti;
	    }
	}

      if (curr->type2 == dpcP || curr->type2 == dpcPS || curr->type2 == dpcPI || curr->type2 == dpcPL )
	{
	  if (! watsoncrick  ||
	      IsRNAComplement(seq[curr->emitk], seq[curr->emitl], TRUE))
	    {
	      ss[curr->emitk] = curr->emitl;
	      ss[curr->emitl] = curr->emitk;
	    }
	}

      if (curr->nxtr) PushTraceknstack(dolist, curr->nxtr);
      if (curr->nxtl) PushTraceknstack(dolist, curr->nxtl);
    }
  FreeTraceknstack(dolist);
  *ret_ss = ss;
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
WriteSeqkn(FILE *outf, ESL_SQ *sq, int *ss, float *ret_cykpairs)
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

  seqlen = strlen(sq->seq);

  strcpy( endstr,"");
  l1 = 0;

  fprintf(outf, "NAM  %s\n", sq->name);

  fprintf(outf, "SEQ\n");

  numline = 1;                /* number seq lines w/ coords  */
  strcpy(endstr, "\n");

  for (i=0, l=0, ibase = 1, lines = 0; i < seqlen; ) {

    if (l1 < 0) 
      l1 = 0;
    else if (l1 == 0) {
      if (numline) 
	fprintf(outf,"%8d ",ibase);
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

    pos[l] = i;
    s[l]   = *(sq->seq+i);
      
    if (ss[i+1] != 0) {
      lss[l]  = ss[i+1];
      cykpairs += 1;
    }
    else 
      lss[l] = 56789;
    
    l++; i++;
    l1++;                 /* don't count spaces for width*/
    if (l1 == width || i == seqlen) {
      s[l]  = '\0';
      lss[l] = 888888;
      
      fprintf(outf, "%s\n", s);
      
      if (numline) 
	fprintf(outf,"         ");
      
      for (j=0; j<tab; j++) 
	fputc(' ',outf);
      
      for (m=0; m<l; m++)
	  if (s[m] != ' '  &&  pos[m] <= 9) 
	    fprintf(outf,"%d   ", *(pos+m));
	  else if (s[m] != ' '  && pos[m] > 9 && pos[m] <= 99) 
	    fprintf(outf,"%d  ", *(pos+m));
	  else if (s[m] != ' ') 
	    fprintf(outf,"%d ", *(pos+m));
      
      fprintf(outf,"\n");
      
      if (numline) 
	fprintf(outf,"         ");
      
      for (j = 0; j < tab; j++) 
	fputc(' ',outf);
      
      for (m = 0; m < l; m++)
	if (s[m] != ' '  && lss[m] <= 9 && lss[m] != 56789)
	  fprintf(outf,"%d   ", *(lss+m));
	else if (s[m] != ' '  && lss[m] > 9 && lss[m] <= 99 && lss[m] != 56789)
	  fprintf(outf,"%d  ", *(lss+m));
	else if (s[m] != ' ' && lss[m] != 56789)
	  fprintf(outf,"%d ", *(lss+m));
	else if (s[m] != ' ') 
	  fprintf(outf, ".   ");
      
      fprintf(outf,"%s\n",endstr);
 
      l = 0; l1 = 0;
      lines++;
      ibase = i+1;
    }
  }

  *ret_cykpairs = cykpairs/2;

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

