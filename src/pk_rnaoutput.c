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

#include "cfg.h"
#include "proto.h"
#include "protovx.h"
#include "protowbx.h"
#include "protowx.h"
#include "squid.h"

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
    Die("malloc failed");
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
    Die("malloc failed");

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
WriteSeqkn(FILE *outf, char *seq, SQINFO *sqinfo, int *ss, int format, float *ret_pairs, float *ret_cykpairs)
{
  int   numline = 0;
  int   lines = 0, spacer = 4, width = 20, tab = 0;
  int   i, j, l, l1, ibase, m;
  char  endstr[10]; 
  char  s[100];			/* buffer for sequence  */
  int   pos[100];		/* buffer for structure */
  int   lss[100];		/* buffer for secondary structure */
  int   checksum = 0;
  int   seqlen;   
  int   dostruc;		/* TRUE to print structure lines*/
  int   pairs = 0;
  int   cykpairs = 0;

  dostruc    = FALSE;		
  seqlen     = (sqinfo->flags & SQINFO_LEN) ? sqinfo->len : strlen(seq);

  strcpy( endstr,"");
  l1 = 0;

  /* 10Nov91: write this out in all possible formats: */
  checksum = GCGchecksum(seq, seqlen);

  fprintf(outf, "NAM  %s\n", sqinfo->name);

  if (sqinfo->flags & (SQINFO_ID | SQINFO_ACC | SQINFO_START | SQINFO_STOP | SQINFO_OLEN))
    fprintf(outf, "SRC  %s %s %d..%d::%d\n",
	    (sqinfo->flags & SQINFO_ID)    ? sqinfo->id     : "-",
	    (sqinfo->flags & SQINFO_ACC)   ? sqinfo->acc    : "-",
	    (sqinfo->flags & SQINFO_START) ? sqinfo->start  : 0,
	    (sqinfo->flags & SQINFO_STOP)  ? sqinfo->stop   : 0,
	    (sqinfo->flags & SQINFO_OLEN)  ? sqinfo->olen   : 0);

  if (sqinfo->flags & SQINFO_DESC)
    fprintf(outf, "DES  %s\n", sqinfo->desc);

  if (sqinfo->flags & SQINFO_SS) {
    fprintf(outf, "SEQ  +SS\n");
    dostruc = TRUE;	/* print structure lines too */
  }
  else
    fprintf(outf, "SEQ\n");

  numline = 1;                /* number seq lines w/ coords  */
  strcpy(endstr, "\n");

  if(format == kSquid || format == kSelex)
    for (i = 0; i < seqlen; i++) 
      if (sqinfo->ss[i] != '.') 
	pairs += 1;
  
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
    s[l]   = *(seq+i);
      
    if (ss[i] != -1) {
      lss[l]  = ss[i];
      cykpairs += 1;
    }
    else 
      lss[l] = 56789;
    
    l++; i++;
    l1++;                 /* don't count spaces for width*/
    if (l1 == width || i == seqlen) {
      s[l]  = '\0';
      lss[l] = 888888;
      
      if (dostruc) {
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
      }
      else {
	if (i == seqlen) fprintf(outf,"%s%s\n",s,endstr);
	else fprintf(outf,"%s\n",s);
      }
      l = 0; l1 = 0;
      lines++;
      ibase = i+1;
    }
  }

  *ret_pairs    = pairs/2;
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

void
CalculatePairs(SQINFO sqinfo, int *ss, int format, float *ret_pairs, float *ret_cykpairs)
{
  int i;
  float pairs = 0.;
  float cykpairs = 0.;

  for (i = 0; i < sqinfo.len; i++) {
    if ((format == kSquid || format == kSelex ) && sqinfo.ss[i] != '.') 
      pairs += 1.;
    if (ss[i] != -1) 
      cykpairs += 1.;
  }

  *ret_pairs    = pairs/2.;
  *ret_cykpairs = cykpairs/2.;
}



void
PrintCtSeq(FILE *ofp, SQINFO *sqinfo, char *seq, int start, int L, char *ss)
{
  int   i, pos;
  int   line = 50;
  int   nlines = 0;

  pos = start;
  fprintf (ofp, "\n");
  
  while (L/(pos+1)) {
     fprintf (ofp, "\n\t");
     for (i = 0; i < line; i++)
       if (i == 0)         fprintf (ofp, "%15.15s  %c", "SS", ss[pos+i]);
       else if (pos+i < L) fprintf (ofp, "%c", ss[pos+i]);
     
     fprintf (ofp, "\n\t");
     for (i = 0; i < line; i++)
       if (i == 0)        fprintf (ofp, "%15.15s  %c", sqinfo->name, seq[pos+i]);
       else if (pos+i < L) fprintf (ofp, "%c", seq[pos+i]);
     
     nlines++;
     pos += line;
     fprintf (ofp, "\n");
   }
   
   fprintf (ofp, "\n");
}

/* Function: CompareRNAStrutures()
 * 
 * ER, Tue Oct  9 15:19:38 CDT 2001 
 *
 * Purpose:  given 2 secondary structures, compare them
 *
 *
 * C coefficient: Matthews correction factor.
 *
 *   Matthews, BW (1975) 
 *   Comparison of the predicted and observed secondary structure of T4 phage Lysozyme.
 *   Biochem. Biophys. Acta, 405, 442-451.
 *
 *   Pt = true  positives = agree.
 *   Pf = false positives = pairs - agree.
 *
 *   Nt = ture  negatives = total unpaired in both sequences.
 *   Nf = false negatives = pairs_true - agree;
 *
 *
 *   C =( Pt*Nt - Pf*Nf) / sqrt[ (Nt+Nf) * (Nt+Pf) * (Pt+Nf) * (Pt+Pf) ]
 *
 *   under the approximations Nf/Nt -> 0 and Pf/Nt ->0 for N -> infty
 *                            Pt > 0 with at least Pt \sim Pf or Pt sim Nf
 *
 *   C_app = sqrt[ Pt/(Pt+Nf) * Pt/(Pt+Pf) ] (geometric mean of sensitivity and specifity)
 *
 */
void
CompareRNAStructures(FILE *ofp, int start, int L, char *ss_true, int *cc)
{
  int  *cc_true;
  int   i;
  float sen, spe;
  float c;                  /* Matthews coefficient             */
  float c_ap;               /* approximate Matthews coefficient */
  float agree = 0.;         /* agree      = Pt                  */
  float pairs = 0.;         /* pairs      = Pt + Pf             */
  float pairs_true = 0.;    /* pairs_true = Pt + Nf             */    
  float Pf, Pt;
  float Nf;
  float Nt = 0.;
  int   verbose = FALSE;
 
  KHS2ct(start+ss_true, L, TRUE, &cc_true);

  if (verbose) 
    for (i = 0; i < L; i++) printf(" %d %d %d\n", i, cc_true[i], cc[i]);
  
  for (i = 0; i < L; i++) {
    if (cc_true[i] != -1) pairs_true += 1.; else Nt += 1.;
    if (cc[i]      != -1) pairs      += 1.; else Nt += 1.;

    if (cc_true[i] != -1 && cc_true[i] == cc[i]) agree += 1.;
  }

  Pt = agree;
  Pf = pairs      - agree;
  Nf = pairs_true - agree;

  sen = 100.*Pt/(Pt+Nf);
  spe = 100.*Pt/(Pt+Pf);

  c    = 100.*( Pt*Nt - Pf*Nf) / sqrt( (Nt+Nf) * (Nt+Pf) * (Pt+Nf) * (Pt+Pf) ); /* not  very useful */
  c_ap = sqrt(sen * spe);

  fprintf(ofp, "true pairs %.1f\t found pairs %.1f\t agree %.1f (sen=%.2f spe=%.2f -- C_ap = %.2f)\n", 
	  pairs_true/2.0, pairs/2.0, agree/2.0, 
	  sen, spe, c_ap);

  if (((int)pairs%2 != 0) || ((int)pairs_true%2 != 0) || ((int)agree%2 != 0)) 
    Die("Error in CompareRNAStrutures(); odd number of paired nucleotides\n");

  if ( ((int)Nt%2 != 0) ) 
    Die("Error in CompareRNAStrutures(); odd number of total unpaired nucleotides\n");


  free(cc_true);
}






















