/* misc.c
 * Functions with no home.
 * SRE, Fri May 27 15:11:00 1994
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/types.h>
#include <time.h>

#include "cfg.h"
#include "proto.h"
#include "squid.h"

#ifdef MEMDEBUG
#include "dbmalloc.h"
#endif

/* Function: StripDegeneracy()
 * 
 * Purpose:  Convert degenerate nucleotides into a random choice
 *           of ACGU. String is guaranteed to contain only
 *           ACGU when it comes out. (Gaps are removed.)
 *           
 * Args:     seq    - sequence to strip (null-terminated)
 *           
 * Return:   (void)
 */
void
StripDegeneracy(char *seq)
{
  char *wp; /* write pointer (where we're writing seq) */
  char *rp; /* read pointer (where we're reading seq)  */

  for (wp = rp = seq; *rp != '\0'; rp++)
    {
      if (isgap(*rp)) continue;
      if (strchr("ACGU", *rp)) *wp++ = *rp;
      else 
	{
	  /* then it's a degenerate symbol.
	   * According to alphabet, choose a single symbol to represent it.
	   * note the too-clever scheme for random choice: "ABC"[random() % 3]
	   */
	  switch (*rp) {
	  case 'B': *wp++ = "CGU" [random() % 3]; break;
	  case 'D': *wp++ = "AGU" [random() % 3]; break;
	  case 'H': *wp++ = "ACU" [random() % 3]; break;
	  case 'K': *wp++ = "GU"  [random() % 2]; break;
	  case 'M': *wp++ = "AC"  [random() % 2]; break;
	  case 'N': *wp++ = "ACGU"[random() % 4]; break;
	  case 'R': *wp++ = "AG"  [random() % 2]; break;
	  case 'S': *wp++ = "CG"  [random() % 2]; break;
	  case 'T': *wp++ = 'U';                  break;
	  case 'V': *wp++ = "ACG" [random() % 3]; break;
	  case 'W': *wp++ = "AU"  [random() % 2]; break;
	  case 'X': *wp++ = "ACGU"[random() % 4]; break;
	  case 'Y': *wp++ = "CU"  [random() % 2]; break;
	  default: Die("unrecognized character %c in sequence\n", *rp);
	  }
	}
    }
}

/* Function: IntizeSequence()
 * 
 * Purpose:  Convert a sequence of A,C,G,U into a sequence
 *           of integer indices 0,1,2,3
 *           
 * Args:     seq      - sequence (0..N-1) only A,C,G,U allowed
 *           len      - length of seq
 *           ret_iseq - RETURN: integer-ized sequence
 *           
 * Return:   (void)
 *           ret_iseq is alloc'ed here, must be free'd by caller.
 */
void
IntizeSequence(char *seq, int len, int **ret_iseq)
{
  int  i;
  int *iseq;

  if ((iseq = (int *) malloc (len * sizeof(int))) == NULL)
    Die("malloc failed");
  for (i = 0; i < len; i++)
    switch (seq[i]) 
      {
      case 'A': iseq[i] = 0; break;
      case 'C': iseq[i] = 1; break;
      case 'G': iseq[i] = 2; break;
      case 'U': iseq[i] = 3; break;
      default:  iseq[i] = 4; break;
      }
  *ret_iseq = iseq;
}

char *
DupSeq(char *seq, int j, int d)
{
  char *new;
  int   mid;

  if ((new = (char *) malloc ((d+1) * sizeof(char))) == NULL)
    Die("malloc failed");
  for (mid = 0; mid <= d; mid++)
    new[mid] = seq[j-d+mid];
  return new;
}

/* Function: ShuffleSequence()
 * 
 * Purpose:  Convert a sequence of A,C,G,U into a shuffled sequence
 *           
 *           
 * Args:     seq      - sequence (0..N-1) only A,C,G,U allowed
 *           len      - length of seq
 *           
 * Return:   (void)
 */
void
ShuffleSequence(char *seq, int len, int endpos, int verbose)
{
  int    d, s;
  int    i;
  int    intpos, pos, pos_s;
  int    count = 0;  /* number of times we shuffle */
  int    seed;
  char   seq_s;
  
  seed = (int) time ((time_t *) NULL);

  sre_srandom(seed); /* reinit sre_random each time you shuffle a sequence */

  intpos = endpos - len + 1;

  if (verbose) {
    printf("\n\n");
    for (d = 0; d < len; d++) printf("%c ", seq[intpos+d]);
    printf("\n\n");
  }
  
  while (count < 5000) {
    for (d = 0; d < len; d++) {
      s = (int)(sre_random()*(len-1));

      pos   = intpos + d;
      pos_s = intpos + s;

      seq_s = seq[pos_s];
      seq[pos_s] = seq[pos];
      seq[pos] = seq_s;
    }
    
    count++;
  }

  if (verbose) {
    printf("\n\n");
    for (d = 0; d < len; d++) printf("%c ", seq[intpos+d]);
    printf("\n\n"); 
  }
}







