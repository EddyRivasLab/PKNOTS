/************************************************************
 * QRNA - Comparative analysis of biological sequences 
 *        with pair hidden Markov models and pair stochastic context-free grammars
 * Copyright (C) 2000-Howard Hughes Medical Institute
 * All Rights Reserved
 * 
 *     This source code is distributed under the terms of the
 *     GNU General Public License. See the files COPYING and LICENSE
 *     for details.
 ************************************************************/

/* rnaio.c
 *
 * E. Rivas [St. Louis]
 * 
 * 9 april 1999.
 *
 *  I/O for RNA secondary structure
 * 
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

/* Function: KHS2ct()
 * 
 * Purpose:  Convert a secondary structure string to an array of integers
 *           representing what position each position is base-paired 
 *           to (0..len-1), or -1 if none. This is off-by-one from a
 *           Zuker .ct file representation.
 *           
 *           The .ct representation can accomodate pseudoknots but the 
 *           secondary structure string cannot easily; the string contains
 *           "Aa", "Bb", etc. pairs as a limited representation of
 *           pseudoknots. The string contains "><" for base pairs.
 *           Other symbols are ignored. If allow_pseudoknots is FALSE,
 *           the pseudoknot symbols will be ignored and these positions
 *           will be treated as single stranded.
 *           
 * Return:   ret_ct is allocated here and must be free'd by caller.
 *           Returns 1 on success, 0 if ss is somehow inconsistent.
 */
int 
KHS2ct(char *ss, int len, int allow_pseudoknots, int **ret_ct)
{
  struct intstack_s *dolist[27];
  int *ct;
  int  i;
  int  pos, pair;
  int  status = 1;              /* success or failure return status */
  int  verbose = FALSE;

  if (verbose) {
    for (pos = 0; ss[pos] != '\0'; pos++) printf("%c ", ss[pos]);
    printf("\n");
  }
    
  for (i = 0; i < 27; i++)
    dolist[i] = InitIntStack();

  if ((ct = (int *) malloc (len * sizeof(int))) == NULL)
    Die("malloc failed");
  for (pos = 0; pos < len; pos++)
    ct[pos] = -1;

  for (pos = 0; ss[pos] != '\0'; pos++)
    {
      if (ss[pos] > 127) status = 0; /* bulletproof against SGI buggy ctype.h */

      else if (ss[pos] == '>')  /* left side of a pair: push onto stack 0 */
        PushIntStack(dolist[0], pos);
      else if (ss[pos] == '<')  /* right side of a pair; resolve pair */
        {
          if (! PopIntStack(dolist[0], &pair))
            { status = 0; }
          else
            {
              ct[pos]  = pair;
              ct[pair] = pos;
            }
        }
                                /* same stuff for pseudoknots */
      else if (allow_pseudoknots && isupper((int) ss[pos]))
        PushIntStack(dolist[ss[pos] - 'A' + 1], pos);
      else if (allow_pseudoknots && islower((int) ss[pos]))
        {
          if (! PopIntStack(dolist[ss[pos] - 'a' + 1], &pair))
            { status = 0; }
          else
            {
              ct[pos]  = pair;
              ct[pair] = pos;
            }
        }
      else if (allow_pseudoknots && !isgap(ss[pos])) status = 0; /* bad character */
    }

  for (i = 0; i < 27; i++)
    if ( FreeIntStack(dolist[i]) > 0)
      status = 0;

  *ret_ct = ct;
  return status;
}



/* Function: VerifyKHS()
 * 
 * Purpose:  Examine a possibly bad structure string, and print out diagnostics 
 *           about it if wordy is TRUE.
 *
 * Return:   1 if string is OK, 0 if string is bad.
 */
int
VerifyKHS(char *name, char *ss, int wordy)
{
  int symcount[27];             /* 0 is normal pairs. 1-26 for pseudoknots */
  int i;
  int pos;
  int status = 1;

  for (i = 0; i < 27; i++)
    symcount[i] = 0;

  for (pos = 0; ss[pos] != '\0'; pos++)
    {
      if (!sre_isascii(ss[pos]))        /* SGI ctype.h non-ANSI compliant */
        {
          status = 0;
          if (wordy)
            Warn("VerifyKHS: Sequence %s no good. structure has garbagesymbol (val %d) at position %d", 
		 name, (int)ss[pos], pos);
        }
      else if (ss[pos] == '>')
        symcount[0] ++;
      else if (ss[pos] == '<')
        symcount[0] --;
      else if (isupper((int) ss[pos]))
        symcount[ss[pos] - 'A' + 1] ++; /* pseudoknot-left  */
      else if (islower((int) ss[pos])) 
        symcount[ss[pos] - 'a' + 1] --; /* pseudoknot-right */
      else if (!isgap(ss[pos]))
        {
          status = 0;
          if (wordy)
            Warn("VerifyKHS: Sequence %s no good, structure has invalid symbol %c at position %d", 
		 name, (int)ss[pos], pos);
        }
    }
      
  if (symcount[0] != 0)
    {
      status = 0;
      if (wordy)
        Warn("VerifyKHS: Sequence %s no good, structure has extra paired bases: %d on %s",
             name, abs(symcount[0]), 
             (symcount[0] > 0) ? "left" : "right");
    }

  for (i = 1; i <= 26; i++)
    if (symcount[i] != 0)
      {
        status = 0;
        if (wordy)
          Warn("VerifyKHS: Sequence %s no good, structure has extra paired bases for pseudoknot %c: %d on %s",
               name, (char) (i + 'A' - 1),
               abs(symcount[i]), 
               (symcount[i] > 0) ? "left" : "right");
      }
  return status;
}
