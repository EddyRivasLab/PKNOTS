/* pk_util.c
 ` */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/types.h>
#include <time.h>

#include "pknots.h"
#include "pk_util.h"

#include <easel.h>
#include <esl_sqio.h>
#include <esl_wuss.h>


/* Function: pk_fatal()
 *
 * Purpose:  Print an error message and die. The arguments
 *           are formatted exactly like arguments to printf().
 *
 * Return:   None. Exits the program.
 */
/* VARARGS0 */
void
pk_fatal(char *format, ...)
{
  va_list argp;

  va_start(argp, format);
  vfprintf(stderr, format, argp);
  va_end(argp);
  fprintf(stderr, "\n");
  fflush(stderr);
  exit(1);
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
IntizeSequence(char *seq, int len, int *iseq)
{
  int  i;

 for (i = 0; i < len; i++)
    switch (seq[i]) 
      {
      case 'A': iseq[i] = 0; break;
      case 'C': iseq[i] = 1; break;
      case 'G': iseq[i] = 2; break;
      case 'U': iseq[i] = 3; break;
      default:  iseq[i] = 4; break;
      }
}

