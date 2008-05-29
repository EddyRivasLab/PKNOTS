/* pk_util.c
 ` */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/types.h>
#include <time.h>

#include "cfg.h"
#include "proto.h"

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
