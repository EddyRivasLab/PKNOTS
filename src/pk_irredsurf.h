/* (modified) 
 * proto.h
 * ANSI prototypes for all external functions.
 * 
 */

#include <stdio.h>
#include "pknots.h"

#include <easel.h>
#include <esl_sqio.h>
#include <esl_wuss.h>

extern int F1(ESL_DSQ *s, int len, struct rnapar_2 *rnapar, int **icfg, int j, int d);
extern int F2(ESL_DSQ *s, int len, struct rnapar_2 *rnapar, int **icfg, int j, int d, int d1, int d2);

