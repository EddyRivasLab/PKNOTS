/* (modified) 
 * proto.h
 * ANSI prototypes for all external functions.
 * 
 */

#include <stdio.h>
#include "cfg.h"

#include <easel.h>
#include <esl_sqio.h>
#include <esl_wuss.h>

extern int F1(int *s, int len, struct rnapar_2 *rnapar, int **icfg, int j, int d);
extern int F2(int *s, int len, struct rnapar_2 *rnapar, int **icfg, int j, int d, int d1, int d2);

