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

extern void FillMtx_nested(ESL_DSQ *s, int len, struct rnapar_2 *rnapar,
			   int **icfg, int **wx, int **wbx, int **vx, int *vp, 
			   int allow_coaxials);
extern void FillMtx(ESL_DSQ *s, int len, struct rnapar_2 *rnapar, int **icfg, 
		    int **wx, int **wbx, int **vx, int *vp,
		    int ****whx, int ****vhx, int ****zhx, int ****yhx, int allow_coaxials, int approx);








