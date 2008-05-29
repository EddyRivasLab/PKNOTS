#include <stdio.h>
#include "cfg.h"

#include <easel.h>
#include <esl_sqio.h>
#include <esl_wuss.h>

extern void TraceMtx_nested(FILE *outf, int *s, int len, struct rnapar_2 *rnapar, int **icfg, 
			    int **wx, int **wbx, int **vx, 
			    int j, int d, struct tracekn_s **ret_trace, int traceback);
extern void TraceMtx(FILE *outf, int *s, int len, struct rnapar_2 *rnapar, int **icfg,  
		     int **wx, int **wbx, int **vx, 
		     int ****whx, int ****vhx, int ****zhx, int ****yhx,
		     int j, int d, int d1, int d2, int approx, 
		     struct tracekn_s **ret_trace, int traceback);
