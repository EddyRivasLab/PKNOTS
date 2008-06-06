/* pk_cyk.h 
 * 
 */

#include <stdio.h>
#include "pknots.h"

#include <easel.h>
#include <esl_sq.h>
#include <esl_sqio.h>
#include <esl_wuss.h>

extern void Alloc_Mtx(int len, int ***ret_wx, int ***ret_wbx, int ***ret_vx, int **ret_vp);
extern void Pattern_Mtx(int len, int **wx, int **wbx, int **vx, int *vp);
extern void Free_Mtx(int len, int **wx, int **wbx, int **vx, int *vp);
extern void Alloc_Mgp(int len, int *****ret_whx, int *****ret_vhx,
		      int *****ret_zhx, int *****ret_yhx);
extern void Pattern_Mgp(int len, int ****whx, int ****vhx,
			int ****zhx, int ****yhx);
extern void Free_Mgp(int len, int ****whx, int ****vhx,
		     int ****zhx, int ****yhx);
extern void Print_2DMtx(FILE *fp, int **mtx, int len, struct rnapar_2 *rnapar);
extern void Print_4DMtx(FILE *fp, int ****mtx, int len, struct rnapar_2 *rnapar);
extern int  StructurePredictkn_2IS(FILE *outf, ESL_SQ *sq, int len, struct rnapar_2 *rnapar,
				   int **icfg, int verbose, int traceback, struct tracekn_s **ret_tr, 
				   int *ret_score, int allow_coaxials, int allow_pseudoknots, 
				   int approx);




