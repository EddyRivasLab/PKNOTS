/* (modified) 
 * protovhx.h
 * ANSI prototypes for all external functions in vhxgraphs.c used in filltrvhx.c.
 * 
 */

#include <stdio.h>
#include "cfg.h"
#include "squid.h"

/* from vhxgraph.c:  the diagrams that contribute to vhx
 */

extern int VH1(int *s, int len, struct rnapar_2 *rnapar, int **icfg, int j, int d, int d1, int d2);
extern void trace_VH1(FILE *outf, int ****vhx, int j, int d, int d1, int d2, 
		      struct tracekn_s *curr_tr, struct traceknstack_s *dolist, int traceback);

extern int VH2(int *s, int len, struct rnapar_2 *rnapar, int **icfg, int ****vhx, int j, int d, 
	       int d1, int d2, int mid1, int mid2);
extern void trace_VH2(FILE *outf, int ****vhx, int j, int d, int d1, int d2, int mid1, int mid2, 
		      struct tracekn_s *curr_tr, struct traceknstack_s *dolist, int traceback);

extern int VH3_1(int *s, int len, struct rnapar_2 *rnapar, int **icfg, int ****whx, int j, int d, int d1, int d2);
extern void trace_VH3_1(FILE *outf, int ****whx, int j, int d, int d1, int d2,
			struct tracekn_s *curr_tr, struct traceknstack_s *dolist, int traceback);

extern int VH3_2(int *s, int len, struct rnapar_2 *rnapar, int **icfg, int ****whx, int j, int d, int d1, int d2);
extern void trace_VH3_2(FILE *outf, int ****whx, int j, int d, int d1, int d2,
			struct tracekn_s *curr_tr, struct traceknstack_s *dolist, int traceback);

extern int VH3_3(int *s, int len, struct rnapar_2 *rnapar, int **icfg, int ****whx, int j, int d, int d1, int d2);
extern void trace_VH3_3(FILE *outf, int ****whx, int j, int d, int d1, int d2,
			struct tracekn_s *curr_tr, struct traceknstack_s *dolist, int traceback);

extern int VH3_4(int *s, int len, struct rnapar_2 *rnapar, int **icfg, int ****whx, int j, int d, int d1, int d2);
extern void trace_VH3_4(FILE *outf, int ****whx, int j, int d, int d1, int d2,
			struct tracekn_s *curr_tr, struct traceknstack_s *dolist, int traceback);

extern int VH3_5(int *s, int len, struct rnapar_2 *rnapar, int **icfg, int ****whx, int j, int d, int d1, int d2);
extern void trace_VH3_5(FILE *outf, int ****whx, int j, int d, int d1, int d2,
			struct tracekn_s *curr_tr, struct traceknstack_s *dolist, int traceback);

extern int VH3_6(int *s, int len, struct rnapar_2 *rnapar, int **icfg, int ****whx, int j, int d, int d1, int d2);
extern void trace_VH3_6(FILE *outf, int ****whx, int j, int d, int d1, int d2,
			struct tracekn_s *curr_tr, struct traceknstack_s *dolist, int traceback);

extern int VH3_7(int *s, int len, struct rnapar_2 *rnapar, int **icfg, int ****whx, int j, int d, int d1, int d2);
extern void trace_VH3_7(FILE *outf, int ****whx, int j, int d, int d1, int d2,
			struct tracekn_s *curr_tr, struct traceknstack_s *dolist, int traceback);

extern int VH3_8(int *s, int len, struct rnapar_2 *rnapar, int **icfg, int ****whx, int j, int d, int d1, int d2);
extern void trace_VH3_8(FILE *outf, int ****whx, int j, int d, int d1, int d2,
			struct tracekn_s *curr_tr, struct traceknstack_s *dolist, int traceback);

extern int VH3_9(int *s, int len, struct rnapar_2 *rnapar, int **icfg, int ****whx, int j, int d, int d1, int d2);
extern void trace_VH3_9(FILE *outf, int ****whx, int j, int d, int d1, int d2,
			struct tracekn_s *curr_tr, struct traceknstack_s *dolist, int traceback);

extern int VH3_10(int *s, int len, struct rnapar_2 *rnapar, int **icfg, int ****whx, int j, int d, int d1, int d2);
extern void trace_VH3_10(FILE *outf, int ****whx, int j, int d, int d1, int d2,
			 struct tracekn_s *curr_tr, struct traceknstack_s *dolist, int traceback);

extern int VH3_11(int *s, int len, struct rnapar_2 *rnapar, int **icfg, int ****whx, int j, int d, int d1, int d2);
extern void trace_VH3_11(FILE *outf, int ****whx, int j, int d, int d1, int d2,
			 struct tracekn_s *curr_tr, struct traceknstack_s *dolist, int traceback);

extern int VH3_12(int *s, int len, struct rnapar_2 *rnapar, int **icfg, int ****whx, int j, int d, int d1, int d2);
extern void trace_VH3_12(FILE *outf, int ****whx, int j, int d, int d1, int d2,
			 struct tracekn_s *curr_tr, struct traceknstack_s *dolist, int traceback);

extern int VH3_13(int *s, int len, struct rnapar_2 *rnapar, int **icfg, int ****whx, int j, int d, int d1, int d2);
extern void trace_VH3_13(FILE *outf, int ****whx, int j, int d, int d1, int d2,
			 struct tracekn_s *curr_tr, struct traceknstack_s *dolist, int traceback);

extern int VH3_14(int *s, int len, struct rnapar_2 *rnapar, int **icfg, int ****whx, int j, int d, int d1, int d2);
extern void trace_VH3_14(FILE *outf, int ****whx, int j, int d, int d1, int d2,
			 struct tracekn_s *curr_tr, struct traceknstack_s *dolist, int traceback);

extern int VH3_15(int *s, int len, struct rnapar_2 *rnapar, int **icfg, int ****whx, int j, int d, int d1, int d2);
extern void trace_VH3_15(FILE *outf, int ****whx, int j, int d, int d1, int d2,
			 struct tracekn_s *curr_tr, struct traceknstack_s *dolist, int traceback);

extern int VH3_16(int *s, int len, struct rnapar_2 *rnapar, int **icfg, int ****whx, int j, int d, int d1, int d2);
extern void trace_VH3_16(FILE *outf, int ****whx, int j, int d, int d1, int d2,
			 struct tracekn_s *curr_tr, struct traceknstack_s *dolist, int traceback);










