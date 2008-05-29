/* (modified) 
 * protowx.h
 * ANSI prototypes for all external functions in wxgraphs.c used in filltrwx.c.
 * 
 */

#include <stdio.h>
#include "cfg.h"
#include "squid.h"

/* from wxgraph.c:  the diagrams that contribute to wx
 */

extern int W1(int *s, int len, struct rnapar_2 *rnapar, int **icfg, int **vx, int j, int d);
extern void trace_W1(FILE *outf, int **vx,  int j, int d, 
		     struct tracekn_s *curr_tr, struct traceknstack_s *dolist, int traceback);

extern int W2(int *s, int len, struct rnapar_2 *rnapar, int **icfg, int **vx, int j, int d);
extern void trace_W2(FILE *outf, int **vx,  int j, int d, 
		     struct tracekn_s *curr_tr, struct traceknstack_s *dolist, int traceback);

extern int W3(int *s, int len, struct rnapar_2 *rnapar, int **icfg, int **vx, int j, int d);
extern void trace_W3(FILE *outf, int **vx,  int j, int d,
		     struct tracekn_s *curr_tr, struct traceknstack_s *dolist, int traceback);

extern int W4(int *s, int len, struct rnapar_2 *rnapar, int **icfg, int **vx, int j, int d);
extern void trace_W4(FILE *outf, int **vx,  int j, int d,
		     struct tracekn_s *curr_tr, struct traceknstack_s *dolist, int traceback);

extern int W5(int *s, int len, struct rnapar_2 *rnapar, int **icfg, int **wx, int j, int d);
extern void trace_W5(FILE *outf, int **wx, int j, int d,
		     struct tracekn_s *curr_tr, struct traceknstack_s *dolist, int traceback);

extern int W6 (int *s, int len, struct rnapar_2 *rnapar, int **icfg, int **wx, int j, int d);
extern void trace_W6(FILE *outf, int **wx,  int j, int d,
		     struct tracekn_s *curr_tr, struct traceknstack_s *dolist, int traceback);

extern int W7(int *s, int len, struct rnapar_2 *rnapar, int **icfg, int **wx, int j, int d, int mid);
extern void trace_W7(FILE *outf, int **wx,  int j, int d, int mid,
		     struct tracekn_s *curr_tr, struct traceknstack_s *dolist, int traceback);

extern int W8_1(int *s, int len, struct rnapar_2 *rnapar, int **icfg, int **vx, int j, int d, int mid);
extern int W8_2(int *s, int len, struct rnapar_2 *rnapar, int **icfg, int **vx, int j, int d, int mid);
extern int W8_3(int *s, int len, struct rnapar_2 *rnapar, int **icfg, int **vx, int j, int d, int mid);
extern int W8_4(int *s, int len, struct rnapar_2 *rnapar, int **icfg, int **vx, int j, int d, int mid);
extern int W8_5(int *s, int len, struct rnapar_2 *rnapar, int **icfg, int **vx, int j, int d, int mid);
extern int W8_6(int *s, int len, struct rnapar_2 *rnapar, int **icfg, int **vx, int j, int d, int mid);
extern int W8_7(int *s, int len, struct rnapar_2 *rnapar, int **icfg, int **vx, int j, int d, int mid);
extern int W8_8(int *s, int len, struct rnapar_2 *rnapar, int **icfg, int **vx, int j, int d, int mid);
extern void trace_W8_1(FILE *outf, int **vx,  int j, int d, int mid,
		       struct tracekn_s *curr_tr, struct traceknstack_s *dolist, int traceback);
extern void trace_W8_2(FILE *outf, int **vx,  int j, int d, int mid,
		       struct tracekn_s *curr_tr, struct traceknstack_s *dolist, int traceback);
extern void trace_W8_3(FILE *outf, int **vx,  int j, int d, int mid,
		       struct tracekn_s *curr_tr, struct traceknstack_s *dolist, int traceback);
extern void trace_W8_4(FILE *outf, int **vx,  int j, int d,  int mid,
		       struct tracekn_s *curr_tr, struct traceknstack_s *dolist, int traceback);
extern void trace_W8_5(FILE *outf, int **vx,  int j, int d, int mid,
		       struct tracekn_s *curr_tr, struct traceknstack_s *dolist, int traceback);
extern void trace_W8_6(FILE *outf, int **vx,  int j, int d, int mid,
		       struct tracekn_s *curr_tr, struct traceknstack_s *dolist, int traceback);
extern void trace_W8_7(FILE *outf, int **vx,  int j, int d, int mid,
		       struct tracekn_s *curr_tr, struct traceknstack_s *dolist, int traceback);
extern void trace_W8_8(FILE *outf, int **vx,  int j, int d,  int mid,
		       struct tracekn_s *curr_tr, struct traceknstack_s *dolist, int traceback);

extern int W9(int *s, int len, struct rnapar_2 *rnapar, int **icfg, int ****whx, int j, int d, int mid, int mid1, int mid2);
extern void trace_W9(FILE *outf, int ****whx, int j, int d, int mid, int mid1, int mid2, 
		     struct tracekn_s *curr_tr, struct traceknstack_s *dolist, int traceback);

extern int W10_1(int *s, int len, struct rnapar_2 *rnapar, int **icfg, int ****yhx, 
		 int j, int d, int mid, int mid1, int mid2);
extern int W10_2(int *s, int len, struct rnapar_2 *rnapar, int **icfg, int ****yhx,
		 int j, int d, int mid, int mid1, int mid2);
extern void trace_W10_1(FILE *outf, int ****yhx, int j, int d, int mid, int mid1, int mid2, 
			struct tracekn_s *curr_tr, struct traceknstack_s *dolist, int traceback);
extern void trace_W10_2(FILE *outf, int ****yhx, int j, int d,  int mid, int mid1, int mid2,
			struct tracekn_s *curr_tr, struct traceknstack_s *dolist, int traceback);




