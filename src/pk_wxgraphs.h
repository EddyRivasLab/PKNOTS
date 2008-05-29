/* (modified) 
 * protowx.h
 * ANSI prototypes for all external functions in wxgraphs.c used in filltrwx.c.
 * 
 */

#include <stdio.h>

/* from wxgraph.c:  the diagrams that contribute to wx
 */

extern int W1(ESL_DSQ *s, int len, struct rnapar_2 *rnapar, int **icfg, int **vx, int j, int d);
extern void trace_W1(FILE *outf, int **vx,  int j, int d, 
		     struct tracekn_s *curr_tr, struct traceknstack_s *dolist, int traceback);

extern int W2(ESL_DSQ *s, int len, struct rnapar_2 *rnapar, int **icfg, int **vx, int j, int d);
extern void trace_W2(FILE *outf, int **vx,  int j, int d, 
		     struct tracekn_s *curr_tr, struct traceknstack_s *dolist, int traceback);

extern int W3(ESL_DSQ *s, int len, struct rnapar_2 *rnapar, int **icfg, int **vx, int j, int d);
extern void trace_W3(FILE *outf, int **vx,  int j, int d,
		     struct tracekn_s *curr_tr, struct traceknstack_s *dolist, int traceback);

extern int W4(ESL_DSQ *s, int len, struct rnapar_2 *rnapar, int **icfg, int **vx, int j, int d);
extern void trace_W4(FILE *outf, int **vx,  int j, int d,
		     struct tracekn_s *curr_tr, struct traceknstack_s *dolist, int traceback);

extern int W5(ESL_DSQ *s, int len, struct rnapar_2 *rnapar, int **icfg, int **wx, int j, int d);
extern void trace_W5(FILE *outf, int **wx, int j, int d,
		     struct tracekn_s *curr_tr, struct traceknstack_s *dolist, int traceback);

extern int W6 (ESL_DSQ *s, int len, struct rnapar_2 *rnapar, int **icfg, int **wx, int j, int d);
extern void trace_W6(FILE *outf, int **wx,  int j, int d,
		     struct tracekn_s *curr_tr, struct traceknstack_s *dolist, int traceback);

extern int W7(ESL_DSQ *s, int len, struct rnapar_2 *rnapar, int **icfg, int **wx, int j, int d, int mid);
extern void trace_W7(FILE *outf, int **wx,  int j, int d, int mid,
		     struct tracekn_s *curr_tr, struct traceknstack_s *dolist, int traceback);

extern int W8_1(ESL_DSQ *s, int len, struct rnapar_2 *rnapar, int **icfg, int **vx, int j, int d, int mid);
extern int W8_2(ESL_DSQ *s, int len, struct rnapar_2 *rnapar, int **icfg, int **vx, int j, int d, int mid);
extern int W8_3(ESL_DSQ *s, int len, struct rnapar_2 *rnapar, int **icfg, int **vx, int j, int d, int mid);
extern int W8_4(ESL_DSQ *s, int len, struct rnapar_2 *rnapar, int **icfg, int **vx, int j, int d, int mid);
extern int W8_5(ESL_DSQ *s, int len, struct rnapar_2 *rnapar, int **icfg, int **vx, int j, int d, int mid);
extern int W8_6(ESL_DSQ *s, int len, struct rnapar_2 *rnapar, int **icfg, int **vx, int j, int d, int mid);
extern int W8_7(ESL_DSQ *s, int len, struct rnapar_2 *rnapar, int **icfg, int **vx, int j, int d, int mid);
extern int W8_8(ESL_DSQ *s, int len, struct rnapar_2 *rnapar, int **icfg, int **vx, int j, int d, int mid);
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

extern int W9(ESL_DSQ *s, int len, struct rnapar_2 *rnapar, int **icfg, int ****whx, int j, int d, int mid, int mid1, int mid2);
extern void trace_W9(FILE *outf, int ****whx, int j, int d, int mid, int mid1, int mid2, 
		     struct tracekn_s *curr_tr, struct traceknstack_s *dolist, int traceback);

extern int W10_1(ESL_DSQ *s, int len, struct rnapar_2 *rnapar, int **icfg, int ****yhx, 
		 int j, int d, int mid, int mid1, int mid2);
extern int W10_2(ESL_DSQ *s, int len, struct rnapar_2 *rnapar, int **icfg, int ****yhx,
		 int j, int d, int mid, int mid1, int mid2);
extern void trace_W10_1(FILE *outf, int ****yhx, int j, int d, int mid, int mid1, int mid2, 
			struct tracekn_s *curr_tr, struct traceknstack_s *dolist, int traceback);
extern void trace_W10_2(FILE *outf, int ****yhx, int j, int d,  int mid, int mid1, int mid2,
			struct tracekn_s *curr_tr, struct traceknstack_s *dolist, int traceback);




