/* (modified) 
 * protowbx.h
 * ANSI prototypes for all external functions in wbxgraphs.c used in filltrwbx.c.
 * 
 */

#include <stdio.h>

/* from wbxgraph.c:  the diagrams that contribute to wbx
 */

extern int WB1(ESL_DSQ *s, int len, struct rnapar_2 *rnapar, int **icfg, int **vx, int j, int d);
extern void trace_WB1(FILE *outf, int **vx,  int j, int d, 
		      struct tracekn_s *curr_tr, struct traceknstack_s *dolist, int traceback);

extern int WB2(ESL_DSQ *s, int len, struct rnapar_2 *rnapar, int **icfg, int **vx, int j, int d);
extern void trace_WB2(FILE *outf, int **vx,  int j, int d, 
		      struct tracekn_s *curr_tr, struct traceknstack_s *dolist, int traceback);

extern int WB3(ESL_DSQ *s, int len, struct rnapar_2 *rnapar, int **icfg, int **vx, int j, int d);
extern void trace_WB3(FILE *outf, int **vx,  int j, int d,
		       struct tracekn_s *curr_tr, struct traceknstack_s *dolist, int traceback);

extern int WB4(ESL_DSQ *s, int len, struct rnapar_2 *rnapar, int **icfg, int **vx, int j, int d);
extern void trace_WB4(FILE *outf, int **vx,  int j, int d,
		       struct tracekn_s *curr_tr, struct traceknstack_s *dolist, int traceback);

extern int WB5(ESL_DSQ *s, int len, struct rnapar_2 *rnapar, int **icfg, int **wbx, int j, int d);
extern void trace_WB5(FILE *outf, int **wbx, int j, int d,
		       struct tracekn_s *curr_tr, struct traceknstack_s *dolist, int traceback);

extern int WB6 (ESL_DSQ *s, int len, struct rnapar_2 *rnapar, int **icfg, int **wbx, int j, int d);
extern void trace_WB6(FILE *outf, int **wbx,  int j, int d,
		      struct tracekn_s *curr_tr, struct traceknstack_s *dolist, int traceback);

extern int WB7(ESL_DSQ *s, int len, struct rnapar_2 *rnapar, int **icfg, int **wbx, int j, int d, int mid);
extern void trace_WB7(FILE *outf, int **wbx,  int j, int d, int mid,
		      struct tracekn_s *curr_tr, struct traceknstack_s *dolist, int traceback);

extern int WB8_1(ESL_DSQ *s, int len, struct rnapar_2 *rnapar, int **icfg, int **vx, int j, int d, int mid);
extern int WB8_2(ESL_DSQ *s, int len, struct rnapar_2 *rnapar, int **icfg, int **vx, int j, int d, int mid);
extern int WB8_3(ESL_DSQ *s, int len, struct rnapar_2 *rnapar, int **icfg, int **vx, int j, int d, int mid);
extern int WB8_4(ESL_DSQ *s, int len, struct rnapar_2 *rnapar, int **icfg, int **vx, int j, int d, int mid);
extern int WB8_5(ESL_DSQ *s, int len, struct rnapar_2 *rnapar, int **icfg, int **vx, int j, int d, int mid);
extern int WB8_6(ESL_DSQ *s, int len, struct rnapar_2 *rnapar, int **icfg, int **vx, int j, int d, int mid);
extern int WB8_7(ESL_DSQ *s, int len, struct rnapar_2 *rnapar, int **icfg, int **vx, int j, int d, int mid);
extern int WB8_8(ESL_DSQ *s, int len, struct rnapar_2 *rnapar, int **icfg, int **vx, int j, int d, int mid);
extern void trace_WB8_1(FILE *outf, int **vx,  int j, int d, int mid,
			struct tracekn_s *curr_tr, struct traceknstack_s *dolist, int traceback);
extern void trace_WB8_2(FILE *outf, int **vx,  int j, int d, int mid,
			struct tracekn_s *curr_tr, struct traceknstack_s *dolist, int traceback);
extern void trace_WB8_3(FILE *outf, int **vx,  int j, int d, int mid,
			struct tracekn_s *curr_tr, struct traceknstack_s *dolist, int traceback);
extern void trace_WB8_4(FILE *outf, int **vx,  int j, int d,  int mid,
			struct tracekn_s *curr_tr, struct traceknstack_s *dolist, int traceback);
extern void trace_WB8_5(FILE *outf, int **vx,  int j, int d, int mid,
			struct tracekn_s *curr_tr, struct traceknstack_s *dolist, int traceback);
extern void trace_WB8_6(FILE *outf, int **vx,  int j, int d, int mid,
			struct tracekn_s *curr_tr, struct traceknstack_s *dolist, int traceback);
extern void trace_WB8_7(FILE *outf, int **vx,  int j, int d, int mid,
			struct tracekn_s *curr_tr, struct traceknstack_s *dolist, int traceback);
extern void trace_WB8_8(FILE *outf, int **vx,  int j, int d,  int mid,
			struct tracekn_s *curr_tr, struct traceknstack_s *dolist, int traceback);

extern int WB9(ESL_DSQ *s, int len, struct rnapar_2 *rnapar, int **icfg, int ****whx, int j, int d, int mid, int mid1, int mid2);
extern void trace_WB9(FILE *outf, int ****whx, int j, int d, int mid, int mid1, int mid2, 
		      struct tracekn_s *curr_tr, struct traceknstack_s *dolist, int traceback);

extern int WB10_1(ESL_DSQ *s, int len, struct rnapar_2 *rnapar, int **icfg, int ****yhx, 
		  int j, int d, int mid, int mid1, int mid2);
extern int WB10_2(ESL_DSQ *s, int len, struct rnapar_2 *rnapar, int **icfg, int ****yhx,
		  int j, int d, int mid, int mid1, int mid2);
extern void trace_WB10_1(FILE *outf, int ****yhx, int j, int d, int mid, int mid1, int mid2, 
			 struct tracekn_s *curr_tr, struct traceknstack_s *dolist, int traceback);
extern void trace_WB10_2(FILE *outf, int ****yhx, int j, int d,  int mid, int mid1, int mid2,
			 struct tracekn_s *curr_tr, struct traceknstack_s *dolist, int traceback);


