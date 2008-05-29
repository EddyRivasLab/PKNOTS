/* (modified) 
 * protovx.h
 * ANSI prototypes for all external functions in vxgraphs.c used in filltrvx.c.
 * 
 */

#include <stdio.h>

/* from vxgraph.c:  the diagrams that contribute to vx
 */

extern int V1(ESL_DSQ *s, int len, struct rnapar_2 *rnapar, int **icfg, int j, int d);
extern void trace_V1(FILE *outf, int **vx,  int j, int d, 
		     struct tracekn_s *curr_tr, struct traceknstack_s *dolist, int traceback);

extern int V2(ESL_DSQ *s, int len, struct rnapar_2 *rnapar, int **icfg, int **vx, int j, int d, int mid1, int mid2);
extern void trace_V2(FILE *outf, int **vx, ESL_DSQ *s, int len, struct rnapar_2 *rnapar, int **icfg, 
		     int j, int d, int mid1, int mid2, 
		     struct tracekn_s *curr_tr, struct traceknstack_s *dolist, int traceback);

extern int V3_1(ESL_DSQ *s, int len, struct rnapar_2 *rnapar, int **icfg, int **wbx, int j, int d, int mid);
extern int V3_2(ESL_DSQ *s, int len, struct rnapar_2 *rnapar, int **icfg, int **wbx, int j, int d, int mid);
extern int V3_3(ESL_DSQ *s, int len, struct rnapar_2 *rnapar, int **icfg, int **wbx, int j, int d, int mid);
extern int V3_4(ESL_DSQ *s, int len, struct rnapar_2 *rnapar, int **icfg, int **wbx, int j, int d, int mid);
extern void trace_V3_1(FILE *outf, int **wbx, int j, int d, int mid, 
		       struct tracekn_s *curr_tr, struct traceknstack_s *dolist, int traceback);
extern void trace_V3_2(FILE *outf, int **wbx, int j, int d, int mid,
		       struct tracekn_s *curr_tr, struct traceknstack_s *dolist, int traceback);
extern void trace_V3_3(FILE *outf, int **wbx, int j, int d, int mid,
		       struct tracekn_s *curr_tr, struct traceknstack_s *dolist, int traceback);
extern void trace_V3_4(FILE *outf, int **wbx, int j, int d, int mid, 
		       struct tracekn_s *curr_tr, struct traceknstack_s *dolist, int traceback);

extern int V4_1(ESL_DSQ *s, int len, struct rnapar_2 *rnapar, int **icfg, int **wbx, int **vx, int j, int d, int mid);
extern int V4_2(ESL_DSQ *s, int len, struct rnapar_2 *rnapar, int **icfg, int **wbx, int **vx, int j, int d, int mid);
extern int V4_3(ESL_DSQ *s, int len, struct rnapar_2 *rnapar, int **icfg, int **wbx, int **vx, int j, int d, int mid);
extern int V4_4(ESL_DSQ *s, int len, struct rnapar_2 *rnapar, int **icfg, int **wbx, int **vx, int j, int d, int mid);
extern int V4_5(ESL_DSQ *s, int len, struct rnapar_2 *rnapar, int **icfg, int **wbx, int **vx, int j, int d, int mid);
extern int V4_6(ESL_DSQ *s, int len, struct rnapar_2 *rnapar, int **icfg, int **wbx, int **vx, int j, int d, int mid);
extern int V4_7(ESL_DSQ *s, int len, struct rnapar_2 *rnapar, int **icfg, int **wbx, int **vx, int j, int d, int mid);
extern int V4_8(ESL_DSQ *s, int len, struct rnapar_2 *rnapar, int **icfg, int **wbx, int **vx, int j, int d, int mid);
extern void trace_V4_1(FILE *outf, int **wbx, int **vx, int j, int d, int mid,
		       struct tracekn_s *curr_tr, struct traceknstack_s *dolist, int traceback);
extern void trace_V4_2(FILE *outf, int **wbx, int **vx, int j, int d, int mid, 
		       struct tracekn_s *curr_tr, struct traceknstack_s *dolist, int traceback);
extern void trace_V4_3(FILE *outf, int **wbx, int **vx, int j, int d, int mid,
		       struct tracekn_s *curr_tr, struct traceknstack_s *dolist, int traceback);
extern void trace_V4_4(FILE *outf, int **wbx, int **vx, int j, int d, int mid, 
		       struct tracekn_s *curr_tr, struct traceknstack_s *dolist, int traceback);
extern void trace_V4_5(FILE *outf, int **wbx, int **vx, int j, int d, int mid,
		       struct tracekn_s *curr_tr, struct traceknstack_s *dolist, int traceback);
extern void trace_V4_6(FILE *outf, int **wbx, int **vx, int j, int d, int mid, 
		       struct tracekn_s *curr_tr, struct traceknstack_s *dolist, int traceback);
extern void trace_V4_7(FILE *outf, int **wbx, int **vx, int j, int d, int mid,
		       struct tracekn_s *curr_tr, struct traceknstack_s *dolist, int traceback);
extern void trace_V4_8(FILE *outf, int **wbx, int **vx, int j, int d, int mid, 
		       struct tracekn_s *curr_tr, struct traceknstack_s *dolist, int traceback);

extern int V5_1(ESL_DSQ *s, int len, struct rnapar_2 *rnapar, int **icfg, int **wbx, int **vx, int j, int d, int mid);
extern int V5_2(ESL_DSQ *s, int len, struct rnapar_2 *rnapar, int **icfg, int **wbx, int **vx, int j, int d, int mid);
extern int V5_3(ESL_DSQ *s, int len, struct rnapar_2 *rnapar, int **icfg, int **wbx, int **vx, int j, int d, int mid);
extern int V5_4(ESL_DSQ *s, int len, struct rnapar_2 *rnapar, int **icfg, int **wbx, int **vx, int j, int d, int mid);
extern int V5_5(ESL_DSQ *s, int len, struct rnapar_2 *rnapar, int **icfg, int **wbx, int **vx, int j, int d, int mid);
extern int V5_6(ESL_DSQ *s, int len, struct rnapar_2 *rnapar, int **icfg, int **wbx, int **vx, int j, int d, int mid);
extern int V5_7(ESL_DSQ *s, int len, struct rnapar_2 *rnapar, int **icfg, int **wbx, int **vx, int j, int d, int mid);
extern int V5_8(ESL_DSQ *s, int len, struct rnapar_2 *rnapar, int **icfg, int **wbx, int **vx, int j, int d, int mid);
extern void trace_V5_1(FILE *outf, int **wbx, int **vx, int j, int d, int mid,
		       struct tracekn_s *curr_tr, struct traceknstack_s *dolist, int traceback);
extern void trace_V5_2(FILE *outf, int **wbx, int **vx, int j, int d, int mid, 
		       struct tracekn_s *curr_tr, struct traceknstack_s *dolist, int traceback);
extern void trace_V5_3(FILE *outf, int **wbx, int **vx, int j, int d, int mid, 
		       struct tracekn_s *curr_tr, struct traceknstack_s *dolist, int traceback);
extern void trace_V5_4(FILE *outf, int **wbx, int **vx, int j, int d, int mid, 
		       struct tracekn_s *curr_tr, struct traceknstack_s *dolist, int traceback);
extern void trace_V5_5(FILE *outf, int **wbx, int **vx, int j, int d, int mid,
		       struct tracekn_s *curr_tr, struct traceknstack_s *dolist, int traceback);
extern void trace_V5_6(FILE *outf, int **wbx, int **vx, int j, int d, int mid, 
		       struct tracekn_s *curr_tr, struct traceknstack_s *dolist, int traceback);
extern void trace_V5_7(FILE *outf, int **wbx, int **vx, int j, int d, int mid, 
		       struct tracekn_s *curr_tr, struct traceknstack_s *dolist, int traceback);
extern void trace_V5_8(FILE *outf, int **wbx, int **vx, int j, int d, int mid, 
		       struct tracekn_s *curr_tr, struct traceknstack_s *dolist, int traceback);

extern int V6_1 (ESL_DSQ *s, int len, struct rnapar_2 *rnapar, int **icfg, int **vx, int j, int d, int mid, int mid1, int mid2);
extern int V6_2a(ESL_DSQ *s, int len, struct rnapar_2 *rnapar, int **icfg, int **vx, int j, int d, int mid, int mid1, int mid2);
extern int V6_2b(ESL_DSQ *s, int len, struct rnapar_2 *rnapar, int **icfg, int **vx, int j, int d, int mid, int mid1, int mid2);
extern int V6_2c(ESL_DSQ *s, int len, struct rnapar_2 *rnapar, int **icfg, int **vx, int j, int d, int mid, int mid1, int mid2);
extern int V6_3a(ESL_DSQ *s, int len, struct rnapar_2 *rnapar, int **icfg, int **vx, int j, int d, int mid, int mid1, int mid2);
extern int V6_3b(ESL_DSQ *s, int len, struct rnapar_2 *rnapar, int **icfg, int **vx, int j, int d, int mid, int mid1, int mid2);
extern int V6_3c(ESL_DSQ *s, int len, struct rnapar_2 *rnapar, int **icfg, int **vx, int j, int d, int mid, int mid1, int mid2);
extern int V6_4a(ESL_DSQ *s, int len, struct rnapar_2 *rnapar, int **icfg, int **vx, int j, int d, int mid, int mid1, int mid2);
extern int V6_4b(ESL_DSQ *s, int len, struct rnapar_2 *rnapar, int **icfg, int **vx, int j, int d, int mid, int mid1, int mid2);
extern int V6_4c(ESL_DSQ *s, int len, struct rnapar_2 *rnapar, int **icfg, int **vx, int j, int d, int mid, int mid1, int mid2);
extern int V6_4d(ESL_DSQ *s, int len, struct rnapar_2 *rnapar, int **icfg, int **vx, int j, int d, int mid, int mid1, int mid2);
extern int V6_4e(ESL_DSQ *s, int len, struct rnapar_2 *rnapar, int **icfg, int **vx, int j, int d, int mid, int mid1, int mid2);
extern int V6_4f(ESL_DSQ *s, int len, struct rnapar_2 *rnapar, int **icfg, int **vx, int j, int d, int mid, int mid1, int mid2);
extern int V6_4g(ESL_DSQ *s, int len, struct rnapar_2 *rnapar, int **icfg, int **vx, int j, int d, int mid, int mid1, int mid2);
extern int V6_4h(ESL_DSQ *s, int len, struct rnapar_2 *rnapar, int **icfg, int **vx, int j, int d, int mid, int mid1, int mid2);
extern int V6_4i(ESL_DSQ *s, int len, struct rnapar_2 *rnapar, int **icfg, int **vx, int j, int d, int mid, int mid1, int mid2);
extern void trace_V6_CC(FILE *outf, int **vx, int j, int d, int mid, int mid1, int mid2, 
			struct tracekn_s *curr_tr, struct traceknstack_s *dolist, int traceback);
extern int V6_5 (ESL_DSQ *s, int len, struct rnapar_2 *rnapar, int **icfg, int **vx, int j, int d, int mid, int mid1, int mid2);
extern int V6_6a(ESL_DSQ *s, int len, struct rnapar_2 *rnapar, int **icfg, int **vx, int j, int d, int mid, int mid1, int mid2);
extern int V6_6b(ESL_DSQ *s, int len, struct rnapar_2 *rnapar, int **icfg, int **vx, int j, int d, int mid, int mid1, int mid2);
extern int V6_6c(ESL_DSQ *s, int len, struct rnapar_2 *rnapar, int **icfg, int **vx, int j, int d, int mid, int mid1, int mid2);
extern int V6_7a(ESL_DSQ *s, int len, struct rnapar_2 *rnapar, int **icfg, int **vx, int j, int d, int mid, int mid1, int mid2);
extern int V6_7b(ESL_DSQ *s, int len, struct rnapar_2 *rnapar, int **icfg, int **vx, int j, int d, int mid, int mid1, int mid2);
extern int V6_7c(ESL_DSQ *s, int len, struct rnapar_2 *rnapar, int **icfg, int **vx, int j, int d, int mid, int mid1, int mid2);
extern int V6_8a(ESL_DSQ *s, int len, struct rnapar_2 *rnapar, int **icfg, int **vx, int j, int d, int mid, int mid1, int mid2);
extern int V6_8b(ESL_DSQ *s, int len, struct rnapar_2 *rnapar, int **icfg, int **vx, int j, int d, int mid, int mid1, int mid2);
extern int V6_8c(ESL_DSQ *s, int len, struct rnapar_2 *rnapar, int **icfg, int **vx, int j, int d, int mid, int mid1, int mid2);
extern int V6_8d(ESL_DSQ *s, int len, struct rnapar_2 *rnapar, int **icfg, int **vx, int j, int d, int mid, int mid1, int mid2);
extern int V6_8e(ESL_DSQ *s, int len, struct rnapar_2 *rnapar, int **icfg, int **vx, int j, int d, int mid, int mid1, int mid2);
extern int V6_8f(ESL_DSQ *s, int len, struct rnapar_2 *rnapar, int **icfg, int **vx, int j, int d, int mid, int mid1, int mid2);
extern int V6_8g(ESL_DSQ *s, int len, struct rnapar_2 *rnapar, int **icfg, int **vx, int j, int d, int mid, int mid1, int mid2);
extern int V6_8h(ESL_DSQ *s, int len, struct rnapar_2 *rnapar, int **icfg, int **vx, int j, int d, int mid, int mid1, int mid2);
extern int V6_8i(ESL_DSQ *s, int len, struct rnapar_2 *rnapar, int **icfg, int **vx, int j, int d, int mid, int mid1, int mid2);
extern void trace_V6_NCC(FILE *outf, int **vx, int j, int d, int mid, int mid1, int mid2, 
			 struct tracekn_s *curr_tr, struct traceknstack_s *dolist, int traceback);

extern int V7_1(ESL_DSQ *s, int len, struct rnapar_2 *rnapar, int **icfg, int ****whx, int j, int d, int mid, int mid1, int mid2);
extern int V7_2(ESL_DSQ *s, int len, struct rnapar_2 *rnapar, int **icfg, int ****whx, int j, int d, int mid, int mid1, int mid2);
extern int V7_3(ESL_DSQ *s, int len, struct rnapar_2 *rnapar, int **icfg, int ****whx, int j, int d, int mid, int mid1, int mid2);
extern int V7_4(ESL_DSQ *s, int len, struct rnapar_2 *rnapar, int **icfg, int ****whx, int j, int d, int mid, int mid1, int mid2);
extern void trace_V7_1(FILE *outf, int ****whx, int j, int d, int mid, int mid1, int mid2, 
		       struct tracekn_s *curr_tr, struct traceknstack_s *dolist, int traceback);
extern void trace_V7_2(FILE *outf, int ****whx, int j, int d, int mid, int mid1, int mid2, 
		       struct tracekn_s *curr_tr, struct traceknstack_s *dolist, int traceback);
extern void trace_V7_3(FILE *outf, int ****whx, int j, int d, int mid, int mid1, int mid2, 
		       struct tracekn_s *curr_tr, struct traceknstack_s *dolist, int traceback);
extern void trace_V7_4(FILE *outf, int ****whx, int j, int d, int mid, int mid1, int mid2, 
		       struct tracekn_s *curr_tr, struct traceknstack_s *dolist, int traceback);

extern int V8_1(ESL_DSQ *s, int len, struct rnapar_2 *rnapar, int **icfg, int ****whx, int ****zhx, 
		int j, int d, int mid, int mid1, int mid2);
extern int V8_2(ESL_DSQ *s, int len, struct rnapar_2 *rnapar, int **icfg, int ****whx, int ****zhx, 
		int j, int d, int mid, int mid1, int mid2);
extern int V8_3(ESL_DSQ *s, int len, struct rnapar_2 *rnapar, int **icfg, int ****whx, int ****zhx, 
		int j, int d, int mid, int mid1, int mid2);
extern int V8_4(ESL_DSQ *s, int len, struct rnapar_2 *rnapar, int **icfg, int ****whx, int ****zhx, 
		int j, int d, int mid, int mid1, int mid2);
extern void trace_V8_1(FILE *outf, int ****whx, int ****zhx, int j, int d, int mid, int mid1, 
		       int mid2, struct tracekn_s *curr_tr, struct traceknstack_s *dolist, int traceback);
extern void trace_V8_2(FILE *outf, int ****whx, int ****zhx, int j, int d, int mid, int mid1, 
		       int mid2, struct tracekn_s *curr_tr, struct traceknstack_s *dolist, int traceback);
extern void trace_V8_3(FILE *outf, int ****whx, int ****zhx, int j, int d, int mid, int mid1, 
		       int mid2, struct tracekn_s *curr_tr, struct traceknstack_s *dolist, int traceback);
extern void trace_V8_4(FILE *outf, int ****whx, int ****zhx, int j, int d, int mid, int mid1, 
		       int mid2, struct tracekn_s *curr_tr, struct traceknstack_s *dolist, int traceback);

extern int V9_1(ESL_DSQ *s, int len, struct rnapar_2 *rnapar, int **icfg, int ****whx, int ****zhx,
		int j, int d, int mid, int mid1, int mid2);
extern int V9_2(ESL_DSQ *s, int len, struct rnapar_2 *rnapar, int **icfg, int ****whx, int ****zhx,
		int j, int d, int mid, int mid1, int mid2);
extern int V9_3(ESL_DSQ *s, int len, struct rnapar_2 *rnapar, int **icfg, int ****whx, int ****zhx,
		int j, int d, int mid, int mid1, int mid2);
extern int V9_4(ESL_DSQ *s, int len, struct rnapar_2 *rnapar, int **icfg, int ****whx, int ****zhx,
		int j, int d, int mid, int mid1, int mid2);
extern void trace_V9_1(FILE *outf, int ****whx, int ****zhx, int j, int d, int mid, int mid1, 
		       int mid2, struct tracekn_s *curr_tr, struct traceknstack_s *dolist, int traceback);
extern void trace_V9_2(FILE *outf, int ****whx, int ****zhx, int j, int d, int mid, int mid1, 
		       int mid2, struct tracekn_s *curr_tr, struct traceknstack_s *dolist, int traceback);
extern void trace_V9_3(FILE *outf, int ****whx, int ****zhx, int j, int d, int mid, int mid1, 
		       int mid2, struct tracekn_s *curr_tr, struct traceknstack_s *dolist, int traceback);
extern void trace_V9_4(FILE *outf, int ****whx, int ****zhx, int j, int d, int mid, int mid1, 
		       int mid2, struct tracekn_s *curr_tr, struct traceknstack_s *dolist, int traceback);

extern int V10_1(ESL_DSQ *s, int len, struct rnapar_2 *rnapar, int **icfg, int ****yhx,
		 int j, int d, int mid, int mid1, int mid2);
extern int V10_2(ESL_DSQ *s, int len, struct rnapar_2 *rnapar, int **icfg, int ****yhx,
		 int j, int d, int mid, int mid1, int mid2);
extern int V10_3(ESL_DSQ *s, int len, struct rnapar_2 *rnapar, int **icfg, int ****yhx,
		 int j, int d, int mid, int mid1, int mid2);
extern int V10_4(ESL_DSQ *s, int len, struct rnapar_2 *rnapar, int **icfg, int ****yhx,
		 int j, int d, int mid, int mid1, int mid2);
extern int V10_5(ESL_DSQ *s, int len, struct rnapar_2 *rnapar, int **icfg, int ****yhx,
		 int j, int d, int mid, int mid1, int mid2);
extern int V10_6(ESL_DSQ *s, int len, struct rnapar_2 *rnapar, int **icfg, int ****yhx,
		 int j, int d, int mid, int mid1, int mid2);
extern int V10_7(ESL_DSQ *s, int len, struct rnapar_2 *rnapar, int **icfg, int ****yhx,
		   int j, int d, int mid, int mid1, int mid2);
extern int V10_8(ESL_DSQ *s, int len, struct rnapar_2 *rnapar, int **icfg, int ****yhx,
		 int j, int d, int mid, int mid1, int mid2);
extern void trace_V10_1(FILE *outf, int ****yhx, int j, int d, int mid, int mid1, int mid2, 
			struct tracekn_s *curr_tr, struct traceknstack_s *dolist, int traceback);
extern void trace_V10_2(FILE *outf, int ****yhx, int j, int d, int mid, int mid1, int mid2, 
			struct tracekn_s *curr_tr, struct traceknstack_s *dolist, int traceback);
extern void trace_V10_3(FILE *outf, int ****yhx, int j, int d, int mid, int mid1, int mid2, 
			struct tracekn_s *curr_tr, struct traceknstack_s *dolist, int traceback);
extern void trace_V10_4(FILE *outf, int ****yhx, int j, int d, int mid, int mid1, int mid2, 
			struct tracekn_s *curr_tr, struct traceknstack_s *dolist, int traceback);
extern void trace_V10_5(FILE *outf, int ****yhx, int j, int d, int mid, int mid1, int mid2, 
			struct tracekn_s *curr_tr, struct traceknstack_s *dolist, int traceback);
extern void trace_V10_6(FILE *outf, int ****yhx, int j, int d, int mid, int mid1, int mid2, 
			struct tracekn_s *curr_tr, struct traceknstack_s *dolist, int traceback);
extern void trace_V10_7(FILE *outf, int ****yhx, int j, int d, int mid, int mid1, int mid2, 
			struct tracekn_s *curr_tr, struct traceknstack_s *dolist, int traceback);
extern void trace_V10_8(FILE *outf, int ****yhx, int j, int d, int mid, int mid1, int mid2, 
			struct tracekn_s *curr_tr, struct traceknstack_s *dolist, int traceback);
