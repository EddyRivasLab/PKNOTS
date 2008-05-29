/* (modified) 
 * protovx.h
 * ANSI prototypes for all external functions in vxgraphs.c used in filltrvx.c.
 * 
 */

#include <stdio.h>
#include "cfg.h"
#include "squid.h"

/* from vxgraph.c:  the diagrams that contribute to vx
 */

extern int V1(int *s, int len, struct rnapar_2 *rnapar, int **icfg, int j, int d);
extern void trace_V1(FILE *outf, int **vx,  int j, int d, 
		     struct tracekn_s *curr_tr, struct traceknstack_s *dolist, int traceback);

extern int V2(int *s, int len, struct rnapar_2 *rnapar, int **icfg, int **vx, int j, int d, int mid1, int mid2);
extern void trace_V2(FILE *outf, int **vx, int *s, int len, struct rnapar_2 *rnapar, int **icfg, 
		     int j, int d, int mid1, int mid2, 
		     struct tracekn_s *curr_tr, struct traceknstack_s *dolist, int traceback);

extern int V3_1(int *s, int len, struct rnapar_2 *rnapar, int **icfg, int **wbx, int j, int d, int mid);
extern int V3_2(int *s, int len, struct rnapar_2 *rnapar, int **icfg, int **wbx, int j, int d, int mid);
extern int V3_3(int *s, int len, struct rnapar_2 *rnapar, int **icfg, int **wbx, int j, int d, int mid);
extern int V3_4(int *s, int len, struct rnapar_2 *rnapar, int **icfg, int **wbx, int j, int d, int mid);
extern void trace_V3_1(FILE *outf, int **wbx, int j, int d, int mid, 
		       struct tracekn_s *curr_tr, struct traceknstack_s *dolist, int traceback);
extern void trace_V3_2(FILE *outf, int **wbx, int j, int d, int mid,
		       struct tracekn_s *curr_tr, struct traceknstack_s *dolist, int traceback);
extern void trace_V3_3(FILE *outf, int **wbx, int j, int d, int mid,
		       struct tracekn_s *curr_tr, struct traceknstack_s *dolist, int traceback);
extern void trace_V3_4(FILE *outf, int **wbx, int j, int d, int mid, 
		       struct tracekn_s *curr_tr, struct traceknstack_s *dolist, int traceback);

extern int V4_1(int *s, int len, struct rnapar_2 *rnapar, int **icfg, int **wbx, int **vx, int j, int d, int mid);
extern int V4_2(int *s, int len, struct rnapar_2 *rnapar, int **icfg, int **wbx, int **vx, int j, int d, int mid);
extern int V4_3(int *s, int len, struct rnapar_2 *rnapar, int **icfg, int **wbx, int **vx, int j, int d, int mid);
extern int V4_4(int *s, int len, struct rnapar_2 *rnapar, int **icfg, int **wbx, int **vx, int j, int d, int mid);
extern int V4_5(int *s, int len, struct rnapar_2 *rnapar, int **icfg, int **wbx, int **vx, int j, int d, int mid);
extern int V4_6(int *s, int len, struct rnapar_2 *rnapar, int **icfg, int **wbx, int **vx, int j, int d, int mid);
extern int V4_7(int *s, int len, struct rnapar_2 *rnapar, int **icfg, int **wbx, int **vx, int j, int d, int mid);
extern int V4_8(int *s, int len, struct rnapar_2 *rnapar, int **icfg, int **wbx, int **vx, int j, int d, int mid);
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

extern int V5_1(int *s, int len, struct rnapar_2 *rnapar, int **icfg, int **wbx, int **vx, int j, int d, int mid);
extern int V5_2(int *s, int len, struct rnapar_2 *rnapar, int **icfg, int **wbx, int **vx, int j, int d, int mid);
extern int V5_3(int *s, int len, struct rnapar_2 *rnapar, int **icfg, int **wbx, int **vx, int j, int d, int mid);
extern int V5_4(int *s, int len, struct rnapar_2 *rnapar, int **icfg, int **wbx, int **vx, int j, int d, int mid);
extern int V5_5(int *s, int len, struct rnapar_2 *rnapar, int **icfg, int **wbx, int **vx, int j, int d, int mid);
extern int V5_6(int *s, int len, struct rnapar_2 *rnapar, int **icfg, int **wbx, int **vx, int j, int d, int mid);
extern int V5_7(int *s, int len, struct rnapar_2 *rnapar, int **icfg, int **wbx, int **vx, int j, int d, int mid);
extern int V5_8(int *s, int len, struct rnapar_2 *rnapar, int **icfg, int **wbx, int **vx, int j, int d, int mid);
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

extern int V6_1 (int *s, int len, struct rnapar_2 *rnapar, int **icfg, int **vx, int j, int d, int mid, int mid1, int mid2);
extern int V6_2a(int *s, int len, struct rnapar_2 *rnapar, int **icfg, int **vx, int j, int d, int mid, int mid1, int mid2);
extern int V6_2b(int *s, int len, struct rnapar_2 *rnapar, int **icfg, int **vx, int j, int d, int mid, int mid1, int mid2);
extern int V6_2c(int *s, int len, struct rnapar_2 *rnapar, int **icfg, int **vx, int j, int d, int mid, int mid1, int mid2);
extern int V6_3a(int *s, int len, struct rnapar_2 *rnapar, int **icfg, int **vx, int j, int d, int mid, int mid1, int mid2);
extern int V6_3b(int *s, int len, struct rnapar_2 *rnapar, int **icfg, int **vx, int j, int d, int mid, int mid1, int mid2);
extern int V6_3c(int *s, int len, struct rnapar_2 *rnapar, int **icfg, int **vx, int j, int d, int mid, int mid1, int mid2);
extern int V6_4a(int *s, int len, struct rnapar_2 *rnapar, int **icfg, int **vx, int j, int d, int mid, int mid1, int mid2);
extern int V6_4b(int *s, int len, struct rnapar_2 *rnapar, int **icfg, int **vx, int j, int d, int mid, int mid1, int mid2);
extern int V6_4c(int *s, int len, struct rnapar_2 *rnapar, int **icfg, int **vx, int j, int d, int mid, int mid1, int mid2);
extern int V6_4d(int *s, int len, struct rnapar_2 *rnapar, int **icfg, int **vx, int j, int d, int mid, int mid1, int mid2);
extern int V6_4e(int *s, int len, struct rnapar_2 *rnapar, int **icfg, int **vx, int j, int d, int mid, int mid1, int mid2);
extern int V6_4f(int *s, int len, struct rnapar_2 *rnapar, int **icfg, int **vx, int j, int d, int mid, int mid1, int mid2);
extern int V6_4g(int *s, int len, struct rnapar_2 *rnapar, int **icfg, int **vx, int j, int d, int mid, int mid1, int mid2);
extern int V6_4h(int *s, int len, struct rnapar_2 *rnapar, int **icfg, int **vx, int j, int d, int mid, int mid1, int mid2);
extern int V6_4i(int *s, int len, struct rnapar_2 *rnapar, int **icfg, int **vx, int j, int d, int mid, int mid1, int mid2);
extern void trace_V6_CC(FILE *outf, int **vx, int j, int d, int mid, int mid1, int mid2, 
			struct tracekn_s *curr_tr, struct traceknstack_s *dolist, int traceback);
extern int V6_5 (int *s, int len, struct rnapar_2 *rnapar, int **icfg, int **vx, int j, int d, int mid, int mid1, int mid2);
extern int V6_6a(int *s, int len, struct rnapar_2 *rnapar, int **icfg, int **vx, int j, int d, int mid, int mid1, int mid2);
extern int V6_6b(int *s, int len, struct rnapar_2 *rnapar, int **icfg, int **vx, int j, int d, int mid, int mid1, int mid2);
extern int V6_6c(int *s, int len, struct rnapar_2 *rnapar, int **icfg, int **vx, int j, int d, int mid, int mid1, int mid2);
extern int V6_7a(int *s, int len, struct rnapar_2 *rnapar, int **icfg, int **vx, int j, int d, int mid, int mid1, int mid2);
extern int V6_7b(int *s, int len, struct rnapar_2 *rnapar, int **icfg, int **vx, int j, int d, int mid, int mid1, int mid2);
extern int V6_7c(int *s, int len, struct rnapar_2 *rnapar, int **icfg, int **vx, int j, int d, int mid, int mid1, int mid2);
extern int V6_8a(int *s, int len, struct rnapar_2 *rnapar, int **icfg, int **vx, int j, int d, int mid, int mid1, int mid2);
extern int V6_8b(int *s, int len, struct rnapar_2 *rnapar, int **icfg, int **vx, int j, int d, int mid, int mid1, int mid2);
extern int V6_8c(int *s, int len, struct rnapar_2 *rnapar, int **icfg, int **vx, int j, int d, int mid, int mid1, int mid2);
extern int V6_8d(int *s, int len, struct rnapar_2 *rnapar, int **icfg, int **vx, int j, int d, int mid, int mid1, int mid2);
extern int V6_8e(int *s, int len, struct rnapar_2 *rnapar, int **icfg, int **vx, int j, int d, int mid, int mid1, int mid2);
extern int V6_8f(int *s, int len, struct rnapar_2 *rnapar, int **icfg, int **vx, int j, int d, int mid, int mid1, int mid2);
extern int V6_8g(int *s, int len, struct rnapar_2 *rnapar, int **icfg, int **vx, int j, int d, int mid, int mid1, int mid2);
extern int V6_8h(int *s, int len, struct rnapar_2 *rnapar, int **icfg, int **vx, int j, int d, int mid, int mid1, int mid2);
extern int V6_8i(int *s, int len, struct rnapar_2 *rnapar, int **icfg, int **vx, int j, int d, int mid, int mid1, int mid2);
extern void trace_V6_NCC(FILE *outf, int **vx, int j, int d, int mid, int mid1, int mid2, 
			 struct tracekn_s *curr_tr, struct traceknstack_s *dolist, int traceback);

extern int V7_1(int *s, int len, struct rnapar_2 *rnapar, int **icfg, int ****whx, int j, int d, int mid, int mid1, int mid2);
extern int V7_2(int *s, int len, struct rnapar_2 *rnapar, int **icfg, int ****whx, int j, int d, int mid, int mid1, int mid2);
extern int V7_3(int *s, int len, struct rnapar_2 *rnapar, int **icfg, int ****whx, int j, int d, int mid, int mid1, int mid2);
extern int V7_4(int *s, int len, struct rnapar_2 *rnapar, int **icfg, int ****whx, int j, int d, int mid, int mid1, int mid2);
extern void trace_V7_1(FILE *outf, int ****whx, int j, int d, int mid, int mid1, int mid2, 
		       struct tracekn_s *curr_tr, struct traceknstack_s *dolist, int traceback);
extern void trace_V7_2(FILE *outf, int ****whx, int j, int d, int mid, int mid1, int mid2, 
		       struct tracekn_s *curr_tr, struct traceknstack_s *dolist, int traceback);
extern void trace_V7_3(FILE *outf, int ****whx, int j, int d, int mid, int mid1, int mid2, 
		       struct tracekn_s *curr_tr, struct traceknstack_s *dolist, int traceback);
extern void trace_V7_4(FILE *outf, int ****whx, int j, int d, int mid, int mid1, int mid2, 
		       struct tracekn_s *curr_tr, struct traceknstack_s *dolist, int traceback);

extern int V8_1(int *s, int len, struct rnapar_2 *rnapar, int **icfg, int ****whx, int ****zhx, 
		int j, int d, int mid, int mid1, int mid2);
extern int V8_2(int *s, int len, struct rnapar_2 *rnapar, int **icfg, int ****whx, int ****zhx, 
		int j, int d, int mid, int mid1, int mid2);
extern int V8_3(int *s, int len, struct rnapar_2 *rnapar, int **icfg, int ****whx, int ****zhx, 
		int j, int d, int mid, int mid1, int mid2);
extern int V8_4(int *s, int len, struct rnapar_2 *rnapar, int **icfg, int ****whx, int ****zhx, 
		int j, int d, int mid, int mid1, int mid2);
extern void trace_V8_1(FILE *outf, int ****whx, int ****zhx, int j, int d, int mid, int mid1, 
		       int mid2, struct tracekn_s *curr_tr, struct traceknstack_s *dolist, int traceback);
extern void trace_V8_2(FILE *outf, int ****whx, int ****zhx, int j, int d, int mid, int mid1, 
		       int mid2, struct tracekn_s *curr_tr, struct traceknstack_s *dolist, int traceback);
extern void trace_V8_3(FILE *outf, int ****whx, int ****zhx, int j, int d, int mid, int mid1, 
		       int mid2, struct tracekn_s *curr_tr, struct traceknstack_s *dolist, int traceback);
extern void trace_V8_4(FILE *outf, int ****whx, int ****zhx, int j, int d, int mid, int mid1, 
		       int mid2, struct tracekn_s *curr_tr, struct traceknstack_s *dolist, int traceback);

extern int V9_1(int *s, int len, struct rnapar_2 *rnapar, int **icfg, int ****whx, int ****zhx,
		int j, int d, int mid, int mid1, int mid2);
extern int V9_2(int *s, int len, struct rnapar_2 *rnapar, int **icfg, int ****whx, int ****zhx,
		int j, int d, int mid, int mid1, int mid2);
extern int V9_3(int *s, int len, struct rnapar_2 *rnapar, int **icfg, int ****whx, int ****zhx,
		int j, int d, int mid, int mid1, int mid2);
extern int V9_4(int *s, int len, struct rnapar_2 *rnapar, int **icfg, int ****whx, int ****zhx,
		int j, int d, int mid, int mid1, int mid2);
extern void trace_V9_1(FILE *outf, int ****whx, int ****zhx, int j, int d, int mid, int mid1, 
		       int mid2, struct tracekn_s *curr_tr, struct traceknstack_s *dolist, int traceback);
extern void trace_V9_2(FILE *outf, int ****whx, int ****zhx, int j, int d, int mid, int mid1, 
		       int mid2, struct tracekn_s *curr_tr, struct traceknstack_s *dolist, int traceback);
extern void trace_V9_3(FILE *outf, int ****whx, int ****zhx, int j, int d, int mid, int mid1, 
		       int mid2, struct tracekn_s *curr_tr, struct traceknstack_s *dolist, int traceback);
extern void trace_V9_4(FILE *outf, int ****whx, int ****zhx, int j, int d, int mid, int mid1, 
		       int mid2, struct tracekn_s *curr_tr, struct traceknstack_s *dolist, int traceback);

extern int V10_1(int *s, int len, struct rnapar_2 *rnapar, int **icfg, int ****yhx,
		 int j, int d, int mid, int mid1, int mid2);
extern int V10_2(int *s, int len, struct rnapar_2 *rnapar, int **icfg, int ****yhx,
		 int j, int d, int mid, int mid1, int mid2);
extern int V10_3(int *s, int len, struct rnapar_2 *rnapar, int **icfg, int ****yhx,
		 int j, int d, int mid, int mid1, int mid2);
extern int V10_4(int *s, int len, struct rnapar_2 *rnapar, int **icfg, int ****yhx,
		 int j, int d, int mid, int mid1, int mid2);
extern int V10_5(int *s, int len, struct rnapar_2 *rnapar, int **icfg, int ****yhx,
		 int j, int d, int mid, int mid1, int mid2);
extern int V10_6(int *s, int len, struct rnapar_2 *rnapar, int **icfg, int ****yhx,
		 int j, int d, int mid, int mid1, int mid2);
extern int V10_7(int *s, int len, struct rnapar_2 *rnapar, int **icfg, int ****yhx,
		   int j, int d, int mid, int mid1, int mid2);
extern int V10_8(int *s, int len, struct rnapar_2 *rnapar, int **icfg, int ****yhx,
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
