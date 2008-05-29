/* (modified) 
 * protozhxyhx.h
 * ANSI prototypes for all external functions in vhxgraphs.c used in filltrvhx.c.
 * 
 */

#include <stdio.h>
#include "cfg.h"
#include "squid.h"

/* from zhxyhxgraph.c:  the diagrams that contribute to zhx
 */

extern int ZH1(int *s, int len, struct rnapar_2 *rnapar, int **icfg, int ****vhx, int j, int d, int d1, int d2);
extern void trace_ZH1(FILE *outf, int ****vhx, int j, int d, int d1, int d2, 
		      struct tracekn_s *curr_tr, struct traceknstack_s *dolist, int traceback);

extern int ZH2(int *s, int len, struct rnapar_2 *rnapar, int **icfg, int ****vhx, int j, int d, int d1, int d2);
extern void trace_ZH2(FILE *outf, int ****vhx, int j, int d, int d1, int d2, 
		      struct tracekn_s *curr_tr, struct traceknstack_s *dolist, int traceback);

extern int ZH3(int *s, int len, struct rnapar_2 *rnapar, int **icfg, int ****vhx, int j, int d, int d1, int d2);
extern void trace_ZH3(FILE *outf, int ****vhx, int j, int d, int d1, int d2, 
		      struct tracekn_s *curr_tr, struct traceknstack_s *dolist, int traceback);

extern int ZH4(int *s, int len, struct rnapar_2 *rnapar, int **icfg, int ****vhx, int j, int d, int d1, int d2);
extern void trace_ZH4(FILE *outf, int ****vhx, int j, int d, int d1, int d2, 
		      struct tracekn_s *curr_tr, struct traceknstack_s *dolist, int traceback);

extern int ZH5(int *s, int len, struct rnapar_2 *rnapar, int **icfg, int ****zhx, int j, int d, int d1, int d2);
extern void trace_ZH5(FILE *outf, int ****zhx, int j, int d, int d1, int d2, 
		      struct tracekn_s *curr_tr, struct traceknstack_s *dolist, int traceback);

extern int ZH6(int *s, int len, struct rnapar_2 *rnapar, int **icfg, int ****zhx, int j, int d, int d1, int d2);
extern void trace_ZH6(FILE *outf, int ****zhx, int j, int d, int d1, int d2, 
		      struct tracekn_s *curr_tr, struct traceknstack_s *dolist, int traceback);


extern int ZH7(int *s, int len, struct rnapar_2 *rnapar, int **icfg, int **wbx, int ****zhx, 
	       int j, int d, int d1, int d2, int mid);
extern void trace_ZH7(FILE *outf, int **wbx, int ****zhx, int j, int d, int d1, int d2, int mid,
		      struct tracekn_s *curr_tr, struct traceknstack_s *dolist, int traceback);

extern int ZH8_1(int *s, int len, struct rnapar_2 *rnapar, int **icfg, int **vx, int ****vhx, 
		 int j, int d, int d1, int d2, int mid);
extern void trace_ZH8_1(FILE *outf, int **vx, int ****vhx, int j, int d, int d1, int d2, int mid,
			struct tracekn_s *curr_tr, struct traceknstack_s *dolist, int traceback);

extern int ZH8_2(int *s, int len, struct rnapar_2 *rnapar, int **icfg, int **vx, int ****vhx, 
		 int j, int d, int d1, int d2, int mid);
extern void trace_ZH8_2(FILE *outf, int **vx, int ****vhx, int j, int d, int d1, int d2, int mid,
			struct tracekn_s *curr_tr, struct traceknstack_s *dolist, int traceback);

extern int ZH8_3(int *s, int len, struct rnapar_2 *rnapar, int **icfg, int **vx, int ****vhx, 
		 int j, int d, int d1, int d2, int mid);
extern void trace_ZH8_3(FILE *outf, int **vx, int ****vhx, int j, int d, int d1, int d2, int mid,
			struct tracekn_s *curr_tr, struct traceknstack_s *dolist, int traceback);

extern int ZH8_4(int *s, int len, struct rnapar_2 *rnapar, int **icfg, int **vx, int ****vhx, 
		 int j, int d, int d1, int d2, int mid);
extern void trace_ZH8_4(FILE *outf, int **vx, int ****vhx, int j, int d, int d1, int d2, int mid,
			struct tracekn_s *curr_tr, struct traceknstack_s *dolist, int traceback);

extern int ZH9(int *s, int len, struct rnapar_2 *rnapar, int **icfg, int **wbx, int ****zhx, 
	       int j, int d, int d1, int d2, int mid);
extern void trace_ZH9(FILE *outf, int **wbx, int ****zhx, int j, int d, int d1, int d2, int mid,
		      struct tracekn_s *curr_tr, struct traceknstack_s *dolist, int traceback);

extern int ZH10_1(int *s, int len, struct rnapar_2 *rnapar, int **icfg, int **vx, int ****vhx, 
		  int j, int d, int d1, int d2, int mid);
extern void trace_ZH10_1(FILE *outf, int **vx, int ****vhx, int j, int d, int d1, int d2, int mid,
			 struct tracekn_s *curr_tr, struct traceknstack_s *dolist, int traceback);

extern int ZH10_2(int *s, int len, struct rnapar_2 *rnapar, int **icfg, int **vx, int ****vhx, 
		  int j, int d, int d1, int d2, int mid);
extern void trace_ZH10_2(FILE *outf, int **vx, int ****vhx, int j, int d, int d1, int d2, int mid,
			 struct tracekn_s *curr_tr, struct traceknstack_s *dolist, int traceback);

extern int ZH10_3(int *s, int len, struct rnapar_2 *rnapar, int **icfg, int **vx, int ****vhx, 
		  int j, int d, int d1, int d2, int mid);
extern void trace_ZH10_3(FILE *outf, int **vx, int ****vhx, int j, int d, int d1, int d2, int mid,
			 struct tracekn_s *curr_tr, struct traceknstack_s *dolist, int traceback);

extern int ZH10_4(int *s, int len, struct rnapar_2 *rnapar, int **icfg, int **vx, int ****vhx, 
		  int j, int d, int d1, int d2, int mid);
extern void trace_ZH10_4(FILE *outf, int **vx, int ****vhx, int j, int d, int d1, int d2, int mid,
			 struct tracekn_s *curr_tr, struct traceknstack_s *dolist, int traceback);

extern int ZH11(int *s, int len, struct rnapar_2 *rnapar, int **icfg, int ****zhx, int ****vhx, 
		int j, int d, int d1, int d2, int mid1, int mid2);
extern void trace_ZH11(FILE *outf, int ****zhx, int ****vhx, int j, int d, int d1, int d2, 
		       int mid1, int mid2, struct tracekn_s *curr_tr, 
		       struct traceknstack_s *dolist, int traceback);

extern int ZH12_1(int *s, int len, struct rnapar_2 *rnapar, int **icfg, int ****whx, int j, int d, int d1, int d2);
extern void trace_ZH12_1(FILE *outf, int ****whx, int j, int d, int d1, int d2,
			struct tracekn_s *curr_tr, struct traceknstack_s *dolist, int traceback);
extern int ZH12_2(int *s, int len, struct rnapar_2 *rnapar, int **icfg, int ****whx, int j, int d, int d1, int d2);
extern void trace_ZH12_2(FILE *outf, int ****whx, int j, int d, int d1, int d2,
			struct tracekn_s *curr_tr, struct traceknstack_s *dolist, int traceback);
extern int ZH12_3(int *s, int len, struct rnapar_2 *rnapar, int **icfg, int ****whx, int j, int d, int d1, int d2);
extern void trace_ZH12_3(FILE *outf, int ****whx, int j, int d, int d1, int d2,
			struct tracekn_s *curr_tr, struct traceknstack_s *dolist, int traceback);
extern int ZH12_4(int *s, int len, struct rnapar_2 *rnapar, int **icfg, int ****whx, int j, int d, int d1, int d2);
extern void trace_ZH12_4(FILE *outf, int ****whx, int j, int d, int d1, int d2,
			struct tracekn_s *curr_tr, struct traceknstack_s *dolist, int traceback);


/*  from zhxyhxgraph.c: the diagrams that contribute to yhx
 */

extern int YH1(int *s, int len, struct rnapar_2 *rnapar, int **icfg, int ****vhx, int j, int d, int d1, int d2);
extern void trace_YH1(FILE *outf, int ****vhx, int j, int d, int d1, int d2, 
		      struct tracekn_s *curr_tr, struct traceknstack_s *dolist, int traceback);

extern int YH2(int *s, int len, struct rnapar_2 *rnapar, int **icfg, int ****vhx, int j, int d, int d1, int d2);
extern void trace_YH2(FILE *outf, int ****vhx, int j, int d, int d1, int d2, 
		      struct tracekn_s *curr_tr, struct traceknstack_s *dolist, int traceback);

extern int YH3(int *s, int len, struct rnapar_2 *rnapar, int **icfg, int ****vhx, int j, int d, int d1, int d2);
extern void trace_YH3(FILE *outf, int ****vhx, int j, int d, int d1, int d2, 
		      struct tracekn_s *curr_tr, struct traceknstack_s *dolist, int traceback);

extern int YH4(int *s, int len, struct rnapar_2 *rnapar, int **icfg, int ****vhx, int j, int d, int d1, int d2);
extern void trace_YH4(FILE *outf, int ****vhx, int j, int d, int d1, int d2, 
		      struct tracekn_s *curr_tr, struct traceknstack_s *dolist, int traceback);

extern int YH5(int *s, int len, struct rnapar_2 *rnapar, int **icfg, int ****yhx, int j, int d, int d1, int d2);
extern void trace_YH5(FILE *outf, int ****yhx, int j, int d, int d1, int d2, 
		      struct tracekn_s *curr_tr, struct traceknstack_s *dolist, int traceback);

extern int YH6(int *s, int len, struct rnapar_2 *rnapar, int **icfg, int ****yhx, int j, int d, int d1, int d2);
extern void trace_YH6(FILE *outf, int ****yhx, int j, int d, int d1, int d2, 
		      struct tracekn_s *curr_tr, struct traceknstack_s *dolist, int traceback);

extern int YH7(int *s, int len, struct rnapar_2 *rnapar, int **icfg, int **wbx, int ****yhx, 
	       int j, int d, int d1, int d2, int mid);
extern void trace_YH7(FILE *outf, int **wbx, int ****yhx, int j, int d, int d1, int d2, int mid,
		      struct tracekn_s *curr_tr, struct traceknstack_s *dolist, int traceback);

extern int YH8_1(int *s, int len, struct rnapar_2 *rnapar, int **icfg, int **vx, int ****vhx, 
		 int j, int d, int d1, int d2, int mid);
extern void trace_YH8_1(FILE *outf, int **vx, int ****vhx, int j, int d, int d1, int d2, int mid,
			struct tracekn_s *curr_tr, struct traceknstack_s *dolist, int traceback);

extern int YH8_2(int *s, int len, struct rnapar_2 *rnapar, int **icfg, int **vx, int ****vhx, 
		 int j, int d, int d1, int d2, int mid);
extern void trace_YH8_2(FILE *outf, int **vx, int ****vhx, int j, int d, int d1, int d2, int mid,
			struct tracekn_s *curr_tr, struct traceknstack_s *dolist, int traceback);

extern int YH8_3(int *s, int len, struct rnapar_2 *rnapar, int **icfg, int **vx, int ****vhx, 
		 int j, int d, int d1, int d2, int mid);
extern void trace_YH8_3(FILE *outf, int **vx, int ****vhx, int j, int d, int d1, int d2, int mid,
			struct tracekn_s *curr_tr, struct traceknstack_s *dolist, int traceback);

extern int YH8_4(int *s, int len, struct rnapar_2 *rnapar, int **icfg, int **vx, int ****vhx, 
		 int j, int d, int d1, int d2, int mid);
extern void trace_YH8_4(FILE *outf, int **vx, int ****vhx, int j, int d, int d1, int d2, int mid,
			struct tracekn_s *curr_tr, struct traceknstack_s *dolist, int traceback);

extern int YH9(int *s, int len, struct rnapar_2 *rnapar, int **icfg, int **wbx, int ****yhx, 
	       int j, int d, int d1, int d2, int mid);
extern void trace_YH9(FILE *outf, int **wbx, int ****yhx, int j, int d, int d1, int d2, int mid,
		      struct tracekn_s *curr_tr, struct traceknstack_s *dolist, int traceback);

extern int YH10_1(int *s, int len, struct rnapar_2 *rnapar, int **icfg, int **vx, int ****vhx, 
		  int j, int d, int d1, int d2, int mid);
extern void trace_YH10_1(FILE *outf, int **vx, int ****vhx, int j, int d, int d1, int d2, int mid,
			 struct tracekn_s *curr_tr, struct traceknstack_s *dolist, int traceback);

extern int YH10_2(int *s, int len, struct rnapar_2 *rnapar, int **icfg, int **vx, int ****vhx, 
		  int j, int d, int d1, int d2, int mid);
extern void trace_YH10_2(FILE *outf, int **vx, int ****vhx, int j, int d, int d1, int d2, int mid,
			 struct tracekn_s *curr_tr, struct traceknstack_s *dolist, int traceback);

extern int YH10_3(int *s, int len, struct rnapar_2 *rnapar, int **icfg, int **vx, int ****vhx, 
		  int j, int d, int d1, int d2, int mid);
extern void trace_YH10_3(FILE *outf, int **vx, int ****vhx, int j, int d, int d1, int d2, int mid,
			 struct tracekn_s *curr_tr, struct traceknstack_s *dolist, int traceback);

extern int YH10_4(int *s, int len, struct rnapar_2 *rnapar, int **icfg, int **vx, int ****vhx, 
		  int j, int d, int d1, int d2, int mid);
extern void trace_YH10_4(FILE *outf, int **vx, int ****vhx, int j, int d, int d1, int d2, int mid,
			 struct tracekn_s *curr_tr, struct traceknstack_s *dolist, int traceback);

extern int YH11(int *s, int len, struct rnapar_2 *rnapar, int **icfg, int ****yhx, int ****vhx, 
		int j, int d, int d1, int d2, int mid1, int mid2);
extern void trace_YH11(FILE *outf, int ****yhx, int ****vhx, int j, int d, int d1, int d2, 
		       int mid1, int mid2, struct tracekn_s *curr_tr,
		       struct traceknstack_s *dolist, int traceback);

extern int YH12_1(int *s, int len, struct rnapar_2 *rnapar, int **icfg, int ****whx, int j, int d, int d1, int d2);
extern void trace_YH12_1(FILE *outf, int ****whx, int j, int d, int d1, int d2,
			struct tracekn_s *curr_tr, struct traceknstack_s *dolist, int traceback);
extern int YH12_2(int *s, int len, struct rnapar_2 *rnapar, int **icfg, int ****whx, int j, int d, int d1, int d2);
extern void trace_YH12_2(FILE *outf, int ****whx, int j, int d, int d1, int d2,
			struct tracekn_s *curr_tr, struct traceknstack_s *dolist, int traceback);
extern int YH12_3(int *s, int len, struct rnapar_2 *rnapar, int **icfg, int ****whx, int j, int d, int d1, int d2);
extern void trace_YH12_3(FILE *outf, int ****whx, int j, int d, int d1, int d2,
			struct tracekn_s *curr_tr, struct traceknstack_s *dolist, int traceback);
extern int YH12_4(int *s, int len, struct rnapar_2 *rnapar, int **icfg, int ****whx, int j, int d, int d1, int d2);
extern void trace_YH12_4(FILE *outf, int ****whx, int j, int d, int d1, int d2,
			struct tracekn_s *curr_tr, struct traceknstack_s *dolist, int traceback);

