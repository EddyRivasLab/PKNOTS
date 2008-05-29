/* (modified) 
 * protovhx.h
 * ANSI prototypes for all external functions in whxgraphs.c used in filltrwhx.c.
 * 
 */

#include <stdio.h>

/* from whxgraph.c:  the diagrams that contribute to whx
 */

extern int WH1_1(ESL_DSQ *s, int len, struct rnapar_2 *rnapar, int **icfg, int ****vhx, int j, int d, int d1, int d2);
extern void trace_WH1_1(FILE *outf, int ****vhx, int j, int d, int d1, int d2, 
		        struct tracekn_s *curr_tr, struct traceknstack_s *dolist, int traceback);
extern int WH1_2(ESL_DSQ *s, int len, struct rnapar_2 *rnapar, int **icfg, int ****vhx, int j, int d, int d1, int d2);
extern void trace_WH1_2(FILE *outf, int ****vhx, int j, int d, int d1, int d2, 
		        struct tracekn_s *curr_tr, struct traceknstack_s *dolist, int traceback);

extern int WH1_3(ESL_DSQ *s, int len, struct rnapar_2 *rnapar, int **icfg, int ****vhx, int j, int d, int d1, int d2);
extern void trace_WH1_3(FILE *outf, int ****vhx, int j, int d, int d1, int d2, 
		        struct tracekn_s *curr_tr, struct traceknstack_s *dolist, int traceback);

extern int WH1_4(ESL_DSQ *s, int len, struct rnapar_2 *rnapar, int **icfg, int ****vhx, int j, int d, int d1, int d2);
extern void trace_WH1_4(FILE *outf, int ****vhx, int j, int d, int d1, int d2, 
		        struct tracekn_s *curr_tr, struct traceknstack_s *dolist, int traceback);

extern int WH1_5(ESL_DSQ *s, int len, struct rnapar_2 *rnapar, int **icfg, int ****vhx, int j, int d, int d1, int d2);
extern void trace_WH1_5(FILE *outf, int ****vhx, int j, int d, int d1, int d2, 
		        struct tracekn_s *curr_tr, struct traceknstack_s *dolist, int traceback);

extern int WH1_6(ESL_DSQ *s, int len, struct rnapar_2 *rnapar, int **icfg, int ****vhx, int j, int d, int d1, int d2);
extern void trace_WH1_6(FILE *outf, int ****vhx, int j, int d, int d1, int d2, 
		        struct tracekn_s *curr_tr, struct traceknstack_s *dolist, int traceback);

extern int WH1_7(ESL_DSQ *s, int len, struct rnapar_2 *rnapar, int **icfg, int ****vhx, int j, int d, int d1, int d2);
extern void trace_WH1_7(FILE *outf, int ****vhx, int j, int d, int d1, int d2, 
		        struct tracekn_s *curr_tr, struct traceknstack_s *dolist, int traceback);

extern int WH1_8(ESL_DSQ *s, int len, struct rnapar_2 *rnapar, int **icfg, int ****vhx, int j, int d, int d1, int d2);
extern void trace_WH1_8(FILE *outf, int ****vhx, int j, int d, int d1, int d2, 
		        struct tracekn_s *curr_tr, struct traceknstack_s *dolist, int traceback);

extern int WH1_9(ESL_DSQ *s, int len, struct rnapar_2 *rnapar, int **icfg, int ****vhx, int j, int d, int d1, int d2);
extern void trace_WH1_9(FILE *outf, int ****vhx, int j, int d, int d1, int d2, 
		        struct tracekn_s *curr_tr, struct traceknstack_s *dolist, int traceback);

extern int WH1_10(ESL_DSQ *s, int len, struct rnapar_2 *rnapar, int **icfg, int ****vhx, int j, int d, int d1, int d2);
extern void trace_WH1_10(FILE *outf, int ****vhx, int j, int d, int d1, int d2, 
			 struct tracekn_s *curr_tr, struct traceknstack_s *dolist, int traceback);

extern int WH2_1(ESL_DSQ *s, int len, struct rnapar_2 *rnapar, int **icfg, int ****zhx, int j, int d, int d1, int d2);
extern void trace_WH2_1(FILE *outf, int ****zhx, int j, int d, int d1, int d2, 
		        struct tracekn_s *curr_tr, struct traceknstack_s *dolist, int traceback);
extern int WH2_2(ESL_DSQ *s, int len, struct rnapar_2 *rnapar, int **icfg, int ****zhx, int j, int d, int d1, int d2);
extern void trace_WH2_2(FILE *outf, int ****zhx, int j, int d, int d1, int d2, 
		        struct tracekn_s *curr_tr, struct traceknstack_s *dolist, int traceback);
extern int WH2_3(ESL_DSQ *s, int len, struct rnapar_2 *rnapar, int **icfg, int ****zhx, int j, int d, int d1, int d2);
extern void trace_WH2_3(FILE *outf, int ****zhx, int j, int d, int d1, int d2, 
		        struct tracekn_s *curr_tr, struct traceknstack_s *dolist, int traceback);
extern int WH2_4(ESL_DSQ *s, int len, struct rnapar_2 *rnapar, int **icfg, int ****zhx, int j, int d, int d1, int d2);
extern void trace_WH2_4(FILE *outf, int ****zhx, int j, int d, int d1, int d2, 
		        struct tracekn_s *curr_tr, struct traceknstack_s *dolist, int traceback);


extern int WH3_1(ESL_DSQ *s, int len, struct rnapar_2 *rnapar, int **icfg, int ****yhx, int j, int d, int d1, int d2);
extern void trace_WH3_1(FILE *outf, int ****yhx, int j, int d, int d1, int d2, 
		        struct tracekn_s *curr_tr, struct traceknstack_s *dolist, int traceback);
extern int WH3_2(ESL_DSQ *s, int len, struct rnapar_2 *rnapar, int **icfg, int ****yhx, int j, int d, int d1, int d2);
extern void trace_WH3_2(FILE *outf, int ****yhx, int j, int d, int d1, int d2, 
		        struct tracekn_s *curr_tr, struct traceknstack_s *dolist, int traceback);
extern int WH3_3(ESL_DSQ *s, int len, struct rnapar_2 *rnapar, int **icfg, int ****yhx, int j, int d, int d1, int d2);
extern void trace_WH3_3(FILE *outf, int ****yhx, int j, int d, int d1, int d2, 
		        struct tracekn_s *curr_tr, struct traceknstack_s *dolist, int traceback);
extern int WH3_4(ESL_DSQ *s, int len, struct rnapar_2 *rnapar, int **icfg, int ****yhx, int j, int d, int d1, int d2);
extern void trace_WH3_4(FILE *outf, int ****yhx, int j, int d, int d1, int d2, 
		        struct tracekn_s *curr_tr, struct traceknstack_s *dolist, int traceback);

extern int WH4_1(ESL_DSQ *s, int len, struct rnapar_2 *rnapar, int **icfg, int ****whx, int j, int d, int d1, int d2);
extern void trace_WH4_1(FILE *outf, int ****whx, int j, int d, int d1, int d2, 
		        struct tracekn_s *curr_tr, struct traceknstack_s *dolist, int traceback);
extern int WH4_2(ESL_DSQ *s, int len, struct rnapar_2 *rnapar, int **icfg, int ****whx, int j, int d, int d1, int d2);
extern void trace_WH4_2(FILE *outf, int ****whx, int j, int d, int d1, int d2, 
		        struct tracekn_s *curr_tr, struct traceknstack_s *dolist, int traceback);
extern int WH4_3(ESL_DSQ *s, int len, struct rnapar_2 *rnapar, int **icfg, int ****whx, int j, int d, int d1, int d2);
extern void trace_WH4_3(FILE *outf, int ****whx, int j, int d, int d1, int d2, 
		        struct tracekn_s *curr_tr, struct traceknstack_s *dolist, int traceback);
extern int WH4_4(ESL_DSQ *s, int len, struct rnapar_2 *rnapar, int **icfg, int ****whx, int j, int d, int d1, int d2);
extern void trace_WH4_4(FILE *outf, int ****whx, int j, int d, int d1, int d2, 
		        struct tracekn_s *curr_tr, struct traceknstack_s *dolist, int traceback);

extern int WH5(ESL_DSQ *s, int len, struct rnapar_2 *rnapar, int **icfg, int **wbx, int j, int d, int d1, int d2);
extern void trace_WH5(FILE *outf, int **wbx, int j, int d, int d1, int d2, 
		        struct tracekn_s *curr_tr, struct traceknstack_s *dolist, int traceback);

extern int WH6(ESL_DSQ *s, int len, struct rnapar_2 *rnapar, int **icfg, int **wbx, int ****whx, 
	       int j, int d, int d1, int d2, int mid);
extern void trace_WH6(FILE *outf, int **wbx, int ****whx, int j, int d, int d1, int d2, 
		      int mid, struct tracekn_s *curr_tr, struct traceknstack_s *dolist, int traceback);

extern int WH7_1(ESL_DSQ *s, int len, struct rnapar_2 *rnapar, int **icfg, int **vx, int ****zhx, 
		 int j, int d, int d1, int d2, int mid);
extern void trace_WH7_1(FILE *outf, int **vx, int ****zhx, int j, int d, int d1, int d2, 
			int mid, struct tracekn_s *curr_tr, struct traceknstack_s *dolist, int traceback);
extern int WH7_2(ESL_DSQ *s, int len, struct rnapar_2 *rnapar, int **icfg, int **vx, int ****zhx, 
		 int j, int d, int d1, int d2, int mid);
extern void trace_WH7_2(FILE *outf, int **vx, int ****zhx, int j, int d, int d1, int d2, 
			int mid, struct tracekn_s *curr_tr, struct traceknstack_s *dolist, int traceback);
extern int WH7_3(ESL_DSQ *s, int len, struct rnapar_2 *rnapar, int **icfg, int **vx, int ****zhx, 
		 int j, int d, int d1, int d2, int mid);
extern void trace_WH7_3(FILE *outf, int **vx, int ****zhx, int j, int d, int d1, int d2, 
			int mid, struct tracekn_s *curr_tr, struct traceknstack_s *dolist, int traceback);
extern int WH7_4(ESL_DSQ *s, int len, struct rnapar_2 *rnapar, int **icfg, int **vx, int ****zhx, 
		 int j, int d, int d1, int d2, int mid);
extern void trace_WH7_4(FILE *outf, int **vx, int ****zhx, int j, int d, int d1, int d2, 
			int mid, struct tracekn_s *curr_tr, struct traceknstack_s *dolist, int traceback);

extern int WH8(ESL_DSQ *s, int len, struct rnapar_2 *rnapar, int **icfg, int **wbx, int ****whx, 
	       int j, int d, int d1, int d2, int mid);
extern void trace_WH8(FILE *outf, int **wbx, int ****whx, int j, int d, int d1, int d2, 
		      int mid, struct tracekn_s *curr_tr, struct traceknstack_s *dolist, int traceback);

extern int WH9_1(ESL_DSQ *s, int len, struct rnapar_2 *rnapar, int **icfg, int **vx, int ****yhx, 
		 int j, int d, int d1, int d2, int mid);
extern void trace_WH9_1(FILE *outf, int **vx, int ****yhx, int j, int d, int d1, int d2, 
			int mid, struct tracekn_s *curr_tr, struct traceknstack_s *dolist, int traceback);
extern int WH9_2(ESL_DSQ *s, int len, struct rnapar_2 *rnapar, int **icfg, int **vx, int ****yhx, 
		 int j, int d, int d1, int d2, int mid);
extern void trace_WH9_2(FILE *outf, int **vx, int ****yhx, int j, int d, int d1, int d2, 
			int mid, struct tracekn_s *curr_tr, struct traceknstack_s *dolist, int traceback);
extern int WH9_3(ESL_DSQ *s, int len, struct rnapar_2 *rnapar, int **icfg, int **vx, int ****yhx, 
		 int j, int d, int d1, int d2, int mid);
extern void trace_WH9_3(FILE *outf, int **vx, int ****yhx, int j, int d, int d1, int d2, 
			int mid, struct tracekn_s *curr_tr, struct traceknstack_s *dolist, int traceback);
extern int WH9_4(ESL_DSQ *s, int len, struct rnapar_2 *rnapar, int **icfg, int **vx, int ****yhx, 
		 int j, int d, int d1, int d2, int mid);
extern void trace_WH9_4(FILE *outf, int **vx, int ****yhx, int j, int d, int d1, int d2, 
			int mid, struct tracekn_s *curr_tr, struct traceknstack_s *dolist, int traceback);

extern int WH10(ESL_DSQ *s, int len, struct rnapar_2 *rnapar, int **icfg, int **wbx, int ****whx, 
		int j, int d, int d1, int d2, int mid);
extern void trace_WH10(FILE *outf, int **wbx, int ****whx, int j, int d, int d1, int d2, 
		       int mid, struct tracekn_s *curr_tr, 
		       struct traceknstack_s *dolist, int traceback);

extern int WH11_1(ESL_DSQ *s, int len, struct rnapar_2 *rnapar, int **icfg, int **vx, int ****zhx, 
		  int j, int d, int d1, int d2, int mid);
extern void trace_WH11_1(FILE *outf, int **vx, int ****zhx, int j, int d, int d1, int d2, 
			 int mid, struct tracekn_s *curr_tr, 
			 struct traceknstack_s *dolist, int traceback);
extern int WH11_2(ESL_DSQ *s, int len, struct rnapar_2 *rnapar, int **icfg, int **vx, int ****zhx, 
		  int j, int d, int d1, int d2, int mid);
extern void trace_WH11_2(FILE *outf, int **vx, int ****zhx, int j, int d, int d1, int d2, 
			 int mid, struct tracekn_s *curr_tr, 
			 struct traceknstack_s *dolist, int traceback);
extern int WH11_3(ESL_DSQ *s, int len, struct rnapar_2 *rnapar, int **icfg, int **vx, int ****zhx, 
		  int j, int d, int d1, int d2, int mid);
extern void trace_WH11_3(FILE *outf, int **vx, int ****zhx, int j, int d, int d1, int d2, 
			 int mid, struct tracekn_s *curr_tr, 
			 struct traceknstack_s *dolist, int traceback);
extern int WH11_4(ESL_DSQ *s, int len, struct rnapar_2 *rnapar, int **icfg, int **vx, int ****zhx, 
		  int j, int d, int d1, int d2, int mid);
extern void trace_WH11_4(FILE *outf, int **vx, int ****zhx, int j, int d, int d1, int d2, 
			 int mid, struct tracekn_s *curr_tr, 
			 struct traceknstack_s *dolist, int traceback);

extern int WH12(ESL_DSQ *s, int len, struct rnapar_2 *rnapar, int **icfg, int **wbx, int ****whx, 
		int j, int d, int d1, int d2, int mid);
extern void trace_WH12(FILE *outf, int **wbx, int ****whx, int j, int d, int d1, int d2, 
		       int mid, struct tracekn_s *curr_tr, struct traceknstack_s *dolist, int traceback);

extern int WH13_1(ESL_DSQ *s, int len, struct rnapar_2 *rnapar, int **icfg, int **vx, int ****yhx, 
		  int j, int d, int d1, int d2, int mid);
extern void trace_WH13_1(FILE *outf, int **vx, int ****yhx, int j, int d, int d1, int d2, 
			 int mid, struct tracekn_s *curr_tr, struct traceknstack_s *dolist, int traceback);
extern int WH13_2(ESL_DSQ *s, int len, struct rnapar_2 *rnapar, int **icfg, int **vx, int ****yhx, 
		  int j, int d, int d1, int d2, int mid);
extern void trace_WH13_2(FILE *outf, int **vx, int ****yhx, int j, int d, int d1, int d2, 
			 int mid, struct tracekn_s *curr_tr, struct traceknstack_s *dolist, int traceback);
extern int WH13_3(ESL_DSQ *s, int len, struct rnapar_2 *rnapar, int **icfg, int **vx, int ****yhx, 
		  int j, int d, int d1, int d2, int mid);
extern void trace_WH13_3(FILE *outf, int **vx, int ****yhx, int j, int d, int d1, int d2, 
			 int mid, struct tracekn_s *curr_tr, 
			 struct traceknstack_s *dolist, int traceback);
extern int WH13_4(ESL_DSQ *s, int len, struct rnapar_2 *rnapar, int **icfg, int **vx, int ****yhx, 
		  int j, int d, int d1, int d2, int mid);
extern void trace_WH13_4(FILE *outf, int **vx, int ****yhx, int j, int d, int d1, int d2, 
			 int mid, struct tracekn_s *curr_tr, 
			 struct traceknstack_s *dolist, int traceback);

extern int WH14(ESL_DSQ *s, int len, struct rnapar_2 *rnapar, int **icfg, int ****zhx, int ****yhx, 
		int j, int d, int d1, int d2, int mid1, int mid2);
extern void trace_WH14(FILE *outf, int ****zhx, int ****yhx, int j, int d, int d1, int d2, 
		       int mid1, int mid2, struct tracekn_s *curr_tr, 
		       struct traceknstack_s *dolist, int traceback);

extern int WH15(ESL_DSQ *s, int len, struct rnapar_2 *rnapar, int **icfg, int ****whx,  
		int j, int d, int d1, int d2, int mid1, int mid2);
extern void trace_WH15(FILE *outf, int ****whx, int j, int d, int d1, int d2, 
		       int mid1, int mid2, struct tracekn_s *curr_tr, 
		       struct traceknstack_s *dolist, int traceback);

extern int WH16(ESL_DSQ *s, int len, struct rnapar_2 *rnapar, int **icfg, int ****whx,  
		int j, int d, int d1, int d2, int mid1, int mid2);
extern void trace_WH16(FILE *outf, int ****whx, int j, int d, int d1, int d2, 
		       int mid1, int mid2, struct tracekn_s *curr_tr, 
		       struct traceknstack_s *dolist, int traceback);

extern int WH17(ESL_DSQ *s, int len, struct rnapar_2 *rnapar, int **icfg, int ****whx, 
		int j, int d, int d1, int d2, int mid1, int mid2);
extern void trace_WH17(FILE *outf, int ****whx, int j, int d, int d1, int d2, 
		       int mid1, int mid2, struct tracekn_s *curr_tr, 
		       struct traceknstack_s *dolist, int traceback);

extern int WH18_1(ESL_DSQ *s, int len, struct rnapar_2 *rnapar, int **icfg, int ****zhx, int ****yhx, 
		  int j, int d, int d1, int d2, int mid1, int mid2);
extern void trace_WH18_1(FILE *outf, int ****zhx, int ****yhx, int j, int d, int d1, int d2, 
			 int mid1, int mid2, struct tracekn_s *curr_tr, 
			 struct traceknstack_s *dolist, int traceback);
extern int WH18_2(ESL_DSQ *s, int len, struct rnapar_2 *rnapar, int **icfg, int ****zhx, int ****yhx, 
		  int j, int d, int d1, int d2, int mid1, int mid2);
extern void trace_WH18_2(FILE *outf, int ****zhx, int ****yhx, int j, int d, int d1, int d2, 
			 int mid1, int mid2, struct tracekn_s *curr_tr, 
			 struct traceknstack_s *dolist, int traceback);
extern int WH18_3(ESL_DSQ *s, int len, struct rnapar_2 *rnapar, int **icfg, int ****zhx, int ****yhx, 
		  int j, int d, int d1, int d2, int mid1, int mid2);
extern void trace_WH18_3(FILE *outf, int ****zhx, int ****yhx, int j, int d, int d1, int d2, 
			 int mid1, int mid2, struct tracekn_s *curr_tr, 
			 struct traceknstack_s *dolist, int traceback);
extern int WH18_4(ESL_DSQ *s, int len, struct rnapar_2 *rnapar, int **icfg, int ****zhx, int ****yhx, 
		  int j, int d, int d1, int d2, int mid1, int mid2);
extern void trace_WH18_4(FILE *outf, int ****zhx, int ****yhx, int j, int d, int d1, int d2, 
			 int mid1, int mid2, struct tracekn_s *curr_tr, 
			 struct traceknstack_s *dolist, int traceback);

extern int WH19_1(ESL_DSQ *s, int len, struct rnapar_2 *rnapar, int **icfg, int ****yhx, 
		  int j, int d, int d1, int d2, int mid1, int mid2);
extern void trace_WH19_1(FILE *outf, int ****yhx, int j, int d, int d1, int d2, 
			 int mid1, int mid2, struct tracekn_s *curr_tr, 
			 struct traceknstack_s *dolist, int traceback);
extern int WH19_2(ESL_DSQ *s, int len, struct rnapar_2 *rnapar, int **icfg, int ****yhx, 
		  int j, int d, int d1, int d2, int mid1, int mid2);
extern void trace_WH19_2(FILE *outf, int ****yhx, int j, int d, int d1, int d2, 
			 int mid1, int mid2, struct tracekn_s *curr_tr, 
			 struct traceknstack_s *dolist, int traceback);
extern int WH19_3(ESL_DSQ *s, int len, struct rnapar_2 *rnapar, int **icfg, int ****yhx, 
		  int j, int d, int d1, int d2, int mid1, int mid2);
extern void trace_WH19_3(FILE *outf, int ****yhx, int j, int d, int d1, int d2, 
			 int mid1, int mid2, struct tracekn_s *curr_tr, 
			 struct traceknstack_s *dolist, int traceback);
extern int WH19_4(ESL_DSQ *s, int len, struct rnapar_2 *rnapar, int **icfg, int ****yhx, 
		  int j, int d, int d1, int d2, int mid1, int mid2);
extern void trace_WH19_4(FILE *outf, int ****yhx, int j, int d, int d1, int d2, 
			 int mid1, int mid2, struct tracekn_s *curr_tr, 
			 struct traceknstack_s *dolist, int traceback);

extern int WH20(ESL_DSQ *s, int len, struct rnapar_2 *rnapar, int **icfg, int ****whx, 
		int j, int d, int d1, int d2, int mid1, int mid2);
extern void trace_WH20(FILE *outf, int ****whx, int j, int d, int d1, int d2, 
		       int mid1, int mid2, struct tracekn_s *curr_tr, 
		       struct traceknstack_s *dolist, int traceback);

extern int WH21_1(ESL_DSQ *s, int len, struct rnapar_2 *rnapar, int **icfg, int ****zhx, int ****yhx, 
		  int j, int d, int d1, int d2, int mid1, int mid2);
extern void trace_WH21_1(FILE *outf, int ****zhx, int ****yhx, int j, int d, int d1, int d2, 
			 int mid1, int mid2, struct tracekn_s *curr_tr, 
			 struct traceknstack_s *dolist, int traceback);
extern int WH21_2(ESL_DSQ *s, int len, struct rnapar_2 *rnapar, int **icfg, int ****zhx, int ****yhx, 
		  int j, int d, int d1, int d2, int mid1, int mid2);
extern void trace_WH21_2(FILE *outf, int ****zhx, int ****yhx, int j, int d, int d1, int d2, 
			 int mid1, int mid2, struct tracekn_s *curr_tr, 
			 struct traceknstack_s *dolist, int traceback);
extern int WH21_3(ESL_DSQ *s, int len, struct rnapar_2 *rnapar, int **icfg, int ****zhx, int ****yhx, 
		  int j, int d, int d1, int d2, int mid1, int mid2);
extern void trace_WH21_3(FILE *outf, int ****zhx, int ****yhx, int j, int d, int d1, int d2, 
			 int mid1, int mid2, struct tracekn_s *curr_tr, 
			 struct traceknstack_s *dolist, int traceback);
extern int WH21_4(ESL_DSQ *s, int len, struct rnapar_2 *rnapar, int **icfg, int ****zhx, int ****yhx, 
		  int j, int d, int d1, int d2, int mid1, int mid2);
extern void trace_WH21_4(FILE *outf, int ****zhx, int ****yhx, int j, int d, int d1, int d2, 
			 int mid1, int mid2, struct tracekn_s *curr_tr, 
			 struct traceknstack_s *dolist, int traceback);

extern int WH22_1(ESL_DSQ *s, int len, struct rnapar_2 *rnapar, int **icfg, int ****yhx, 
		  int j, int d, int d1, int d2, int mid1, int mid2);
extern void trace_WH22_1(FILE *outf, int ****yhx, int j, int d, int d1, int d2, 
			 int mid1, int mid2, struct tracekn_s *curr_tr, 
			 struct traceknstack_s *dolist, int traceback);
extern int WH22_2(ESL_DSQ *s, int len, struct rnapar_2 *rnapar, int **icfg, int ****yhx, 
		  int j, int d, int d1, int d2, int mid1, int mid2);
extern void trace_WH22_2(FILE *outf, int ****yhx, int j, int d, int d1, int d2, 
			 int mid1, int mid2, struct tracekn_s *curr_tr, 
			 struct traceknstack_s *dolist, int traceback);
extern int WH22_3(ESL_DSQ *s, int len, struct rnapar_2 *rnapar, int **icfg, int ****yhx, 
		  int j, int d, int d1, int d2, int mid1, int mid2);
extern void trace_WH22_3(FILE *outf, int ****yhx, int j, int d, int d1, int d2, 
			 int mid1, int mid2, struct tracekn_s *curr_tr, 
			 struct traceknstack_s *dolist, int traceback);
extern int WH22_4(ESL_DSQ *s, int len, struct rnapar_2 *rnapar, int **icfg, int ****yhx, 
		  int j, int d, int d1, int d2, int mid1, int mid2);
extern void trace_WH22_4(FILE *outf, int ****yhx, int j, int d, int d1, int d2, 
			 int mid1, int mid2, struct tracekn_s *curr_tr, 
			 struct traceknstack_s *dolist, int traceback);







