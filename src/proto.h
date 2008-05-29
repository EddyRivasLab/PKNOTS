/* (modified) 
 * proto.h
 * ANSI prototypes for all external functions.
 * 
 */

#include <stdio.h>
#include "cfg.h"

#include <easel.h>
#include <esl_sqio.h>
#include <esl_wuss.h>

/* from cfgio.c: I/O of models to/from files
 */
extern int  SaveSCFG(FILE *ofp, float **cfg);
extern int  ReadSCFG(FILE *ofp, float ***ret_cfg);
extern void WriteRdbSCFG(FILE *ofp, float **cfg);
extern void WriteRdbISCFG(FILE *ofp, int **icfg);
extern void WriteRdbSummary(FILE *ofp, float **cfg);

/* from dewachter.c:
 * Parsing De Wachter rRNA databases
 */
extern int ParseBelgianRNA(char *rna, int **helnum, int nhel,
			   char **ret_aseq, char **ret_ss, 
			   int *ret_alen, int *ret_rlen,
			   int *ret_nstem, int *ret_nseg, int *ret_nassign,
			   char errdiag[128]);
extern void MakeHelnum(char *dewach, int alen, int **ret_helnum);

/* from fillmtx.c: fill matrices
 */
extern void FillMtx_nested(int *s, int len, struct rnapar_2 *rnapar,
			   int **icfg, int **wx, int **wbx, int **vx, int *vp, 
			   int allow_coaxials);
extern void FillMtx(int *s, int len, struct rnapar_2 *rnapar, int **icfg, 
		    int **wx, int **wbx, int **vx, int *vp,
		    int ****whx, int ****vhx, int ****zhx, int ****yhx, int allow_coaxials, int approx);

/* from filltrvx.c: fill and traceback vx matrix
 */
extern void FillVP(int *s, int len, struct rnapar_2 *rnapar, int **icfg, 
		   int **vx, int *vp, int j, int d);
extern void FillVX_nested(int *s, int len, struct rnapar_2 *rnapar, int **icfg, 
			  int **wbx, int **vx, int *vp, int j, int d, 
			  int allow_coaxials);
extern void TraceVX_nested(FILE *outf, int *s, int len, struct rnapar_2 *rnapar, int **icfg, 
			   int **wbx, int **vx, int j, int d, int *flag,
			   struct tracekn_s *curr_tr, struct traceknstack_s *dolist, int traceback);
extern void  FillVX(int *s, int len, struct rnapar_2 *rnapar, int **icfg, 
		    int **wbx, int **vx, int *vp,
		    int ****whx, int ****zhx, int ****yhx, int j, int d, int allow_coaxials, int approx);
extern void TraceVX(FILE *outf, int *s, int len, struct rnapar_2 *rnapar, 
		    int **icfg, int **wbx, int **vx,  
		    int ****whx, int ****zhx, int ****yhx, int j, int d, int approx,
		    struct tracekn_s *curr_tr, struct traceknstack_s *dolist, int traceback);

/* from filltrwx.c: fill and traceback wx, wbx matrices
 */
extern void  FillWX_nested(int *s, int len, struct rnapar_2 *rnapar, 
			   int **icfg, int **wx, int **vx, int j, int d);
extern void TraceWX_nested(FILE *outf, int *s, int len, struct rnapar_2 *rnapar, int **icfg, 
			   int **wx, int **vx, int j, int d, int *flag,  
			   struct tracekn_s *curr_tr, struct traceknstack_s *dolist, int traceback);
extern void  FillWX(int *s, int len, struct rnapar_2 *rnapar, 
		    int **icfg, int **wx, int **vx, 
		    int ****whx, int ****yhx, int j, int d);
extern void TraceWX(FILE *outf, int *s, int len, struct rnapar_2 *rnapar, 
		    int **icfg, int **wx, int **vx, 
		    int ****whx, int ****yhx, int j, int d, 
		    struct tracekn_s *curr_tr, struct traceknstack_s *dolist, int traceback);

/* from filltrwxwbx.c: fill and traceback wx, wbx matrices
 */
extern void FillWBX_nested(int *s, int len, struct rnapar_2 *rnapar, 
			   int **icfg, int **wbx, int **vx, int j, int d);
extern void TraceWBX_nested(FILE *outf, int *s, int len, struct rnapar_2 *rnapar, int **icfg, 
			    int **wbx, int **vx, int j, int d, int *flag,  
			    struct tracekn_s *curr_tr, struct traceknstack_s *dolist, int traceback);
extern void  FillWBX(int *s, int len, struct rnapar_2 *rnapar, int **icfg, 
		     int **wbx, int **vx, 
		     int ****whx, int ****yhx, int j, int d, int approx);
extern void TraceWBX(FILE *outf, int *s, int len, struct rnapar_2 *rnapar, 
		     int **icfg, int **wbx, int **vx, 
		     int ****whx, int ****yhx, int j, int d, int approx, 
		     struct tracekn_s *curr_tr, struct traceknstack_s *dolist, int traceback);

/* from filltrvhx.c: fill and traceback vhx matrix
 */
extern void FillVHX(int *s, int len, struct rnapar_2 *rnapar, 
		    int **icfg, int ****whx, int ****vhx, 
		    int j, int d, int d1, int d2, int min_kn);
extern void TraceVHX(FILE *outf, int *s, int len, struct rnapar_2 *rnapar, 
		     int **icfg, int ****whx, int ****vhx, 
		     int j, int d, int d1, int d2, 
		     struct tracekn_s *curr_tr, struct traceknstack_s *dolist, int traceback);

/* from filltrzhxyhx.c: fill and traceback zhx, yhx  matrices
 */
extern void FillZHX(int *s, int len, struct rnapar_2 *rnapar, 
		    int **icfg,int **wbx, int **vx, 
		    int ****whx, int ****vhx, int ****zhx, 
		    int j, int d, int d1, int d2);
extern void TraceZHX(FILE *outf, int *s, int len, struct rnapar_2 *rnapar, 
		     int **icfg, int **wbx, int **vx, 
		     int ****whx, int ****vhx, int ****zhx,
		     int j, int d, int d1, int d2, 
		     struct tracekn_s *curr_tr, struct traceknstack_s *dolist, int traceback);
extern void FillYHX(int *s, int len, struct rnapar_2 *rnapar,
		    int **icfg, int **wbx, int **vx, 
		    int ****whx, int ****vhx, int ****yhx,
		    int j, int d, int d1, int d2, int min_kn);
extern void TraceYHX(FILE *outf, int *s, int len, struct rnapar_2 *rnapar, 
		     int **icfg, int **wbx, int **vx, 
		     int ****whx, int ****vhx, int ****yhx,
		     int j, int d, int d1, int d2, 
		     struct tracekn_s *curr_tr, struct traceknstack_s *dolist, int traceback);

/* from filltrwhx.c: fill and traceback whx matrix
 */
extern void FillWHX(int *s, int len, struct rnapar_2 *rnapar, 
		    int **icfg, int **wbx, int **vx, 
		    int ****whx, int ****vhx, int ****zhx, int ****yhx,
		    int j, int d, int d1, int d2);
extern void TraceWHX(FILE *outf, int *s, int len, struct rnapar_2 *rnapar, 
		     int **icfg, int **wbx, int **vx, 
		     int ****whx, int ****vhx, int ****zhx, int ****yhx,
		     int j, int d, int d1, int d2,
		     struct tracekn_s *curr_tr, struct traceknstack_s *dolist, int traceback);

/* from fillvxscan.c: fill and traceback vx matrix
 */
extern void  FillVPScan(int *s, int len, int win, struct rnapar_2 *rnapar, int **icfg, 
			int **vx, int *vp, int j, int jmod, int d);
extern void  FillVX_nestedScan(int *s, int len, int win, struct rnapar_2 *rnapar, int **icfg, 
			       int **wbx, int **vx, int *vp, int j, int jmod, int d);

/* from fillwbxscan.c: fill and traceback wx, wbx matrices
 */
extern void FillWBX_nestedScan(int *s, int len, int win, struct rnapar_2 *rnapar, 
			       int **icfg, int **wbx, int **vx,  
			       int j, int jmod, int d);

/* from fillwxscan.c: fill and traceback wx, wbx matrices 
 */
extern void  FillWX_nestedScan(int *s, int len, int win, struct rnapar_2 *rnapar, 
			       int **icfg, int **wx, int **vx,  
			       int j, int jmod, int d);

/* from irredsurf.c: calculates F1 and F2
 */
extern int F1(int *s, int len, struct rnapar_2 *rnapar, int **icfg, int j, int d);
extern int F2(int *s, int len, struct rnapar_2 *rnapar, int **icfg, int j, int d, int d1, int d2);

/* from misc.c: miscellaneous functions with no home
 */
extern void  StripDegeneracy(char *seq);
extern void  IntizeSequence(char *seq, int len, int **ret_iseq);
extern void  ShuffleSequence(char *seq, int len, int endpos, int verbose);
extern char *DupSeq(char *seq, int j, int d);

/* from model.c: stuff to do with handling a CFG as an object
 */
extern float **AllocSCFG(void);
extern int   **AllocIntSCFG(void);
extern void    AllocBaseFreq(float **ret_basefreq);
extern void    BaseFreq(char *s, int j, int d, float **ret_basefreq);
extern float **DupSCFG(float **cfg);
extern void    FreeSCFG(float **cfg);
extern void    FreeIntSCFG(int **icfg);
extern int   **LogifySCFG(float **cfg);
extern int   **LogoddsifySCFG(float **cfg);
extern void    RandomSCFG(float **cfg);
extern void    NormalizeSCFG(float **cfg);
extern int   **NussinovIntSCFG(void);
extern void    ProbifySCFG(float **cfg);
extern int     Stype(int node, int symi, int symj, int size, int asy, int tlp);

/* from rnaio.c:
 */
extern int KHS2ct(char *ss, int len, int allow_pseudoknots, int **ret_ct);
extern int VerifyKHS(char *name, char *ss, int wordy);

/* from rnaoutput.c: conversion of internal data structures into
 *                   more user-friendly forms
 */
extern void  Tracekn(struct tracekn_s *tr, char *seq, int rlen, int watsoncrick, char **ret_ss);     
extern void  Traceintkn(struct tracekn_s *tr, char *seq, int rlen, int watsoncrick, int **ret_ss);
extern int   WriteSeqkn(FILE *outf, char *seq, SQINFO *sqinfo, int *ss, int format, 
			float *ret_pairs, float *cykret_pairs);
extern int   IsRNAComplement(char sym1, char sym2, int allow_gu);
extern void  CalculatePairs(SQINFO sqinfo, int *ss, int format, float *ret_pairs, float *ret_cykpairs);
extern void  PrintCtSeq(FILE *ofp, SQINFO *sqinfo, char *s, int start, int L, char *ss);
extern void  CompareRNAStructures(FILE *ofp, int start, int L, char *ss_true, int *cc);


/* from rnaparam.c: read thermodinamic order 2  parameters
 */
extern void    Parameters2_Zkn(struct rnapar_2 *rnapar);
extern int   **ParamIntSCFG(struct rnapar_2 *rnapar);
extern int     IntizeScale (double val);

/* from tracemtx.c: trace back best score matrices
 */
extern void TraceMtx_nested(FILE *outf, int *s, int len, struct rnapar_2 *rnapar, int **icfg, 
			    int **wx, int **wbx, int **vx, 
			    int j, int d, struct tracekn_s **ret_trace, int traceback);
extern void TraceMtx(FILE *outf, int *s, int len, struct rnapar_2 *rnapar, int **icfg,  
		     int **wx, int **wbx, int **vx, 
		     int ****whx, int ****vhx, int ****zhx, int ****yhx,
		     int j, int d, int d1, int d2, int approx, 
		     struct tracekn_s **ret_trace, int traceback);

/* 
 * from trace.c
 */
extern struct trace_s *InitTrace(void);
extern struct tracekn_s *InitTracekn(void);
extern struct trace_s *AttachTrace(struct trace_s *parent, 
				   int emitr, int emitl, int type);
extern struct tracekn_s *AttachTracekn(struct tracekn_s *parent, 
				       int emiti, int emitj, int emitk, int emitl, 
                                       int type1, int type2);
extern void   FreeTrace(struct trace_s *tr);
extern void   FreeTracekn(struct tracekn_s *tr);
extern void   DeleteTracenode(struct trace_s *oldtr);
extern void   DeleteTraceknnode(struct tracekn_s *oldtr);
extern struct tracestack_s *InitTracestack(void);
extern struct traceknstack_s *InitTraceknstack(void);
extern void   PushTracestack(struct tracestack_s *stack, struct trace_s *node);
extern void   PushTraceknstack(struct traceknstack_s *stack, struct tracekn_s *node);
extern struct trace_s *PopTracestack(struct tracestack_s *stack);
extern struct tracekn_s *PopTraceknstack(struct traceknstack_s *stack);
extern void   FreeTracestack(struct tracestack_s *stack);
extern void   FreeTraceknstack(struct traceknstack_s *stack);
extern void   TraceCount(char *seq, int len, float wgt, 
			 struct trace_s *tr, float **cfg);
extern void   TraceknCount(char *seq, int len, float wgt, 
			   struct tracekn_s *tr, float **cfg);
extern float  TraceScore(char *seq, int len, struct trace_s *tr, int **icfg);
extern float  TraceknScore(char *seq, int len, struct tracekn_s *tr, int **icfg);
extern void   GraphicTrace(FILE *fp, struct trace_s *tr, char *seq, 
			   int len, int **icfg);
extern void   GraphicTracekn(FILE *fp, struct tracekn_s *tr, char *seq, 
			     int len, int **icfg);


/* from tying.c: knocking down the number of free parameters
 */
extern int   *EnslaveStates(void);
extern void   TieCounts(float **cfg, int *enslave);
extern void   CountFreeParameters(int *enslave, int *ret_free, int *ret_nonzero);
extern float *CountsPerState(float **cfg);
extern void   BootstrapConfidence(float **cfg, float *counts, int *enslave, 
				  int nboot, float conf,
				  float ***ret_high, float ***ret_low);


/* from viterbi_nuss.c:  the dynamic programming alignment algorithm
 */
extern void StructurePredictkn_nuss(char *seq, int len, int **icfg, int verbose, 
				    struct tracekn_s **ret_tr, float *ret_score);
extern int *fill_F1(int *s, int len, int **icfg, int j, int d);
extern int *fill_F2(int *s, int len, int **icfg, int j, int d, int d1, int d2);

/* from viterbi.c:  the dynamic programming alignment algorithm
 */
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
extern void StructurePredictkn_2IS(FILE *outf, char *seq, int len, struct rnapar_2 *rnapar,
				   int **icfg, int verbose, int traceback, struct tracekn_s **ret_tr, 
				   int *ret_score, int allow_coaxials, int allow_pseudoknots, 
				   int approx);

/* from zscores.c:  functions to calculate Zscores
 */
extern void CalculateZscores(FILE *outf, char *seq, int len, int **icfg, int verbose, float score, 
			     float *ret_zscore, int allow_coaxials, int allow_pseudoknots, int approx);


/* from zscorescan.c:  functions to calculate Zscores
 */
extern void CalculateZscoreScan(FILE *outf, char *seq, int seqlen, int len, int pos_end, 
				int win, int off, int **icfg, int verbose, int allow_coaxials, 
				int allow_pseudoknots, int approx);
extern void IdxWindow(int win, int off, int j, int *jmax, int *dmax);










