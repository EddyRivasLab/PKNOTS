extern void FillVP(ESL_DSQ *s, int len, struct rnapar_2 *rnapar, int **icfg, 
		   int **vx, int *vp, int j, int d);
extern void FillVX_nested(ESL_DSQ *s, int len, struct rnapar_2 *rnapar, int **icfg, 
			  int **wbx, int **vx, int *vp, int j, int d, 
			  int allow_coaxials);
extern void TraceVX_nested(FILE *outf, ESL_DSQ *s, int len, struct rnapar_2 *rnapar, int **icfg, 
			   int **wbx, int **vx, int j, int d, int *flag,
			   struct tracekn_s *curr_tr, struct traceknstack_s *dolist, int traceback);
extern void  FillVX(ESL_DSQ *s, int len, struct rnapar_2 *rnapar, int **icfg, 
		    int **wbx, int **vx, int *vp,
		    int ****whx, int ****zhx, int ****yhx, int j, int d, int allow_coaxials, int approx);
extern void TraceVX(FILE *outf, ESL_DSQ *s, int len, struct rnapar_2 *rnapar, 
		    int **icfg, int **wbx, int **vx,  
		    int ****whx, int ****zhx, int ****yhx, int j, int d, int approx,
		    struct tracekn_s *curr_tr, struct traceknstack_s *dolist, int traceback);

