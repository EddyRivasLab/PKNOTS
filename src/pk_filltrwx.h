extern void  FillWX_nested(ESL_DSQ *s, int len, struct rnapar_2 *rnapar, 
			   int **icfg, int **wx, int **vx, int j, int d);
extern void TraceWX_nested(FILE *outf, ESL_DSQ *s, int len, struct rnapar_2 *rnapar, int **icfg, 
			   int **wx, int **vx, int j, int d, int *flag,  
			   struct tracekn_s *curr_tr, struct traceknstack_s *dolist, int traceback);
extern void  FillWX(ESL_DSQ *s, int len, struct rnapar_2 *rnapar, 
		    int **icfg, int **wx, int **vx, 
		    int ****whx, int ****yhx, int j, int d);
extern void TraceWX(FILE *outf, ESL_DSQ *s, int len, struct rnapar_2 *rnapar, 
		    int **icfg, int **wx, int **vx, 
		    int ****whx, int ****yhx, int j, int d, 
		    struct tracekn_s *curr_tr, struct traceknstack_s *dolist, int traceback);

