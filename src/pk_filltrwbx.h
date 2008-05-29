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

