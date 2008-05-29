extern void FillWHX(ESL_DSQ *s, int len, struct rnapar_2 *rnapar, 
		    int **icfg, int **wbx, int **vx, 
		    int ****whx, int ****vhx, int ****zhx, int ****yhx,
		    int j, int d, int d1, int d2);
extern void TraceWHX(FILE *outf, ESL_DSQ *s, int len, struct rnapar_2 *rnapar, 
		     int **icfg, int **wbx, int **vx, 
		     int ****whx, int ****vhx, int ****zhx, int ****yhx,
		     int j, int d, int d1, int d2,
		     struct tracekn_s *curr_tr, struct traceknstack_s *dolist, int traceback);

