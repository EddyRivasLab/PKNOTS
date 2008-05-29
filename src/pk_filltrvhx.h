extern void FillVHX(int *s, int len, struct rnapar_2 *rnapar, 
		    int **icfg, int ****whx, int ****vhx, 
		    int j, int d, int d1, int d2, int min_kn);
extern void TraceVHX(FILE *outf, int *s, int len, struct rnapar_2 *rnapar, 
		     int **icfg, int ****whx, int ****vhx, 
		     int j, int d, int d1, int d2, 
		     struct tracekn_s *curr_tr, struct traceknstack_s *dolist, int traceback);



