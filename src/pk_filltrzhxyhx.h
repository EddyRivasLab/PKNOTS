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

