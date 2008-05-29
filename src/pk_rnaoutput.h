extern void  Tracekn(struct tracekn_s *tr, char *seq, int rlen, int watsoncrick, char **ret_ss);     
extern void  Traceintkn(struct tracekn_s *tr, char *seq, int rlen, int watsoncrick, int **ret_ss);
extern int   IsRNAComplement(char sym1, char sym2, int allow_gu);
extern int   WriteSeqkn(FILE *outf, ESL_SQ *sq, int *ss, float *ret_cykpairs);

