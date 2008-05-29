extern void  Tracekn(struct tracekn_s *tr, char *seq, int rlen, int watsoncrick, char **ret_ss);     
extern void  Traceintkn(struct tracekn_s *tr, char *seq, int rlen, int watsoncrick, int **ret_ss);
extern int   WriteSeqkn(FILE *outf, char *seq, SQINFO *sqinfo, int *ss, int format, 
			float *ret_pairs, float *cykret_pairs);
extern int   IsRNAComplement(char sym1, char sym2, int allow_gu);
extern void  CalculatePairs(SQINFO sqinfo, int *ss, int format, float *ret_pairs, float *ret_cykpairs);
extern void  PrintCtSeq(FILE *ofp, SQINFO *sqinfo, char *s, int start, int L, char *ss);
extern void  CompareRNAStructures(FILE *ofp, int start, int L, char *ss_true, int *cc);


