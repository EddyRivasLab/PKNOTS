extern int   Tracekn(struct tracekn_s *tr, ESL_SQ *sq, int watsoncrick);     
extern int   Traceintkn(struct tracekn_s *tr, ESL_SQ *sq, int watsoncrick, int **ret_ct);
extern int   IsRNAComplement(char sym1, char sym2, int allow_gu);
extern int   IsRNAComplementDigital(int sym1, int sym2, int allow_gu);
extern int   WriteSeqkn(FILE *outf, ESL_ALPHABET *abc, ESL_SQ *sq, struct tracekn_s *tr, int ctoutput, int stoutput, 
			struct rnapar_2 *zkn_param, int format, int shuffleseq, 
			int allow_pseudoknots, int approx, CYKVAL sc);
