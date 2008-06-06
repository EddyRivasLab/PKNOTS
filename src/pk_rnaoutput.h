extern void  Tracekn(struct tracekn_s *tr, ESL_DSQ *dsq, int rlen, int watsoncrick, char *ss);     
extern int   IsRNAComplement(char sym1, char sym2, int allow_gu);
extern int   IsRNAComplementDigital(int sym1, int sym2, int allow_gu);
extern int   WriteSeqkn(FILE *outf, ESL_ALPHABET *abc, ESL_SQ *sq, int *ct, int ctoutput, 
			struct rnapar_2 zkn_param, int format, int shuffleseq, 
			int allow_pseudoknots, int approx, CYKVAL sc);
