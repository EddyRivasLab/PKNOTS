extern struct trace_s *InitTrace(void);
extern struct tracekn_s *InitTracekn(void);
extern struct trace_s *AttachTrace(struct trace_s *parent, 
				   int emitr, int emitl, int type);
extern struct tracekn_s *AttachTracekn(struct tracekn_s *parent, 
				       int emiti, int emitj, int emitk, int emitl, 
                                       int type1, int type2);
extern void   FreeTrace(struct trace_s *tr);
extern void   FreeTracekn(struct tracekn_s *tr);
extern void   DeleteTracenode(struct trace_s *oldtr);
extern void   DeleteTraceknnode(struct tracekn_s *oldtr);
extern struct tracestack_s *InitTracestack(void);
extern struct traceknstack_s *InitTraceknstack(void);
extern void   PushTracestack(struct tracestack_s *stack, struct trace_s *node);
extern void   PushTraceknstack(struct traceknstack_s *stack, struct tracekn_s *node);
extern struct trace_s *PopTracestack(struct tracestack_s *stack);
extern struct tracekn_s *PopTraceknstack(struct traceknstack_s *stack);
extern void   FreeTracestack(struct tracestack_s *stack);
extern void   FreeTraceknstack(struct traceknstack_s *stack);
extern void   TraceCount(char *seq, int len, float wgt, 
			 struct trace_s *tr, float **cfg);
extern void   TraceknCount(char *seq, int len, float wgt, 
			   struct tracekn_s *tr, float **cfg);
extern float  TraceScore(char *seq, int len, struct trace_s *tr, int **icfg);
extern float  TraceknScore(char *seq, int len, struct tracekn_s *tr, int **icfg);
extern void   GraphicTrace(FILE *fp, struct trace_s *tr, char *seq, 
			   int len, int **icfg);
extern void   GraphicTracekn(FILE *fp, struct tracekn_s *tr, char *seq, 
			     int len, int **icfg);





