extern int  SaveSCFG(FILE *ofp, float **cfg);
extern int  ReadSCFG(FILE *ofp, float ***ret_cfg);
extern void WriteRdbSCFG(FILE *ofp, float **cfg);
extern void WriteRdbISCFG(FILE *ofp, int **icfg);
extern void WriteRdbSummary(FILE *ofp, float **cfg);

