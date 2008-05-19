/* globals.c
 * 
 * 
 * 
 */

#include <stdio.h>
#include "cfg.h"

#ifdef MEMDEBUG
#include "dbmalloc.h"
#endif

/* Connects defines which node types connect to each other.
 * The order must be the same as in cfg.h, as this is indexed
 * by dpcP, dpcL, etc. #define's.
 */
int Connects[NNODES][NNODES] = {
/* from:   to:  dpcP dpcL dpcR dpcB dpcS dpcML dpcMR dpcPS dpcPL dpcPI dpcX  dpcV  dpcW  dpcWA dpcTL dpcE  */
/* dpcP */   {   1,   0,   0,   1,   0,   1,    1,    1,    0,    0,    1,     1,    1,   0,    1,    0 },
/* dpcL */   {   1,   1,   1,   1,   0,   0,    0,    0,    0,    0,    0,     0,    0,   0,    0,    1 },
/* dpcR */   {   1,   0,   1,   1,   0,   0,    0,    0,    0,    0,    0,     0,    0,   0,    0,    0 },
/* dpcB */   {   0,   0,   0,   0,   1,   1,    1,    0,    0,    0,     0,    0,   0,    0,    0 },
/* dpcS */   {   1,   1,   1,   1,   0,   0,    0,    1,    1,    0,     1,    0,   0,    1,    1 },
/* dpcML */  {   1,   0,   0,   1,   0,   1,    1,    0,    0,    0,     0,    0,   0,    0,    0 },
/* dpcMR */  {   1,   0,   0,   1,   0,   1,    1,    0,    0,    0,     0,    0,   0,    0,    0 },
/* dpcPL*/   {   0,   0,   0,   0,   0,   0,    0,    0,    0,    1,     0,    0,   0,    0,    0 },

/* dpcPS*/   {   1,   0,   0,   1,   0,   0,    0,    1,    1,    1,     1,    1,    1,   0,    1,    0 },

/* dpcPI*/   {   0,   0,   0,   0,   0,   0,    0,    0,    0,    0,     0,    1,   0,    0,    0 },
/* dpcX */   {   0,   0,   0,   0,   0,   0,    0,    0,    0,    0,     0,    0,   0,    0,    1 },
/* dpcV */   {   1,   1,   1,   1,   0,   0,    0,    1,    1,    0,     0,    0,   0,    0,    1 },
/* dpcW */   {   1,   1,   1,   1,   0,   0,    0,    1,    1,    0,     0,    0,   0,    0,    1 },
/* dpcWA*/   {   1,   1,   1,   1,   0,   0,    0,    1,    1,    0,     0,    0,   0,    0,    1 },
/* dpcTL*/   {   0,   0,   0,   0,   0,   0,    0,    0,    0,    0,     0,    0,   0,    0,    1 },
/* dpcE */   {   0,   0,   0,   0,   0,   0,    0,    0,    0,    0,     0,    0,   0,    0,    0 },
};

/* Statenum tells how many states there are for each node type
 */
int Statenum[NNODES] = {
/*   dpcP dpcL dpcR dpcB dpcS dpcML dpcMR dpcPS dpcPL dpcPI dpcX dpcV dpcW dpcWA dpcTL dpcE  */
      16,  4,   4,   1,   1,   4,    4,    16,   16,   16,   31,  31,  31,  4,    30,   1 
};

/* Startidx tells where the state indices begin for a given node type.
 * This has to match up with idxP, etc. in cfg.h
 */
int Startidx[NNODES] = {
/*   dpcP dpcL dpcR dpcB dpcS dpcML dpcMR dpcPS dpcPL dpcPI dpcX  dpcV  dpcW  dpcWA dpcTL dpcE  */
      0,   16,  20,  24,  25,  26,   30,   34,   50,   66,   82,   113,  144,  175,  179,  4275  
};


/* Ntype maps state types to node types. The order
 * in this array is tied to definitions in cfg.h.
 */
int Ntype[NSTATES] = {
  dpcP, dpcP, dpcP, dpcP, dpcP, dpcP, dpcP, dpcP,
  dpcP, dpcP, dpcP, dpcP, dpcP, dpcP, dpcP, dpcP,
  dpcL, dpcL, dpcL, dpcL,
  dpcR, dpcR, dpcR, dpcR,
  dpcB,
  dpcS,
  dpcPL, dpcPL, dpcPL, dpcPL, dpcPL, dpcPL, dpcPL, dpcPL,
  dpcPL, dpcPL, dpcPL, dpcPL, dpcPL, dpcPL, dpcPL, dpcPL,
  dpcPI, dpcPI, dpcPI, dpcPI, dpcPI, dpcPI, dpcPI, dpcPI,
  dpcPI, dpcPI, dpcPI, dpcPI, dpcPI, dpcPI, dpcPI, dpcPI,
  dpcX, dpcX, dpcX, dpcX, dpcX, dpcX, dpcX, dpcX, dpcX, dpcX,
  dpcX, dpcX, dpcX, dpcX, dpcX, dpcX, dpcX, dpcX, dpcX, dpcX,
  dpcX, dpcX, dpcX, dpcX, dpcX, dpcX, dpcX, dpcX, dpcX, dpcX, dpcX, 
  dpcV, dpcV, dpcV, dpcV, dpcV, dpcV, dpcV, dpcV, dpcV, dpcV,
  dpcV, dpcV, dpcV, dpcV, dpcV, dpcV, dpcV, dpcV, dpcV, dpcV,
  dpcV, dpcV, dpcV, dpcV, dpcV, dpcV, dpcV, dpcV, dpcV, dpcV, dpcV,
  dpcW, dpcW, dpcW, dpcW, dpcW, dpcW, dpcW, dpcW, dpcW, dpcW,
  dpcW, dpcW, dpcW, dpcW, dpcW, dpcW, dpcW, dpcW, dpcW, dpcW,
  dpcW, dpcW, dpcW, dpcW, dpcW, dpcW, dpcW, dpcW, dpcW, dpcW, dpcW, 
  dpcWA, dpcWA, dpcWA, dpcWA,
  dpcTL, dpcTL, dpcTL, dpcTL, dpcTL, dpcTL, dpcTL, dpcTL, dpcTL, dpcTL,  
  dpcTL, dpcTL, dpcTL, dpcTL, dpcTL, dpcTL, dpcTL, dpcTL, dpcTL, dpcTL,
  dpcTL, dpcTL, dpcTL, dpcTL, dpcTL, dpcTL, dpcTL, dpcTL, dpcTL, dpcTL,
  dpcE
};

/* The two NAME arrays are used for debugging purposes,
 * to convert integer index values to strings
 */
char *dpcNAME[NNODES] = { 
  "pair", 
  "ssL",     "ssR", 
  "bifurc", 
  "start", 
  "bulge",  
  "int",   
  "hairpn",
  "bulge",  
  "int",
  "int asy",
  "tetraloops",
  "end" 
};

char *stNAME[NSTATES] = {
  "P.aa", "P.ac", "P.ag", "P.au",
  "P.ca", "P.cc", "P.cg", "P.cu",
  "P.ga", "P.gc", "P.gg", "P.gu",
  "P.ua", "P.uc", "P.ug", "P.uu",
  "L.a", "L.c", "L.g", "L.u",
  "R.a", "R.c", "R.g", "R.u",
  "B",
  "S",
  "PL.aa", "PL.ac", "PL.ag", "PL.au",
  "PL.ca", "PL.cc", "PL.cg", "PL.cu",
  "PL.ga", "PL.gc", "PL.gg", "PL.gu",
  "PL.ua", "PL.uc", "PL.ug", "PL.uu",
  "PR.aa", "PR.ac", "PR.ag", "PR.au",
  "PR.ca", "PR.cc", "PR.cg", "PR.cu",
  "PR.ga", "PR.gc", "PR.gg", "PR.gu",
  "PR.ua", "PR.uc", "PR.ug", "PR.uu",
  "X.0",  "X.1",  "X.2",  "X.3",  "X.4",  "X.5",  "X.6",  "X.7",  "X.8",  "X.9",  
  "X.10", "X.11", "X.12", "X.13", "X.14", "X.15", "X.16", "X.17", "X.18", "X.19",  
  "X.20", "X.21", "X.22", "X.23", "X.24", "X.25", "X.26", "X.27", "X.28", "X.29", "X.30",  
  "V.0",  "V.1",  "V.2",  "V.3",  "V.4",  "V.5",  "V.6",  "V.7",  "V.8",  "V.9",  
  "V.10", "V.11", "V.12", "V.13", "V.14", "V.15", "V.16", "V.17", "V.18", "V.19",  
  "V.20", "V.21", "V.22", "V.23", "V.24", "V.25", "V.26", "V.27", "V.28", "V.29", "V.30",  
  "W.0",  "W.1",  "W.2",  "W.3",  "W.4",  "W.5",  "W.6",  "W.7",  "W.8",  "W.9",  
  "W.10", "W.11", "W.12", "W.13", "W.14", "W.15", "W.16", "W.17", "W.18", "W.19",  
  "W.20", "W.21", "W.22", "W.23", "W.24", "W.25", "W.26", "W.27", "W.28", "W.29", "W.30",  
  "WA.1",  "WA.2",  "WA.3",  "WA.4",
  "TL.GGGGAC", "TL.GGUGAC", "TL.CGAAAG", "TL.GGAGAC", "TL.CGCAAG", "TL.GGAAAC", "TL.CGGAAG", "TL.CUUCGG", "TL.CGUGAG", "TL.CGAAGG", 
  "TL.CUACGG", "TL.GGCAAC", "TL.CGCGAG", "TL.UGAGAG", "TL.CGAGAG", "TL.AGAAAU", "TL.CGUAAG", "TL.CUAACG", "TL.UGAAAG", "TL.GGAAGC", 
  "TL.GGGAAC", "TL.UGAAAA", "TL.AGCAAU", "TL.AGUAAU", "TL.CGGGAG", "TL.AGUGAU", "TL.GGCGAC", "TL.GGGAGC", "TL.GUGAAC", "TL.UGGAAA", 
  "E",
};





