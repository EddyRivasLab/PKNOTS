#ifndef CFG_H_INCLUDED
#define CFG_H_INCLUDED

/* (modified)
 * cfg.h
 * SRE, Wed May 25 09:17:47 1994. (original)
 * Thu Apr 18 10:30:58 1996. (with more extensive loop model)
 * 
 * Master header file.
 * Defines structures and macros for the RNA folding SCFG.
 */

/* Dimensions for malloc of 2D and 4D matrices 
 */
/* Dim2len(x) = \sum_{j=0}^{x-1}\sum_{d=0}^{j} 1 = x(x+1)/2
 */
#define Dim2len(x) ((x) * ((x)+1) / 2)

/* Dim2j(x) = \sum_{d=0}^{x}\sum_{d1=0}^{d} 1 = (x+1)(x+2)/2
 */
#define Dim2j(x)   (((x)+1) * ((x)+2) / 2)

/* Dim2d(x,y) = \sum_{d1=0}^{x}\sum_{d2=0}^{y-d1} 1 = (x+1)(2y-x+2)/2
 */
#define Dim2d(x,y) (((x)+1) * ((y)*2-(x)+2) / 2)  

/* Dim3(x) = \sum_{j=0}^{x-1}\sum_{d=0}^{j}\sum_{d1=0}^{d} 1 = x(x^2+3x+2)/6
 */
#define Dim3(x)    ((x) * ((x)*(x)+3*(x)+2) / 6)

/* Dim3j(x) = \sum_{d=0}^{x}\sum_{d1=0}^{d}\sum_{d2=0}^{d-d1} 1 = (x+1)(x^2+5x+6)/6
 */
#define Dim3j(x)   (((x)+1) * ((x)*(x)+5*(x)+6)  / 6)

/* Dim4(x) = \sum_{j=0}^{x-1}\sum_{d=0}^{j}\sum_{d1=0}^{d}\sum_{d2=0}^{d-d1} 1 = x(x^3+6x^2+11x+6)/24
 */
#define Dim4(x)    ((x) * ((x)*(x)*(x)+6*(x)*(x)+11*(x)+6) /24)


/* Cutoffs  
 */
#define MIN_SCORE   0.0   /* minimum significant score */
#define MIN_LEN     20    /* minimum length of significant structures */

typedef int CYKVAL;

/* The dynamic programming matrix keeps different cells
 * for the different main types of non-terminals in the grammar.
 * dpcE (END) must always be last, to allow it to be
 * conveniently ignored by DP algorithms.
 */
#define dpcP     0		/* pairwise non-terminal                */
#define dpcL     1		/* leftwise non-terminal                */
#define dpcR     2		/* rightwise non-terminal               */
#define dpcB     3		/* bifurcation non-terminal             */
#define dpcS     4		/* start non-terminal                   */
#define dpcML    5		/* leftwise non-terminal                */
#define dpcMR    6		/* rightwise non-terminal               */
#define dpcPL    7		/* hairpin pairwise non-terminal        */
#define dpcPS    8		/* stem pairwise non-terminal           */
#define dpcPI    9		/* interior loop pairwise non-terminal  */
#define dpcX     10             /* hairpin loop                 */
#define dpcV     11             /* bulge loop                   */
#define dpcW     12		/* interior  loop               */
#define dpcWA    13		/* interior  loop               */
#define dpcTL    14		/* tetra-hairpin loop           */
#define dpcE     15             /* end, not kept by DP          */
#define NCELLS   16		/* # of DP cells kept           */
#define NNODES   17             /* total # of node types        */

/* These globals allow us to name the non-terminals (useful for debugging)
 * and to define which transitions are allowed and disallowed.
 * See globals.c for their definitions.
 */
extern char *dpcNAME[NNODES];          /* ASCII cell names  */
extern int   Connects[NNODES][NNODES]; /* allowed transitions in the model */
extern int   Statenum[NNODES];	       /* number of states for this node (globals.c)   */
extern int   Startidx[NNODES];	       /* starting point in state indexing (globals.c) */
extern int   Tying[NNODES][NNODES];    /* tying of nodes for parameter reduction */

/* These macros define the layout of the state transition matrix.
 *   0..15  : pairwise states. AA, AC, .. CA, CC, .. UG, UU.
 *  16..19  : singlet-left states a,c,g,u
 *  20..23  : singlet-right states a,c,g,u
 *  24      : bifurcation
 *  25      : begin
 *  26..29  : singlet-left states in multiloop 
 *  30..33  : singlet-right states in multiloop 
 *  34..49  : pairwise state at the end of a hairpin loop
 *  50..65  : pairwise state at stems
 *  66..81  : pairwise state at the end of an interior loop
 *  82..112 : hairpin from 0 up tp size 30 (size= #interior nucleotides)
 * 113..143 : bulge  from 0 up to size 30
 * 144..174 : interior loop from 0 up to size 30 
 * 175..178 : interior loop asymmetry (4 states)
 * 179..434 : ultrastable tetraloops (13 of a total of 256)
 * 435      : end
 */
#define idxP(symi,symj)    (symi*4+symj)          /* 16 pairwise non-terminals           */
#define idxL(symi)         (16+symi)              /*  4 leftwise non-terminals           */
#define idxR(symj)         (20+symj)              /*  4 rightwise non-terminals          */
#define idxB               24	                  /*  bifurcation non-terminal           */
#define idxS               25                     /*  start non-terminal                 */
#define idxML(symi)        (26+symi)              /*  4 leftwise non-terminals           */
#define idxMR(symj)        (30+symj)              /*  4 rightwise non-terminals          */
#define idxPL(symi,symj)   (34+symi*4+symj)       /* 16 pairwise hairpin loop            */
#define idxPS(symi,symj)   (50+symi*4+symj)       /* 16 pairwise stem                    */
#define idxPI(symi,symj)   (66+symi*4+symj)       /* 16 pairwise internal loop           */

#define idxX(size)         (82+size)              /*  31 hairpin  nonterminal            */
#define idxV(size)         (113+size)             /*  31 bulge nonterminal               */
#define idxW(size)         (144+size)             /*  31 interior nonterminal            */
#define idxWA(asy)         (175+asy)              /*   4 interior assymetry nonterminal  */
#define idxTL(tlp)         (179+tlp)              /* 256 ultrastable tetraloops          */

#define idxE               4275                   /*  end non-terminal                   */
#define NSTATES            4276

/* Maps from state index to meaningful names or node indices
 * see globals.c for definitions
 */
extern char *stNAME[NSTATES];	/* ASCII state names */
extern int   Ntype[NSTATES];	/* maps states to node types        */


/* Useful macros
 */
#define LN2         0.69314718056
#define LN2INV      1.44269504089
#define LOG2(x)     ((x) == 0.0 ? -HUGE_VAL : log(x) * LN2INV)
#define FLOATSCALE  10.0
#define INTSCALE    10
#define BIGINT      9999999	/* prohibition in viterbi.c without overflow */
#define BIGFLOAT    99999.99	/* prohibition in viterbi.c without overflow */

/* Structure: trace_s
 * 
 * Binary tree structure for storing an RNA structure,
 * derived as a traceback of an alignment of SCFG to sequence.
 */
struct trace_s {
  int emitl;			/* i position (0..N-1) or -1 if nothing   */
  int emitr;			/* j position (0..N-1) or -1 if nothing   */

  int type;                     /* (i,j) state */
    
  struct trace_s *nxtl;		/* ptr to left (or only) branch, or NULL for end    */
  struct trace_s *nxtr;		/* ptr to right branch, for o(mx)^2 only, else NULL */
  struct trace_s *prv;          /* ptr to parent                                    */
};

struct tracekn_s {
  int emiti;			/* i position (0..N-1) or -1 if nothing   */
  int emitj;			/* j position (0..N-1) or -1 if nothing   */
  int emitk;			/* k position (0..N-1) or -1 if nothing   */
  int emitl;			/* l position (0..N-1) or -1 if nothing   */

  int type1;                    /* (i,j) status : dpcP (paired), dcpL (unparied) */
  int type2;                    /* (k,l) status :          dpcS (unknown)        */
    
  struct tracekn_s *nxtl;	/* ptr to left (or only) branch, or NULL for end    */
  struct tracekn_s *nxtr;	/* ptr to right branch, for o(mx)^2 only, else NULL */
  struct tracekn_s *prv;        /* ptr to parent                                    */
};

/* Structure: tracestack_s
 *
 * A pushdown stack used for traversing a binary tree of trace_s structures.
 */
struct tracestack_s {
  struct trace_s      *node;
  struct tracestack_s *nxt;
};

/* Structure: traceknstack_s
 *
 * A pushdown stack used for traversing a binary tree of trace_s structures.
 */
struct traceknstack_s {
  struct tracekn_s      *node;
  struct traceknstack_s *nxt;
};

/* type of matrices */
#define VHX 0
#define VX  1
#define ZHX 2
#define YHX 3
#define WX  4
#define WBX 5
#define WHX 6

  /* Freier/Turner/Tinoco Rules
   * Thermodynamic parameters for RNA structure
   */
#define MAXRNALOOP  30
#define PROHIBIT    50.0
#define MAXPEN      -1 * IntizeScale(3.0)  /* maximun correction for asymmetric internal loop*/
  
#define PRELOG     1.079        /* internal(>32), bulge(>32), hairpin (>31)            */   
#define EPARAM1    0.0          /* extra term for stem                                 */
#define EPARAM2    0.0          /* extra term for bulge                                */
#define EPARAM3    3.0          /* extra term for interior loop                        */
#define EPARAM4    0.0          /* extra term for hairpin loop                         */
#define EPARAM5    4.6          /* extra term for multiloop (IS>2) inside vx           */
#define EPARAM5P   7.0          /* extra term for multiloop (IS>2) inside vhx          */
#define EPARAM6    0.4          /* extra term for dangling (inside vx)                 */
#define EPARAM6P   0.2          /* extra term for dangling (in pseudoknots)            */
#define EPARAM7    MAXRNALOOP   /* max size of an interior loop (l1+l2-2)              */
#define EPARAM8    MAXRNALOOP   /* max loopsidednes of an interior loop |l1-l2|        */
#define EPARAM9   -500*INTSCALE /* bonus energy used to force base pairs               */
#define EPARAM10   0.1          /* extra term for base pair in regular multiloops      */
#define EPARAM10P  0.1          /* extra term for base pair in pseudoknot multiloops   */
#define EPARAM11   7.0          /* penalty for using hole structures under wx          */
#define EPARAM12   6.0          /* penalty for using hole structures under whx         */
#define EPARAM13  16.0          /* penalty for using hole structures under vx and wbx  */

#define wsf 100   /* weight given to no-knot structures */  
#define wkn 83    /* weight given to knots structures */ 

#define lng    4          /* allow hairpinloops of length (j-i) > lng    */
#define lng_kn lng + 0    /* allow        knots of length (l-k) > lng_kn */



struct rnapar_2 {

  int P1;
  int P2;
  int P3;
  int P4;
  int P5;
  int P5P;
  int P6;
  int P6P;
  int P9;
  int P10;
  int P10P;
  int P11;
  int P12;
  int P13;

  float  stack[4][4][4][4];      /* stem          5'-AC...GU-3' is [A][U][C][G]  */
  float tstckh[4][4][4][4];      /* loop terminal 5'-AX...YU-3' is [A][U][X][Y]  */
  float tstcki[4][4][4][4];      /* internal-loop terminal 5'-AX...G  C...YU-3' is [A][U][X][Y] */
                                 /*                        5'-A...XG  CY...U-3' is [C][G][Y][X] */
  float dangle3[4][4][4];        /* dangle3: 5'-A...UX-3' is [X][A][U], 5'-AX...U-3' is [X][U][A] */
  float dangle5[4][4][4];        /* dangle5: 5'-XA...U-3' is [X][A][U], 5'-A...XU-3' is [X][U][A] */

  float hairpin[MAXRNALOOP+1];   /* length contribution to hairpin score */
  float   bulge[MAXRNALOOP+1];   /* length contribution to bulge score   */
  float   inter[MAXRNALOOP+1];   /* total-length contribution to internal-loop score      */
  float             poppen[4];   /* asymmetric-length contribution to internal-loop score */

  char   tloop[30][32];          /* ultrastable tetraloops */
  float  tetraloop[4096];        /* tetraloops' scores     */
};

#define sre_isascii(c)   (!((c) & ~0177))


#endif /* CFG_H_INCLUDED */


















