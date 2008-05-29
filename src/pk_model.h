/* (modified) 
 * proto.h
 * ANSI prototypes for all external functions.
 * 
 */

#include <stdio.h>

#include <easel.h>
#include <esl_sqio.h>
#include <esl_wuss.h>

extern float **AllocSCFG(void);
extern int   **AllocIntSCFG(void);
extern void    AllocBaseFreq(float **ret_basefreq);
extern void    BaseFreq(char *s, int j, int d, float **ret_basefreq);
extern float **DupSCFG(float **cfg);
extern void    FreeSCFG(float **cfg);
extern void    FreeIntSCFG(int **icfg);
extern int   **LogifySCFG(float **cfg);
extern int   **LogoddsifySCFG(float **cfg);
extern void    NormalizeSCFG(float **cfg);
extern int   **NussinovIntSCFG(void);
extern void    ProbifySCFG(float **cfg);
extern int     Stype(int node, int symi, int symj, int size, int asy, int tlp);

