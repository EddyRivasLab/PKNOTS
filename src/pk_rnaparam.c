/* rnaparam.c
 * rna parameters at order 2 in IS,
 * mostly zuker style + pseudoknots (P11, P12, P13, wkn)
 *
 * 
 */
          
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include <easel.h>
#include <esl_alphabet.h>
#include <esl_sq.h>
#include <esl_sqio.h>

#include "pknots.h"
#include "pk_model.h"
#include "pk_rnaparam.h"
#include "pk_trace.h"
#include "pk_util.h"

/* Function: Parameters2_Zkn()
 * 
 * Purpose:  Input of thermodynamic RNA folding parameters.
 *           with Zuker sign convention, adding knots
 *
 *           ranpar  - RETURN: structure containing data.   
 *           
 * Return:   (void)
 *           rnapar is filled.          
 */

void
Parameters2_Zkn(struct rnapar_2 **ret_rnapar)
{
  struct rnapar_2  *rnapar = NULL;
  int               i,j;
  int               k,l;
  int             **iseq = NULL;
  int               tlp, ntlp, isc;
  int               status;
  
  ESL_ALLOC(rnapar, sizeof(struct rnapar_2));

  /* Freier/Turner/Tinoco Rules
   * Thermodynamic parameters for RNA structure
   */

  /* stacking energies */
  double stack03[4][4] = {{+0.4, +0.4, +0.4, -0.9},
			 {+0.4, +0.4, -2.1, +0.4},
			 {+0.4, -1.7, +0.4, -0.5},
			 {-0.9, +0.4, -1.0, +0.4}};

  double stack30[4][4] = {{+0.4, +0.4, +0.4, -1.1},
			 {+0.4, +0.4, -2.3, +0.4},
			 {+0.4, -1.8, +0.4, -0.8},
			 {-0.9, +0.4, -1.1, +0.4}};

  double stack12[4][4] = {{+0.4, +0.4, +0.4, -1.8},
			 {+0.4, +0.4, -2.9, +0.4},
			 {+0.4, -2.0, +0.4, -1.2},
			 {-1.7, +0.4, -1.9, +0.4}};

  double stack21[4][4] = {{+0.4, +0.4, +0.4, -2.3},
			 {+0.4, +0.4, -3.4, +0.4},
			 {+0.4, -2.9, +0.4, -1.4},
			 {-2.1, +0.4, -2.1, +0.4}};

  double stack23[4][4] = {{+0.4, +0.4, +0.4, -1.1},
			 {+0.4, +0.4, -2.1, +0.4},
			 {+0.4, -1.9, +0.4, -0.4},
			 {-1.0, +0.4, +1.5, +0.4}};

  double stack32[4][4] = {{+0.4, +0.4, +0.4, -0.8},
			 {+0.4, +0.4, -1.4, +0.4},
			 {+0.4, -1.2, +0.4, -0.2},
			 {-0.5, +0.4, -0.4, +0.4}};

  /* hairpin terminal energies */

  double tstckh03[4][4] = {{-0.8, -1.0, -1.7, -1.0},
			  {-0.7, -0.7, -0.7, -0.7},
			  {-1.5, -1.0, -1.0, -1.0},
			  {-0.8, -0.8, -0.8, -0.8}};

  double tstckh30[4][4] = {{-1.0, -0.8, -1.8, -0.9},
			  {-0.7, -0.6, -0.3, -0.5},
			  {-1.8, -0.9, -1.2, -0.9},
			  {-0.3, -0.6, -0.3, -0.5}};

  double tstckh12[4][4] = {{-1.4, -2.0, -2.1, -1.9},
			  {-1.0, -1.1, -1.0, -0.8},
			  {-2.1, -1.9, -1.4, -1.9},
			  {-1.4, -1.5, -1.4, -1.2}};

  double tstckh21[4][4] = {{-1.1, -1.3, -2.0, -1.3},
			  {-1.1, -0.6, -0.6, -0.5},
			  {-2.3, -1.5, -1.4, -1.5},
			  {-0.8, -0.8, -0.8, -0.7}};

  double tstckh23[4][4] = {{-0.8, -1.0, -1.7, -1.0},
			  {-0.7, -0.7, -0.7, -0.7},
			  {-1.5, -1.0, -1.0, -1.0},
			  {-0.8, -0.8, -0.8, -0.8}};

  double tstckh32[4][4] = {{-1.2, -1.4, -2.0, -1.4},
			  {-0.9, -0.9, -0.7, -0.7},
			  {-2.0, -1.4, -1.3, -1.4},
			  {-0.9, -1.1, -0.9, -0.9}};


  /* interior loop terminal energies */

  double tstcki03[4][4] = {{-1.0, -1.0, -2.2, -0.5},
			  {-1.0, -1.0, -0.2, -1.0},
			  {-2.2, -0.5, -1.0, -0.5},
			  {-0.3, -1.0, -0.3, -2.0}};
  
  double tstcki30[4][4] = {{-1.0, -1.0, -2.2, -0.4},
			  {-1.0, -1.0, +0.2, -1.0},
			  {-2.2, -0.4, -1.0, -0.4},
			  {+0.2, -1.0, +0.2, -2.0}};
  
  double tstcki12[4][4] = {{-1.5, -1.5, -2.7, -1.9},
			  {-1.5, -1.5, -1.0, -1.5},
			  {-2.7, -1.9, -1.5, -1.9},
			  {-1.4, -1.5, -1.4, -2.5}};
  
  double tstcki21[4][4] = {{-1.5, -1.5, -2.7, -1.3},
			  {-1.5, -1.5, -0.6, -1.5},
			  {-2.7, -1.5, -1.5, -1.5},
			  {-0.8, -1.5, -0.8, -2.5}};
  
  double tstcki23[4][4] = {{-1.5, -1.5, -2.7, -1.3},
			  {-1.5, -1.5, -0.6, -1.5},
			  {-2.7, -1.5, -1.5, -1.5},
			  {-0.8, -1.5, -0.8, -2.5}};
  
  double tstcki32[4][4] = {{-1.0, -1.0, -2.2, -0.4},
			  {-1.0, -1.0, +0.2, -1.0},
			  {-2.2, -0.4, -1.0, -0.4},
			  {+0.2, -1.0, +0.2, -2.0}};
  
  /* dangling energies */
  
  double dangle50[4][4] =  {{BIGFLOAT, BIGFLOAT, BIGFLOAT,     -0.3},
			   {BIGFLOAT, BIGFLOAT,     -0.5, BIGFLOAT},
			   {BIGFLOAT,    -0.2, BIGFLOAT,      -0.2},
			   {    -0.3, BIGFLOAT,    -0.2, BIGFLOAT}};
  
  double dangle51[4][4] = {{BIGFLOAT, BIGFLOAT, BIGFLOAT,     -0.3},
			  {BIGFLOAT, BIGFLOAT,     -0.3, BIGFLOAT},
			  {BIGFLOAT,     -0.3, BIGFLOAT,     -0.2},
			  {    -0.1, BIGFLOAT,     -0.2, BIGFLOAT}};
  
  double dangle52[4][4] = {{BIGFLOAT, BIGFLOAT, BIGFLOAT,     -0.4},
			  {BIGFLOAT, BIGFLOAT,     -0.2, BIGFLOAT},
			  {BIGFLOAT,     +0.0, BIGFLOAT,     -0.2},
			  {    -0.2, BIGFLOAT,     -0.2, BIGFLOAT}};
  
  double dangle53[4][4] =  {{BIGFLOAT, BIGFLOAT, BIGFLOAT,     -0.2},
			   {BIGFLOAT, BIGFLOAT,     -0.1, BIGFLOAT},
			   {BIGFLOAT,     +0.0, BIGFLOAT,     -0.2},
			   {    -0.2, BIGFLOAT,     -0.2, BIGFLOAT}};
  
  
  double dangle30[4][4] =  {{BIGFLOAT, BIGFLOAT, BIGFLOAT,     -0.7},
			   {BIGFLOAT, BIGFLOAT,     -1.1, BIGFLOAT},
			   {BIGFLOAT,     -1.7, BIGFLOAT,     -1.2},
			   {    -0.8, BIGFLOAT,     -0.8, BIGFLOAT}};
  
  double dangle31[4][4] = {{BIGFLOAT, BIGFLOAT, BIGFLOAT,     -0.1},
			  {BIGFLOAT, BIGFLOAT,     -0.4, BIGFLOAT},
			  {BIGFLOAT,     -0.8, BIGFLOAT,     -0.5},
			  {    -0.5, BIGFLOAT,     -0.5, BIGFLOAT}};
  
  double dangle32[4][4] = {{BIGFLOAT, BIGFLOAT, BIGFLOAT,     -0.7},
			  {BIGFLOAT, BIGFLOAT,     -1.3, BIGFLOAT},
			  {BIGFLOAT,     -1.7, BIGFLOAT,     -1.2},
			  {    -0.8, BIGFLOAT,     -0.8, BIGFLOAT}};
  
  double dangle33[4][4] = {{BIGFLOAT, BIGFLOAT, BIGFLOAT,     -0.1},
			  {BIGFLOAT, BIGFLOAT,     -0.6, BIGFLOAT},
			  {BIGFLOAT,     -1.2, BIGFLOAT,     -0.7},
			  {    -0.6, BIGFLOAT,     -0.6, BIGFLOAT}};
  
  /* DATE: Thu Jun  6 10:50:13 CDT 2002
   *
   * ultrastable tetraloops have changes substantially in MFOLD version 3.1
   * 
   * of the original 13, only 9 have survived, and 5 more have been added.
   * the energy for tetraloops now depends on the neighbouring base pair.
   * 
   * At the end, including the closing base pair, only 30 cases get
   * a negative energy
   */
  char tloop[30][6] = {"GGGGAC",
		       "GGUGAC",
		       "CGAAAG",
		       "GGAGAC",
		       "CGCAAG",
		       "GGAAAC",
		       "CGGAAG",
		       "CUUCGG",
		       "CGUGAG",
		       "CGAAGG",
		       "CUACGG",
		       "GGCAAC",
		       "CGCGAG",
		       "UGAGAG",
		       "CGAGAG",
		       "AGAAAU",
		       "CGUAAG",
		       "CUAACG",
		       "UGAAAG",
		       "GGAAGC",
		       "GGGAAC",
		       "UGAAAA",
		       "AGCAAU",
		       "AGUAAU",
		       "CGGGAG",
		       "AGUGAU",
		       "GGCGAC",
		       "GGGAGC",
		       "GUGAAC",
		       "UGGAAA"
  };

  /* loop energies */
  
  rnapar->inter[0]  = BIGFLOAT;
  rnapar->inter[1]  = BIGFLOAT;
  rnapar->inter[2]  = 4.1;
  rnapar->inter[3]  = 5.1;
  rnapar->inter[4]  = 4.9;
  rnapar->inter[5]  = 5.3;
  rnapar->inter[6]  = 5.7;
  rnapar->inter[7]  = 5.9;
  rnapar->inter[8]  = 6.0;
  rnapar->inter[9]  = 6.1;
  rnapar->inter[10] = 6.3;
  rnapar->inter[11] = 6.4;
  rnapar->inter[12] = 6.4;
  rnapar->inter[13] = 6.5;
  rnapar->inter[14] = 6.6;
  rnapar->inter[15] = 6.7;
  rnapar->inter[16] = 6.8;
  rnapar->inter[17] = 6.8;
  rnapar->inter[18] = 6.9;
  rnapar->inter[19] = 6.9;
  rnapar->inter[20] = 7.0;
  rnapar->inter[21] = 7.1;
  rnapar->inter[22] = 7.1;
  rnapar->inter[23] = 7.1;
  rnapar->inter[24] = 7.2;
  rnapar->inter[25] = 7.2;
  rnapar->inter[26] = 7.3;
  rnapar->inter[27] = 7.3;
  rnapar->inter[28] = 7.4;
  rnapar->inter[29] = 7.4;
  rnapar->inter[30] = 7.4;

  rnapar->bulge[0]  = BIGFLOAT;
  rnapar->bulge[1]  = 3.9;
  rnapar->bulge[2]  = 3.1;
  rnapar->bulge[3]  = 3.5;
  rnapar->bulge[4]  = 4.2;
  rnapar->bulge[5]  = 4.8;
  rnapar->bulge[6]  = 5.0;
  rnapar->bulge[7]  = 5.2;
  rnapar->bulge[8]  = 5.3;
  rnapar->bulge[9]  = 5.4;
  rnapar->bulge[10] = 5.5;
  rnapar->bulge[11] = 5.7;
  rnapar->bulge[12] = 5.7;
  rnapar->bulge[13] = 5.8;
  rnapar->bulge[14] = 5.9;
  rnapar->bulge[15] = 6.0;
  rnapar->bulge[16] = 6.1;
  rnapar->bulge[17] = 6.1;
  rnapar->bulge[18] = 6.2;
  rnapar->bulge[19] = 6.2;
  rnapar->bulge[20] = 6.3;
  rnapar->bulge[21] = 6.3;
  rnapar->bulge[22] = 6.4;
  rnapar->bulge[23] = 6.4;
  rnapar->bulge[24] = 6.5;
  rnapar->bulge[25] = 6.5;
  rnapar->bulge[26] = 6.5;
  rnapar->bulge[27] = 6.6;
  rnapar->bulge[28] = 6.7;
  rnapar->bulge[29] = 6.7;
  rnapar->bulge[30] = 6.7;

  rnapar->hairpin[0]  = BIGFLOAT;
  rnapar->hairpin[1]  = BIGFLOAT;
  rnapar->hairpin[2]  = BIGFLOAT;
  rnapar->hairpin[3]  = 4.1;
  rnapar->hairpin[4]  = 4.9;
  rnapar->hairpin[5]  = 4.4;
  rnapar->hairpin[6]  = 4.7;
  rnapar->hairpin[7]  = 5.0;
  rnapar->hairpin[8]  = 5.1;
  rnapar->hairpin[9]  = 5.2;
  rnapar->hairpin[10] = 5.3;
  rnapar->hairpin[11] = 5.4;
  rnapar->hairpin[12] = 5.5;
  rnapar->hairpin[13] = 5.6;
  rnapar->hairpin[14] = 5.7;
  rnapar->hairpin[15] = 5.8;
  rnapar->hairpin[16] = 5.8;
  rnapar->hairpin[17] = 5.9;
  rnapar->hairpin[18] = 5.9;
  rnapar->hairpin[19] = 6.0;
  rnapar->hairpin[20] = 6.1;
  rnapar->hairpin[21] = 6.1;
  rnapar->hairpin[22] = 6.2;
  rnapar->hairpin[23] = 6.2;
  rnapar->hairpin[24] = 6.3;
  rnapar->hairpin[25] = 6.3;
  rnapar->hairpin[26] = 6.3;
  rnapar->hairpin[27] = 6.4;
  rnapar->hairpin[28] = 6.4;
  rnapar->hairpin[29] = 6.5;
  rnapar->hairpin[30] = 6.5;

  /* internal loops asymmetry factor */

  rnapar->poppen[0] = 0.3;
  rnapar->poppen[1] = 0.3;
  rnapar->poppen[2] = 0.3;
  rnapar->poppen[3] = 0.3;


  /* Prohibit any  danglings and stackings to start with */
  for (i = 0; i < 4; i++)
    for (j = 0; j < 4; j++)
      for (k = 0; k < 4; k++)
	{ 
	 rnapar->dangle5[i][j][k] = BIGFLOAT;
	 rnapar->dangle3[i][j][k] = BIGFLOAT;

	 for (l = 0; l < 4; l++){
	   rnapar->stack[i][j][k][l]  = BIGFLOAT;
	   rnapar->tstckh[i][j][k][l] = BIGFLOAT;
	   rnapar->tstcki[i][j][k][l] = BIGFLOAT;}
	}

  /* introduce the allowed danglings and stackings */
  for (i = 0; i < 4; i++)
    for (j = 0; j < 4; j++)
      {
	rnapar->stack[0][3][i][j] = stack03[i][j];
	rnapar->stack[3][0][i][j] = stack30[i][j];
	rnapar->stack[1][2][i][j] = stack12[i][j];
	rnapar->stack[2][1][i][j] = stack21[i][j];
	rnapar->stack[2][3][i][j] = stack23[i][j];
	rnapar->stack[3][2][i][j] = stack32[i][j];

	rnapar->tstckh[0][3][i][j] = tstckh03[i][j];
	rnapar->tstckh[3][0][i][j] = tstckh30[i][j];
	rnapar->tstckh[1][2][i][j] = tstckh12[i][j];
	rnapar->tstckh[2][1][i][j] = tstckh21[i][j];
	rnapar->tstckh[2][3][i][j] = tstckh23[i][j];
	rnapar->tstckh[3][2][i][j] = tstckh32[i][j];

	rnapar->tstcki[0][3][i][j] = tstcki03[i][j];
	rnapar->tstcki[3][0][i][j] = tstcki30[i][j];
	rnapar->tstcki[1][2][i][j] = tstcki12[i][j];
	rnapar->tstcki[2][1][i][j] = tstcki21[i][j];
	rnapar->tstcki[2][3][i][j] = tstcki23[i][j];
	rnapar->tstcki[3][2][i][j] = tstcki32[i][j];

	rnapar->dangle5[0][i][j] = dangle50[i][j];
	rnapar->dangle5[1][i][j] = dangle51[i][j];
	rnapar->dangle5[2][i][j] = dangle52[i][j];
	rnapar->dangle5[3][i][j] = dangle53[i][j];

	rnapar->dangle3[0][i][j] = dangle30[i][j];
	rnapar->dangle3[1][i][j] = dangle31[i][j];
	rnapar->dangle3[2][i][j] = dangle32[i][j];
	rnapar->dangle3[3][i][j] = dangle33[i][j];
      }

  /* ultrastable tetrallops.
   * initialize all  tetraloop scores to zero 
   */
   for (tlp = 0; tlp < 4096; tlp++)
    rnapar->tetraloop[tlp]  = 0.0;

   ESL_ALLOC(iseq, sizeof(int *) * 30);
   for (ntlp = 0; ntlp < 30; ntlp++)
     ESL_ALLOC(iseq[ntlp], sizeof(int) * 6);
   
   /* assign the sequences of the ultrastable tetraloops */
   for (ntlp = 0; ntlp < 30; ntlp++)
     for (i = 0; i < 6; i++)
       rnapar->tloop[ntlp][i] = tloop[ntlp][i];
   
                 /* assign the scores of the ultrastable tetraloops */
  for (ntlp = 0; ntlp < 30; ntlp++)
    {
      /* convert ultrastable sequences to 0..3 */
      IntizeSequence(rnapar->tloop[ntlp], 6, iseq[ntlp]);
 
      isc = iseq[ntlp][0] + iseq[ntlp][1]*4 
	+ iseq[ntlp][2]*4*4 + iseq[ntlp][3]*4*4*4 
	+ iseq[ntlp][4]*4*4*4*4 + iseq[ntlp][5]*4*4*4*4*4;		

      if      (ntlp <  9) rnapar->tetraloop[isc]  = -3.0;	
      else if (ntlp < 14) rnapar->tetraloop[isc]  = -2.5;	
      else if (ntlp < 19) rnapar->tetraloop[isc]  = -2.0;	
      else if (ntlp < 30) rnapar->tetraloop[isc]  = -1.5;	
    }

   for (ntlp = 0; ntlp < 30; ntlp++)
     if (iseq[ntlp] != NULL) free(iseq[ntlp]);
  if (iseq != NULL) free(iseq);

  *ret_rnapar = rnapar;
  return;

 ERROR:
  printf("bad allocation\n");
  exit(1);
}

/* Function: ParamIntSCFG()
 * 

 * Purpose:  Converts rna parameters at second order to
 *           integers, and changes sign.
 *     
 *           
 * Args:     rnapar - rna parameters
 *           cfg   - model (already allocated)                  
 */
int **
ParamIntSCFG(struct rnapar_2 *rnapar)
{
  int **icfg;
  int i, j, ip, jp;
  int size;
  int tlp;
  int asy; 

  rnapar->P1   = -wsf*IntizeScale(EPARAM1);
  rnapar->P2   = -wsf*IntizeScale(EPARAM2);
  rnapar->P3   = -wsf*IntizeScale(EPARAM3);
  rnapar->P4   = -wsf*IntizeScale(EPARAM4);
  rnapar->P5   = -wsf*IntizeScale(EPARAM5);
  rnapar->P5P  = -wsf*IntizeScale(EPARAM5P);
  rnapar->P6   = -wsf*IntizeScale(EPARAM6);
  rnapar->P6P  = -wsf*IntizeScale(EPARAM6P);
  rnapar->P9   = -wsf*EPARAM9*INTSCALE;
  rnapar->P10  = -wsf*IntizeScale(EPARAM10);
  rnapar->P10P = -wkn*IntizeScale(EPARAM10P);
  rnapar->P11  = -wsf*IntizeScale(EPARAM11);
  rnapar->P12  = -wsf*IntizeScale(EPARAM12);
  rnapar->P13  = -wsf*IntizeScale(EPARAM13);

  icfg = AllocIntSCFG();

  /* Set everything flat to zero.
   */
  for (i = 0; i < NSTATES; i++)
    for (j = 0; j < NSTATES; j++)
      icfg[i][j] = 0;

  /* Use negative stacking energies as
   * scores for P-> transitions for 16 pair states
   */
  for (i = 0; i < 4; i++)
    for (j = 0; j < 4; j++)	 
      {
       if ((i+j == 3) || (i+j == 5)) /* determines base pairs (0,3),(3,0),(1,2),(2,1),(2,3),(3,2) */
	icfg[idxS][idxP(i,j)] = 1 * INTSCALE;   /* Watson-Crick */
       else
	icfg[idxS][idxP(i,j)] = -1 * BIGINT; /* PROHIBIT */
       
       for (ip = 0; ip < 4; ip++)
	 {
           icfg[idxL(ip)][idxP(i,j)] = -1 * IntizeScale(rnapar->dangle5[ip][i][j]);
	   icfg[idxR(ip)][idxP(i,j)] = -1 * IntizeScale(rnapar->dangle3[ip][i][j]);

	   for (jp = 0; jp < 4; jp++)
	     {
	       icfg[idxPS(i,j)][idxPS(ip,jp)]   = 
		                                -1 * IntizeScale(rnapar->stack[i][j][ip][jp]);
	       icfg[idxPL(i,j)][idxPL(ip,jp)] = 
                                                -1 * IntizeScale(rnapar->tstckh[i][j][ip][jp]);
	       icfg[idxPI(i,j)][idxPI(ip,jp)] = 
                                                -1 * IntizeScale(rnapar->tstcki[i][j][ip][jp]);
	     }
	 }
      }

  for (size = 0; size <= MAXRNALOOP; size++)
    {
      icfg[idxS][idxX(size)] =  -1 * IntizeScale(rnapar->hairpin[size]);
      icfg[idxS][idxV(size)] =  -1 * IntizeScale(rnapar->bulge[size]);
      icfg[idxS][idxW(size)] =  -1 * IntizeScale(rnapar->inter[size]);
    }
	
  for (asy = 0; asy < 4; asy++)
    icfg[idxS][idxWA(asy)] =  -1 * IntizeScale(rnapar->poppen[asy]);

  for (tlp = 0; tlp < 4096; tlp++)
    icfg[idxS][idxTL(tlp)] =  -1 * IntizeScale(rnapar->tetraloop[tlp]);

  return icfg;
}



int
IntizeScale (double val) {

  int    ival;
  double scaleval;
  double precs;
  int    verbose = FALSE;

  scaleval = val * FLOATSCALE;

  ival = (int)(scaleval);

  if (fabs(scaleval) > 0.0) {
    precs = scaleval/ival;
    
    if (precs > 1.01) ival += (val > 0.0)? (int)precs : -(int)precs;
  }
  
  if (verbose) printf("VAL %f %d\n", val, ival);

  return ival;
}
























