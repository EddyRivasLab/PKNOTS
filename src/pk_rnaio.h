#include <stdio.h>
#include "cfg.h"

#include <easel.h>
#include <esl_sqio.h>
#include <esl_wuss.h>

extern int KHS2ct(char *ss, int len, int allow_pseudoknots, int **ret_ct);
extern int VerifyKHS(char *name, char *ss, int wordy);

