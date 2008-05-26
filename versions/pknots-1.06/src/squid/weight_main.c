/* SQUID - A C function library for biological sequence analysis
 * Copyright (C) 1992-1996 Sean R. Eddy	
 *
 *    This source code is distributed under terms of the
 *    GNU General Public License. See the files COPYING 
 *    and GNULICENSE for further details.
 *
 */

/* weight_main.c
 * SRE, Thu Mar  3 13:43:39 1994
 * 
 * Calculate weights for a sequence alignment.
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

#include "squid.h"

#ifdef MEMDEBUG
#include "dbmalloc.h"
#endif

static struct opt_s OPTIONS[] = {
  { "-b", TRUE, sqdARG_FLOAT  },
  { "-f", TRUE, sqdARG_FLOAT  }, 
  { "-h", TRUE, sqdARG_NONE   },
  { "-o", TRUE, sqdARG_STRING },
  { "-s", TRUE, sqdARG_INT    }, 
  { "-v", TRUE, sqdARG_NONE   },
};
#define NOPTIONS (sizeof(OPTIONS) / sizeof(struct opt_s))

char usage[] = "\
Usage: weight [-options] <aligned seqfile>\n\
   Calculate sequence weights, compensate for biased representation.\n\
   Available options:\n\
     -b <f>    : use BLOSUM weighting scheme at <f> fractional identity\n\
     -f <f>    : filter out seqs w/ fractional ident > <x> [0-1]\n\
     -h        : help; print version and usage info\n\
     -o <file> : save weight-annotated alignment in <outfile>\n\
     -s <n>    : sample <n> sequences at random into a new alignment\n\
     -v        : use Voronoi weight scheme\n";

int
main(int argc, char **argv)
{
  char  *seqfile;               /* file containing aligned seqs */
  char **aseqs;                 /* aligned sequences (flushed)  */
  AINFO  ainfo;                 /* associated alignment info    */
  int    idx;
  int    fmt;
  
  char  *outfile;               /* output file for weighted alignment */
  FILE  *ofp;                   /* open outfile                       */
  int    do_voronoi;            /* use Sibbald/Argos Voronoi scheme   */
  int    do_blosum;		/* use BLOSUM weighting scheme        */
  int    do_filter;		/* use filtering scheme               */
  float  idlevel;		/* identity level to fiter at, [0-1]  */
  int    samplesize;		/* if >0, don't weight, random sample */

  char *optname;		/* name of option found by Getopt() */
  char *optarg;			/* argument found by Getopt()       */
  int   optind;		        /* index in argv[]                  */

  /***********************************************
   * Parse command line
   ***********************************************/

  outfile    = NULL;
  do_voronoi = FALSE;
  do_filter  = FALSE;
  do_blosum  = FALSE; 
  samplesize = 0;

  while (Getopt(argc, argv, OPTIONS, NOPTIONS, usage,
		&optind, &optname, &optarg))
    {
      
      if     (strcmp(optname, "-b")  == 0)  
	{ do_blosum = TRUE; idlevel = atof(optarg); }
      else if (strcmp(optname, "-f")  == 0)  
	{ do_filter = TRUE; idlevel = atof(optarg); }
      else if (strcmp(optname, "-o")  == 0) outfile    = optarg;
      else if (strcmp(optname, "-s")  == 0) samplesize = atoi(optarg);
      else if (strcmp(optname, "-v")  == 0) do_voronoi = TRUE;
      else if (strcmp(optname, "-h")  == 0)
	{
	  printf("weight %s, %s\n%s\n", squid_version, squid_date, usage);
	  exit(EXIT_SUCCESS);
	}
    }

  if (argc -optind != 1)
    Die("Wrong number of arguments specified on command line\n%s\n", usage);

  seqfile = argv[optind];

  if (outfile == NULL)
    ofp = stdout;
  else if ((ofp = fopen(outfile, "w")) == NULL)
    Die("Failed to open alignment output file %s", outfile);

  if (do_voronoi || samplesize > 0)
    sre_srandom(time(0));

  /***********************************************
   * Read sequences
   ***********************************************/

  if (! SeqfileFormat(seqfile, &fmt, NULL))
    Die("Failed to determine format of file %s\n", seqfile);

  if (! ReadAlignment(seqfile, fmt, &aseqs, &ainfo))
    Die("Failed to read sequences from file %s", seqfile);

  for (idx = 0; idx < ainfo.nseq; idx++)
    s2upper(aseqs[idx]);

  /***********************************************
   * Main routine
   ***********************************************/
  printf("weight %s, %s\n", squid_version, squid_date);

  if (do_filter || samplesize > 0)
    {
      char **anew;
      int    nnew;
      AINFO *newinfo;

      if (do_filter)
	FilterAlignment(aseqs, ainfo.nseq, &ainfo, idlevel,
			&anew, &nnew, &newinfo);
      else if (samplesize > 0)
	SampleAlignment(aseqs, ainfo.nseq, &ainfo, samplesize,
			&anew, &nnew, &newinfo);

      WriteSELEX(ofp, anew, newinfo, 60);      
      FreeAlignment(aseqs, &ainfo);
      FreeAlignment(anew, newinfo);
    }
  else				
    {
      if (do_voronoi)     VoronoiWeights(aseqs, &ainfo);
      else if (do_blosum) BlosumWeights(aseqs, &ainfo, idlevel);
      else                GSCWeights(aseqs, &ainfo);

      WriteSELEX(ofp, aseqs, &ainfo, 60);    
      FreeAlignment(aseqs, &ainfo);
    }


  return EXIT_SUCCESS;
}

