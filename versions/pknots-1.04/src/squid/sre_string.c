/* SQUID - A C function library for biological sequence analysis
 * Copyright (C) 1992-1996 Sean R. Eddy	
 *
 *    This source code is distributed under terms of the
 *    GNU General Public License. See the files COPYING 
 *    and GNULICENSE for further details.
 *
 */

/* sre_string.c
 * 
 * my library of extra string functions. Some for portability
 * across UNIXes
 */

#include <stdio.h>
#include <string.h>
#include <stdarg.h>
#include "gnuregex.h"
#include "squid.h"


#ifdef MEMDEBUG
#include "dbmalloc.h"
#endif

/* global sqd_parse[] are managed by Strparse().
 */
char *sqd_parse[10];



/* Obsolete. Anyone who doesn't have strstr() is
 * not ANSI-compliant and must die.
 */
#ifdef NOSTR
char  *strstr(char *s, char *subs)
{
  int  i;

  for ( ; *s != 0; s++) {
    if (*s == *subs) {
      for (i = 1; subs[i] != 0 && subs[i] == s[i]; i++) ;
      if (subs[i] == 0) return(s);
      }
    }
  return (NULL);
}
#endif /* NOSTR */


char *
Strdup(char *s)
{
  char *new;
  if ((new = (char *) malloc (strlen(s) +1)) == NULL) return NULL;
  strcpy(new, s);
  return new;
}

int
Strinsert(char  *s1,            /* string to insert a char into  */
	  char   c,		/* char to insert                */
	  int    pos)		/* position in s1 to insert c at */
{
  char    oldc;
  char   *s;

  for (s = s1 + pos; c; s++)
    {
				/* swap current char for inserted one */
      oldc = *s;		/* pick up current */
      *s   = c;   		/* put down inserted one    */
      c    = oldc;		/* old becomes next to insert */
    }
  *s = '\0';

  return 1;
}


int
Strdelete(char *s1,             /* string to delete a char from       */
	  int   pos)		/* position of char to delete 0..n-1  */
{
  char *s;                      

  for (s = s1 + pos; *s; s++)
    *s = *(s + 1);

  return 1;
}

void
s2lower(char *s)
{
  for (; *s != '\0'; s++)
    *s = sre_tolower((int) *s);
}

void
s2upper(char *s)
{
  for (; *s != '\0'; s++)
    *s = sre_toupper((int) *s);
}


void *
MallocOrDie(size_t size)
{
  void *ptr;

  if ((ptr = malloc (size)) == NULL)
    Die("malloc failed");
  return ptr;
}

void *
ReallocOrDie(void *p, size_t size)
{
  void *ptr;

  if ((ptr = realloc(p, size)) == NULL)
    Die("realloc failed");
  return ptr;
}


/* Function: Strparse()
 * 
 * Purpose:  Match a regexp to a string.
 *           Return 0 if it matches, REG_NOMATCH if it doesn't.
 *
 *           Much like Perl, Strparse() makes copies of the matching
 *           substrings available via globals, sqd_parse[].
 *           sqd_parse[0] contains a copy of the complete matched
 *           text. sqd_parse[1-9] contain copies of up to nine
 *           different substrings matched within parentheses.
 *           The memory for these strings is internally managed and
 *           volatile; the next call to Strparse() may destroy them.
 *           If the caller needs the matched substrings to persist
 *           beyond a new Strparse() call, it must make its own 
 *           copies.
 *           
 *           A minor drawback of the memory management is that
 *           there will be a small amount of unfree'd memory being
 *           managed by Strparse() when a program exits; this may
 *           confuse memory debugging (Purify, dbmalloc).
 * 
 *           Uses an extended POSIX regular expression interface.
 *           A copylefted GNU implementation is included in the squid
 *           implementation (gnuregex.c) for use on non-POSIX compliant
 *           systems. POSIX 1003.2-compliant systems (all UNIX,
 *           some WinNT, I believe) can omit the GNU code if necessary.
 *           
 *           I built this for ease of use, not speed nor efficiency.
 *
 * Example:  Strparse("foo-...-baz", "foo-bar-baz")  returns 0
 *           Strparse("foo-(...)-baz", "foo-bar-baz")
 *              returns 0; sqd_parse[0] is "foo-bar-baz";
 *              sqd_parse[1] is "bar".
 *              
 * Args:     rexp  - regular expression, extended POSIX form
 *           s     - string to match against
 *           ntok  - number of () substrings we will save (maximum 9)
 *                   
 * Return:   0 on match, REG_NOMATCH on failure to match
 */
int
Strparse(char *rexp, char *s, int ntok)
{
  regex_t     pat;
  int         code;
  regmatch_t *pmatch;
  int         len;
  int         i;

				/* sanity check */
  if (ntok > 9)  Die("Strparse(): ntok must be <= 9"); 

  /* Free previous global substring buffers
   */
  for (i = 0; i <= ntok; i++)
    if (sqd_parse[i] != NULL) 
      { 
	free(sqd_parse[i]);
	sqd_parse[i] = NULL;
      }

  /* Compile and match the pattern
   */
  if ((code = regcomp(&pat, rexp, REG_EXTENDED)) != 0)
    Die("POSIX regex compilation failed.");
  pmatch = (regmatch_t *) MallocOrDie (sizeof(regmatch_t) * 10);
  code = regexec(&pat, s, 10, pmatch, 0);

  /* Fill the global substring buffers
   */
  if (code == 0) 
    for (i = 0; i <= ntok; i++)
      {
	len = pmatch[i].rm_eo - pmatch[i].rm_so;
	sqd_parse[i] = (char *) MallocOrDie(sizeof(char) * (len+1));
	strncpy(sqd_parse[i], s+pmatch[i].rm_so, len);
	sqd_parse[i][len] = '\0';
      }

  free(pmatch);
  regfree(&pat);
  return code;
}


/* Function: StrShuffle()
 * 
 * Purpose:  Returns a shuffled version of s2, in s1.
 *  
 * Args:     s1 - allocated space for shuffled string.
 *           s2 - string to shuffle.
 *           
 * Return:   void
 */
void
StrShuffle(char *s1, char *s2)
{
  int  len;
  int  pos;
  char c;
  
  strcpy(s1, s2);
  for (len = strlen(s1); len > 1; len--)
    {				
      pos       = CHOOSE(len);
      c         = s1[pos];
      s1[pos]   = s1[len-1];
      s1[len-1] = c;
    }
}
  
/* Function: RandomSequence()
 * 
 * Purpose:  Generate an iid symbol sequence according
 *           to some alphabet, alphabet_size, probability
 *           distribution, and length. Return the
 *           sequence.
 *           
 * Args:     alphabet  - e.g. "ACGT"
 *           p         - probability distribution [0..n-1]
 *           n         - number of symbols in alphabet
 *           len       - length of generated sequence 
 *           
 * Return:   ptr to random sequence, or NULL on failure.
 */
char *
RandomSequence(char *alphabet, float *p, int n, int len)
{
  char *s;
  int   x;

  s = (char *) MallocOrDie (sizeof(char) * (len+1));
  for (x = 0; x < len; x++)
    s[x] = alphabet[FChoose(p,n)];
  s[x] = '\0';
  return s;
}

