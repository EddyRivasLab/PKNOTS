/* (modified)
 * trace.c  -  support for traceback tree structure 
 *     (RNA structure representation)
 *
 * cove 1.0: Mon May 17 09:38:14 1993
 * moved to cove 2.0, Mon Sep  6 13:34:55 1993
 * modified for yarn, Sun Aug 28 10:01:42 1994
 * 
 * The traceback of an SCFG-sequence alignment is an RNA structure,
 * which is kept as a tree of trace_s structures. Here
 * we provide support for the traceback data structures: the
 * tree itself, and a pushdown stack used for traversing the
 * tree.
 * 
 * The trace tree structure has a dummy node at its beginning
 * and dpcE END states at the leaves. Unlike COVE, these ends are
 * created by the functions that create traces, not maintained
 * automatically. Non-BIFURC states have a NULL right branch. 
 * 
 * The pushdown stack structure has a dummy begin node, and the
 * end is signified by a final NULL ptr.
 */


#include <stdio.h>
#include <string.h>
#include <stdlib.h>

#include <easel.h>
#include <esl_sqio.h>

#include "pknots.h"		/* struct trace_s and struct tracestack_s */
#include "pk_trace.h"
#include "pk_util.h"

#ifdef MEMDEBUG
#include "dbmalloc.h"
#endif
#ifdef ASSERT
#include <assert.h>
#endif

/* Function: InitTrace()
 * 
 * Purpose:  Initialize a traceback tree structure.
 *
 * Return:   ptr to the new tree.
 */          
struct trace_s *
InitTrace(void)
{
  struct trace_s *new;

  if ((new = (struct trace_s *) malloc (sizeof(struct trace_s))) == NULL)
    pk_fatal("Memory allocation failure at %s line %d", __FILE__, __LINE__);
  
  new->emitl = -1;
  new->emitr = -1;
  new->type  = -1;
  new->nxtl  = NULL;
  new->nxtr  = NULL;
  new->prv   = NULL;
  return new;
}

/* Function: InitTracekn()
 * 
 * Purpose:  Initialize a traceback tree structure.
 *
 * Return:   ptr to the new tree.
 */          
struct tracekn_s *
InitTracekn(void)
{
  struct tracekn_s *new;

  if ((new = (struct tracekn_s *) malloc (sizeof(struct tracekn_s))) == NULL)
    pk_fatal("Memory allocation failure at %s line %d", __FILE__, __LINE__);
  
  new->emiti = -1;
  new->emitj = -1;
  new->emitk = -1;
  new->emitl = -1;
  new->type1 = -1;
  new->type2 = -1;
  new->nxtl  = NULL;
  new->nxtr  = NULL;
  new->prv   = NULL;
  return new;
}

/* Function: AttachTrace()
 * 
 * Purpose:  attach a new node to a tracetree node.
 *           There are no dummy END nodes, in contrast to COVE.
 *           
 *           Because of the mechanics of tracebacks through a Viterbi matrix,
 *           we have to make sure that BIFURC children are attached
 *           right first, left second.
 *           
 * Returns:  ptr to the new node, or NULL on failure.
 */          
struct trace_s *
AttachTrace(struct trace_s *parent,
	        int         emitl,
	        int         emitr,
                int         type)
{
  struct trace_s *new;

  if (parent->nxtr != NULL)
    pk_fatal("That trace node is already full, fool.");

  /* If left branch is already connected to something, swap it over to the
   * right (thus enforcing the necessary rule that BIFURCS attach to the right
   * branch first), and attach a new dummy end to the left branch. 
   */
  if (parent->nxtl != NULL)
    {
      parent->nxtr = parent->nxtl;
      parent->nxtl = NULL;
    }

  if ((new = (struct trace_s *) malloc (sizeof(struct trace_s))) == NULL)
    pk_fatal("Memory allocation failure at %s line %d", __FILE__, __LINE__);
  new->nxtr    = NULL;
  new->nxtl    = NULL;
  new->prv     = parent;
  new->emitl   = emitl;
  new->emitr   = emitr;
  new->type    = type;
  parent->nxtl = new;

  return new;
}

/* Function: AttachTracekn()
 * 
 * Purpose:  attach a new node to a tracetree node.
 *           There are no dummy END nodes, in contrast to COVE.
 *           
 *           Because of the mechanics of tracebacks through a Viterbi matrix,
 *           we have to make sure that BIFURC children are attached
 *           right first, left second.
 *           
 * Returns:  ptr to the new node, or NULL on failure.
 */          
struct tracekn_s *
AttachTracekn(struct tracekn_s *parent,
	      int                emiti,
	      int                emitj,
	      int                emitk,
	      int                emitl,
	      int                type1,
	      int                type2)
{
  struct tracekn_s *new;

  if (parent->nxtr != NULL)
    pk_fatal("That trace node is already full, fool.");

  /* If left branch is already connected to something, swap it over to the
   * right (thus enforcing the necessary rule that BIFURCS attach to the right
   * branch first), and attach a new dummy end to the left branch. 
   */
  if (parent->nxtl != NULL)
    {
      parent->nxtr = parent->nxtl;
      parent->nxtl = NULL;
    }

  if ((new = (struct tracekn_s *) malloc (sizeof(struct tracekn_s))) == NULL)
    pk_fatal("Memory allocation failure at %s line %d", __FILE__, __LINE__);
  new->nxtr    = NULL;
  new->nxtl    = NULL;
  new->prv     = parent;
  new->emiti   = emiti;
  new->emitj   = emitj;
  new->emitk   = emitk;
  new->emitl   = emitl;
  new->type1   = type1;
  new->type2   = type2;
  parent->nxtl = new;

  return new;
}


void
FreeTrace(struct trace_s *tr)
{
  struct tracestack_s *stack;
  struct trace_s      *currtr;

  stack = InitTracestack();
  PushTracestack(stack, tr);

  while ((currtr = PopTracestack(stack)) != NULL)
    {
      if (currtr->nxtr != NULL)
	PushTracestack(stack, currtr->nxtr);
      if (currtr->nxtl != NULL)
	PushTracestack(stack, currtr->nxtl);
      free(currtr);
    }
  FreeTracestack(stack);
}
 void
FreeTracekn(struct tracekn_s *tr)
{
  struct traceknstack_s *stack;
  struct tracekn_s      *currtr;

  stack = InitTraceknstack();
  PushTraceknstack(stack, tr);

  while ((currtr = PopTraceknstack(stack)) != NULL)
    {
      if (currtr->nxtr != NULL)
	PushTraceknstack(stack, currtr->nxtr);
      if (currtr->nxtl != NULL)
	PushTraceknstack(stack, currtr->nxtl);
      free(currtr);
    }
  FreeTraceknstack(stack);
}

void
DeleteTraceknnode(struct tracekn_s *oldtr)
{
  struct tracekn_s *parent;

  parent = oldtr->prv;

  parent->nxtl = oldtr->nxtl;
  parent->nxtr = oldtr->nxtr;
  oldtr->nxtl->prv = parent;
  if (oldtr->nxtr) oldtr->nxtr->prv = parent;
  free(oldtr);
}

/* Function : InitTracestack()
 *            
 * Purpose:   Implementation of the pushdown stack for
 *            traversing traceback trees.
 */            
struct tracestack_s *
InitTracestack(void)
{
  struct tracestack_s *stack;

  if ((stack = (struct tracestack_s *) malloc (sizeof(struct tracestack_s))) == NULL)
    pk_fatal("Memory allocation failure at %s line %d", __FILE__, __LINE__);
  stack->nxt = NULL;
  return stack;
}

/* Functions: InitTraceknstack()
 *            
 * Purpose:   Implementation of the pushdown stack for
 *            traversing traceback trees.
 */            
struct traceknstack_s *
InitTraceknstack(void)
{
  struct traceknstack_s *stack;

  if ((stack = (struct traceknstack_s *) malloc (sizeof(struct traceknstack_s))) == NULL)
    pk_fatal("Memory allocation failure at %s line %d", __FILE__, __LINE__);
  stack->nxt = NULL;
  return stack;
}

void 
PushTracestack(struct tracestack_s *stack,
	         struct trace_s      *tracenode)
{
  struct tracestack_s *new;

  if ((new = (struct tracestack_s *) malloc (sizeof(struct tracestack_s))) == NULL)
    pk_fatal("Memory allocation failure at %s line %d", __FILE__, __LINE__);
  new->node = tracenode;

  new->nxt = stack->nxt;
  stack->nxt = new;
}

void 
PushTraceknstack(struct traceknstack_s *stack,
	         struct tracekn_s      *tracenode)
{
  struct traceknstack_s *new;

  if ((new = (struct traceknstack_s *) malloc (sizeof(struct traceknstack_s))) == NULL)
    pk_fatal("Memory allocation failure at %s line %d", __FILE__, __LINE__);
  new->node = tracenode;

  new->nxt = stack->nxt;
  stack->nxt = new;
}

struct trace_s *
PopTracestack(struct tracestack_s *stack)
{
  struct trace_s *node;
  struct tracestack_s *old;

  if (stack->nxt == NULL)
    return NULL;

  old = stack->nxt;
  stack->nxt = old->nxt;

  node = old->node;
  free(old);
  return node;
}

struct tracekn_s *
PopTraceknstack(struct traceknstack_s *stack)
{
  struct tracekn_s *node;
  struct traceknstack_s *old;

  if (stack->nxt == NULL)
    return NULL;

  old = stack->nxt;
  stack->nxt = old->nxt;

  node = old->node;
  free(old);
  return node;
}

void
FreeTracestack(struct tracestack_s *stack)
{
  while (PopTracestack(stack) != NULL)
    ;
  free(stack);
}

void
FreeTraceknstack(struct traceknstack_s *stack)
{
  while (PopTraceknstack(stack) != NULL)
    ;
  free(stack);
}




















