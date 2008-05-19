# Makefile for RNAfold_knots: Development version
#
#

VERSION = 1.0
VDATE   = FEB 1999

## where you want things installed
BINDIR  = $(HOME)/bin
MANDIR  = $(HOME)/man

## where are the squid library, include, and source files?
SQUIDDIR    = squid
SQUIDSRCDIR = squid

## your compiler
#CC = gcc
CC = cc    # for SGI Origin200 compiler#

## any special compiler flags you want
#  -pedantic clashes with -DMEMDEBUG??
#
#CFLAGS = -g -Wall -pedantic
CFLAGS = -O
#CFLAGS = -g
#CFLAGS = -TARG:exc_max -64 -Ofast=ip27   #specific SGI Origin200  compiler#

## machine specific definitions
#  None should be needed
#
# Development:
# -DASSERT        Diagnostics level 1. Assertions. 
# -DVERBOSE       Diagnostics level 2. Verbosity.  (unused at the moment)
# -DMEMDEBUG      Diagnostics level 3. Memory checking. (maybe slow!)
#                                      (must link w/ -ldbmalloc) 
# -DLINTING       Suppress overzealous warnings from gcc -Wall and lint   
#
#MDEFS = -DLINTING -DASSERT -DMEMDEBUG
#MDEFS = -DASSERT 
 
## how to install the man pages 
## cp is fine
INSTMAN   = cp
MANSUFFIX = 1

## Your compression program
COMPRESS = gzip               # GNU compression (better, smaller)
#COMPRESS = compress          # generic UNIX compress (portable)

#######
## should not need to modify below this line
#######
SHELL  = /bin/sh
LIBS   = -lsquid -lm
#LIBS = -ldbmalloc -lsquid -lm 

READMES = README 
#INSTALL COPYING DEMO GNULICENSE FILES

PROGS = pknots 

MANSRC= 

DOCS =  

HDRS =  cfg.h proto.h protovhx.h protovx.h protowhx.h protowbx.h protowx.h protozhxyhx.h 

MAIN =  knotf2_main.c\

SQUIDSRC = alignio.c sqerror.c sqio.c iupac.c msf.c revcomp.c\
	   selex.c sre_ctype.c sre_string.c types.c

SQUIDHDRS = squid.h sqfuncs.h

SRC =   cfgio.c  fillmtx.c\
	filltrvhx.c filltrvx.c filltrwbx.c filltrwhx.c filltrwx.c filltrzhxyhx.c \
	globals.c irredsurf.c misc.c model.c rnaio.c rnaoutput.c \
	rnaparam.c trace.c tracemtx.c viterbi.c \
        vhxgraphs.c vxgraphs.c wbxgraphs.c whxgraphs.c wxgraphs.c zhxyhxgraphs.c


OBJ =   cfgio.o  fillmtx.o\
	filltrvhx.o filltrvx.o filltrwbx.o filltrwhx.o filltrwx.o filltrzhxyhx.o \
	globals.o irredsurf.o misc.o model.o rnaio.o rnaoutput.o \
	rnaparam.o trace.o tracemtx.o viterbi.o \
        vhxgraphs.o vxgraphs.o wbxgraphs.o whxgraphs.o wxgraphs.o zhxyhxgraphs.o

DISTFILES = $(MAIN) $(SRC) $(HDRS) $(MANSRC) $(READMES) $(DOCS)

TESTSUITE = 

DEMOS = 

all: 	$(PROGS)

pknots:  $(OBJ) knotf2_main.o
	 $(CC) $(CFLAGS) -o $@ -L$(SQUIDDIR) knotf2_main.o $(OBJ) $(LIBS)

install: $(PROGS)
	cp $(PROGS) $(BINDIR)/

clean:
	-rm -f *.o *~ *.Addrs *.Counts *.pixie Makefile.bak core $(PROGS) TAGS

tags:
	etags -t $(DISTFILES)

dist:
	@if test -d knot-$(VERSION);        then rm -rf knot-$(VERSION);    fi
	@if test -f knot-$(VERSION).tar.gz; then rm knot-$(VERSION).tar.gz; fi
	@if test -f knot-$(VERSION).tar.Z;  then rm knot-$(VERSION).tar.Z;  fi
	mkdir yarn-$(VERSION)/
	cp $(DISTFILES) knot-$(VERSION)/
	for squidfile in $(SQUIDSRC); do \
	  cp $(SQUIDSRCDIR)/$$squidfile knot-$(VERSION)/;\
	done
	for squidfile in $(SQUIDHDRS); do \
	  cp $(SQUIDSRCDIR)/$$squidfile knot-$(VERSION)/;\
	done
	cp Makefile.dist knot-$(VERSION)/Makefile
	mkdir knot-$(VERSION)/Testsuite
	for testfile in $(TESTSUITE); do \
	  cp Testsuite/$$testfile knot-$(VERSION)/Testsuite/;\
	done
	mkdir knot-$(VERSION)/Demos
	for demofile in $(DEMOS); do\
	  cp Demos/$$demofile knot-$(VERSION)/Demos/;\
	done
	tar cvf knot-$(VERSION).tar knot-$(VERSION)
	$(COMPRESS) knot-$(VERSION).tar

test:   
	(cd Testsuite; testsuite.sh)

lint:
	lint -h -I$(SQUIDDIR) $(RFLAGS) $(MAIN) $(SRC) $(SQUIDDIR)/llib-lsquid.ln /usr/lib/lint/llib-lm.ln

.c.o:
	$(CC) $(CFLAGS) $(MDEFS) -I$(SQUIDDIR) -c $<		

.m.o:
	$(CC) $(CFLAGS) $(MDEFS) -I$(SQUIDDIR) -c $<







