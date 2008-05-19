# squid - a library of functions for biological sequence analysis
# Makefile for distribution version 
# SRE, Fri Oct 27 09:27:35 1995

## Where you want things installed
#
#  I assume you're installing this for yourself.
#  For a system-wide installation, change $(HOME) to something like /usr/local
#
BINDIR     = $(HOME)/bin
SCRIPTDIR  = $(HOME)/scripts
MANDIR     = $(HOME)/man

## Your compiler
#
#  Use an ANSI compliant C compiler.
#  Usually this will be cc (your vendor's compiler) or gcc (Gnu C).
#  On pre-Solaris SunOS, you must use gcc. SunOS cc is not ANSI-compliant.
#  On an Alliant, use fxc.
#
CC = cc
#CC = gcc
#CC = fxc

## Any special compiler flags you want.
#
#  Note on optimization: The -O flag is not likely to be the best
#  you can do; but it's portable. 
#  If you feel the need for speed, try:
#     gcc (any platform) : -O3
#     SGI IRIX 5.x:        -O2 
#     DEC Alpha OSF/1:     -O4 -migrate
#     Linux (cc == gcc):   -O3
#     Sun Solaris 5.x:     ?? 
#
#  Some optimizers may be buggy. I've had trouble with Solaris,
#  for instance. If the code builds cleanly but seg faults, try
#  rebuilding without optimization.
#
CFLAGS = -O
#CFLAGS = -g 
#CFLAGS = -TARG:exc_max -64 -Ofast=ip27    #specific SGI Origin200 compiler#

## The archiver command and flags (probably correct for any system)
#
#  On Silicon Graphics (or other SYSV systems) ranlib may not exist.
#  Don't worry. The Makefile tests for ranlib existence and will work
#  fine on SGIs even if RANLIB is set.
#
AR     = ar rcv
RANLIB = ranlib

## Machine specific definitions
#
#  You shouldn't need any. The specific #define's in squid are historical.
#  However, if you do run into compilation trouble, you might need:
# -DNO_STRDUP    no strdup() function (lots of BSDish Unix)
# -DNOSTR        no strstr() function (rarer)
# 
MDEFS =


## Instructions for installing man pages
#
#  since there aren't any man pages to speak of yet,
#  don't worry about this.
#
INSTMAN  = cp
MANSUFFIX = 1


#######
## should not need to modify below this line
#######
SHELL     = /bin/sh
LIBTARGET = libsquid.a
LIBS      = -lm

PROGS = alistat compalign compstruct getseq reformat revcomp seqstat\
	shuffle translate weight

SCRIPTS = fasta2gsi.pl genbank2gsi.pl genpept2gsi.pl pir2gsi.pl swiss2gsi.pl

HDRS =  rk.h sqfuncs.h squid.h gnuregex.h

READMES = COPYING FILES GNULICENSE INSTALL README

MANS =  dbio io	selexio

SRC =   aligneval.c alignio.c cluster.c dayhoff.c file.c getopt.c\
	interleaved.c sqerror.c\
	gnuregex.c sqio.c iupac.c msf.c revcomp.c rk.c selex.c\
	seqencode.c sre_ctype.c sre_math.c\
	sre_string.c stack.c translate.c types.c weight.c \
	compstruct_main.c getseq_main.c reformat_main.c revcomp_main.c\
	seqstat_main.c shuffle_main.c weight_main.c


OBJ =   aligneval.o alignio.o cluster.o dayhoff.o file.o getopt.o\
	interleaved.o sqerror.o\
	gnuregex.o sqio.o iupac.o msf.o revcomp.o rk.o selex.o\
	seqencode.o sre_ctype.o sre_math.o\
	sre_string.o stack.o translate.o types.o weight.o

DIST = $(SRC) $(HDRS) $(READMES) Makefile

all: $(LIBTARGET)  

alistat: $(OBJ) alistat_main.o
	$(CC) $(CFLAGS) $(MDEFS) -o alistat $(OBJ) alistat_main.o $(LIBS)

compalign: $(OBJ) compalign_main.o
	$(CC) $(CFLAGS) $(MDEFS) -o compalign $(OBJ) compalign_main.o $(LIBS)

compstruct: $(OBJ) compstruct_main.o
	$(CC) $(CFLAGS) $(MDEFS) -o compstruct $(OBJ) compstruct_main.o $(LIBS)

getseq: $(OBJ) getseq_main.o
	$(CC) $(CFLAGS) $(MDEFS) -o getseq $(OBJ) getseq_main.o $(LIBS)

reformat: $(OBJ) reformat_main.o
	$(CC) $(CFLAGS) $(MDEFS) -o reformat $(OBJ) reformat_main.o $(LIBS)

revcomp: $(OBJ) revcomp_main.o
	$(CC) $(CFLAGS) $(MDEFS) -o revcomp $(OBJ) revcomp_main.o $(LIBS)

seqstat: $(OBJ) seqstat_main.o
	$(CC) $(CFLAGS) $(MDEFS) -o seqstat $(OBJ) seqstat_main.o $(LIBS)

shuffle: $(OBJ) shuffle_main.o
	$(CC) $(CFLAGS) $(MDEFS) -o shuffle $(OBJ) shuffle_main.o $(LIBS)

translate:  $(OBJ) translate_main.o
	$(CC) $(CFLAGS) $(MDEFS) -o translate $(OBJ) translate_main.o $(LIBS)

weight: $(OBJ) weight_main.o
	$(CC) $(CFLAGS) $(MDEFS) -o weight $(OBJ) weight_main.o $(LIBS)


$(LIBTARGET): $(OBJ)
	$(AR) $(LIBTARGET) $(OBJ)
	-if test -x /bin/ranlib; then /bin/ranlib $(LIBTARGET); else exit 0; fi
	-if test -x /usr/bin/ranlib; then /usr/bin/ranlib $(LIBTARGET); else exit 0; fi
	-if test -x $(RANLIB); then $(RANLIB) $(LIBTARGET); else exit 0; fi
	-chmod 644 $(LIBTARGET)

install: $(LIBTARGET) $(PROGS)
	test -d $(SQUIDHOME) || mkdir -p $(SQUIDHOME)
	test -d $(BINDIR)    || mkdir -p $(BINDIR)
	test -d $(SCRIPTDIR) || mkdir -p $(SCRIPTDIR)
	test -d $(MANDIR)/man$(MANSUFFIX) || mkdir -p $(MANDIR)/man$(MANSUFFIX)
	cp $(LIBTARGET) $(SQUIDHOME)/
	cp $(HDRS)  $(SQUIDHOME)/
	cp $(PROGS) $(BINDIR)/
	for scriptfile in $(SCRIPTS); do \
	  cp Scripts/$$scriptfile $(SCRIPTDIR)/;\
	done
#	@for manpage in $(MANS); do\
#	  $(INSTMAN) $$manpage.man $(MANDIR)/man$(MANSUFFIX)/$$manpage.$(MANSUFFIX);\
#	done

depend:
	makedepend $(SRC)

clean:
	-rm -f *.o *~ core TAGS llib-lsquid.ln $(PROGS) $(LIBTARGET)

lint:
	lint -u $(SRC)

lintlib:
# Note that other systems will use "-o squid" instead to produce a lint library
#	lint -Csquid $(SRC)
	lint -x -u -o squid $(SRC)
	cp llib-lsquid.ln $(SQUIDHOME)/

tags:
	etags -t $(SRC) $(HDRS)

.c.o:
	$(CC) $(CFLAGS) $(MDEFS) -c $<		
