# Makefile for pknots
#
#
MYLIBDIR = -L $(EASELDIR)
MYINCDIR = -I $(EASELDIR)

LIBS     = -lm
MYLIBS   = -leasel
LFLAGS   = -L ~/src/easel
IFLAGS   = -I ~/src/easel

CC       = gcc
CFLAGS   = -g -Wall


PROGS = pknots 

HDRS =  cfg.h proto.h protovhx.h protovx.h protowhx.h protowbx.h protowx.h protozhxyhx.h 

MAIN =  pknots.c\

OBJS =  cfgio.o  fillmtx.o\
 	filltrvhx.o filltrvx.o filltrwbx.o filltrwhx.o filltrwx.o filltrzhxyhx.o \
	globals.o irredsurf.o pk_util.o model.o rnaio.o rnaoutput.o \
	rnaparam.o trace.o tracemtx.o viterbi.o \
        vhxgraphs.o vxgraphs.o wbxgraphs.o whxgraphs.o wxgraphs.o zhxyhxgraphs.o

all: ${PROGS} ~/src/easel/libeasel.a

${PROGS}: %: %.o ${OBJS}
	${CC} ${MYINCDIR} ${CFLAGS} ${IFLAGS} ${LFLAGS} -o $@ $@.o ${OBJS} ${LIBS} ${MYLIBS} 

.c.o:
	${CC} ${MYINCDIR} ${IFLAGS} ${CFLAGS} -c $<

clean:
	rm *.o *~ ${PROGS}

tags:
	etags Makefile *.c *.h
