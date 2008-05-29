# Makefile for pknots
#
#
MYLIBDIR = -L $(EASELDIR)
MYINCDIR = -I $(EASELDIR)

LIBS     = -lm
MYLIBS   = -leasel
LFLAGS   = -L $(EASELDIR)
IFLAGS   = -I $(EASELDIR)

CC       = gcc
CFLAGS   = -g -Wall

PROGS =	pknots 

OBJS =	pk_cfgio.o\
	pk_cyk.o\
	pk_fillmtx.o\
	pk_filltrvhx.o\
	pk_filltrvx.o\
	pk_filltrwbx.o\
	pk_filltrwhx.o\
	pk_filltrwx.o\
	pk_filltrzhxyhx.o\
	pk_globals.o\
	pk_irredsurf.o\
	pk_util.o\
	pk_model.o\
	pk_rnaoutput.o\
	pk_rnaparam.o\
	pk_trace.o\
	pk_tracemtx.o\
	pk_vhxgraphs.o\
	pk_vxgraphs.o\
	pk_wbxgraphs.o\
	pk_whxgraphs.o\
	pk_wxgraphs.o\
	pk_zhxyhxgraphs.o

HDRS =	pknots.h\
	pk_cfgio.h\
	pk_cyk.h\
	pk_fillmtx.h\
	pk_filltrvhx.h\
	pk_filltrvx.h\
	pk_filltrwbx.h\
	pk_filltrwhx.h\
	pk_filltrwx.h\
	pk_filltrzhxyhx.h\
	pk_irredsurf.h\
	pk_util.h\
	pk_model.h\
	pk_rnaoutput.h\
	pk_rnaparam.h\
	pk_trace.h\
	pk_tracemtx.h\
	pk_vhxgraphs.h\
	pk_vxgraphs.h\
	pk_wbxgraphs.h\
	pk_whxgraphs.h\
	pk_wxgraphs.h\
	pk_zhxyhxgraphs.h


all: ${PROGS} $(EASELDIR)/libeasel.a

${PROGS}: %: ${OBJS} %.o 
	${CC} ${MYINCDIR} ${CFLAGS} ${IFLAGS} ${LFLAGS} -o $@ $@.o ${OBJS} ${LIBS} ${MYLIBS} 

.c.o:
	${CC} ${MYINCDIR} ${IFLAGS} ${CFLAGS} -c $<

clean:
	rm *.o *~ ${PROGS}

tags:
	etags Makefile *.c *.h
