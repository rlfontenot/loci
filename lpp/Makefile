ACCENT = /home/mz33/accent-0.7/accent/accent
ART    = /home/mz33/accent-0.7/art/art.o
LEX    = flex
CC     = gcc

default: cins

cins:	FRC
	$(ACCENT) parser.acc
	$(LEX) parser.lex
	$(CC) -o lpp -g memory.c hash.c proc.c yygrammar.c lex.yy.c auxil.c $(ART)

FRC:

clean:
	rm -fr lpp lex.yy.c yygrammar.c yygrammar.h

clean1:
	rm -fr lpp lex.yy.c yygrammar.c yygrammar.h *~
