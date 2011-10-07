
CC 	= gcc
CFLAGS	= -g


all: CellNOptODEs


CellNOptODEs: CellNOptODEs.o

CellNOptODEs.o : CellNOptODEs.c

clean:
	rm -f CellNOptODEs CellNOptODEs.o

#	# $@ evaluate the current target
#	$(CC) $(LDFLAGS) -o $@ $^
#
#
Rlib:
	#create a so lib for R
	R CMD SHLIB  CellNOptODEs.c



