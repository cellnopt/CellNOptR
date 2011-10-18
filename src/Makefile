
CC 	= gcc
CFLAGS	= -g
INCLS = -I./CVODES/include/cvodes/ -I./CVODES/include/nvector -I./CVODES/include/sundials
LIBS = -L. -lm  -L./CVODES -lcvodes 


OBJS = findStates.o getNumBits.o getTruthTables.o printAdjMat.o printNminiTerms.o rhsODE.o decimal2binary.o getAdjacencyMatrix.o  getNumInputs.o  normHill.o        printInterMat.o  printTruthTables.o  simulateODE.o


all: CellNOptODEs



$(OBJS): CNOStructure.h
	$(CC) $(CFLAGS) $(INCLS) $(LIBS) -c $*.c




CellNOptODEs: $(OBJS)
	make -C CVODES 
	$(CC) -o CellNOptODEs $(CFLAGS) $(INCLS)  $(OBJS) $(LIBS) CellNOptODEs.c



clean:
	rm -f CellNOptODEs CellNOptODEs.o $(OBJS)

#	# $@ evaluate the current target
#	$(CC) $(LDFLAGS) -o $@ $^
#
#
Rlib:
	#create a so lib for R
	R CMD SHLIB  CellNOptODEs.c



