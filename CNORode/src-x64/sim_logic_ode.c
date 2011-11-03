/*
 ============================================================================
 Name        : CellNOptODEs.c
 Author      : DH, TC 
 Version     :
 Copyright   : Your copyright notice
 Description : Hello World in C, Ansi-style
 ============================================================================
 */

#include <stdio.h>
#include <stdlib.h>
#include "CNOStructure.h"

#include <math.h>
#include <string.h>


#include <R.h>
#include <Rinternals.h>
#include <R_ext/Print.h>

SEXP sim_logic_ode
(
		SEXP interMat_in,			SEXP notMat_in,				SEXP nRows_in,
		SEXP nCols_in,				SEXP nPars_in,				SEXP timeSignals_in,
		SEXP valueInhibitors_in,	SEXP valueSignals_in,		SEXP valueStimuli_in,
		SEXP nTimes_in,				SEXP nExperiments_in,		SEXP nSignals_in,
		SEXP indexSignals_in,		SEXP nStimuli_in,			SEXP indexStimuli_in,
		SEXP nInhibitors_in,		SEXP indexInhibitors_in,	SEXP odeParameters_in,
		SEXP verbose_in
)
{

	SEXP arg,ans;
	 double *x;
	 int countProtected=0;
	 int i,j,k,nStates,counter;
	 int *indexSig;
	 int *indexStim;
	 int *indexInh;
	 double* timeSig;
	 int **interMAT;
	 int **notMAT;
	 double *odePARAMETERS;
	 double **valueSIGNALS;
	 double **valueINHIBITORS;
	 double **valueSTIMULI;
	 CNOStructure tempData;
	 CNOStructure *data;
	 double* state_array;
	 double*** simResults;
	 double* inhibitor_array;

	 int nRows = INTEGER(nRows_in)[0];
	 int nCols=INTEGER(nCols_in)[0];
	 int nSignals=INTEGER(nSignals_in)[0];
	 int nStimuli=INTEGER(nStimuli_in)[0];
	 int nInhibitors=INTEGER(nInhibitors_in)[0];
	 int nPars=INTEGER(nPars_in)[0];
	 int nExperiments=INTEGER(nExperiments_in)[0];
	 int nTimes=INTEGER(nTimes_in)[0];
	 int verbose=INTEGER(verbose_in)[0];

	 counter=0;
	 indexSig=(int*)malloc(nSignals*sizeof(int));
	 for (i = 0; i < nSignals; i++)
	 {
		 indexSig[i] = INTEGER(indexSignals_in)[counter++];
	 }

	 counter=0;
	 indexStim=(int*)malloc(nStimuli*sizeof(int));
	 for (i = 0; i < nStimuli; i++)
	 {
		 indexStim[i] = INTEGER(indexStimuli_in)[counter++];
	 }

	 counter=0;
	 indexInh=(int*)malloc(nInhibitors*sizeof(int));
	 for (i = 0; i < nInhibitors; i++)
	 {
		 indexInh[i] = INTEGER(indexInhibitors_in)[counter++];
	 }

	 odePARAMETERS=(double*)malloc(nPars*sizeof(double));
	 for (i = 0; i < nPars; i++)
	 {
		 odePARAMETERS[i]=REAL(odeParameters_in)[i];
	 }

	 counter=0;
	 interMAT = (int**) malloc(nRows * sizeof(int*));
	  for (i = 0; i < nRows; i++)
	  {
	    interMAT[i] = (int*) malloc(nCols * sizeof(int));
	    for (j = 0; j < nCols; j++)
	    {
	      interMAT[i][j] = INTEGER(interMat_in)[counter++];
	    }
	  }

	  counter=0;
	  notMAT = (int**)malloc(nRows * sizeof(int*));
	  for (i = 0; i < nRows; i++)
	  {
		  notMAT[i] = (int*)malloc(nCols*sizeof(int));
		  for (j = 0; j < nCols; j++)
		  {
			  notMAT[i][j]=INTEGER(notMat_in)[counter++];
		  }
	  }

	  counter=0;
	  valueSIGNALS = (double**)malloc(nExperiments * sizeof(double*));
	  for (i = 0; i < nExperiments; i++)
	  {
		  valueSIGNALS[i] = (double*)malloc(nSignals*sizeof(double));
		  for (j = 0; j < nSignals; j++)
		  {
			  valueSIGNALS[i][j]= REAL(valueSignals_in)[counter++];
		  }
	  }
	  counter=0;
	  valueINHIBITORS = (double**)malloc(nExperiments * sizeof(double*));
	  for (i = 0; i < nExperiments; i++)
	  {
	    valueINHIBITORS[i] = (double*)malloc(nInhibitors*sizeof(double));
	    for (j = 0; j < nInhibitors; j++)
	    {
	      valueINHIBITORS[i][j]=REAL(valueInhibitors_in)[counter++];
	    }
	  }

	  counter=0;
	  valueSTIMULI = (double**)malloc(nExperiments * sizeof(double*));
	  for (i = 0; i < nExperiments; i++)
	  {
		  valueSTIMULI[i] = (double*)malloc(nStimuli*sizeof(double));
		  for (j = 0; j < nStimuli; j++)
		  {
			  valueSTIMULI[i][j]=REAL(valueStimuli_in)[counter++];
		  }
	  }

	  timeSig= (double*)malloc(nTimes*sizeof(double));
	  for (i = 0; i < nTimes; i++)
	  {
		  timeSig[i]=REAL(timeSignals_in)[i];
	  }


	  // Fill the CNOStructure
	  tempData.interMat=interMAT;
	  tempData.notMat=notMAT;
	  tempData.valueSignals=valueSIGNALS;
	  tempData.valueInhibitors=valueINHIBITORS;
	  tempData.valueStimuli=valueSTIMULI;
	  tempData.indexSignals=indexSig;
	  tempData.indexStimuli=indexStim;
	  tempData.indexInhibitors=indexInh;
	  tempData.timeSignals=timeSig;
	  tempData.odeParameters=odePARAMETERS;
	  tempData.nPars=nPars;
	  tempData.nRows=nRows;
	  tempData.nCols=nCols;
	  tempData.nStimuli=nStimuli;
	  tempData.nInhibitors=nInhibitors;
	  tempData.nSignals=nSignals;
	  tempData.nTimes=nTimes;

	  tempData.adjacencyMatrix = getAdjacencyMatrix(tempData.interMat,tempData.nRows,tempData.nCols);
	  tempData.numInputs = getNumInputs(tempData.adjacencyMatrix,tempData.nRows);

	  tempData.numBits = getNumBits(tempData.numInputs,tempData.nRows);
	  tempData.isState = findStates(tempData.adjacencyMatrix,tempData.nRows);


	  tempData.truthTables = getTruthTables(tempData.adjacencyMatrix,tempData.interMat,
	  tempData.notMat,tempData.isState,tempData.numInputs,tempData.numBits,tempData.nRows,tempData.nCols);

	  state_array= (double*)malloc(tempData.nRows*sizeof(double));
	  inhibitor_array=(double*)malloc((tempData.nRows)*sizeof(double));

	  tempData.state_index=(int*)getStateIndex(tempData.adjacencyMatrix,tempData.nRows);
	  tempData.inhibitor_array=inhibitor_array;
	  tempData.state_array=state_array;

	  counter=0;
	  for (i = 0; i < nRows; ++i)
		  if(tempData.isState[i]) counter++;

	  tempData.nStates = counter;
	  //simResults=(double**)malloc(nExperiments*sizeof(double*));
	  simResults=(double***)malloc(nExperiments*sizeof(double**));
	  for (i = 0; i <nExperiments; ++i)
	  {
		  simResults[i]=(double**)malloc(nTimes*sizeof(double*));
		  //  simResults=(double**)malloc(nTimes*sizeof(double*));
		  for (j = 0; j <nTimes; ++j)
		  {
			  simResults[i][j]=(double*)malloc(tempData.nStates*sizeof(double));
		  }
	  }

	  tempData.sim_results=simResults;

	  data=malloc(sizeof(tempData));

	  *data=tempData;

	  for (i = 0; i <nExperiments; ++i)
	  {
		  simulateODE(data,i,verbose);
	  }

	  //Put the data into a LIST
	  PROTECT(ans = allocVector(VECSXP, nTimes));
	  countProtected++;
	  for (i = 0; i <nTimes; ++i)
	  {
		  PROTECT(arg = allocMatrix(REALSXP, nExperiments,nRows));
		  countProtected++;
		  counter=0;
			  for (k = 0; k < nRows; ++k)
			  {
				  for (j = 0; j < nExperiments; ++j)
				  {
					  if(tempData.isState[k])
					  	  REAL(arg)[counter++]=simResults[j][i][tempData.state_index[k]];
				  	  else  REAL(arg)[counter++]=NA_REAL;
				  }
		  }
		  SET_VECTOR_ELT(ans,i,arg);
		  UNPROTECT(1);
	  }
	  UNPROTECT(1);
	  //ALL THE DATA IS NOW INSIDE ANS

	  for (i = 0; i <nExperiments; ++i)
	  {
		  for (j = 0; j <nTimes; ++j)
		  {
			  free(simResults[i][j]);
		  }
		  free(simResults[i]);
	  }
	  free(simResults);

	  free(data);
	  free(indexSig);
	  free(indexStim);
	  free(indexInh);
	  free(odePARAMETERS);

	  for (i = 0; i < nRows; i++)
	  {
		free(tempData.adjacencyMatrix[i]);
	    free(tempData.truthTables[i]);
	  }
	  free(tempData.truthTables);
	  free(tempData.adjacencyMatrix);

	  for (i = 0; i < nExperiments; i++)
		  free(valueSIGNALS[i]);
	  free(valueSIGNALS);

	  for (i = 0; i < nExperiments; i++)
		  free(valueSTIMULI[i]);
	  free(valueSTIMULI);

	  for (i = 0; i < nExperiments; i++)
		  free(valueINHIBITORS[i]);
	  free(valueINHIBITORS);

	  for (i = 0; i < nRows; i++)
		  free(notMAT[i]);
	  free(notMAT);

	  for (i = 0; i < nRows; i++)
		  free(interMAT[i]);
	  free(interMAT);

	  free(state_array);
	  free(inhibitor_array);
	  free(tempData.state_index);
	  free(tempData.numBits);
	  free(tempData.numInputs);

	  return(ans);
}

