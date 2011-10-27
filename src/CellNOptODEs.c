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
#include "Toy_Model_MMB_Feedback_OriginalPars_C_test.h"

int main(void)
{
  int i,j,counter;
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
  double** simResults;
  double* inhibitor_array;;

  indexSig=(int*)malloc(nSignals*sizeof(int));
  for (i = 0; i < nSignals; i++)
  {
    indexSig[i] = indexSignals[i];
  }

  indexStim=(int*)malloc(nStimuli*sizeof(int));
  for (i = 0; i < nStimuli; i++)
  {
    indexStim[i] = indexStimuli[i];
  }

  indexInh=(int*)malloc(nInhibitors*sizeof(int));
  for (i = 0; i < nInhibitors; i++)
  {
    indexInh[i] = indexInhibitors[i];
  }

  odePARAMETERS=(double*)malloc(nPars*sizeof(double));
  for (i = 0; i < nPars; i++)
  {
    odePARAMETERS[i]=odeParameters[i];
  }

  interMAT = (int**) malloc(nRows * sizeof(int*));
  for (i = 0; i < nRows; i++)
  {
    interMAT[i] = (int*) malloc(nCols * sizeof(int));
    for (j = 0; j < nCols; j++)
    {
      interMAT[i][j] = (int)interMat[i][j];
    }
  }

  valueSIGNALS = (double**)malloc(nExperiments * sizeof(double*));
  for (i = 0; i < nExperiments; i++)
  {
    valueSIGNALS[i] = (double*)malloc(nSignals*sizeof(double));
    for (j = 0; j < nSignals; j++)
    {
      valueSIGNALS[i][j]= valueSignals[i][j];
    }
  }

  valueINHIBITORS = (double**)malloc(nExperiments * sizeof(double*));
  for (i = 0; i < nExperiments; i++)
  {
    valueINHIBITORS[i] = (double*)malloc(nInhibitors*sizeof(double));
    for (j = 0; j < nInhibitors; j++)
    {
      valueINHIBITORS[i][j]=valueInhibitors[i][j];
    }
  }

  valueSTIMULI = (double**)malloc(nExperiments * sizeof(double*));
  for (i = 0; i < nExperiments; i++)
  {
    valueSTIMULI[i] = (double*)malloc(nStimuli*sizeof(double));
    for (j = 0; j < nStimuli; j++)
    {
      valueSTIMULI[i][j]=(double)valueStimuli[i][j];
    }
  }

  timeSig= (double*)malloc(nTimes*sizeof(double*));
   for (i = 0; i < nTimes; i++)
   {
	   timeSig[i]=timeSignals[i];
   }

  notMAT = (int**)malloc(nRows * sizeof(int*));
  for (i = 0; i < nRows; i++)
  {
    notMAT[i] = (int*)malloc(nCols*sizeof(int));
    for (j = 0; j < nCols; j++)
    {
      notMAT[i][j]=notMat[i][j];
    }
  }

  /* Fill the CNOStructure */
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
  tempData.numInputs = getNumInputs(tempData.adjacencyMatrix,tempData.nCols);
  tempData.numBits = getNumBits(tempData.numInputs,tempData.nRows);
  tempData.isState = findStates(tempData.adjacencyMatrix,tempData.nRows);
  tempData.truthTables = getTruthTables(tempData.adjacencyMatrix,tempData.interMat,
		  tempData.notMat,tempData.isState,tempData.numInputs,tempData.numBits,tempData.nRows,tempData.nCols);

  state_array= malloc((tempData.nRows)*sizeof(double));
  simResults= malloc(tempData.nTimes*sizeof(double*));
  inhibitor_array= malloc((tempData.nRows)*sizeof(double));

  tempData.state_index=getStateIndex(tempData.adjacencyMatrix,tempData.nRows);
  tempData.inhibitor_array=inhibitor_array;
  tempData.state_array=state_array;

  counter=0;
  for (i = 0; i < nRows; ++i)
	  if(tempData.isState[i]) counter++;

  tempData.nStates = counter;

  for (i = 0; i <tempData.nTimes; ++i)
  {
	  simResults[i]=malloc(tempData.nStates*sizeof(double));
  }

  tempData.sim_results=simResults;

  data=malloc(sizeof(tempData));

  *data=tempData;

  simulateODE(data);

  free(data);

  free(indexSig);
  free(indexStim);
  free(indexInh);
  free(odePARAMETERS);


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

   for (i = 0; i <tempData.nTimes; ++i)
	   free(simResults[i]);
   free(simResults);

   free(state_array);
   free(malloc(tempData.nTimes*sizeof(double*)));
   free(inhibitor_array);

return EXIT_SUCCESS;
}


