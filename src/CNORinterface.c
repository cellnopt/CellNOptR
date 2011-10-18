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
#include <R.h>


void hello(int *in)
{
    int i, n=in[0];
 
    for(i=0; i < n; i++) {
        Rprintf("Hello, world!\n");
    }
}

void hello2(double *xIn, int *nrowsIn, int *ncolsIn)
{
    int nRows = nrowsIn[0];
    int nCols = ncolsIn[0];

    int i, j;
 
    for(i=0; i < nRows; i++) {
        for(j=0; j < nCols; j++) {
        Rprintf("%f ", xIn[i*nRows + j]);
    }
        Rprintf("\n");
    }
}


 
int CNORinterfaceTest(int *interMat, int *notMat, int *nRowsIn, int *nColsIn, 
            double * valueInhibitors, int * indexInhibitors, int *nInhibitorsIn,
            double * valueSignals,    int * indexSignals,    int *nSignalsIn,
            double * valueStimuli,    int * indexStimuli,    int *nStimuliIn,
            double *  timeSignals,     int *nTimesIn,
            int *nExperimentsIn,
            double *odeParameters, int *nParsIn) 

{
 /*Reading input*/
 int nSignals = *nSignalsIn;
 int nStimuli = *nStimuliIn;
 int nTimes = *nTimesIn;
 int nExperiments = *nExperimentsIn;
 int nInhibitors = *nInhibitorsIn;
 int nPars = *nParsIn;
 int nRows = *nRowsIn;
 int nCols = *nColsIn;





  int i,j,nStates,counter;
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


  indexSig = (int*) malloc(nSignals*sizeof(int));
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
      interMAT[i][j] = 2 * (int)interMat[i * nRows + j];
      interMat[i * nRows + j]*=2;
    }
  }

  valueSIGNALS = (double**) malloc(nExperiments * sizeof(double*));
  for (i = 0; i < nExperiments; i++)
  {
    valueSIGNALS[i] = (double*) malloc(nSignals * sizeof(double));
    for (j = 0; j < nSignals; j++)
    {
      valueSIGNALS[i][j]= valueSignals[i * nExperiments + j];
    }
  }

  valueINHIBITORS = (double**)malloc(nExperiments * sizeof(double*));
  for (i = 0; i < nExperiments; i++)
  {
    valueINHIBITORS[i] = (double*)malloc(nInhibitors*sizeof(double));
    for (j = 0; j < nInhibitors; j++)
    {
      valueINHIBITORS[i][j]=valueInhibitors[i * nExperiments + j];
    }
  }

  valueSTIMULI = (double**)malloc(nExperiments * sizeof(double*));
  for (i = 0; i < nExperiments; i++)
  {
    valueSTIMULI[i] = (double*)malloc(nStimuli*sizeof(double));
    for (j = 0; j < nStimuli; j++)
    {
      valueSTIMULI[i][j]=(double)valueStimuli[i * nExperiments + j];
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
      notMAT[i][j]=notMat[i * nRows  + j];
    }
  }



  /* Fill the CNOStructure */
  /*tempData.interMat=interMAT;
  tempData.notMat=notMAT;
  tempData.valueSignals=valueSIGNALS;
  tempData.valueInhibitors=valueINHIBITORS;
  tempData.valueStimuli=valueSTIMULI;
  tempData.indexSignals=indexSig;
  tempData.indexStimuli=indexStim;
  tempData.indexInhibitors=indexInh;
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

  counter=0;
  for (i = 0; i < nRows; ++i){
    if(tempData.isState[i]){
        counter++;
    }
  }

  tempData.nStates = counter;

  data=malloc(sizeof(tempData));

  *data=tempData;

  //simulateODE(data);

  free(data);
*/

  /* simple pointers first */
  free(indexSig);
  free(indexStim);
  free(indexInh);
  free(odePARAMETERS);

  /* ** pointer then */
  for (i = 0; i < nExperiments; i++)
    free(valueSIGNALS[i]);
  free(valueSIGNALS);
  valueSignals[0] = 10;


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

}

