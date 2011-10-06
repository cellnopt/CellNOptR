/*
 ============================================================================
 Name        : CellNOptODEs.c
 Author      : 
 Version     :
 Copyright   : Your copyright notice
 Description : Hello World in C, Ansi-style
 ============================================================================
 */

#include <stdio.h>
#include <stdlib.h>
#include "CNOStructure.h"

static int interMat [4][3]=
{
	{1,0,1,0},
	{1,0,0,1},
	{1,0,0,1}
};

static int notMat[4][3]=
{
		{1,0,1,0},
		{1,0,0,1},
		{1,0,0,1}
};

static double valueSignals[4][3]=
{
	{1e-04,0,1,0},
	{1,0,0,1},
	{1,0,0,1}
};

static double valueInhibitors[4][3]=
{
	{1e-04,0,1,0},
	{1,0,0,1},
	{1,0,0,1}
};

static double valueStimuli[4][3]=
{
	{1e-04,0,1,0},
	{1,0,0,1},
	{1,0,0,1}
};

static int indexSignals[3]={1,2,3};
static int indexStimuli[3]={1,2,3};
static int indexInhibitors[3]={1,2,3};
static double odeParameters[5]={1,2,2,3,5};
static int nPars=3;
static int nRows=4;
static int nCols=3;
static int nStimuli=3;
static int nInhibitors=3;
static int nSignals=3;
static int nTimes=1;
static int nExperiments;

int main(void)
{
	int i,j;
	int *indexSig;
	int *indexStim;
	int *indexInh;
	int **interMAT;
	int **notMAT;
	double *odePARAMETERS;
	double **valueSIGNALS;
	double **valueINHIBITORS;
	double **valueSTIMULI;
    CNOStructure tempData;
	CNOStructure *data;

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

    for (i = 0; i < nRows; i++)
    {
    	interMAT[i] =malloc(nCols*sizeof(int));
        for (j = 0; j < nCols; j++)
        {
            interMAT[i][j]=(int)interMat[j];
        }
    }

    for (i = 0; i < nExperiments; i++)
    {
    	valueSIGNALS[i] =malloc(nSignals*sizeof(int));
        for (j = 0; j < nSignals; j++)
        {
            valueSignals[i][j]=(int)interMat[j];
        }
    }

    for (i = 0; i < nExperiments; i++)
    {
    	valueINHIBITORS[i] =malloc(nInhibitors*sizeof(int));
        for (j = 0; j < nInhibitors; j++)
        {
            valueINHIBITORS[i][j]=(int)valueInhibitors[j];
        }
    }

    for (i = 0; i < nExperiments; i++)
    {
    	valueSTIMULI[i] =malloc(nStimuli*sizeof(int));
        for (j = 0; j < nStimuli; j++)
        {
            valueSTIMULI[i][j]=(int)valueStimuli[j];
        }
    }

    for (i = 0; i < nRows; i++)
    {
    	notMAT[i] =malloc(nCols*sizeof(int));
    	for (j = 0; j < nCols; j++)
    	{
    		notMAT[i][j]=(int)interMat[j];
    	}
    }

	tempData.interMat=interMAT;
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
	tempData.nExperiments=nExperiments;

    data=malloc(sizeof(tempData));
	*data=tempData;

	printf("haha",interMat[0][0]);
	printf("%f",valueSignals[0][0]);
	puts("!!!Hello World!!!"); /* prints !!!Hello World!!! */
	return EXIT_SUCCESS;
}
