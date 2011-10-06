typedef struct
{
	int** interMat;
	int** notMat;
	double** valueSignals;
	double** valueInhibitors;
	double** valueStimuli;
	int* indexSignals;
	int* indexStimuli;
	int* indexInhibitors;
	double* odeParameters;
	int nPars;
	int nRows;
	int nCols;
	int nStimuli;
	int nInhibitors;
	int nSignals;
	int nTimes;
	int nExperiments;

}CNOStructure;
