/* Prototype for a CNOstruct */


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
	double* timeSignals;
	int* isState;
	int* isInput;
	int** adjacencyMatrix;
	int** truthTables;
	int* numInputs;
	int* numBits;
	double* odeParameters;
	int nPars;
	int nRows;
	int nCols;
	int nStimuli;
	int nInhibitors;
	int nSignals;
	int nTimes;
	int nExperiments;
	int nStates;
	//Use this only inside simulations
	double* state_array;
	int* state_index;

}CNOStructure;
