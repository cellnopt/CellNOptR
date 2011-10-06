static double interMat[12][12]=
{
	{-1,-1,0,0,0,0,0,0,0,0,0,1,},
	{0,0,0,0,0,0,0,0,0,0,1,-1,},
	{0,1,0,0,0,0,0,0,0,0,0,0,},
	{0,0,-1,0,0,0,0,0,1,1,0,0,},
	{0,0,0,-1,0,0,0,1,0,0,0,0,},
	{0,0,0,0,0,1,1,0,0,0,0,0,},
	{0,0,1,0,0,0,0,0,0,0,0,-1,},
	{0,0,0,1,0,0,0,-1,0,0,0,0,},
	{0,0,0,0,-1,-1,0,-1,-1,0,0,0,},
	{0,0,0,0,1,0,0,0,0,0,0,0,},
	{0,0,0,0,0,0,0,0,0,-1,-1,0,},
	{1,0,0,0,0,0,-1,0,0,0,0,0,}
};
static double notMat[12][12]=
{
	{1,1,1,1,1,1,1,1,1,1,1,1,},
	{1,1,1,1,1,1,1,1,1,1,1,1,},
	{1,1,1,1,1,1,1,1,1,1,1,1,},
	{1,1,1,1,1,1,1,1,1,1,1,1,},
	{1,1,1,1,1,1,1,1,1,1,1,1,},
	{1,1,1,1,1,1,1,1,1,1,1,1,},
	{1,1,1,1,1,1,1,1,1,1,1,0,},
	{1,1,1,1,1,1,1,0,1,1,1,1,},
	{1,1,1,1,1,1,1,1,1,1,1,1,},
	{1,1,1,1,1,1,1,1,1,1,1,1,},
	{1,1,1,1,1,1,1,1,1,1,1,1,},
	{1,1,1,1,1,1,1,1,1,1,1,1,}
};
static double valueSignals[9][7]=
{
	{0,0,0,0,0,0,0,},
	{0,0,0,0,0,0,0,},
	{0,0,0,0,0,0,0,},
	{0,0,0,0,0,0,0,},
	{0,0,0,0,0,0,0,},
	{0,0,0,0,0,0,0,},
	{0,0,0,0,0,0,0,},
	{0,0,0,0,0,0,0,},
	{0,0,0,0,0,0,0,}
};
static double valueStimuli[9][2]=
{
	{1,0,},
	{0,1,},
	{1,1,},
	{1,0,},
	{0,1,},
	{1,1,},
	{1,0,},
	{0,1,},
	{1,1,}
};
static double valueInhibitors[9][2]=
{
	{0,0,},
	{0,0,},
	{0,0,},
	{1,0,},
	{1,0,},
	{1,0,},
	{0,1,},
	{0,1,},
	{0,1,}
};
static int indexSignals[7]=
{
	7,	6,	10,	12,	3,	5,	8
};
static int indexInhibitors[2]=
{
	2,	4
};
static int indexStimuli[2]=
{
	11,	9
};
static double timeSignals [2]=
{
	0,	10
};
static double odeParameters[38]=
{
	2.203934e+000,	6.262312e-001,	1.359371e+000,	8.331686e-001,	8.413469e-001,	1.088145e+000,	5.515802e-001,	4.646723e-001,	1.968970e+000,	6.239473e-001,	8.709632e-001,	4.088349e+000,	4.171849e-001,	4.874472e+000,	8.243882e-001,	1.261108e-001,	2.448862e+000,	6.793143e-001,	2.056699e+000,	4.625627e-001,	5.942217e-001,	4.564655e+000,	4.420256e-001,	3.202982e+000,	1.512126e-001,	1.608986e-001,	2.505269e+000,	3.879579e-001,	3.688188e-001,	1.786170e+000,	7.341170e-001,	2.064854e-001,	4.098042e+000,	2.482593e-001,	9.666349e-001,	4.161417e+000,	4.353775e-001,	2.235785e-001
};
static int nPars=38;
static int nRows=12;
static int nCols=12;
static int nStimuli=2;
static int nTimes=2;
static int nSignals=7;
static int nInhibitors=2;
static int nTimes=2;
static int nExperiments=9;
