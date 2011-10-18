#include <math.h>
#include <stdio.h>

/* Sundials Header Files */

#include "CVODES/include/cvodes/cvodes.h"           /* prototypes for CVODES fcts. and consts. */
#include "CVODES/include/cvodes/cvodes_dense.h"     /* prototype for CVDense */
#include "CVODES/include/cvodes/cvodes_band.h"
#include "CVODES/include/cvodes/cvodes_diag.h"
#include "CVODES/include/nvector/nvector_serial.h"  /* serial N_Vector types, fcts., and macros */
#include "CVODES/include/sundials/sundials_dense.h" /* definitions DenseMat and DENSE_ELEM */
#include "CVODES/include/sundials/sundials_types.h" /* definition of type realtype */
#include "CNOStructure.h"

#define Ith(v,i) ( NV_DATA_S(v)[i] )
#define IJth(A,i,j) DENSE_ELEM(A,i,j)

typedef int rhs_func(realtype t, N_Vector y, N_Vector ydot, void *f_data);

int rhsODE(realtype t, N_Vector y, N_Vector ydot, void *data);

static int check_flag(void *flagvalue, char *funcname, int opt);

int simulateODE(CNOStructure* data)
{
	printf("\n\n Test %d",(*data).nStates);
	int i,j,neq,index_signal,counter,flag, flagr, iout;
	int* state_array=malloc((*data).nRows*sizeof(int));
	realtype t, tout, ti, tf;
	N_Vector y, abstol;
	void *cvode_mem;
	int exp_num=4;
	rhs_func *rhs;
	double maxStepSize=0.0;
	int maxNumSteps=10000;
	 realtype *atol2;
	 realtype atol=0.000000001;
	 realtype reltol=0.00000001;
	 double** simResults=malloc((*data).nTimes*sizeof(double*));

	int verbose=1;

	cvode_mem = NULL;
	abstol = NULL;
	y = NULL;

	neq=(*data).nStates;

    y = N_VNew_Serial(neq);
	if (check_flag((void *)y, "N_VNew_Serial", 0))
    {
		if(verbose)printf("\nSolver failed. . .\n");
		return;
	}

    abstol = N_VNew_Serial(neq);
	if (check_flag((void *)abstol, "N_VNew_Serial", 0))
    {
		if(verbose)printf("\nSolver failed. . .\n");
		return;
	}


    /* Initialize y */
	for(i=0; i<(*data).nRows; i++)
			state_array[i] = 0;

	for(i=0; i<(*data).nRows; ++i)
	{
		if((*data).isState[i])
		{
			for (j = 0; j < (*data).nSignals; j++)
			{
			//	The passed indexes are from 1 to N intead from 0 to N-1
				if((*data).indexSignals[j]==i+1)
				{
					state_array[i] = (*data).valueSignals[exp_num][(*data).indexSignals[j]-1];
				}
			}
		}
		else
		{
			for (j = 0; j < (*data).nStimuli; j++)
			{
				//	The passed indexes are from 1 to N intead from 0 to N-1
				if((*data).indexStimuli[j]==i+1)
				{
					state_array[i] = (*data).valueStimuli[exp_num][(*data).indexStimuli[j]-1];
				}
			}
		}
	}

	for (i = 0; i <(*data).nStates; ++i)
	{
		simResults[i]=malloc((*data).nStates*sizeof(double));
	}

	for(i=0; i<(*data).nRows; i++)
	{
		if((*data).isState[i])
		{
			Ith(y,i) = state_array[i];
			simResults[0][i]=state_array[i];

		}
	}

	//Dont use this field out of simulation
	(*data).state_array=state_array;

	cvode_mem = CVodeCreate(CV_BDF, CV_NEWTON);
	if (check_flag((void *)cvode_mem, "CVodeCreate", 0)) {
		if(verbose)printf("\nSolver failed. . .\n");
		return(0);
	}

	if(verbose)printf("Solver Memory Allocated\n");

	/* Set f_data */
	 flag = CVodeSetFdata(cvode_mem, data);
	 if(check_flag(&flag, "CVodeSetFdata", 1))
	 {
		 if(verbose)printf("\nSolver failed. . .\n");
		 return(0);
	 }

	 atol2 = &atol;

	 rhs =*rhsODE;

	 ti=(*data).timeSignals[0];
	 flag = CVodeMalloc(cvode_mem, rhs, ti, y, CV_SS, reltol, atol2);

	 if (check_flag(&flag, "CVodeMalloc", 1))
	 {
		if(verbose)printf("\nSolver failed. . .\n");
	 	return(0);
	 }

	 flag = CVDense(cvode_mem, neq);
	 if (check_flag(&flag, "CVDense", 1))
	 {
		 if(verbose)printf("\nSolver failed. . .\n");
		 return(0);
	 }
	 if(verbose)printf("CVDENSE Solver Initiated\n");

	 /* Set maxnumsteps */
	 flag = CVodeSetMaxNumSteps(cvode_mem, maxNumSteps);
	 if(check_flag(&flag, "CVodeSetMaxNumSteps", 1))
	 {
		if(verbose)printf("\nSolver failed. . .\n");
		return(0);
	 }
	 if(verbose)printf("Max number of steps: %i\n", maxNumSteps);

	  /* Set maxstepsize */
	  flag = CVodeSetMaxStep(cvode_mem, maxStepSize);
	  if(check_flag(&flag, "CVodeSetMaxStep", 1))
	  {
		  if(verbose)printf("\nSolver failed. . .\n");
		  return;
	  }
	  if(verbose)printf("Max step size: %f\n", maxStepSize);

	  for (i = 1; i < (*data).nTimes; ++i)
	  {
		  tout=(*data).timeSignals[i];
		  flag = CVode(cvode_mem, tout, y, &t, CV_NORMAL);

		  for (j = 0; j < (*data).nStates; ++j)
		  {
			  simResults[i][j]= (double) Ith(y,i);
		  }

		  if (check_flag(&flag, "CVode", 1))
		  {
			  if(verbose)printf("\nSolver failed. . .\n");
			  break;
		  }
	  }

	return(0);
}

static int check_flag(void *flagvalue, char *funcname, int opt)
{
  int *errflag;

  /* Check if SUNDIALS function returned NULL pointer - no memory allocated */
  if (opt == 0 && flagvalue == NULL) {
    fprintf(stderr, "\nSUNDIALS_ERROR: %s() failed - returned NULL pointer\n\n",
	    funcname);
    return(1); }

  /* Check if flag < 0 */
  else if (opt == 1) {
    errflag = (int *) flagvalue;
    if (*errflag < 0) {
      fprintf(stderr, "\nSUNDIALS_ERROR: %s() failed with flag = %d\n\n",
	      funcname, *errflag);
      return(1); }}

  /* Check if function returned NULL pointer - no memory allocated */
  else if (opt == 2 && flagvalue == NULL) {
    fprintf(stderr, "\nMEMORY_ERROR: %s() failed - returned NULL pointer\n\n",
	    funcname);
    return(1); }

  return(0);
}



