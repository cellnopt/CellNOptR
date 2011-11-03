#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <R.h>

/* Sundials Header Files */

#include <cvodes/cvodes.h>
#include <nvector/nvector_serial.h>
#include <sundials/sundials_types.h>
#include <sundials/sundials_math.h>


#include "CNOStructure.h"

#define Ith(v,i) ( NV_DATA_S(v)[i] )
#define IJth(A,i,j) DENSE_ELEM(A,i,j)

typedef int rhs_func(realtype t, N_Vector y, N_Vector ydot, void *f_data);

int rhsODE(realtype t, N_Vector y, N_Vector ydot, void *data);

static int check_flag(void *flagvalue, char *funcname, int opt);

int simulateODE(CNOStructure* data,int exp_num, int verbose)
{
	int i,j,neq,counter,flag, flagr, iout;
	realtype t, tout, ti, tf,abstol;
	N_Vector y;
	void *cvode_mem;
	double maxStepSize=0.0;
	int maxNumSteps=10000;
	realtype atol=1e-5;
	realtype reltol=1e-3;

	cvode_mem = NULL;
	y = NULL;

	neq=(*data).nStates;

    y = N_VNew_Serial(neq);
	if (check_flag((void *)y, "N_VNew_Serial", 0))
    {
		if(verbose)printf("\nSolver failed. . .\n");
		return(0);
	}

    /* Initialize y */
	for(i=0; i<(*data).nRows; i++)
	{
		(*data).state_array[i] = 0;
		(*data).inhibitor_array[i]=0;
	}

	for(i=0; i<(*data).nRows; ++i)
	{
		if((*data).isState[i])
		{
			for (j = 0; j < (*data).nSignals; j++)
			{
			//	The passed indexes are from 1 to N intead from 0 to N-1
				if((*data).indexSignals[j]==i+1)
				{
					(*data).state_array[i] = (*data).valueSignals[exp_num][j];
				}
			}

			for (j = 0; j < (*data).nInhibitors; j++)
			{
				//	The passed indexes are from 1 to N intead from 0 to N-1
				if((*data).indexInhibitors[j]==i+1)
				{
					(*data).inhibitor_array[i] = (*data).valueInhibitors[exp_num][j];
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
					(*data).state_array[i] = (*data).valueStimuli[exp_num][j];
				}
			}
		}
	}

	counter=0;
	for(i=0; i<(*data).nRows; i++)
	{
		if((*data).isState[i])
		{
			Ith(y,(*data).state_index[i]) =
					(*data).state_array[i];
			(*data).sim_results[exp_num][0][(*data).state_index[i]]=
					(*data).state_array[i];
			if(verbose)
			{
				printf("species number initial value %d",i);
				printf("\t%f\n",Ith(y,(*data).state_index[i]));
			}
		}
	}

	cvode_mem = CVodeCreate(CV_BDF, CV_NEWTON);
	if (check_flag((void *)cvode_mem, "CVodeCreate", 0)) {
		if(verbose)printf("\nSolver failed. . .\n");
		N_VDestroy_Serial(y);
		return(0);
	}

	 ti=(*data).timeSignals[0];
	 tf=(*data).timeSignals[(*data).nTimes-1];
	  flag = CVodeInit(cvode_mem,*rhsODE, ti, y);
	 if (check_flag(&flag, "CVodeMalloc", 1))
	 {
		if(verbose)printf("\nSolver failed. . .\n");
		  N_VDestroy_Serial(y);
		  /* Free integrator memory */
		  CVodeFree(&cvode_mem);
	 	return(0);
	 }

	if(verbose)printf("Solver Memory Allocated\n");

	/* Set f_data */
	 flag = CVodeSetUserData(cvode_mem, data);
	 if(check_flag(&flag, "CVodeSetFdata", 1))
	 {
		 if(verbose)printf("\nSolver failed. . .\n");
		  N_VDestroy_Serial(y);
		  /* Free integrator memory */
		  CVodeFree(&cvode_mem);
		 return(0);
	 }

	 flag = CVodeSStolerances(cvode_mem, reltol, atol);
	  if(check_flag(&flag, "CVodeSStolerances", 1)) return(1);


	  flag = CVDense(cvode_mem, neq);
	  if (check_flag(&flag, "CVDense", 1)) return(1);
	 if(verbose)printf("CVDENSE Solver Initiated\n");

	 /* Set maxnumsteps */
	 flag = CVodeSetMaxNumSteps(cvode_mem, maxNumSteps);
	 if(check_flag(&flag, "CVodeSetMaxNumSteps", 1))
	 {
		if(verbose)printf("\nSolver failed. . .\n");
		return(0);
		  N_VDestroy_Serial(y);
		  /* Free integrator memory */
		  CVodeFree(&cvode_mem);
	 }
	 if(verbose)printf("Max number of steps: %i\n", maxNumSteps);

	  /* Set maxstepsize */
	  flag = CVodeSetMaxStep(cvode_mem, maxStepSize);
	  if(check_flag(&flag, "CVodeSetMaxStep", 1))
	  {
		  if(verbose)printf("\nSolver failed. . .\n");
		  N_VDestroy_Serial(y);
		  /* Free integrator memory */
		  CVodeFree(&cvode_mem);
		  return(0);
	  }
	  if(verbose)printf("Max step size: %f\n", maxStepSize);

	  for (i = 1; i < (*data).nTimes; ++i)
	  {
		  tout=(*data).timeSignals[i];
		  flag = CVode(cvode_mem, tout, y, &tf, CV_NORMAL);

		  for (j = 0; j < (*data).nStates; j++)
		  {
			//  (*data).sim_results[exp_num][i][j]= (double) Ith(y,j);
			  (*data).sim_results[exp_num][i][j]= (double)Ith(y,j);
			if(verbose)printf("%f\t",Ith(y,j));
		  }

		  if (check_flag(&flag, "CVode", 1))
		  {
			  if(verbose)fprintf(stderr,"\nSolver failed. . .\n");
			  break;
		  }
		  if(verbose)printf("\n");
	  }
	  if(verbose)printf("\n");

	  N_VDestroy_Serial(y);
	  /* Free integrator memory */
	  CVodeFree(&cvode_mem);

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



