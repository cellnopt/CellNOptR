#include <stdio.h>
#include <math.h>
#include <R.h>
#include <Rinternals.h>
#include <R_ext/Print.h>

#include "include/ida/ida.h"
#include "include/ida/ida_dense.h"
#include "include/nvector/nvector_serial.h"
#include "include/sundials/sundials_types.h"
#include "include/sundials/sundials_math.h"

/* Need number of equations and number or roots to be global */

realtype neq;

/* Prototypes of functions called by IDA */

typedef int resrob_func(realtype tres, N_Vector yy, N_Vector yp, 
           N_Vector resval, void *rdata);

typedef int jacrob_func(long int Neq, realtype tt, N_Vector yy, N_Vector yp,
           N_Vector resvec, realtype cj, void *jdata, DenseMat JJ,
           N_Vector tempv1, N_Vector tempv2, N_Vector tempv3);
		   
typedef int root_func(realtype t, N_Vector yy, N_Vector yp, realtype *gout,
                void *g_data);

/* Prototypes of private functions */
static void PrintHeader(int neq);
static void PrintOutput(void *mem, realtype t, N_Vector y, int neq);
static void PrintRootInfo(int *roots);
static void RprintfinalStats(void *mem);
static int check_flag(void *flagvalue, char *funcname, int opt);

/*
 *--------------------------------------------------------------------
 * Main Program
 *--------------------------------------------------------------------
 */
 
SEXP ida(SEXP r_y, SEXP r_yp, SEXP r_times, SEXP r_resfunc, SEXP r_jacfunc, SEXP r_rootfunc, SEXP r_data,
		 SEXP r_numroots, SEXP r_rtol, SEXP r_atol, SEXP r_maxnumsteps, SEXP r_maxstep, SEXP r_verbose, SEXP r_lasttime) {
		 
	/* Badthing */
	SEXP badthing;
	PROTECT (badthing = allocVector(REALSXP, 1));
	REAL(badthing)[0] = 0;
	UNPROTECT(1);
	
	
	/* IDA Data */
	void *mem;
	double *data; data = NULL; data = (double *) malloc(LENGTH(r_data) * sizeof (double));
	N_Vector yy, yp, avtol;
	realtype rtol, *yval, *ypval, *atol;
	realtype ti, tf, tout1, tout, tret;
	int retval, retvalr;
	int i = 0,  j = 0, tcount = 0;
	int verbose = INTEGER(r_verbose)[0];
	int lasttime = INTEGER(r_lasttime)[0];
	// Root info
	int numroots = INTEGER(r_numroots)[0];
	if (numroots == 0) numroots = 1;	// Otherwise next line will complain
	int rootsfound[numroots];
	
	// User residual function
	resrob_func *resrob = (resrob_func *) R_ExternalPtrAddr(r_resfunc);
	

	if (verbose) Rprintf("SUNDIALS IDADENSE Linear Solver\n");
	
	/* Problem Constants */
	neq = RCONST(length(r_y));
	ti = RCONST(REAL(r_times)[0]); tcount++; tout = RCONST(REAL(r_times)[tcount]);
	
	if (verbose) {
		Rprintf("Number of Equations: %lg \n", neq);
		Rprintf("Integration Limits: %lg to ", ti); Rprintf("%lg \n", REAL(r_times)[LENGTH(r_times) - 1]);
	}
	
	
	mem = NULL;
	yy = yp = avtol = NULL;
	yval = ypval = NULL;


	/* Allocate N-vectors. */
	yy = N_VNew_Serial(neq);
	if(check_flag((void *)yy, "N_VNew_Serial", 0)) {
		Rprintf("\nSolver failed. . .\n");
		return(badthing);
	}
	yp = N_VNew_Serial(neq);
	if(check_flag((void *)yp, "N_VNew_Serial", 0)) {
		Rprintf("\nSolver failed. . .\n");
		return(badthing);
	}
	avtol = N_VNew_Serial(neq);
	if(check_flag((void *)avtol, "N_VNew_Serial", 0)) {
		Rprintf("\nSolver failed. . .\n");
		return(badthing);
	}
	
	
	/* Create and initialize  y, y' */
	yval  = NV_DATA_S(yy);
	for (i=0; i<neq; i++) { yval[i] = RCONST(REAL(r_y)[i]); }
	
	ypval = NV_DATA_S(yp);
	for (i=0; i<neq; i++) { ypval[i] = RCONST(REAL(r_yp)[i]); }


	mem = IDACreate();
	if(check_flag((void *)mem, "IDACreate", 0)) {
		Rprintf("\nSolver failed. . .\n");
		return(badthing);
	}
	
	/* Set f_data */
	if (!isNull(r_data)) {
		for (i=0; i<LENGTH(r_data); i++) data[i] = REAL(r_data)[i];
		retval = IDASetRdata(mem, data);
		if(check_flag(&retval, "IDASetRdata", 1)) {
			Rprintf("\nSolver failed. . .\n");
			return(badthing);
		}
	}
	
	if (verbose) Rprintf("Solver Memory Allocated\n");
	
	/* Set maxnumsteps */
	retval = IDASetMaxNumSteps(mem, INTEGER(r_maxnumsteps)[0]);
	if(check_flag(&retval, "IDASetMaxNumSteps", 1)) {
		Rprintf("\nSolver failed. . .\n");
		return(badthing);
	} else if (verbose) Rprintf("Max Number of Steps: %i\n", INTEGER(r_maxnumsteps)[0]);
	
	/* Set maxstep */
	retval = IDASetMaxStep(mem, INTEGER(r_maxstep)[0]);
	if(check_flag(&retval, "IDASetMaxStep", 1)) {
		Rprintf("\nSolver failed. . .\n");
		return(badthing);
	} else if (verbose) Rprintf("Max step size: %i\n", INTEGER(r_maxstep)[0]);
	
	
	/* Initialize relative and absolute tolerances */
	rtol = RCONST(REAL(r_rtol)[0]);
	
	if (verbose) Rprintf("Relative Tolerance: %lg \n", rtol);

	if (LENGTH(r_atol) == 1) {
		realtype atol;
		atol = RCONST(REAL(r_atol)[0]);
		realtype *atol2 = &atol;	
		retval = IDAMalloc(mem, resrob, ti, yy, yp, IDA_SS, rtol, atol2);
		
		if(check_flag(&retval, "IDAMalloc", 1)) {
			Rprintf("\nSolver failed. . .\n");
			return(badthing);
		}
		
		if (verbose) Rprintf("Absolute Tolerance: %lg\n",atol);
	} else {
		atol = NV_DATA_S(avtol);
		for (i=0; i<neq; i++) atol[i] = REAL(r_atol)[i]; 
		retval = IDAMalloc(mem, resrob, ti, yy, yp, IDA_SV, rtol, avtol);
		
		if(check_flag(&retval, "IDAMalloc", 1)) {
			Rprintf("\nSolver failed. . .\n");
			return(badthing);
		}
		
		if (verbose) {
			Rprintf("Absolute Tolerances: ");
			for (i=0; i<LENGTH(r_atol); i++) Rprintf("%lg ",atol[i]);
			Rprintf("\n");
		}
		
		N_VDestroy_Serial(avtol);
	}

	/* Call IDARootInit to specify the root function grob with 2 components */
	if (!isNull(r_rootfunc)) {
		root_func *rfunc = (root_func *) R_ExternalPtrAddr(r_rootfunc);
		retval = IDARootInit(mem, numroots, rfunc, data);
		if(check_flag(&retval, "IDARootInit", 1)) {
			Rprintf("\nSolver failed. . .\n");
			return(badthing);
		}
		if (verbose) Rprintf("Root Function Initialized for %i functions\n", numroots);
	}

	/* Call IDADense and set up the linear solver. */
	retval = IDADense(mem, neq);
	if(check_flag(&retval, "IDADense", 1)) {
		Rprintf("\nSolver failed. . .\n");
		return(badthing);
	}
	
	/* Set up the Jacobian */
	if (!isNull(r_jacfunc)) {
		jacrob_func *jacrob = (jacrob_func *) R_ExternalPtrAddr(r_jacfunc);
		retval = IDADenseSetJacFn(mem, jacrob, data);
		if(check_flag(&retval, "IDADenseSetJacFn", 1)) {
			Rprintf("\nSolver failed. . .\n");
			return(badthing);
		}
		if (verbose) Rprintf("Jacobian Function Initialized\n");
	}
	
	if (verbose) Rprintf("IDADense Linear Solver Initialized\n");
	if (verbose) PrintHeader(neq);
	
	/* Set up an array for lasttime */
	int x = LENGTH(r_times);
	int y = neq;
	double yvals[x][y + 1];
	yvals[0][0] = REAL(r_times)[0];
	for (i=0; i<neq; i++) yvals[0][i + 1] = REAL(r_y)[i];
	
	/* In loop, call IDASolve, print results, and test for error.
	   Break out of loop when NOUT preset output times have been reached. */
	while(tcount < LENGTH(r_times)) {
		retval = IDASolve(mem, tout, &tret, yy, yp, IDA_NORMAL);

		if (verbose) PrintOutput(mem,tret,yy,neq);
		if (!lasttime) {
				realtype *yval;
				yval = NV_DATA_S(yy);
				//yvals[tcount][0] = REAL(r_times)[tcount];
				for (i=0; i<neq; i++) yvals[tcount][i] = yval[i];
		}

		if(check_flag(&retval, "IDASolve", 1)) {
			Rprintf("\nSolver failed. . .\n");
			return(badthing);
		}

		if (retval == IDA_ROOT_RETURN) {
			retvalr = IDAGetRootInfo(mem, rootsfound);
			if(check_flag(&retval, "IDASolve", 1)) {
				Rprintf("\nSolver failed. . .\n");
				return(badthing);
			}
			PrintRootInfo(rootsfound);
		}
	
		if (retval == IDA_SUCCESS) {
			tcount++;
			tout = REAL(r_times)[tcount];
		}
  }
  

  if (verbose) RprintfinalStats(mem);
  
  /* Save solutions */
  SEXP solutions;
  if (!lasttime) {
	PROTECT (solutions = allocMatrix(REALSXP, LENGTH(r_times), neq));
	for (i=0; i<LENGTH(r_times); i++) {
		for (j=0; j<neq; j++) {
			REAL(solutions)[i + (LENGTH(r_times))*j] = yvals[i][j];
		}
	}
	UNPROTECT(1);
  } else {
	realtype *y;
	y = NV_DATA_S(yy);
	
	PROTECT (solutions = allocVector(REALSXP, neq));
	for (i=0; i<neq; i++) REAL(solutions)[i] =y[i];
	UNPROTECT(1);
  }

  /* Free memory */
  IDAFree(&mem);
  N_VDestroy_Serial(yy);
  N_VDestroy_Serial(yp);

  return (solutions);
}


/*
 *--------------------------------------------------------------------
 * Private functions
 *--------------------------------------------------------------------
 */

/* 
 * Print first lines of output (problem description)
 */

static void PrintHeader(int neq)
{
  int i = 0;
  Rprintf("-----------------------------------------------------------------------\n");
  Rprintf("  t             ");
  for (i=0; i<neq; i++) {
	Rprintf("y%d",i);
	if (i<neq - 1) Rprintf("           "); else Rprintf("      ");
  }
  Rprintf("| nst  k      h\n");
  Rprintf("-----------------------------------------------------------------------\n");
}


/*
 * Save Times
 */
void SaveTimes(void *mem, realtype t, N_Vector y) {
	
	
}

/*
 * Print Output
 */
static void PrintOutput(void *mem, realtype t, N_Vector y, int neq)
{
  realtype *yval;
  int retval, kused;
  long int nst;
  realtype hused;
  int i = 0;

  yval  = NV_DATA_S(y);

  retval = IDAGetLastOrder(mem, &kused);
  check_flag(&retval, "IDAGetLastOrder", 1);
  retval = IDAGetNumSteps(mem, &nst);
  check_flag(&retval, "IDAGetNumSteps", 1);
  retval = IDAGetLastStep(mem, &hused);
  check_flag(&retval, "IDAGetLastStep", 1);
  
#if defined(SUNDIALS_EXTENDED_PRECISION)
  Rprintf("%10.4le ", t);
  for (i=0; i<neq; i++) Rprintf("%12.4le ", yval[i]);
  Rprintf(" | %3ld  %1d %12.4le\n", nst, kused, hused);
#elif defined(SUNDIALS_DOUBLE_PRECISION)
  Rprintf("%10.4le ", t);
  for (i=0; i<neq; i++) Rprintf("%12.4le ", yval[i]);
  Rprintf(" | %3ld  %1d %12.4le\n", nst, kused, hused);
#else
  Rprintf("%10.4le ", t);
  for (i=0; i<neq; i++) Rprintf("%12.4le ", yval[i]);
  Rprintf(" | %3ld  %1d %12.4le\n", nst, kused, hused);
#endif
}

static void PrintRootInfo(int *roots)
{
  int count = 0;
  while (roots[count] == 0 || roots[count] == 1) {
	if (roots[count] == 1) Rprintf("    A root was found for function %i\n",count);
	count++;
  }
  
  return;
}

/*
 * Print final integrator statistics
 */

static void RprintfinalStats(void *mem)
{
  int retval;
  long int nst, nni, nje, nre, nreLS, netf, ncfn, nge;

  retval = IDAGetNumSteps(mem, &nst);
  check_flag(&retval, "IDAGetNumSteps", 1);
  retval = IDAGetNumResEvals(mem, &nre);
  check_flag(&retval, "IDAGetNumResEvals", 1);
  retval = IDADenseGetNumJacEvals(mem, &nje);
  check_flag(&retval, "IDADenseGetNumJacEvals", 1);
  retval = IDAGetNumNonlinSolvIters(mem, &nni);
  check_flag(&retval, "IDAGetNumNonlinSolvIters", 1);
  retval = IDAGetNumErrTestFails(mem, &netf);
  check_flag(&retval, "IDAGetNumErrTestFails", 1);
  retval = IDAGetNumNonlinSolvConvFails(mem, &ncfn);
  check_flag(&retval, "IDAGetNumNonlinSolvConvFails", 1);
  retval = IDADenseGetNumResEvals(mem, &nreLS);
  check_flag(&retval, "IDADenseGetNumResEvals", 1);
  retval = IDAGetNumGEvals(mem, &nge);
  check_flag(&retval, "IDAGetNumGEvals", 1);

  Rprintf("\nFinal Run Statistics: \n\n");
  Rprintf("Number of steps                    = %ld\n", nst);
  Rprintf("Number of residual evaluations     = %ld\n", nre+nreLS);
  Rprintf("Number of Jacobian evaluations     = %ld\n", nje);
  Rprintf("Number of nonlinear iterations     = %ld\n", nni);
  Rprintf("Number of error test failures      = %ld\n", netf);
  Rprintf("Number of nonlinear conv. failures = %ld\n", ncfn);
  Rprintf("Number of root fn. evaluations     = %ld\n", nge);
}

/*
 * Check function return value...
 *   opt == 0 means SUNDIALS function allocates memory so check if
 *            returned NULL pointer
 *   opt == 1 means SUNDIALS function returns a flag so check if
 *            flag >= 0
 *   opt == 2 means function allocates memory so check if returned
 *            NULL pointer 
 */

static int check_flag(void *flagvalue, char *funcname, int opt)
{
  int *errflag;
  /* Check if SUNDIALS function returned NULL pointer - no memory allocated */
  if (opt == 0 && flagvalue == NULL) {
    Rprintf("\nSUNDIALS_ERROR: %s() failed - returned NULL pointer\n\n", funcname);
    return(1);
  } else if (opt == 1) {
    /* Check if flag < 0 */
    errflag = (int *) flagvalue;
    if (*errflag < 0) {
      Rprintf("\nSUNDIALS_ERROR: %s() failed with flag = %d\n\n", funcname, *errflag);
      return(1); 
    }
  } else if (opt == 2 && flagvalue == NULL) {
    /* Check if function returned NULL pointer - no memory allocated */
    Rprintf("\nMEMORY_ERROR: %s() failed - returned NULL pointer\n\n",funcname);
    return(1);
  }

  return(0);
}
