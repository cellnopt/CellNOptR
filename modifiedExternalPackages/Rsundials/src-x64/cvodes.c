#include <stdio.h>
#include <math.h>
#include <string.h>

/* R Header Files */


#include <R.h>
#include <Rinternals.h>
#include <R_ext/Print.h>

/* Sundials Header Files */

#include "include/cvodes/cvodes.h"           /* prototypes for CVODES fcts. and consts. */
#include "include/cvodes/cvodes_dense.h"     /* prototype for CVDense */
#include "include/cvodes/cvodes_band.h"
#include "include/cvodes/cvodes_diag.h"
#include "include/nvector/nvector_serial.h"  /* serial N_Vector types, fcts., and macros */
#include "include/sundials/sundials_dense.h" /* definitions DenseMat and DENSE_ELEM */
#include "include/sundials/sundials_types.h" /* definition of type realtype */

/* User-defined vector and matrix accessor macros: Ith, IJth */

#define Ith(v,i) ( NV_DATA_S(v)[i] )
#define IJth(A,i,j) DENSE_ELEM(A,i,j)

/* Functions Called by the Solver */

typedef int rhs_func(realtype t, N_Vector y, N_Vector ydot, void *f_data);

typedef int root_func(realtype t, N_Vector y, realtype *gout, void *g_data);

typedef int dense_jac_func(long int N, DenseMat J, realtype t,
               N_Vector y, N_Vector fy, void *jac_data,
               N_Vector tmp1, N_Vector tmp2, N_Vector tmp3);
			   
typedef int band_jac_func(long int N, long int mu, long int ml, BandMat J,
               realtype t, N_Vector u, N_Vector fu, void *jac_data,
               N_Vector tmp1, N_Vector tmp2, N_Vector tmp3);

/* Private functions to output results */

static void PrintOutput(realtype t, N_Vector y, int neq);
static void PrintRootInfo(int *roots);
static void PrintOutputS(N_Vector *uS, int num);

/* Private function to print final statistics */

static void PrintFinalStats(void *cvode_mem);

/* Private function to check function return values */

static int check_flag(void *flagvalue, char *funcname, int opt);


/*
 *-------------------------------
 * Main Program
 *-------------------------------
 */

SEXP cvodes(SEXP r_y, SEXP r_times, SEXP r_rhs, SEXP r_data, SEXP r_jacfunc, SEXP r_rootfunc, SEXP r_numroots, SEXP r_solver, SEXP r_rtol, SEXP r_atol, SEXP r_maxnumsteps, SEXP r_maxstep, 
			SEXP r_verbose, SEXP r_lasttime) {		
		 
	/* Badthing */
	SEXP badthing;
	PROTECT (badthing = allocVector(REALSXP, 1));
	REAL(badthing)[0] = 0;
	UNPROTECT(1);
	
	
	/* CVODES Data */
	realtype reltol, t, tout, ti, tf;
	realtype *atol;
	N_Vector y, abstol;
	void *cvode_mem;
	double *data; data = NULL; data = (double *) malloc(LENGTH(r_data) * sizeof (double));
	int flag, flagr, iout, neq;
	int i = 0, j = 0, tcount = 0;
	// Booleans
	int verbose = INTEGER(r_verbose)[0];
	int lasttime = INTEGER(r_lasttime)[0];
	// Root info
	int numroots = INTEGER(r_numroots)[0];
	if (numroots == 0) numroots = 1;	// Otherwise next line will complain
	int rootsfound[numroots];

	// User RHS Function
	rhs_func *rhs = (rhs_func *) R_ExternalPtrAddr(r_rhs);
	
	y = NULL;
	cvode_mem = NULL;
	abstol = NULL;
	
	if (verbose) Rprintf("SUNDIALS CVODES Linear Solver\n");

	/* Problem Constants */
	neq = LENGTH(r_y);
	//if (REAL(r_times)[tcount] == 0) tcount++; ti = RCONST(REAL(r_times)[tcount]);
	ti = REAL(r_times)[tcount]; tcount++; 	tout = REAL(r_times)[tcount];
	
	if (verbose) {
		Rprintf("Number of Equations: %i \n", neq);
		Rprintf("Integration Limits: %lg to ", ti); Rprintf("%lg \n", REAL(r_times)[LENGTH(r_times) - 1]);
	}

	/* Create serial vector of length NEQ for I.C. and abstol */
	y = N_VNew_Serial(neq);
	if (check_flag((void *)y, "N_VNew_Serial", 0)) {
		Rprintf("\nSolver failed. . .\n");
		return(badthing);
	}
  
	abstol = N_VNew_Serial(neq); 
	if (check_flag((void *)abstol, "N_VNew_Serial", 0)) {
		Rprintf("\nSolver failed. . .\n");
		return(badthing);
	}

	/* Initialize y */
	for(i=0; i<neq; i++) Ith(y,i) = REAL(r_y)[i];
	
	if (verbose) {
		Rprintf("Initial Values:  ");
		for(i=0; i<neq; i++) Rprintf("y%i = %lg   ",i,Ith(y,i));
		Rprintf("\n");
	}
	
  /* 
     Call CVodeCreate to create the solver memory:
     
     CV_BDF     specifies the Backward Differentiation Formula (or CV_ADAMS)
     CV_NEWTON  specifies a Newton iteration (or CV_FUNCTIONAL)

     A pointer to the integrator problem memory is returned and stored in cvode_mem.
  */	
	
	cvode_mem = CVodeCreate(CV_BDF, CV_NEWTON);
	if (check_flag((void *)cvode_mem, "CVodeCreate", 0)) {
		Rprintf("\nSolver failed. . .\n");
		return(badthing);
	}
	if (verbose) Rprintf("Solver Memory Allocated\n");
	
	/* Set f_data */
	if (!isNull(r_data)) {
		for (i=0; i<LENGTH(r_data); i++) data[i] = REAL(r_data)[i];
		flag = CVodeSetFdata(cvode_mem, data);
		if(check_flag(&flag, "CVodeSetFdata", 1)) {
			Rprintf("\nSolver failed. . .\n");
			return(badthing);
		}
	}
	
	// Set and initialize error tolerances
	reltol = RCONST(REAL(r_rtol)[0]);

	if (verbose) Rprintf("Relative Tolerance: %lg \n", reltol);
	
	if (LENGTH(r_atol) == 1) {			// Scalar absolute tolerance
		realtype atol;
		atol = RCONST(REAL(r_atol)[0]);
		realtype *atol2 = &atol;	
		flag = CVodeMalloc(cvode_mem, rhs, ti, y, CV_SS, reltol, atol2);
		
		if (check_flag(&flag, "CVodeMalloc", 1)) {
			Rprintf("\nSolver failed. . .\n");
			return(badthing);
		}
		
		if (verbose) Rprintf("Absolute Tolerance: %lg\n",atol);
	} else {							// Vector absolute tolerance
		for(i=0; i<neq; i++) Ith(abstol,i) = REAL(r_atol)[i];
		flag = CVodeMalloc(cvode_mem, rhs, ti, y, CV_SV, reltol, abstol);
		if (check_flag(&flag, "CVodeMalloc", 1)) {
			Rprintf("\nSolver failed. . .\n");
			return(badthing);
		}
		if (verbose) {
			Rprintf("Absolute Tolerances: ");
			for (i=0; i<LENGTH(r_atol); i++) Rprintf("%lg   ",Ith(abstol,i));
			Rprintf("\n");
		}
		N_VDestroy_Serial(abstol);
	}
	
  
  /* 
     Call CVodeMalloc to initialize the integrator memory: 
     
     cvode_mem is the pointer to the integrator memory returned by CVodeCreate
     f         is the user's right hand side function in y'=f(t,y)
     T0        is the initial time
     y         is the initial dependent variable vector
     CV_SV     specifies scalar relative and vector absolute tolerances
     &reltol   is a pointer to the scalar relative tolerance
     abstol    is the absolute tolerance vector
  */


  /* Call CVodeRootInit to specify the root function g with 2 components */
	if (!isNull(r_rootfunc)) {
		root_func *rfunc = (root_func*) R_ExternalPtrAddr(r_rootfunc);
		flag = CVodeRootInit(cvode_mem, numroots, rfunc, data);
		if(check_flag(&flag, "CVodeRootInit", 1)) {
			Rprintf("\nSolver failed. . .\n");
			return(badthing);
		}
		if (verbose) Rprintf("Root Function Initialized for %i functions\n", numroots);
	}
	
  /* Specify the CVDENSE dense linear solver */
	  if (INTEGER(r_solver)[0] == 1) {
		flag = CVDense(cvode_mem, neq);
		if (check_flag(&flag, "CVDense", 1)) {
			Rprintf("\nSolver failed. . .\n");
			return(badthing);
		}
		
	    
	   // Set the Jacobian routine to Jac (user-supplied)
	    if (!isNull(r_jacfunc)) {
			dense_jac_func *jacrob = (dense_jac_func *) R_ExternalPtrAddr(r_jacfunc);
			flag = CVDenseSetJacFn(cvode_mem, jacrob, data);
			if (check_flag(&flag, "CVDenseSetJacFn", 1)) {
				Rprintf("\nSolver failed. . .\n");
				return(badthing);
			}
		}
		if (verbose) Rprintf("CVDENSE Solver Initiated\n");
	  } /* else if (INTEGER(r_solver)[0] == 2) {
		flag = CVBand(cvode_mem, neq, INTEGER(r_mlower)[0], INTEGER(r_mupper)[0]);
		if (check_flag(&flag, "CVBand", 1)) {
			Rprintf("\nSolver failed. . .\n");
			return(badthing);
		}	

	   // Set the Jacobian routine to Jac (user-supplied)
		if (!isNull(r_jacfunc)) {
			band_jac_func *jacrob = (band_jac_func *) R_ExternalPtrAddr(r_jacfunc);
			flag = CVBandSetJacFn(cvode_mem, jacrob, NULL);
			if (check_flag(&flag, "CVDenseSetJacFn", 1)) {
				Rprintf("\nSolver failed. . .\n");
				return(badthing);
			}
		}
		if (verbose) Rprintf("CVBAND Solver Initiated\n");
	} else if (INTEGER(r_solver)[0] == 3) {
		flag = CVDiag(cvode_mem);
		if (check_flag(&flag, "CVDiag", 1)) {
			Rprintf("\nSolver failed. . .\n");
			return(badthing);
		}
		if (verbose) Rprintf("CVDIAG Solver Initiated\n");
	}
	*/
	
	
	/* Set maxnumsteps */
	flag = CVodeSetMaxNumSteps(cvode_mem, INTEGER(r_maxnumsteps)[0]);
	if(check_flag(&flag, "CVodeSetMaxNumSteps", 1)) {
		Rprintf("\nSolver failed. . .\n");
		return(badthing);
	} else if (verbose) Rprintf("Max number of steps: %i\n", INTEGER(r_maxnumsteps)[0]);
	
	/* Set maxstepsize */
	flag = CVodeSetMaxStep(cvode_mem, INTEGER(r_maxstep)[0]);
	if(check_flag(&flag, "CVodeSetMaxStep", 1)) {
		Rprintf("\nSolver failed. . .\n");
		return(badthing);
	} else if (verbose) Rprintf("Max step size: %i\n", INTEGER(r_maxstep)[0]);
	
	/* Set up an array for lasttime */
	int a = LENGTH(r_times);
	int b = neq;
	double yvals[a][b + 1];
	yvals[0][0] = REAL(r_times)[0];
	for (i=0; i<neq; i++) yvals[0][i + 1] = REAL(r_y)[i];
	if (verbose) Rprintf("Requesting data for all time points.\n");
	
	
  /* In loop, call CVode, print results, and test for error.
     Break out of loop when NOUT preset output times have been reached.  */
	/*Rprintf("\n");*/
	while(tcount < LENGTH(r_times)) {
		flag = CVode(cvode_mem, tout, y, &t, CV_NORMAL);
		if (verbose) PrintOutput(t, y, neq);
		
		if (!lasttime) {
			//yvals[tcount][0] = REAL(r_times)[tcount]; // this is for adding timestamps which is a bad idea
			for (i=0; i<neq; i++) yvals[tcount][i] = (double) Ith(y,i);
		}
	
		
		

		if (flag == CV_ROOT_RETURN) {
			flagr = CVodeGetRootInfo(cvode_mem, rootsfound);
			if (check_flag(&flagr, "CVodeGetRootInfo", 1)) {
				Rprintf("\nSolver failed. . .\n");
				return(badthing);
			}
			PrintRootInfo(rootsfound);
		}

		if (check_flag(&flag, "CVode", 1)) {
			Rprintf("\nSolver failed. . .\n");
			break;
		}
			
		if (flag == CV_SUCCESS) {
			tcount++;
			tout = REAL(r_times)[tcount];
		}
	}

	/* Print some final statistics */
	if (verbose) PrintFinalStats(cvode_mem);

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
		PROTECT (solutions = allocVector(REALSXP, neq));
		for (i=0; i<neq; i++) REAL(solutions)[i] = Ith(y,i);
		UNPROTECT(1);
  }

	/* Free y vector */
	N_VDestroy_Serial(y);
	/* Free integrator memory */
	CVodeFree(&cvode_mem);
	/* free userdata */
	free(data);

	

	return(solutions);
}


/*
 *-------------------------------
 * Private helper functions
 *-------------------------------
 */

static void PrintOutput(realtype t, N_Vector y, int neq)
{
 int i = 0;

#if defined(SUNDIALS_EXTENDED_PRECISION)
  Rprintf("At t = %0.4le      ", t);
  for (i=0; i<neq; i++) Rprintf("%14.6le  ", Ith(y,i));
  Rprintf("\n");
#elif defined(SUNDIALS_DOUBLE_PRECISION)
  Rprintf("At t = %0.4le      ", t);
  for (i=0; i<neq; i++) Rprintf("%14.6le  ", Ith(y,i));
  Rprintf("\n");
#else
  Rprintf("At t = %0.4le      ", t);
  for (i=0; i<neq; i++) Rprintf("%14.6le  ", Ith(y,i));
  Rprintf("\n");
#endif

  return;
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
 * Get and print some final statistics
 */

static void PrintFinalStats(void *cvode_mem)
{
  long int nst, nfe, nsetups, nje, nfeLS, nni, ncfn, netf, nge;
  int flag;

  flag = CVodeGetNumSteps(cvode_mem, &nst);
  check_flag(&flag, "CVodeGetNumSteps", 1);
  flag = CVodeGetNumRhsEvals(cvode_mem, &nfe);
  check_flag(&flag, "CVodeGetNumRhsEvals", 1);
  flag = CVodeGetNumLinSolvSetups(cvode_mem, &nsetups);
  check_flag(&flag, "CVodeGetNumLinSolvSetups", 1);
  flag = CVodeGetNumErrTestFails(cvode_mem, &netf);
  check_flag(&flag, "CVodeGetNumErrTestFails", 1);
  flag = CVodeGetNumNonlinSolvIters(cvode_mem, &nni);
  check_flag(&flag, "CVodeGetNumNonlinSolvIters", 1);
  flag = CVodeGetNumNonlinSolvConvFails(cvode_mem, &ncfn);
  check_flag(&flag, "CVodeGetNumNonlinSolvConvFails", 1);

  //flag = CVDenseGetNumJacEvals(cvode_mem, &nje);
  //check_flag(&flag, "CVDenseGetNumJacEvals", 1);
  //flag = CVDenseGetNumRhsEvals(cvode_mem, &nfeLS);
  //check_flag(&flag, "CVDenseGetNumRhsEvals", 1);

  flag = CVodeGetNumGEvals(cvode_mem, &nge);
  check_flag(&flag, "CVodeGetNumGEvals", 1);

  Rprintf("\nFinal Run Statistics: \n\n");
  Rprintf("Number of steps                    = %ld\n", nst);
  Rprintf("Number of RHS evaluations          = %ld\n", nfe);
  Rprintf("Number of linear solver setups     = %ld\n", nsetups);
  Rprintf("Number of nonlinear iterations     = %ld\n", nni);
  Rprintf("Number of error test failures      = %ld\n", netf);
  Rprintf("Number of nonlinear conv. failures = %ld\n", ncfn);
  Rprintf("Number of root fn. evaluations     = %ld\n", nge);
}

static void PrintOutputS(N_Vector *uS, int num) {
  realtype *sdata;
  int i = 0; int j = 0;
	
  for (i=0; i<num; i++) {
	sdata = NV_DATA_S(uS[i]);
	Rprintf("                  Sensitivity %i  ",i);

	#if defined(SUNDIALS_EXTENDED_PRECISION)
	 for(j=0; j<num;j++) Rprintf("%12.4le ",sdata[j]);
	 Rprintf("\n"); 
	#elif defined(SUNDIALS_DOUBLE_PRECISION)
	 for(j=0; j<num;j++) Rprintf("%12.4le ",sdata[j]);
	 Rprintf("\n"); 
	#else
	 for(j=0; j<num;j++) Rprintf("%12.4le ",sdata[j]);
	 Rprintf("\n"); 
	#endif
	}
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
