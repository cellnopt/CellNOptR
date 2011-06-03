\name{ida} 
\alias{ida} 
\title{SUNDIALS IDA Solver for Differential-Algebraic Equation Systems} 
\description{ 
	This function allows for evaluation of DAE systems using the SUNDIALS IDA solver hard coded in C.
}
\usage{ 
	ida(y, yp, times, package, resfunc, jacfunc = NULL, rootfunc = NULL, fndata = NULL,
		numroots = 0, rtol = 1e-6, atol = 1e-6, maxnumsteps = 500, maxstep = 0, verbose = FALSE, 
		lasttime = FALSE)
} 
\arguments{ 
	\item{y}{a numeric vector of the initial values of y of the DAE system.}
	\item{yp}{a numeric vector of the initial values of y' of the DAE system.}
	\item{times}{a numeric vector of times at which explicit estimates for y are desired. The first element should be the initial time (often 0). The IDA solver takes steps until the next requested time step t is overtaken. It then computes y(t) by interpolation.

	To easily create geometric sequences, the function gseq can be used.}
	\item{package}{a character vector of the package name containing the user-defined residual function. If the residual function is contained within the file "usrfcns.c" for example, the command R CMD SHLIB usrfcns.c will create the package "usrfcns" which can be loaded using the command dyn.load("usrfcns.so") or dyn.load("usrfcns.dll") depending on your system.}
	\item{resfunc}{a character vector defining name of the user-defined residual function contained within the package specified by the \code{package} variable. This need not match the name of the package itself.}
	\item{jacfunc}{a character vector defining the name of the user-defined Jacobian (if desired). By default, the IDA solver will approximate the Jacobian by difference quotients.}
	\item{rootfunc}{a character vector defining the name of the user-defined method (in C) for computing the roots of some set of functions. The number of functions is set by /code{numroots}.}
	\item{fndata}{a numeric vector of data to be used by the user-defined functions.}
	\item{numroots}{the number of functions defined by rootfunc.}
	\item{rtol}{a scalar that defines the relative error tolerance. The relative error tolerance is not recommended to exceed 1.0E-3. On the other hand, \code{reltol} should not be so small that it is comparable to the unit roundoff of the machine arithmetic (generally around 1.0E-15).}
	\item{atol}{a scalar or vector that defines the absolute error tolerance(s). The absolute tolerances need to be set to control absolute errors when any components of the solution vector \code{y} may be so small that pure relative error control is meaningless.}
	\item{maxnumsteps}{a scalar that sets the maximum number of steps per timestep. The default is 500.}
	\item{maxstep}{a scalar that sets the maximum absolute step size. The default is 0 = infinity.}
	\item{verbose}{a logical value that, when TRUE, triggers verbose output from the solver.}
	\item{lasttime}{by default, with lasttime = FALSE, a data frame is returned that includes all of the values y at each time point. If lasttime is set to TRUE,  only the last values of y will be returned (i.e. at the last time point specified).}
}
\examples{
### See Manual for description and examples
}
\keyword{file} 
