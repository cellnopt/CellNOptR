\name{getLBodeSimFunction}

\alias{getLBodeSimFunction}

\title{
	Get a function to simulate a logic based ODE model.
}

\description{
	This function is internally used by CNORode to configure the simulation function with default arguments.
}

\usage{
	getLBodeSimFunction(cnolist1, model1, adjMatrix1, indices1, odeParameters1,
	time1 = 1, verbose1 = 0, transfer_function1 = 3, reltol1 = 1e-04, atol1 = 0.001,
	maxStepSize1 = Inf, maxNumSteps1 = 1e+05, maxErrTestsFails1 = 50)
}

%- maybe also 'usage' for other objects documented here.
\arguments{
 	\item{cnolist1}{
		A list containing the experimental design and data.
	}
	\item{model1}{
		The logic model to be simulated.
	}
 	\item{adjMatrix1}{
		An adjacency matrix from the model.
	}
	\item{indices1}{
		Indices to map data in the model. Obtained with indexFinder function from CellNOptR. 
	}
	\item{odeParameters1}{
		A list with the ODEs parameter information. Obtained with \code{\link{createLBodeContPars}}.
	}	
	\item{time1}{
		An integer with the index of the time point to start the simulation. Default is 1.
	}
	\item{verbose1}{
		A logical value that triggers a set of comments.
	}
	\item{transfer_function1}{
		The type of used transfer. Use 1 for no transfer function, 2 for Hill function and
		3 for normalized Hill function.
	}
	\item{reltol1}{
		Relative Tolerance for numerical integration.
	}
	\item{atol1}{
		Absolute tolerance for numerical integration.
	}
	\item{maxStepSize1}{
		The maximum step size allowed to ODE solver.
	}
	\item{maxNumSteps1}{
		The maximum number of internal steps between two points being sampled before
		the solver fails.
	}
	\item{maxErrTestsFails1}{
		Specifies the maximum number of error test failures permitted in attempting
		one step.
	}
}

\value{
	A function that returns a simulated model.
}

\author{
	David Henriques, Thomas Cokelaer
}

\note{
	This function is for CNORode internal use.
}

\seealso{
	\code{\link{CellNOptR}}
	\code{\link{CNORode}}
}

% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{CVODES}
