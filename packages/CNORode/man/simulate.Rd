\name{simulate}

\alias{simulate}

\title{
	Simulate value signals a CNO list With Logic-Based ODEs.
}

\description{
	This function receives a set of inputs, namely the cnolist and the model and
	returns a list with the same size of the cnolist$valueSignals. 
}


\usage{
	simulate(cnolist, model, ode_parameters=NULL, indices=NULL,
	adjMatrix=NULL, time=1, verbose=0, transfer_function=3,
	reltol=1e-04, atol=0.001, maxStepSize=Inf, maxNumSteps=1e+05, 
	maxErrTestsFails=50)
}

	\arguments{
  	\item{cnolist}{
  		A list containing the experimental design and data.
	}
 	 \item{model}{
	A list with the ODEs parameter information. Obtained with \code{\link{createLBodeContPars}}.
	}
  	\item{ode_parameters}{
		A list with the ODEs parameter information. Obtained with makeParameterList function.
	}
	  \item{indices}{
		Indices to map data in the model. Obtained with indexFinder function from CellNOptR. 
	}
	\item{adjMatrix}{
		The adjacency matrix. Recomputed if not provided
	}
  	\item{time}{
		An integer with the index of the time point to start the simulation. Default is 1.
	}
	\item{verbose}{
		A logical value that triggers a set of comments.
	}
	\item{transfer_function}{
		The type of used transfer. Use 1 for no transfer function, 2 for Hill function and
		3 for normalized Hill function.
	}
	\item{reltol}{
		Relative Tolerance for numerical integration.
	}
 	\item{atol}{
		Absolute tolerance for numerical integration.
	}
	\item{maxStepSize}{
		The maximum step size allowed to ODE solver.
	}
	\item{maxNumSteps}{
		The maximum number of internal steps between two points being sampled before
		the solver fails.
	}
  	\item{maxErrTestsFails}{
		Specifies the maximum number of error test failures permitted in attempting
		one step.
	}
}

\details{
	Check \code{\link{CellNOptR}} for details about the cnolist and the model format.
	For more details in the configuration of the ODE solver check the CVODES manual. 
}

\value{
	Returns a list with simulated data that has the same structure as the cnolist$valueSignals.
	One matrix for each time-point.
}

\author{
	David Henriques, Thomas Cokelaer
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
	\code{\link{CellNOptR}}
	\code{\link{parEstimationLBode}}
	\code{\link{parEstimationLBodeSSm}}
}

\examples{
library(CNORode)
data("ToyCNOlist",package="CNORode");
data("ToyModel",package="CNORode");
data("ToyIndices",package="CNORode");
dataSimulation = simulate(cnolistCNORodeExample, model,indices=indices);
}

