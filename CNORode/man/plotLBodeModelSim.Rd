\name{plotLBodeModelSim}

\alias{plotLBodeModelSim}
%- Also NEED an '\alias' for EACH other topic documented here.

\title{
 	Simulate the model and plot the obtained with the different experimental conditions.
}

\description{
	Plots the simulated values of the logic based ODE model.
	Only dynamic states are plotted, i.e. those that are not inputs. 
	a blue line. Additionally this functions returns the the simulated
	values.
}

\usage{
	plotLBodeModelSim(cnolist, model, ode_parameters = NULL, indices = NULL, 
	adjMatrix = NULL, timeSignals=NULL, time = 1, verbose = 0, transfer_function = 3,
	reltol = 1e-04, atol = 0.001, maxStepSize = Inf, maxNumSteps = 1e+05,
	maxErrTestsFails = 50, large = FALSE, nsplit = 4)
}
%- maybe also 'usage' for other objects documented here.

\arguments{
	\item{cnolist}{
		A list containing the experimental design and data.
	}
	\item{model}{
		The logic model to be simulated.
	}
  \item{ode_parameters}{
	A list with the ODEs parameter information. Obtained with \code{\link{createLBodeContPars}}.
	}
	\item{indices}{
		Indices to map data in the model. Obtained with indexFinder function from CellNOptR. 
	}
	\item{adjMatrix}{
		Model representation in the form of an adjacency matrix. When not provided
		will be automatically computed based in the model.
	}
	\item{timeSignals}{
		An array containing a different timeSignals. If you use this argument, it will also
		modify the dimensions from valueSignals.
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
	\item{large}{
		Boolean variable defining if the plot should split into several subplots.
	}
 	\item{nsplit}{
		In case the large plot options is selected define how many subplots will exist. Default is 4.
	}
}

\value{
	Returns a list with simulated Model values. One matrix of size number of species by number of
	 experimental conditions for each time-point.
}

\author{
	David Henriques, Thomas Cokelaer
}

\seealso{
	\code{\link{CellNOptR}}
	\code{\link{createLBodeContPars}}
}

\examples{
    	library(CNORode)
	data("ToyCNOlist",package="CNORode");
	data("ToyModel",package="CNORode");
	data("ToyIndices",package="CNORode");
	modelSimulation=plotLBodeModelSim(cnolistCNORodeExample, model,indices=indices);
}

% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
