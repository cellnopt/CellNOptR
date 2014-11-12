\name{getLBodeContObjFunction}

\alias{getLBodeContObjFunction}

\title{
	Returns the objective function to perform parameter estimation.
}

\description{
	This function configures returns the objective function that can be used to evaluate
	the fitness of a logic based ODE model using a particular set of parameters.
	This function can be particularly useful if you are planing to couple a nonlinear
	optimization solver. The returned value of the objective function corresponds to 
	the mean squared value normalized by the number of data points.
}

\usage{
	getLBodeContObjFunction(cnolist, model, ode_parameters, indices=NULL, time = 1, 
	verbose = 0, transfer_function = 3, reltol = 1e-04, atol = 0.001, maxStepSize = Inf, 
	maxNumSteps = 1e+05, maxErrTestsFails = 50, nan_fac = 1, useVariances = F,initial_state=0.1)
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
  	\item{nan_fac}{
		A penalty for each data point the model is not able to simulate. We recommend higher than
		0 and smaller that 1.
	}

    \item{useVariances}{
if True, use the variance in the fitness
}

    \item{initial_state}{
initial state of the dynamic nodes (non-measured) (Defaults to 0.1)
}

}

\details{
	Check \code{\link{CellNOptR}} for details about the cnolist and the model format.
	For more details in the configuration of the ODE solver check the CVODES manual. 
}

\value{
	Returns a function to evaluate the model fitness. This function receives a vector
	containing both continuous parameters and integer values representing which reactions
	should be kept in the model.
}

\author{
	David Henriques, Thomas Cokelaer
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
	\code{\link{CellNOptR}}
	\code{\link{createLBodeContPars}}
}

\examples{
 	library(CNORode)
	data("ToyCNOlist",package="CNORode");
	data("ToyModel",package="CNORode");
	data("ToyIndices",package="CNORode");
	
	ode_parameters=createLBodeContPars(model,random=TRUE);
	minlp_obj_function=getLBodeContObjFunction(cnolistCNORodeExample, model,ode_parameters,indices);
	
	x=ode_parameters$parValues;
	
	f=minlp_obj_function(x);
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.

