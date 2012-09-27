\name{parEstimationLBodeGA}

\alias{parEstimationLBodeGA}
%- Also NEED an '\alias' for EACH other topic documented here.

\title{
	Perform parameter estimation using a genetic algorithm (package genalg).
}

\description{
	This function uses a genetic algorithm (package genalg) to perform parameter estimation. The objective function
	is the same as the one provided by \code{\link{getLBodeContObjFunction}}.
}

\usage{
	parEstimationLBodeGA(cnolist, model, ode_parameters = NULL, indices = NULL, mutationChance = NA, popSize = 200, iters = 100, 
		elitism = NA, time = 1, monitor = TRUE, verbose = 0, transfer_function = 3, reltol = 1e-04, 
		atol = 0.001, maxStepSize = Inf, maxNumSteps = 1e+05, maxErrTestsFails = 50, nan_fac = 1)
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
  	\item{mutationChance}{
		the chance that a gene in the chromosome mutates. By default 1/(size+1). It affects
		the convergence rate and the probing of search space: a low chance results
		in quicker convergence, while a high chance increases the span of the search
		space.
	}	
	\item{popSize}{
		the population size.	
	}
	\item{iters}{
		the number of iterations.
	}
	\item{elitism}{
		the number of chromosomes that are kept into the next generation. By default is
		about 20\% of the population size
	}
	\item{time}{
		An integer with the index of the time point to start the simulation. Default is 1.
	}
	 \item{monitor}{
		 If TRUE a plot will be generated to monitor the objective function
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
}

\value{
	\item{LB_n}{
		A numeric value to be used as lower bound for all parameters of type n.
	}
	\item{LB_k}{
		A numeric value to be used as lower bound for all parameters of type k.
	}
	\item{LB_tau}{
		A numeric value to be used as lower bound for all parameters of type tau.
	}
	\item{UB_n}{
		A numeric value to be used as upper bound for all parameters of type n.
	}
 	\item{UB_k}{
		A numeric value to be used as upper bound for all parameters of type k.
	}
  	\item{UB_tau}{
		A numeric value to be used as upper bound for all parameters of type tau. 
	}
  	\item{default_n}{
		The default parameter to be used for every parameter of type n. 
	}
 	 \item{default_k}{
		The default parameter to be used for every parameter of type k. 
	}
	\item{default_tau}{
		The default parameter to be used for every parameter of type tau. 
	}
	\item{LB_in}{
		An array with the the same length as ode_parameters$parValues with lower bounds
		for each specific parameter. 
	}
	\item{UB_in}{
		An array with the the same length as ode_parameters$parValues with upper bounds
		for each specific parameter.
	}
	\item{opt_n}{
		Add all parameter n to the index of parameters to be fitted.
	}
	\item{opt_k}{
		Add all parameter k to the index of parameters to be fitted.
	}
	\item{opt_tau}{
		Add all parameter tau to the index of parameters to be fitted.
	}
	\item{random}{
		A logical value that determines that a random solution is for the parameters 
		to be optimized.
	}
	\item{res}{
		A list containing the information provided by the nonlinear optimization solver (genalg).
	}
}


\author{
	David Henriques, Thomas Cokelaer
}
%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
	\code{\link{CellNOptR}}
	\code{\link{createLBodeContPars}}
	\code{\link{rbga}}
}

\examples{
data("ToyCNOlist",package="CNORode");
data("ToyModel",package="CNORode");
data("ToyIndices",package="CNORode");
	
ode_parameters=createLBodeContPars(model,random=TRUE);
#Visualize intial simulation
#simulatedData=plotLBodeFitness(cnolistCNORodeExample, model,ode_parameters,indices=indices)

ode_parameters=parEstimationLBodeGA(cnolistCNORodeExample,model,ode_parameters=ode_parameters,
indices=indices,maxStepSize=1,atol=1e-3,reltol=1e-5,transfer_function=2,popSize=10,iter=40);

#Visual solution after optimization
simulatedData=plotLBodeFitness(cnolistCNORodeExample, model,indices=indices,ode_parameters=ode_parameters);
}

% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{genetic}
\keyword{algorithm}% __ONLY ONE__ keyword per line
\keyword{CellNOptR}
\keyword{logic}
\keyword{model}
