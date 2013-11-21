\name{plotLBodeFitness}

\alias{plotLBodeFitness}
%- Also NEED an '\alias' for EACH other topic documented here.

\title{
	Plot data against simulated values. 
}

\description{
	Plots the simulated values with the logic-based ODE against the the data contained
	contained the data contained in the cnolist. The data values are represented with
	a black line and the simulated values with a blue line. Additionally this functions
	returns the the simulated values.
}

\usage{
	plotLBodeFitness(cnolist, model, ode_parameters = NULL, indices = NULL,
 		adjMatrix = NULL, time = 1, verbose = 0, transfer_function = 3, reltol = 1e-04,
		atol = 0.001, maxStepSize = Inf, maxNumSteps = 1e+05, maxErrTestsFails = 50,
   		plot_index_signals = NULL, plot_index_experiments = NULL,
plot_index_cues = NULL, colormap="heat", plotParams=list(margin=0.1, width=15, height=12,
                  cmap_scale=1, cex=1.6, ymin=NULL)
  

)
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
  	\item{plot_index_signals}{
		In case you only want to plot some signals, provide an integer vector with  the indexes.
	}
	\item{plot_index_experiments}{
		In case you only want to plot some experiments, provide an integer vector with  the indexes.
	}
  	\item{plot_index_cues}{
		In case you only want to plot some cues, provide an integer vector with  the indexes.
	}
    \item{colormap}{
    Uses the same colormap as in CellNOptR by default. If set to "green", it uses the deprecated colormap.
    }

    \item{plotParams}{
    additional parameters to refine the ploggin. See plotOptimResultsPan
function in CellNOptR for more details.
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
	dataSimulation=plotLBodeFitness(cnolistCNORodeExample, model,indices=indices);
}

% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
