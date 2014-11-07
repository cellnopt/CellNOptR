\name{simdata2cnolist}

\alias{simdata2cnolist}

\title{
	converts output of getLBodeModelSim to cnolist
}

\description{
	This function converts the simulated data returned by getLBodeModelSim into
    a valid CNOlist data structure.
}

\usage{
	simdata2cnolist(sim_data, cnolist, model)
}

\arguments{
    \item{sim_data}{
        structure returned by getLBodeModelSim
    }
	\item{cnolist}{
		A list containing the experimental design and data. 
	}
	\item{model}{
		The logic model to be simulated.
	}
}


\value{
	a CNOlist
}

\author{
	Thomas Cokelaer
}
%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
	\code{\link{CellNOptR}}
	\code{\link{createLBodeContPars}}
}

\examples{
data('ToyCNOlist',package='CNORode');
data('ToyModel',package='CNORode');
data('ToyIndices',package='CNORode');
simdata = getLBodeModelSim(cnolistCNORodeExample, model,indices=indices)
cnolist = simdata2cnolist(simdata, cnolistCNORodeExample, model)

cnolist = simdata2cnolist(simdata, cnolistCNORodeExample, model)
}

