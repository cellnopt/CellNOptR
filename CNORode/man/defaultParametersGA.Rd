\name{defaultParametersGA}

\alias{defaultParametersGA}

\title{
	Create default options to perform parameter estimation with a genetic algorithm.
}

\description{
	This function returns a list with several arguments for performing parameter estimation with the genetic 
	algorithm from the package genalg.
}
\usage{
	defaultParametersGA()
}
%- maybe also 'usage' for other objects documented here.

\value{
	\item{mutationChance}{NA}
	\item{popSize}{200}
	\item{iters}{100}
	\item{elitism}{NA}
	\item{time}{1}
	\item{monitor}{TRUE}
	\item{verbose}{0}
	\item{transfer_function}{3}
	\item{reltol}{1e-04}
	\item{atol}{0.001}
	\item{maxStepSize}{Inf}
	\item{maxNumSteps}{1e+05}
	\item{maxErrTestsFails}{50}
	\item{nan_fac = 1}{0}
}

\author{
	David Henriques, Thomas Cokelaer
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
	\code{\link{CellNOptR}}
	\code{\link{parEstimationLBode}}
	\code{\link{parEstimationLBodeGA}}
}

% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.

\keyword{default }
\keyword{parameters}
\keyword{genetic}
\keyword{algorithm}