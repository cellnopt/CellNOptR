\name{getStates}

\alias{getStates}
%- Also NEED an '\alias' for EACH other topic documented here.

\title{
	Find which species in the model are states.
}

\description{
	Receives an adjacency matrix (model$interMat from CellNoptR) and finds which species are states (i.e. not inputs).
}

\usage{
	getStates(adjacency)
}
%- maybe also 'usage' for other objects documented here.

\arguments{
 	\item{adjacency}{
		An adjacency matrix from the model.
	}
}

\value{
	A numeric vector with 0's for positions which are states and 1's  for positions which are.
}

\note{
	For internal use of CNORode.
}

\author{
	David Henriques, Thomas Cokelaer
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
	\code{\link{incidence2Adjacency}}
}

% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{states}
\keyword{adjacency}
