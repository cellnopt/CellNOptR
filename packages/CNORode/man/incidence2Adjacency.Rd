\name{incidence2Adjacency}

\alias{incidence2Adjacency}

\title{
	Convert an incidence matrix into an adjacency matrix.
}

\description{
	Convert the incidence matrix (model representation of CellNoptR) into an adjacency matrix. Denotes the inputs/output relationships.
}

\usage{
	incidence2Adjacency(model)
}

\arguments{
 	\item{model}{
		Model from CellNoptR.
	}
}

\value{
	Directed Adjacency matrix of size n_species by n_species. 
}

\note{
	For internal use of CNORode.
}

\author{
	David Henriques, Thomas Cokelaer
}

\seealso{
	\code{\link{CellNOptR}}
}

% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{incidence}
\keyword{adjacency}
\keyword{matrix}
