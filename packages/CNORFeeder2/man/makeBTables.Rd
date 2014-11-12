\name{makeBTables}
\alias{makeBTables}

\title{
  Make Boolean tables
}

\description{
  This function uses data (CNOlist) to infer a Boolean table for each measured protein, codifying if a particular stimulus inhibitor combination affects the protein. A stimulus or an inhibitor significantly affects an output protein if it is able to modify its activity level of a quantity that exceeds the uncertainty associated with its measurement.
}

\usage{
makeBTables(CNOlist, k=2, measErr=c(0.1, 0), timePoint=NA)
}


\arguments{
  \item{CNOlist}{a CNOlist structure, as produced by \link{makeCNOlist}}
  \item{k}{a parameter that determine the threshold of significancy of the effect of stimuli and inhibitors, default to 2}
  \item{measErr}{a 2 value vector (err1, err2) defining the error model of the data as sd^2 = err1^2 + (err2*data)^2, default to c(0.1, 0)}
  \item{timePoint}{the time point to be considered for the inference of the Boolean tables (i.e. "t1" or "t2"), if not specified all time points are consideres}
}
\details{
This function computes the fist step of FEED to reverse engineer the network strictly from data, i.e. the inference of Boolean tables, as described in (Eduati et al., PLoS ONE, 2010). For each protein, a Boolean table is inferred having one columns for each stimulus and one row for each inhibitor. If a stimulus produces a significant effect on the activity level of the protein this is codified with a 1 in the corresponding column, if also the inhibitor affects the protein there is a 2 in the corresponding cell. The sign of the regulation is coded in separate tables.
}

\value{
this function returns a list with fields:
\item{namesSignals}{a vector of names of signals}
\item{tables}{a list with one Boolean table for each protein codifying the effect of stimuli (columns) and inhibitors (rows), 1 if the stimulus affect the protein, 2 if also the inhibior does}
\item{NotMatStim}{has the same format as tables but just contains a 1 if the regulation has a negative effect, and 0 otherwise}
\item{NotMatInhib}{has the same format as tables but just contains a 1 if the regulation has a negative effect, and 0 otherwise}

}
\references{
F. Eduati, A. Corradin, B. Di Camillo, G. Toffolo. A Boolean approach to linear prediction for signaling network modeling. PLoS ONE; 5(9): e12789.
}
\author{
    F.Eduati
}

\seealso{
\link{makeCNOlist}, \link{mapBTables2model}
}

\examples{

data(CNOlistDREAM,package="CellNOptR")
BTable <- makeBTables(CNOlist=CNOlistDREAM, k=2, measErr=c(0.1, 0))
    
}

