\name{linksRanking}
\alias{linksRanking}

\title{
  Ranking of links inferred from data
}

\description{
  This function uses data (CNOlist) to rank links based on measurement error model as used by FEED method to reverse-engineer the network.
}

\usage{
linksRanking(CNOlist, measErr=c(0.1, 0), savefile=FALSE)
}


\arguments{
  \item{CNOlist}{a CNOlist structure, as produced by \link{makeCNOlist}}
  \item{measErr}{a 2 value vector (err1, err2) defining the error model of the data as sd^2 = err1^2 + (err2*data)^2, default to c(0.1, 0)} 
  \item{savefile}{TRUE to save the file in txt format, FALSE not. Default is FALSE.}
  }
\details{
This function is similar to the fist step of FEED to reverse engineer the network strictly from data, i.e. the inference of Boolean tables, as described in (Eduati et al., PLoS ONE, 2010) and implemented in \link{makeBTables}. Links are ranked according to the upper limit value of parameterk allowing the presence of the link, where k is the parameter which is multiplied by the measurement error in order to assess the relevance of a link. The function returs link in decreasing order of importance and associate to each link a value (maximum value of k allowing the presence of the link) quantifying its relevance.  
}

\value{
this function returns a list with fields:
\item{Lrank}{a matrix in which each link is associated with a numerical value, links are ordered in decreasing order of reliability)}
}
\references{
F. Eduati, A. Corradin, B. Di Camillo, G. Toffolo. A Boolean approach to linear prediction for signaling network modeling. PLoS ONE; 5(9): e12789.
}
\author{
    F.Eduati
}

\seealso{
\link{makeCNOlist}, \link{makeBTables}
}

\examples{

data(CNOlistDREAM,package="CellNOptR")
Lrank <- linksRanking(CNOlist=CNOlistDREAM, measErr=c(0.1, 0))
    
}

