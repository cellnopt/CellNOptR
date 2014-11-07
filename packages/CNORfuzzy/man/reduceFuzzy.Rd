\name{reduceFuzzy}

\alias{reduceFuzzy}

\title{Remove unnecessary interactions from cFL model}

\description{Determine if interactions in cFL model are necessary to fit the data}

\usage{reduceFuzzy(firstCutOff, CNOlist, model,  res, params )}

\arguments{
  \item{firstCutOff}{Threshold for removing or replacing and edge.  If the score doesn't get any worse than this, its removed or replaced}
  \item{CNOlist}{a CNOlist on which the score is based (based on
      \code{valueSignals[[2]]}, i.e. data at t1)}
  \item{model}{Model input into \link{gaDiscreteT1} (uncut)}
  \item{res}{Optimum returned by gaDiscrete}
  \item{params}{parameters (as generated from \link{defaultParametersFuzzy} (only Fuzzy parameters required)}
}

\details{After gaDiscrete removes interactions that are inconsistent with the
data, this function asks if they are necessary to fit the data.  For all
interactions, it asks if they can be removed without making the score worse by a
value of the firstCutOff.  For AND interactions, it asks if they can be replaced
with an OR gate.}

\value{
  \item{redModel}{}
  \item{redSimList}{}
  \item{bitString}{}
  \item{intString}{}
  \item{MSE}{}
}

\author{M.K. Morris}

\seealso{
    \link{gaDiscreteT1}
}


