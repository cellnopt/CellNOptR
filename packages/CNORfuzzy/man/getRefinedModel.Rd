\name{getRefinedModel}

\alias{getRefinedModel}

\title{Refinement of Parameters of cFL model}

\description{Performs refinement of cFL model parameters}

\usage{getRefinedModel(res,CNOlist,cutModel, cutSimList, refParams)}

\arguments{
  \item{res}{Optimum returned by \link{gaDiscreteT1}}
  \item{CNOlist}{a CNOlist on which the score is based (based on valueSignals[[2]], i.e. data at T1)}
  \item{cutModel}{Model (with unnesssary edges cut in \link{reduceFuzzy} or \link{interpretDiscreteGA}}
  \item{cutSimList}{Fields for simulation based on cutModel (again, cut with reduceFuzzy or interpretDiscreteGA)}
  \item{refParams}{parameter list object as returned by \link{defaultParametersFuzzy}.}
}

\details{After the discreteGA chooses transfer functions from a discrete set of transfer functions and removing interactions inconsistent with the data, this function 'refines' the parameters by using optim to go to the local minimum of error to data}

\value{
    \item{refModel}{A refined model}
    \item{finalSet}{Set of final fuzzy parameter}
    \item{MSE}{The MSE value}
}

\author{M.K. Morris}

\seealso{
    \link{gaDiscreteT1}, \link{reduceFuzzy}, \link{interpretDiscreteGA}
}

