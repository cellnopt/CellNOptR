\name{defaultParametersFuzzy}
\alias{defaultParametersFuzzy}
\title{
Create a list of default parameters
}
\description{
Parameters are required at different levels in the Fuzzy optimisation. This
function provides a list with all parameters that are necessary.}
\usage{
defaultParametersFuzzy(data=NA, model=NA, nTF=7)
}

\arguments{
  \item{data}{
   the CNOlist that contains the data that you will use
  }
  \item{model}{
    the model that you want to optimise
  }
  \item{nTF}{number of discrete values to be used for each transfer function
parameter.
  }
}


\details{
The list contains 3 types of parameters: 

    \itemize{
        \item{ Fuzzy parameters (e.g, Type1Funs, Type2Funs, RedThresh, DoRefinement)}
        \item{ GA parameters similar to those used in CellNOptR package (see
            \link{gaDiscreteT1} or \link{defaultParametersFuzzy}
            }

        \item{ optimisation parameters related to the refinement step.
            \enumerate{
                \item{algorithm='NLOPT_LN_SBPLX'}
                \item{xtol_abs=0.001}
                \item{maxEval=1000}
                \item{maxTime=5*60}
            }
        }
        

        \item{In addition, you can set Model and Data (CNOlist).  }
    }
}
\value{
\item{params}{a list of default parameters. }
}
\author{
 T. Cokelaer
}
\examples{
    data(ToyModel,package="CellNOptR")
    data(CNOlistToy,package="CellNOptR")
    params = defaultParametersFuzzy(CNOlistToy, ToyModel)
}
