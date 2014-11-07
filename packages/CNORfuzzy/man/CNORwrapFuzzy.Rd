\name{CNORwrapFuzzy}
\alias{CNORwrapFuzzy}
\title{CNORfuzzy analysis wrapper}

\description{This function is a wrapper around the whole CNOR Fuzzy analysis. 
It performs the following steps:

\enumerate{
    \item Plot the CNOlist
    \item Checks data to model compatibility
    \item Pre-processing steps
    \item Prepare for simulation (see \link{prep4simFuzzy})
    \item Optimisation using Fuzzy transfer function (see \link{gaDiscreteT1})
    \item Refinement and reduction steps (see \link{getRefinedModel} and
        \link{reduceFuzzy}). 
    }
}

\usage{
CNORwrapFuzzy(data, model, paramsList=NULL, verbose=TRUE)
}

\arguments{

    \item{data}{
    a CNOlist structure (as created by \code{\link{makeCNOlist}}) that contains the data that you will use (see \link{readMIDAS}
    and \link{readSIF} from CellNOptR).
    }

    \item{model}{
    the model that you want to optimise
    }

    \item{paramsList}{
    Use defaultParametersFuzzy function to create a template. Entries
    are 3-types: (i) GA algorithm parameters for the optimisation, (ii) Fuzzy
    parameters for the transfer functions and (iii) internal optimisation
    parameters within the refinement step.

    See \link{defaultParametersFuzzy} function for details on the parameters. 
    }

    \item{verbose}{

    }

}

\details{
If you do not provide a parameter list argument, \link{defaultParametersFuzzy} is called
internally to populate the paramsList argument. 
}
\value{
This function returns an object containing the results that can be used by other
functions such as compileRes.
}

\seealso{\link{defaultParametersFuzzy}, \link{compileMultiRes}}

\references{CNORWrap function from \link{CellNOptR} package}

\author{M.K. Morris}


\examples{

    # Load some data
    data(CNOlistToy, package="CellNOptR")
    data(ToyModel, package="CellNOptR")
    # Get some default parameters to play with, limiting the duration of the GA
    # algorithm and optimisation step
    paramsList = defaultParametersFuzzy() 
    paramsList$maxTime = 20
    paramsList$optimisation$maxtime = 10
    results = CNORwrapFuzzy(CNOlistToy, ToyModel, paramsList)


}

