\name{compileMultiRes}

\alias{compileMultiRes}

\title{Compiles results from multiple runs and produces graph for choosing Post Refinement Threshold}

\description{This function takes a list of objects returned by CNORwrapfuzzy (run using
identical parameters, models, and data) and packages them together so they can
be compared with plotMeanFuzzyFit and writeFuzzyNetwork.  Because almost all
training of cFL models are underdetermined problems, analyzing multiple runs
together is essential.}

\usage{
	compileMultiRes(allRes, tag=NULL, show=TRUE)
}

\arguments{
  \item{allRes}{list of objects returned by the CNORwrapFuzzy function. }
  \item{tag}{If provided, save the results in 3 files. Each file starts
with the string "filename" that is provided. (<tag>_allRes.RData,
<tag>_allFinalMSEs.RData and <tag>_allFinalNumParams.RData)}
  \item{show}{plot the MSE and mean number of parameters versus threshold. Can
be switch off if \code{show=FALSE}}
}


\author{M.K. Morris, T. Cokelaer}


\examples{


    data(ToyModel, package="CellNOptR")
    data(CNOlistToy,package="CellNOptR")
    paramsList = defaultParametersFuzzy(CNOlistToy, ToyModel)
    N = 10
    allRes = list()
    \dontrun{
        for (i in 1:N){
            Res = CNORwrapFuzzy(CNOlistToy, ToyModel, paramsList)
            allRes[[i]] = Res
        }

        summary = compileMultiRes(allRes)
        summary$allFinalMSEs
        summary$allFinalNumParams
    
        # You can save the resuls in files using the tag argument
        compileMultiRes(allRes, "output")
    }

}

