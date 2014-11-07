\name{plotMeanFuzzyFit}

\alias{plotMeanFuzzyFit}

\title{Simulates models returned from multiple cFL runs and plots mean fit to data}

\description{Uses post refinement threshold (selection threshold) to choose reduced refined model resulting from each run.  Simulates model and plots result and fit to data}

\usage{
  plotMeanFuzzyFit(postRefThresh, allFinalMSEs, allRes, plotPDF=FALSE, tag=NULL,
show=TRUE, plotParams=list(cex=0.8, cmap_scale=1))
}

\arguments{
  \item{postRefThresh}{Post refinement threshold (selection threshold) chosen
from plot produced by compileMultiRes}.
  \item{allFinalMSEs}{matrix containing MSEs produced by compileMultiRes}
  \item{allRes}{list containing results of several \link{CNORwrapFuzzy} runs}
  \item{plotPDF}{TRUE or FALSE depending on if a PDF file should be saved}
  \item{tag}{String to include in filename of PDF plot}
   \item{show}{If the plot should be displayed}
   \item{plotParams}{
    a list of option related to the PDF and plotting outputs. (1) cex is the
font size of the header. (2) cmap_scale below 1 allows to put more emphasizes on low
errors (default 1 means all colors have the same weight). See
\code{plotOptimResultsPan} from CellNOptR  for other fields.
}

}


\value{
    This function does not have any output, it just plots and saves results if
applicable.
}

\author{M.K. Morris}

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
    plotMeanFuzzyFit(0.1, summary$allFinalMSEs, allRes)

    }
}





