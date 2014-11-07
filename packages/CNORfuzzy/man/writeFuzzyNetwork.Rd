\name{writeFuzzyNetwork}
\alias{writeFuzzyNetwork}
\title{
Despict the network results of training a cFL model to data in multiple runs.
}
\description{
Fuzzy network results output.
}
\usage{
writeFuzzyNetwork(postRefThresh, allFinalMSEs, allRes, tag=NULL,verbose=FALSE)
}
\arguments{
  \item{postRefThresh}{Post refinement threshold (selection threshold) chosen from plot produced by compileMultiRes}
  \item{allFinalMSEs}{matrix containing MSEs produced by compileMultiRes}
  \item{allRes}{list containing all results produced by compileMultiRes}
  \item{tag}{String to include in filename of pdf plot}
	\item{verbose}{If extra warnings should be displayed}
}

\details{
The weights of the edges are computed as the mean across models using post refinement threshold (selection threshold) to choose reduced refined model resulting from each run. 

As with writeNetwork, this function maps back the edges weights from the optimised (expanded and compressed) model to the original model. The mapping back only works if the path has length 2 at most (i.e. you have node1-comp1-comp2-node2, where comp refer to nodes that have been compressed).
}
\value{
This function does not have any output, it just writes a SIF file, an edge attribute file, and a node attribute file
}
\author{
M.K. Morris based on code by C. Terfve
}
\note{
The mapback of this function is still an open question, even in the Matlab version.  Future developments will include more robust versions of the mapping back algorithm, probably as a separate mapback function.
}
\seealso{
writeNetwork
}
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
    
        writeFuzzyNetwork(postRefThresh, summary$allFinalMSEs, allRes)
    }

}
