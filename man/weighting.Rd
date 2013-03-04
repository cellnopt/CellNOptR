\name{weighting}
\alias{weighting}

\title{
  Weight integrated links.
}

\description{
  This function weights links integrated in the model using additional penalty and/or information from protien-protein interactions networks (PINs).
}

\usage{
weighting(modelIntegr,PKNmodel,CNOlist,integrFac,UniprotID,PPI)
}


\arguments{
  \item{modelIntegr}{the integrated model as created by \link{mapDDN2model} or \link{mapBTables2model}}
  \item{PKNmodel}{the model of the original prior-knowledge network}
  \item{CNOlist}{a CNOlisi, as created by \link{makeCNOlist}}
  \item{integrFac}{a number indicating the penalty for integrated links}
  \item{UniprotID}{a list with the Uniprot identifiers of proteins in the PKN}
  \item{PPI}{an igraph of the PIN to be used, if no network is provided (=NULL) this information is not used. Default is NULL.}
 
}

\details{
Integrated links are less reliable than links from the PKN, thus should be penalized in the optimization process. This function allows to include a panalty for integrated links (integrFact). Furthermore links can be differently prioritized based on information derived from pritein interaction networks (PIN): the basic idea is that if, for a directed link A -> B integrated in the PKN, there is a corresponding path in the PIN, it is more plausible that there is a molecular pathway A -> B. Because shorter paths are more feasible, as a first approximation the shortest path length between A and B in the PIN can be used as a reliability score for the integrated link. Since the optimization is performed on a compressed version of the PKN, one link integrated in the compressed network generally corresponds to multiple possible links integrated in the PKN and the shortes path of all. The weight for each integrated link in the compressed network is thus computed as (1 + the inverse of the sum of the inverse of the corresponding PKN of the shortest paths in the PIN). A high quality network of known human physical protein-protein interaction assembled from multiple databases is provided with the package: interactions were included only if validated by at least one binary experimental method in a published paper and the number of experimental evidences was reported for each interaction.
}

\value{
\item{modelIntegr}{the input modelIntegr with an additional field: a vector with the weights of the integrated links}
}

\author{
    F.Eduati
}


\seealso{
    \link{mapDDN2model}, \link{mapBTables2model}, \link{gaBinaryT1W} 
}


\examples{

data(CNOlistDREAM,package="CellNOptR")
data(DreamModel,package="CellNOptR")
data(UniprotIDdream,package="CNORfeeder")

model<-preprocessing(data=CNOlistDREAM, model=DreamModel)

BTable <- makeBTables(CNOlist=CNOlistDREAM, k=2, measErr=c(0.1, 0))
modelIntegr <- mapBTables2model(BTable=BTable,model=model,allInter=TRUE)

modelIntegrWeight <- weighting(modelIntegr=modelIntegr, PKNmodel=DreamModel,
					 CNOlist=CNOlistDREAM, integrFac=10)

# weighting using PPI might take some minutes
\dontrun{
data(UniprotIDdream,package="CNORfeeder")
data(PPINigraph,package="CNORfeeder")
modelIntegrWeight2 <- weighting(modelIntegr=modelIntegr, PKNmodel=DreamModel,
					  CNOlist=CNOlistDREAM, integrFac=10, UniprotID=UniprotIDdream,
					  PPI=PPINigraph)
}


}

