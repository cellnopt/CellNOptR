\name{PPIweight}
\alias{PPIweight}

\title{
  Weight links using protien-protein interactions
}

\description{
  This function weights links integrated in the model using information derived from protein-protein interaction networks (PINs).
}

\usage{
PPIweight(modelIntegr,PKNmodel,CNOlist,UniprotID,PPINigraph)
}


\arguments{
  \item{modelIntegr}{the integrated model as created by \link{MapDDN2Model} or \link{MapBTables2Model}}
  \item{PKNmodel}{the model of the original prior-knowledge network}
  \item{CNOlist}{a CNOlisi, as created by \link{makeCNOlist}}
  \item{UniprotID}{a list with the Uniprot identifiers of proteins in the PKN}
  \item{PPINigraph}{the igraph (igraph package) of the PIN}
 
}

\details{
The basic idea is that if, for a directed link A → B integrated in the PKN,
there is a corresponding path in the PIN, it is more plausible that there is a
molecular pathway A -> B. Because shorter paths are more feasible, as a first approximation the shortest path length between A and B in the PIN can be used as a reliability score for the integrated link. Since the optimization is performed on a compressed version of the PKN, one link integrated in the compressed network generally corresponds to multiple possible links integrated in the PKN and the shortes path of all. The weight for each integrated link in the compressed network is thus computed as (1 + the inverse of the sum of the inverse of the corresponding PKN of the shortest paths in the PIN). A high quality network of known human physical protein-protein interaction assembled from multiple databases is provided with the package: interactions were included only if validated by at least one ‘binary’ experimental method in a published paper and the number of experimental evidences was reported for each interaction.
}

\value{
This function returns a list with elements:
\item{modelIntegr}{the input modelIntegr with an additional field: a vector with the weights of the integrated links}
\item{PPINigraph}{the input PPINigraph with added colour attributes for edges and nodes}
\item{ListPaths}{the list with, for each integrated link, the list of the shortes paths in the PIN corresponding to the paths in the PKN}
\item{saveShortestPath}{has the same structure of ListPaths but contains only the length of the paths in terms of number of edges}
}

\author{
    F.Eduati
}


\seealso{
    \link{MapDDN2Model}, \link{MapBTables2Model}, \link{gaBinaryT1int} 
}


\examples{

library(CellNOptR)
data(CNOlistDREAM,package="CellNOptR")
data(DreamModel,package="CellNOptR")
data(UniprotIDdream,package="CNORfeeder")
data(PPINigraph,package="CNORfeeder")
res<-preprocessing(Data=CNOlistDREAM, Model=DreamModel)
Model<-res$model
    
BTable <- makeBTables(CNOlist=CNOlistDREAM, k=2, measErr=c(0.1, 0))
modelIntegr <- MapBTables2Model(BTable=BTable,Model=Model,allInter=TRUE)

# the followig step may take a while
\dontrun{
resPPIweight <- PPIweight(modelIntegr=modelIntegr,PKNmodel=DreamModel,
                CNOlist=CNOlistDREAM,UniprotID=UniprotIDdream,
                PPINigraph=PPINigraph)}


}

