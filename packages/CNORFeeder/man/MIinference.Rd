\name{MIinference}
\alias{MIinference}

\title{
  Mutual information based network inference
}

\description{
  This function uses data (CNOlist) to infer a data-driven network using the mutual information based appoaches ARACNe and CLR as implemented in the minet package.
}

\usage{
MIinference(CNOlist, method="ARACNE", PKNgraph=NULL,
            filename="ARACNE")
}


\arguments{
  \item{CNOlist}{a CNOlist structure, as produced by \link{makeCNOlist}}
  \item{method}{a character, the name of the method to be used: ARACNE or CLR. Default, ARACNE}
  \item{PKNgraph}{a network to be used for comparison to assess the directionality of some links. Default is NULL.}
  \item{filename}{name of the sif file saved, default ARACNE}
}
\details{
This function transforms the data in a format compatible with minet package, infers the network using aracne or clr as implemented in the minet package and returns the network in the sif format. It is important to notice that mutual information approaches do not allow for determining the directionality of the links thus both directions are considered. The function allows to give as input a network in graph format (graph package, see \link{sif2graph} to convert from sif to graph format) to be used as comparison to assess the directionality of some links, e.g. PKN.
}

\value{
\item{sif}{the inferred data-driven network in sif format}
}
\references{
P. E. Meyer, F. Lafitte and G. Bontempi (2008). MINET: An open source R/Bioconductor Package for Mutual Information based Network Inference. BMC Bioinformatics, 9(1), 2008
}
\author{
    F.Eduati
}

\seealso{
\link{mapDDN2model}, \link{sif2graph}, \link{model2sif}
}

\examples{

data(CNOlistDREAM,package="CellNOptR")
data(DreamModel,package="CellNOptR")
PKNgraph<-sif2graph(model2sif(DreamModel))

method="ARACNE"
#method="CLR"
DDN<-MIinference(CNOlist=CNOlistDREAM, method=method,
                 PKNgraph=PKNgraph, filename=method)

    
}

