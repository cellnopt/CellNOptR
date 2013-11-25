\name{mapBTables2model}
\alias{mapBTables2model}

\title{
  Integrate Boolean tables with the model
}

\description{
  This function infers the network from the Boolean tables and integrates it with the network encoded in the model (generally derived from prior knowledge), adding links that are missing. 
}

\usage{
mapBTables2model(BTable,model,optimRes=NA,allInter=TRUE,compressed=TRUE)
}


\arguments{
  \item{BTable}{a BTable list, as created by \link{makeBTables}}
  \item{model}{a model list, as created by \link{readSif}}
  \item{optimRes}{a bit string with the reaction of the model to be considered, default considers all reactions}
  \item{allInter}{one new link in the network can correspond to more links in the model, set it to TRUE if you want to add all possible links, FALSE to add only one link, default is TRUE}
  \item{compressed}{this argument is used to decede how to deal with unmeasured and unperturbed nodes (white nodes). As general guideline, it should be set to TRUE if the PKN has been compressed in the preprocessing step, FALSE otherwise. Default is TRUE.}
}

\details{
The function receive as input the Boolean Tables, infers the data-driven network form them (as descibed in (Eduati et al., PLoS ONE, 2010)) and integrates it with the model, returning a new model with the integrated links. If the Model is not given as input (Model=NULL), the data-driven network is returned as model.
}

\value{
a new model with the integrated links and an additional field:
\item{indexIntegr}{a vector with the indexes of the integrated links}
}

\author{
    F.Eduati
}

\references{
F. Eduati, A. Corradin, B. Di Camillo, G. Toffolo. A Boolean approach to linear prediction for signaling network modeling. PLoS ONE; 5(9): e12789.
}

\seealso{
    \link{readSif}, \link{readMIDAS}, \link{makeBTables} 
}


\examples{

data(CNOlistDREAM,package="CellNOptR")
data(DreamModel,package="CellNOptR")
model<-preprocessing(data=CNOlistDREAM, model=DreamModel)
BTable <- makeBTables(CNOlist=CNOlistDREAM, k=2, measErr=c(0.1, 0))
modelIntegr <- mapBTables2model(BTable=BTable,model=model,allInter=TRUE)
# modelIntegr$reacID[modelIntegr$indexIntegr] to see the integrated links
}

