\name{MapBTables2Model}
\alias{MapBTables2Model}

\title{
  Integrate Boolean tables with the model
}

\description{
  This function infers the network from the Boolean tables and integrates it with the network encoded in the model (generally derived from prior knowledge), adding links that are missing. 
}

\usage{
MapBTables2Model(BTable,Model,optimRes=NA,allInter=TRUE)
}


\arguments{
  \item{BTable}{a BTable list, as created by \link{makeBTables}}
  \item{Model}{a Model list, as created by \link{readSif}}
  \item{optimRes}{a bit string with the reaction of the model to be considered, default considers all reactions}
  \item{allInter}{one new link in the network can correspond to more links in the model, set it to TRUE if you want to add all possible links, FALSE to add only one link, default is TRUE}
}

\details{
The function receive as input the Boolean Tables, infers the data-driven network form them (as descibed in (Eduati et al., PLoS ONE, 2010)) and integrates it with the model, returning a new model with the integrated links. If the Model is not given as input (Model=NULL), the data-driven network is returned as model.
}

\value{
a new Model with the integrated links and an additional field:
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

library(CellNOptR)
data(CNOlistDREAM,package="CellNOptR")
data(DreamModel,package="CellNOptR")
res<-preprocessing(Data=CNOlistDREAM, Model=DreamModel)
Model<-res$model
BTable <- makeBTables(CNOlist=CNOlistDREAM, k=2, measErr=c(0.1, 0))
modelIntegr <- MapBTables2Model(BTable=BTable,Model=Model,allInter=TRUE)
# modelIntegr$reacID[modelIntegr$indexIntegr] to see the integrated links
}

