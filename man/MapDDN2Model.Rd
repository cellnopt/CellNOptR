\name{MapDDN2Model}
\alias{MapDDN2Model}

\title{
  Integrate data-drive network with the model
}

\description{
This function integrates the data-driven network (in sif format)  with the network encoded in the model (generally derived from prior knowledge), adding links that are missing. 
}

\usage{
MapDDN2Model(DDN,Model,CNOlist,allInter=TRUE)
}


\arguments{
  \item{DDN}{a sif file encoding a data-driven network, as created by \link{Binference} or \link{MIinference}}
  \item{Model}{a Model list, as created by \link{readSif}}
  \item{CNOlist}{a CNOlist, as created by \link{makeCNOlist}}
  \item{allInter}{one new link in the network can correspond to more links in the model, set it to TRUE if you want to add all possible links, FALSE to add only one link, default is TRUE}
}

\details{
The function receives as input a sif file with the data-driven network, as created by \link{Binference} or \link{MIinference}, and integrates it with the model, returning a new model with the integrated links. 
}

\value{
a new Model with the integrated links and an additional field:
\item{indexIntegr}{a vector with the indexes of the integrated links}
}

\author{
    F.Eduati
}


\seealso{
    \link{readSif}, \link{readMIDAS},
    \link{Binference}, \link{MIinference}
}


\examples{

library(CellNOptR)
data(CNOlistDREAM,package="CellNOptR")
data(DreamModel,package="CellNOptR")
res<-preprocessing(Data=CNOlistDREAM, Model=DreamModel)
Model<-res$model
    
DDN<-Binference(CNOlistDREAM, tempCheckOrders=10, maxIter=100,
                filename="BAYESIAN")

modelIntegr<-MapDDN2Model(DDN=DDN,Model=Model,CNOlist=CNOlistDREAM)
 
}

