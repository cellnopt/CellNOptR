\name{simFuzzyT1}
\alias{simFuzzyT1}
\title{
Simulation of a cFL model
}
\description{
This is the simulator, inspired from CNOfuzzySimEngv23 in the Matlab CellNOpt, to be used on one time point simulations
}
\usage{
simFuzzyT1(CNOlist, model, simList)
}
\arguments{
  \item{CNOlist}{
a CNOlist
}
  \item{model}{
a Model that only contains the reactions to be evaluated
}
  \item{simList}{
a SimList as created by prep4simFuzzy, that has also already been cut to contain only the reactions to be evaluated
}
}
\details{I'm not sure if this applies to the one I wrote, which was based on I think an older version by Camille maybe?
Differences from the BoolSimEngMKM simulator include: 
the valueInhibitors has not been previously flipped;
the function outputs the values across all conditions for all species in the model, instead of only for the signal species. This is because then the  output of this function can be used as initial values for the version of the simulator that works on time point 2 (not implemented in this version).

If you would like to compute the output of a model that contains some of the gates in the model but not all, we suggest that you use the function \code{SimulateT1} and specify in the \code{bStringT1} argument which gates you want to be included.  Indeed, SimulateT1 is a wrapper around simulatorT1 that takes care of cutting the model for you before simulating it.
}
\value{
This function outputs a single matrix of format similar to valueSignals in the CNOlist but that contains an output for each species in the model.  This matrix is the simulated equivalent of valueSignals at time 1, if you consider only the columns given by indexSignals.
}

\author{
M.K. Morris based on function by C. Terfve from \link{CellNOptR} package.
}
\seealso{
\code{SimulateT1} from CellNOptR
}
\examples{
#This computes the output of the full model, which is normally not done on a stand alone basis, but if you have a model and would like to visualise its output compared to your data, then this is what you should do
    data(CNOlistToy,package="CellNOptR")
    data(ToyModel,package="CellNOptR")

    paramsList <- defaultParametersFuzzy(data=CNOlistToy, model=ToyModel)

    indicesToy<-indexFinder(CNOlistToy,ToyModel,verbose=TRUE)
    ToyFields4Sim<-prep4simFuzzy(ToyModel,paramsList)

    Sim<-simFuzzyT1(
    	CNOlist=CNOlistToy,
    	model=ToyModel,
    	simList=ToyFields4Sim)
}

