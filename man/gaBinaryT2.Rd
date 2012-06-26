\name{gaBinaryT2}
\alias{gaBinaryT2}
\title{
Genetic algorithm for time point 2
}
\description{
This is the genetic algorithm for time point 2, that should follow optimisation based on time point 1.
}
\usage{
gaBinaryT2(CNOlist, model, simList, indexList, bStringT1, simResT1, sizeFac = 1e-04, NAFac = 1,
popSize = 50, pMutation = 0.5, maxTime = 60, maxGens = 500,
stallGenMax = 100, selPress = 1.2, elitism = 5, relTol = 0.1, verbose=TRUE,
priorBitString=NULL)
}
\arguments{
  \item{CNOlist}{
a CNOlist on which the score is based (based on valueSignals[[3]], i.e. data at  t2)
}
  \item{model}{
a model structure, as created by \code{readSIF}, normally pre-processed but that
is not a requirement of this function
}
  \item{simList}{a simList structure as created by \code{prep4sim}.
}
  \item{indexList}{
a list of indexes of species stimulated/inhibited/signals, as produced by indexfinder applied on teh model and CNOlist above
}
  \item{bStringT1}{
the optimal bitstring from optimisation at time 1 (i.e. a vector of 0s and 1s)
}
  \item{simResT1}{
the results of optimisation at t1 (i.e. the whole list as output by gabinaryT1)
}
  \item{sizeFac}{
the scaling factor for the size term in the objective function, default to 0.0001
}
  \item{NAFac}{
the scaling factor for the NA term in the objective function, default to 1
}
  \item{popSize}{
the population size for the genetic algorithm, default set to 50
}
  \item{pMutation}{
the mutation probability for the genetic algorithm, default set to 0.5
}
  \item{maxTime}{
the maximum optimisation time in seconds, default set to 60
}
  \item{maxGens}{
the maximum number of generations in the genetic algorithm, default set to 500
}
  \item{stallGenMax}{
the maximum number of stall generations in the genetic algorithm, default to 100
}
  \item{selPress}{
the selective pressure in the genetic algorithm, default set to 1.2
}
  \item{elitism}{
the number of best individuals that are propagated to the next generation in the gen. al. default set to 5
}
  \item{relTol}{
the relative tolerance for the best bitstring reported by the genetic algorithm, i.e.how different from the best solution can solutions be to be reported as well, default set to 0.1
}
  \item{verbose}{
logical (default to TRUE) do you want the statistics of each generation to be printed on the screen
}
\item{priorBitString}{
A bitString of same length at the initial bitstring made of 0, 1 or NA. By
default, this bitstring is set to NULL (equivalent to setting all bits to NA).
If provided, all bitstring in a population will be changed to be in agreement
with the priorBitString list.
}

}
\details{
This function takes in the same input as the T1 ga, but in addition it takes in the bitstring optimised for T1, and does not take an initial bitstring. It also takes a matrix simResT1, which is the results of simulatorT1 ran on the model optimised for T1. Be aware that the bitString that this function returns is one that only includes the bits that it actually looks at, i.e. the bits that were 0 in the bStringT1
}
\value{
This function returns a list with elements:
\item{bString}{the best bitstring}
\item{results}{ a matrix with columns "Generation","Best_score","Best_bitString","Stall_Generation","Avg_Score_Gen","Best_score_Gen","Best_bit_Gen","Iter_time"}
\item{stringsTol}{the bitstrings whose scores are within the tolerance}
\item{stringsTolScores}{the scores of the above-mentionned strings}
}
\author{
C.Terfve
}
\seealso{
\link{getFit}, \link{prep4sim}, \link{indexFinder}, \link{simulatorT1},
\link{simulatorT2}, \link{gaBinaryT2}

}
\examples{
#load data

data(CNOlistToy2,package="CellNOptR")
data(ToyModel2,package="CellNOptR")

#pre-process model

checkSignals(CNOlistToy2,ToyModel2)
results = preprocessing(CNOlistToy2, ToyModel2)
processedModel = results$model
indices = results$indices

#optimise t1
ToyFields4Sim2<-prep4sim(processedModel)
initBstring2<-rep(1,length(processedModel$reacID))
ToyT1opt2<-gaBinaryT1(
	CNOlist=CNOlistToy2,
	model=processedModel,
	simList=ToyFields4Sim2,
	indexList=indices,
	initBstring=initBstring2,
	maxGens=2,
	popSize = 5,
	verbose=FALSE)

#Optimise T2
SimToyT12<-simulateT1(
	CNOlist=CNOlistToy2,
	model=processedModel,
	bStringT1=ToyT1opt2$bString,
	simList=ToyFields4Sim2,
	indexList=indices)

ToyT1opt2T2<-gaBinaryT2(
	CNOlist=CNOlistToy2,
	model=processedModel,
	simList=ToyFields4Sim2,
	indexList=indices,
	bStringT1=ToyT1opt2$bString,
	simResT1=SimToyT12,
	maxGens=2,
	popSize = 5,
	verbose=FALSE)
}
