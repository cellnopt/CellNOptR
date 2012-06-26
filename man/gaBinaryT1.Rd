\name{gaBinaryT1}
\alias{gaBinaryT1}
\title{
Genetic algorithm used to optimise a model
}
\description{
This function is the genetic algorithm to be used to optimise a model by fitting to data containing one time point.
}
\usage{
gaBinaryT1(CNOlist, model, simList, indexList, sizeFac = 1e-04, NAFac = 1,
initBstring, popSize = 50, pMutation = 0.5, maxTime = 60, maxGens = 500,
stallGenMax = 100, selPress = 1.2, elitism = 5, relTol = 0.1, verbose=TRUE,
priorBitString=NULL)
}
\arguments{
  \item{CNOlist}{
a CNOlist on which the score is based (based on valueSignals[[2]], i.e. data at t1)
}
  \item{model}{
a model structure, as created by \code{readSIF}, normally pre-processed but that
is not a requirement of this function
}
  \item{simList}{
a list that contains additional fields for the simulator, as created by prep4sim applied to the model above
}
  \item{indexList}{
a list of indexes of species stimulated/inhibited/signals, as produced by indexfinder applied on the model and CNOlist above
}
  \item{sizeFac}{
the scaling factor for the size term in the objective function, default to 0.0001
}
  \item{NAFac}{
the scaling factor for the NA term in the objective function, default to 1
}
  \item{initBstring}{
an initial bitsrting to be tested, should be of the same size as the number of reactions in the model above
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
the number of best individuals that are propagated to the next generation in the genetic algorithm, default set to 5
}
  \item{relTol}{
the relative tolerance for the best bitstring reported by the genetic algorithm,
i.e., how different from the best solution, default set to 0.1
}
  \item{verbose}{
logical (default to TRUE) do you want the statistics of each generation to be printed on the screen?
}
  \item{priorBitString}{
A bitString of same length at the initial bitstring made of 0, 1 or NA. By
default, this bitstring is set to NULL (equivalent to setting all bits to NA).
If provided, all bitstring in a population will be changed to be in agreement
with the priorBitString list.
}
}
\details{
The whole procedure is described in details in Saez-Rodriguez et al. (2009). The basic principle is that at each generation, the algorithm evaluates a population of models based on excluding or including some gates in the initial pre-processed model (this is encoded in a bitstring with contains 0/1 entries for each gate).  The population is then evolved based on the results of the evaluation of these networks, where the evaluation is obtained by simulating the model (to steady state) under the various conditions present in the data, and then computing the squared deviation from the data, to which a penalty is added for size of the model and for species in the model that do not reach steady state.
}
\value{
This function returns a list with elements:
\item{bString}{the best bitstring}
\item{results}{ a matrix with columns "Generation","Best_score","Best_bitString","Stall_Generation","Avg_Score_Gen","Best_score_Gen","Best_bit_Gen","Iter_time"}
\item{stringsTol}{the bitstrings whose scores are within the tolerance}
\item{stringsTolScores}{the scores of the above-mentioned strings}
}
\references{
J. Saez-Rodriguez, L. G. Alexopoulos, J. Epperlein, R. Samaga, D. A. Lauffenburger, S. Klamt and P. K. Sorger. Discrete logic modeling as a means to link protein signaling networks with functional analysis of
mammalian signal transduction, Molecular Systems Biology, 5:331, 2009.
}
\author{
C. Terfve
}
\seealso{
\link{getFit}, \link{prep4sim}, \link{indexFinder}, \link{simulatorT1}
}
\examples{
data(CNOlistToy,package="CellNOptR")
data(ToyModel,package="CellNOptR")

#pre-process model

results = preprocessing(CNOlistToy, ToyModel)
processedModel = results$model
indices = results$indices

#optimise

ToyFields4Sim<-prep4sim(processedModel)
initBstring<-rep(1,length(processedModel))
ToyT1opt<-gaBinaryT1(
	CNOlist=CNOlistToy,
	model=processedModel,
	simList=ToyFields4Sim,
	indexList=indices,
	initBstring=initBstring,
	maxGens=2, popSize=5, verbose=FALSE)
}
