\name{gaBinaryTimeScale}
\alias{gaBinaryTimeScale}
\title{
Genetic algorithm for optimizing models using multiple time-points with 1 or 2 steady states
}
\description{
The genetic algorithm used to optimize a model by fitting to data containing multiple time points. The data can be fitted in 2 ways: applying a single scaling factor to the boolean simulation or, additionally, dividing the data with a user-provided time-point into 2 time frames (that signify early and late reactions).
}
\usage{
gaBinaryTimeScale(CNOlist, model, simList, indexList, sizeFac = 1e-04, NApenFac = 1,
initBstring, popSize = 50, pMutation = 0.5, maxTime = 60, maxGens = 500, stallGenMax = 100,
selPress = 1.2, elitism = 5, relTol = 0.1, verbose=TRUE, boolUpdates, divTime=NULL,
lowerB=lowerB, upperB=upperB)
}
\arguments{
\item{CNOlist}{
a CNOlist on which the score is based (based on all valueSignals).
}
\item{model}{
a Model list.
}
\item{simList}{
a list that contains additional fields for the simulator, as created by prep4sim applied to the model above.
}
\item{indexList}{
a list of indexes of species stimulated/inhibited/signals, as produced by indexfinder applied on the model and CNOlist above.
}
\item{sizeFac}{
the scaling factor for the size term in the objective function, default to 0.0001.
}
\item{NApenFac}{
the scaling factor for the NA term in the objective function, default to 1.
}
\item{initBstring}{
an initial bitstring to be tested, should be of the same size as the number of reactions in the model above.
}
\item{popSize}{
the population size for the genetic algorithm, default set to 50.
}
\item{pMutation}{
the mutation probability for the genetic algorithm, default set to 0.5.
}
\item{maxTime}{
the maximum optimisation time in seconds, default set to 60.
}
\item{maxGens}{
the maximum number of generations in the genetic algorithm, default set to 500.
}
\item{stallGenMax}{
the maximum number of stall generations in the genetic algorithm, default to 100.
}
\item{selPress}{
the selective pressure in the genetic algorithm, default set to 1.2.
}
\item{elitism}{
the number of best individuals that are propagated to the next generation in the genetic algorithm, default set to 5.
}
\item{relTol}{
the relative tolerance for the best bitstring reported by the genetic algorithm, i.e. how different from the best solution can solutions be to be reported as well, default set to 0.1.
}
\item{verbose}{
logical (default to TRUE), do you want the statistics of each generation to be printed on the screen?
}
\item{boolUpdates}{
the number of synchronous updates performed by the boolean simulator.
}
\item{divTime}{
the experimental time point that divides the early and late steady states, if using the 2-steady state version.
}
\item{lowerB}{
The lower bound for the optimized value of the scaling factor.}
\item{upperB}{
The upper bound for the optimized value of the scaling factor.}
}
\details{
This is the modified version of gaBinaryT1 from the CellNOptR package, which is able to use all data from CNOlist$valueSignals.
}
\value{
This function returns a list with elements:
\item{bString}{the best bitstring.}
\item{Results}{a matrix with columns: "Generation", "Best_score", "Best_bitString", "Stall_Generation", "Avg_Score_Gen", "Best_score_Gen", "Best_bit_Gen", "Iter_time".}
\item{StringsTol}{the bitstrings whose scores are within the tolerance.}
\item{StringsTolScores}{the scores of the above-mentioned strings.}
}
\references{
J. Saez-Rodriguez, L. G. Alexopoulos, J. Epperlein, R. Samaga, D. A. Lauffenburger, S. Klamt and P. K. Sorger. Discrete logic modeling as a means to link protein signaling networks with functional analysis of
mammalian signal transduction, Molecular Systems Biology, 5:331, 2009.
}
\author{
A. MacNamara
}
\seealso{
getFitTimeScale, simulatorTimeScale
}
\examples{
library(CellNOptR)
library(CNORdt)
data(CNOlistPB, package="CNORdt")
data(ToyModelPB, package="CNORdt")

# pre-process model
indexOrig <- indexFinder(CNOlist=CNOlist, model=Model, verbose=TRUE)
indexNONC <- findNONC(model=Model, indexes=indexOrig, verbose=TRUE)
ModelCut <- cutNONC(model=Model, NONCindexes=indexNONC)
indexNONCcut <- indexFinder(CNOlist=CNOlist, model=ModelCut)
ModelCutCompress <- compressModel(model=ModelCut, indexes=indexNONCcut)
indexNONCcutComp <- indexFinder(CNOlist=CNOlist, model=ModelCutCompress)
ModelCutCompressExpand <- expandGates(model=ModelCutCompress)
fields4Sim <- prep4sim(model=ModelCutCompressExpand)
initBstring <- rep(1, length(ModelCutCompressExpand$reacID))

# optimize
opt1 <- gaBinaryTimeScale(CNOlist=CNOlist, model=ModelCutCompressExpand,
simList=fields4Sim, indexList=indexNONCcutComp,
initBstring=initBstring, verbose=TRUE, boolUpdates=30, maxTime=100, lowerB=0.8, upperB=10)

# plot result
cutAndPlotResultsTimeScale(model=ModelCutCompressExpand, bString=opt1$bString, simList=fields4Sim,
CNOlist=CNOlist, indexList=indexNONCcutComp, boolUpdates=30, lowerB=0.8, upperB=10)
}
