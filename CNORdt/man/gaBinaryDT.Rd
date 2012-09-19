\name{gaBinaryDT}
\alias{gaBinaryDT}
\title{
Genetic algorithm for optimizing models using multiple time-points
}
\description{
The genetic algorithm used to optimize a model by fitting to data consisting of multiple time points. The data can be fitted by applying a single scaling factor to the boolean simulation.
}
\usage{
gaBinaryDT(CNOlist, model, initBstring = NULL, sizeFac = 1e-04, NAFac = 1, popSize = 50,
pMutation = 0.5, maxTime = 60, maxGens = 500, stallGenMax = 100, selPress = 1.2,
elitism = 5, relTol = 0.1, verbose = TRUE, priorBitString = NULL, maxSizeHashTable = 5000,
boolUpdates, lowerB = lowerB, upperB = upperB)
}
\arguments{
\item{CNOlist}{
A CNOlist on which the score is based.
}
\item{model}{
A model structure, as created by \code{readSIF}, normally pre-processed but that
is not a requirement of this function.
}
\item{initBstring}{
An initial bitstring to be tested, should be of the same size as the number of
reactions in the model above (model$reacID). The default is all ones.
}
\item{sizeFac}{
The scaling factor for the size term in the objective function, it defaults to 0.0001.
}
\item{NAFac}{
The scaling factor for the NA term in the objective function, it defaults to 1.
}
\item{popSize}{
The population size for the genetic algorithm, it is set to 50.
}
\item{pMutation}{
the mutation probability for the genetic algorithm, default set to 0.5.
}
\item{maxTime}{
the maximum optimisation time in seconds, default set to 60.
}
\item{maxGens}{
The maximum number of generations in the genetic algorithm, default set to 500.
}
\item{stallGenMax}{
The maximum number of stall generations in the genetic algorithm, default to 100.
}
\item{selPress}{
The selective pressure in the genetic algorithm, default set to 1.2.
}
\item{elitism}{
The number of best individuals that are propagated to the next generation in the genetic algorithm, default set to 5.
}
\item{relTol}{
The relative tolerance for the best bitstring reported by the genetic algorithm,
i.e., how different from the best solution, default set to 0.1.
}
\item{verbose}{
Logical (default to TRUE): do you want the statistics of each generation to be printed on the screen?
}
\item{priorBitString}{
At each generation, the GA algorithm creates a population of bitstrings that
will be used to perform the optimisation. If the user knows the values of some
bits, they can be used to overwrite bit values proposed by the GA algorithm. If
provided, the priorBitString must have the same length as the initial bitstring and
be made of 0, 1 or NA (by default, this bitstring is set to NULL, which is
equivalent to setting all bits to NA). Bits that are set
to 0 or 1 are used to replace the bits created by the GA itself (see example).
}
\item{maxSizeHashTable}{
A hash table is used to store bitstrings and related scores. This allows the GA to be very efficient is the case of small models. The
size of the hash table is 5000 by default, which may be too large for large
models. 
}
\item{boolUpdates}{
The number of synchronous updates performed by the boolean simulator.
}
\item{lowerB}{
The lower bound for the optimized value of the scaling factor.
}
\item{upperB}{
The upper bound for the optimized value of the scaling factor.
}
}
\details{
This is the modified version of gaBinaryT1 from the CellNOptR package, which is able to use all data from CNOlist$valueSignals.
}
\value{
This function returns a list with elements:
\item{bString}{The best bitstring.}
\item{Results}{A matrix with columns: "Generation", "Best_score", "Best_bitString", "Stall_Generation", "Avg_Score_Gen", "Best_score_Gen", "Best_bit_Gen", "Iter_time".}
\item{StringsTol}{The bitstrings whose scores are within the tolerance.}
\item{StringsTolScores}{The scores of the above-mentioned strings.}
}
\references{
A. MacNamara, C. Terfve, D. Henriques, B. Penalver Bernabe and J. Saez-Rodriguez, 2012. State-time spectrum of signal transduction logic models. Physical biology, 9(4), p.045003.
}
\author{
A. MacNamara
}
\seealso{
getFitDT, simulatorDT
}
\examples{
library(CellNOptR)
library(CNORdt)
data(CNOlistPB, package="CNORdt")
data(modelPB, package="CNORdt")

# pre-process model
model = preprocessing(CNOlistPB, modelPB)

# optimise
initBstring <- rep(1, length(model$reacID))
opt1 <- gaBinaryDT(CNOlist=CNOlistPB, model=model, initBstring=initBstring,
verbose=TRUE, boolUpdates=10, maxTime=30, lowerB=0.8, upperB=10)
}
