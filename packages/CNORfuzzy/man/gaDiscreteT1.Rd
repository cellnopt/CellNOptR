\name{gaDiscreteT1}
\alias{gaDiscreteT1}
\title{
Genetic algorithm used to optimise a cFL model
}
\description{
This function is the genetic algorithm to be used to optimise a cFL model by fitting to data containing one time point.
}
\usage{
gaDiscreteT1(CNOlist, model, paramsList, initBstring=NULL, sizeFac=0.0001,
    NAFac=1, popSize=50, pMutation=0.5, maxTime=60, maxGens=500,
    stallGenMax=100, selPress=1.2, elitism=5, relTol=0.1,
    verbose=TRUE,maxSizeHashTable = 1000)
}
\arguments{
  \item{CNOlist}{
a CNOlist on which the score is based (based on valueSignals[[2]], i.e. data at t1)
}
  \item{model}{
a Model list
}
  \item{paramsList}{
CellNOptR software parameters (this functions uses transfer functions to choose from)
}
  \item{initBstring}{
an initial bitsring to be tested, should be of the same size as the number of
reactions in the model above.
}
  \item{sizeFac}{
the scaling factor for the size term in the objective function, default to 0.0001
}
  \item{NAFac}{
the scaling factor for the NA term in the objective function, default to 1.
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
the maximum number of generations in the genetic algorithm, default set to 500.
}
  \item{stallGenMax}{
the maximum number of stall generations in the genetic algorithm, default to
100.
}
  \item{selPress}{
the selective pressure in the genetic algorithm, default set to 1.2.
}
  \item{elitism}{
the number of best individuals that are propagated to the next generation in the
genetic algorithm, default set to 5.
}
  \item{relTol}{
the relative tolerance for the best bitstring reported by the genetic algorithm,
i.e.how different from the best solution can solutions be to be reported as
well, default set to 0.1.
}
  \item{verbose}{
logical (default to TRUE) do you want the statistics of each generation to be printed on the screen?
}

\item{maxSizeHashTable}{
a hash table is use to store bitstring and related
score. This allows the GA to be very efficient is the case of small models. The
size of the hash table is 5000 by default, which may be too large for large
models
}

}
\details{
The GA procedure is implemented based on the gaBinaryT1 in CellNOptR (see those
man pages for a basic description). Necessary extensions to optimize a string of
numbers rather than zero and one have been made.  Additionally, since the
scoring function is defined inside the function, it has also been altered for
this function.

The parameters are similar to those used in CellNOptR and the returned list
contains similar results as well. 
}
\value{
This function returns a list with elements:
\item{bString}{The best bitstring}
\item{results}{A matrix with columns "Generation", "Best_score",
"Best_bitString", "Stall_Generation", "Avg_Score_Gen", "Best_score_Gen"
,"Best_bit_Gen", "Iter_time".}
\item{stringsTol}{The bitstrings whose scores are within the tolerance}
\item{stringsTolScores}{The scores of the above-mentioned strings}
}
\author{
M. Morris based on gaBinaryT1 function by C. Terfve (CellNOptR package)
}
\seealso{
\link{prep4simFuzzy}, \link{simFuzzyT1}
}
\examples{

    data(CNOlistToy,package="CellNOptR")
    data(ToyModel,package="CellNOptR")

    #pre-process model
    model <- preprocessing(CNOlistToy, ToyModel, verbose=FALSE)

    #set parameters
    paramsList <- defaultParametersFuzzy(CNOlistToy, ToyModel)


    # the GA algorithm 
    ToyT1opt<-gaDiscreteT1(
    	CNOlist=CNOlistToy,
    	model=model,
    	paramsList=paramsList,
    	maxTime=3,
    	verbose=FALSE)
}
