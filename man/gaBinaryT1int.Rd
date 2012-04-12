\name{gaBinaryT1int}
\alias{gaBinaryT1int}
\title{
Genetic algorithm used to optimise a model
}
\description{
This function is the genetic algorithm to be used to optimise a model by fitting to data containing one time point. It is the function \link{gaBinaryT1} of CellNOptR modified in orter to differently  weights for the integrated links
}
\usage{
gaBinaryT1int(CNOlist, Model, SimList, indexList, sizeFac = 1e-04,
              integrFac=10, NAFac = 1, initBstring, PopSize = 50,
              Pmutation = 0.5, MaxTime = 60, maxGens = 500,
              StallGenMax = 100, SelPress = 1.2, elitism = 5,
              RelTol = 0.1, verbose=TRUE)
}
\arguments{
  \item{CNOlist}{
a CNOlist on which the score is based (based on valueSignals[[2]], i.e. data at t1)
}
  \item{Model}{
a Model list
}
  \item{SimList}{
a list that contains additional fields for the simulator, as created by prep4Sim applied to the model above
}
  \item{indexList}{
a list of indexes of species stimulated/inhibited/signals, as produced by indexfinder applied on the model and CNOlist above
}
  \item{sizeFac}{
the scaling factor for the size term in the objective function, default to 0.0001
}
  \item{integrFac}{
the scaling factor for the integration term in the objective function, default to 10
}

  \item{NAFac}{
the scaling factor for the NA term in the objective function, default to 1
}
  \item{initBstring}{
an initial bitsring to be tested, should be of the same size as the number of reactions in the model above
}
  \item{PopSize}{
the population size for the genetic algorithm, default set to 50
}
  \item{Pmutation}{
the mutation probability for the genetic algorithm, default set to 0.5
}
  \item{MaxTime}{
the maximum optimisation time in seconds, default set to 60
}
  \item{maxGens}{
the maximum number of generations in the genetic algorithm, default set to 500
}
  \item{StallGenMax}{
the maximum number of stall generations in the genetic algorithm, default to 100
}
  \item{SelPress}{
the selective pressure in the genetic algorithm, default set to 1.2
}
  \item{elitism}{
the number of best individuals that are propagated to the next generation in the genetic algorithm, default set to 5
}
  \item{RelTol}{
the relative tolerance for the best bitstring reported by the genetic algorithm,
i.e., how different from the best solution, default set to 0.1
}
  \item{verbose}{
logical (default to TRUE) do you want the statistics of each generation to be printed on the screen?
}
}
\details{
The whole procedure is described in details in Saez-Rodriguez et al. (2009), see \link{gaBinaryT1} for detaied description. The only additional input is integrFac, that is used to penalize more links inferred strictly from data using reverse-engineering approaches and integrated using \link{MapBTables2Model} or \link{MapDDN2Model}. If the input model has the field IntegrPPIscores, as given by \link{PPIweight}, this is used to differently penalize integrated links in the computation of the score for each element of the population.
}
\value{
This function returns a list with elements:
\item{bString}{the best bitstring}
\item{Results}{ a matrix with columns "Generation", "Best_score", "Best_bitString", "Stall_Generation", "Avg_Score_Gen", "Best_score_Gen", "Best_bit_Gen", "Iter_time"}
\item{StringsTol}{the bitstrings whose scores are within the tolerance}
\item{StringsTolScores}{the scores of the above-mentioned strings}
}
\references{
J. Saez-Rodriguez, L. G. Alexopoulos, J. Epperlein, R. Samaga, D. A. Lauffenburger, S. Klamt and P. K. Sorger. Discrete logic modeling as a means to link protein signaling networks with functional analysis of
mammalian signal transduction, Molecular Systems Biology, 5:331, 2009.
}
\author{
    F.Eduati
}
\seealso{
\link{prep4Sim}, \link{indexFinder}, \link{MapDDN2Model}, \link{MapBTables2Model} \link{PPIweight}
}
\examples{
library(CellNOptR)
data(CNOlistDREAM,package="CellNOptR")
data(DreamModel,package="CellNOptR")
data(UniprotIDdream,package="CNORfeeder")
data(PPINigraph,package="CNORfeeder")
res<-preprocessing(Data=CNOlistDREAM, Model=DreamModel)
Model <- res$model
    
BTable <- makeBTables(CNOlist=CNOlistDREAM, k=2, measErr=c(0.1, 0))
modelIntegr <- MapBTables2Model(BTable=BTable,Model=Model,allInter=TRUE)

# the followig step may take a while
\dontrun{
resPPIweight <- PPIweight(modelIntegr=modelIntegr,PKNmodel=DreamModel,
                          CNOlist=CNOlistDREAM,UniprotID=UniprotIDdream,
                          PPINigraph=PPINigraph)
modelIntegr <- resPPIweight$modelIntegr}

DreamFields4Sim <- prep4Sim(modelIntegr)
initBstring <- rep(1,length(modelIntegr$reacID))
DreamT1opt <- gaBinaryT1int(
	CNOlist=CNOlistDREAM,
	Model=modelIntegr,
	SimList=DreamFields4Sim,
	indexList=res$indices,
	initBstring=initBstring,
	maxGens=2,
	PopSize=5,
	verbose=FALSE)
}
