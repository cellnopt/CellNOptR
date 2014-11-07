\name{CNORfuzzy-package}
\alias{CNORfuzzy-package}
\alias{CNORfuzzy}
\docType{package}
\title{
R version of CNOFuzzy, a Constrained Fuzzy Logic Network Optimisation
}
\description{
This package does optimisation of constrained Fuzzy logic networks of signalling
pathways based on a previous knowledge network and a set of data collected upon
perturbation of some of the nodes in the network.
}
\details{
\tabular{ll}{
Package: \tab CNOR\cr
Type: \tab Package\cr
Version: \tab 1.4.0\cr
Date: \tab 2013-08-28\cr
License: \tab GPL-2\cr
LazyLoad: \tab yes\cr
Depends: \tab R (>= 2.15.0), CellNOptR (>= 1.3.29), nloptr (>= 0.8.5) \cr
}
}
\author{
M.K. Morris\cr
Maintainer: T. Cokelaer <cokelaer@ebi.ac.uk>
}
\references{
     \enumerate{
     \item{J. Saez-Rodriguez, L. G. Alexopoulos, J. Epperlein, R. Samaga,
         D. A. Lauffenburger, S. Klamt and P. K. Sorger. Discrete logic modeling as a
        means to link     protein signaling networks with functional analysis of
        mammalian signal transduction, Molecular Systems Biology, 5:331, 2009.}
     \item{Morris MK, Saez-Rodriguez J, Clarke DC, Sorger PK,
        Lauffenburger DA
        (2011). Training Signaling Pathway Maps to Biochemical Data with Constrained
        Fuzzy Logic: Quantitative Analysis of Liver Cell Responses to Inflammatory
        Stimuli. PLoS Comput Biol. 7(3): e1001099.}
        }
}



\keyword{CNORfuzzy, fuzzy logic}
\seealso{
\link{CellNOptR} package.  
}


\examples{

    # Get data from CellNOptR package
    data(CNOlistToy,package="CellNOptR")
    data(ToyModel,package="CellNOptR")

    # Use the default parameters and set Data and Model
    paramsList=defaultParametersFuzzy()
    paramsList$data<-CNOlistToy
    paramsList$model<-ToyModel

\dontrun{
    # Run the simulator
    Res = CNORwrapFuzzy(data=CNOlistToy, model=ToyModel, paramsList=paramsList)
    }


}
