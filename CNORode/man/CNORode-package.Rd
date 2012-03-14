\name{CNORode}
\alias{CNORode}
\docType{package}
\title{
Logic based ODE extension for CellNOptR
}
\description{
This package is used for the simulation and fitting of
logic based ODE models. 
}
\details{
\tabular{ll}{
Package: \tab CNORode\cr
Type: \tab Package\cr
Version: \tab 1.0\cr
Date: \tab 2011-11-03\cr
License: \tab GPL-2\cr
LazyLoad: \tab yes\cr
}

}
\author{
David Henriques, Thomas Cokelaer

Maintainer: David Henriques <dhenriques@ebi.ac.uk>

}

\references{
[1] Maria Rodriguez-Fernandez, Jose A. Egea, and Julio R. Banga. Novel
    metaheuristic for parameter estimation in nonlinear dynamic biological
    systems. BMC Bioinformatics, 7:483+, November 2006.

[2] Dominik Wittmann, Jan Krumsiek, Julio S. Rodriguez, Douglas
    Lauffenburger, Steffen Klamt, and Fabian Theis. Transforming boolean
    models to continuous models: methodology and application to t-cell
    receptor signaling. BMC Systems Biology, 3(1):98+, September 2009.

[3] Jan Krumsiek, Sebastian Polsterl, Dominik Wittmann, and Fabian
    Theis. Odefy - from discrete to continuous models. BMC
    Bioinformatics, 11(1):233+, 2010.
    
[4] Alan C. Hindmarsh, Peter N. Brown, Keith E. Grant, Steven L. Lee,
	Radu Serban, Dan E. Shumaker, and Carol S. Woodward. SUNDIALS:
    Suite of nonlinear and differential/algebraic equation solvers. ACM
    Trans. Math. Softw., 31(3):363-396, September 2005.
}
\keyword{CNORode}
\seealso{
~~ Optional links to other man pages, e.g. ~~
~~ \code{\link{CellNOptR}} ~~
}
\examples{
library(CNORode)
data("ToyCNOlist",package="CNORode");
data("ToyModel",package="CNORode");
data("ToyIndices",package="CNORode");

ode_parameters=createLBodeContPars(model,random=TRUE);
#Visualize intial simulation
#simulatedData=plotLBodeFitness(cnolist, model,ode_parameters,indices=indices)
params = defaultParametersGA()
params$maxtime = 10
params$ndiverse = 10
params$dim_refset = 6
ode_parameters=parEstimationLBode(cnolist,model,ode_parameters=ode_parameters,
	indices=indices,paramsGA=params);
#Visual solution after optimization
#simulatedData=plotLBodeFitness(cnolist, model,indices=indices);
}
