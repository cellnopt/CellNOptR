\name{parEstimationLBode}
\alias{parEstimationLBode}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Perform parameter estimation using a genetic algorithm (package genalg) or ssm
(if package essm available).
}
\description{
	This function is an alias to the parEstimationLBode variants
(\code{\link{parEstimationLBodeGA}} and \code{\link{parEstimationLBodeSSm}})
}
\usage{
parEstimationLBode(cnolist, model, method="ga", ode_parameters = NULL, indices = NULL,
paramsGA=NULL, paramsSSm=NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{cnolist}{
A list containing the experimental design and data.
}
  \item{model}{
The logic model to be simulated.
}

  \item{method}{
Only "ga" or "ssm" arguments are accepted.
}

  \item{ode_parameters}{
A list with the ODEs parameter information. Obtained with \code{\link{createLBodeContPars}}.
}
  \item{indices}{
Indices to map data in the model. Obtained with indexFinder function from CellNOptR. 
}
\item{paramsGA}{
 A list of GA parameters. default is the list returned by \code{defaultParametersGA}.
 } 
\item{paramsSSm}{
 A list of SSm parameters. default is the list returned by\code{defaultParametersSSm}.
 } 
}
\value{
 \item{LB_n}{
A numeric value to be used as lower bound for all parameters of type n.
}
  \item{LB_k}{
A numeric value to be used as lower bound for all parameters of type k.
}

  \item{LB_tau}{
A numeric value to be used as lower bound for all parameters of type tau.
}

  \item{UB_n}{
A numeric value to be used as upper bound for all parameters of type n.
}

  \item{UB_k}{
A numeric value to be used as upper bound for all parameters of type k.
}

  \item{UB_tau}{
A numeric value to be used as upper bound for all parameters of type tau. 
}

  \item{default_n}{
The default parameter to be used for every parameter of type n. 
}
  \item{default_k}{
The default parameter to be used for every parameter of type k. 
}
  \item{default_tau}{
The default parameter to be used for every parameter of type tau. 
}
  \item{LB_in}{
An array with the the same length as ode_parameters$parValues with lower bounds
for each specific parameter. 
}
  \item{UB_in}{
An array with the the same length as ode_parameters$parValues with upper bounds
for each specific parameter.
}
  \item{opt_n}{
Add all parameter n to the index of parameters to be fitted.
}
  \item{opt_k}{
Add all parameter k to the index of parameters to be fitted.
}
  \item{opt_tau}{
Add all parameter tau to the index of parameters to be fitted.
}
  \item{random}{
A logical value that determines that a random solution is for the parameters 
to be optimized.
}

\item{res}{A list containing the information provided by the solver.}

}

\references{
Dominik Wittmann, Jan Krumsiek, Julio S. Rodriguez, Douglas
Lauffenburger, Steffen Klamt, and Fabian Theis. Transforming boolean
models to continuous models: methodology and application to t-cell
receptor signaling. BMC Systems Biology, 3(1):98+, September 2009.

Egea, J.A., Maria, R., Banga, J.R. (2010) An evolutionary method for
complex-process optimization. Computers & Operations Research 37(2):
315-324.
 
Egea, J.A., Balsa-Canto, E., Garcia, M.S.G., Banga, J.R. (2009) Dynamic
optimization of nonlinear processes with an enhanced scatter search
method. Industrial & Engineering Chemistry Research 49(9): 4388-4401. 

Jan Krumsiek, Sebastian Polsterl, Dominik Wittmann, and Fabian
Theis. Odefy - from discrete to continuous models. BMC
Bioinformatics, 11(1):233+, 2010.

R. Serban and A. C. Hindmarsh, "CVODES: the Sensitivity-Enabled ODE Solver in SUNDIALS,"
Proceedings of IDETC/CIE 2005, Sept. 2005, Long Beach, CA. Also available as LLNL
technical report UCRL-JP-200039.     
}


\author{
David Henriques, Thomas Cokelaer
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link{CellNOptR}}
\code{\link{createLBodeContPars}}
\code{\link{rbga}}
}

\examples{
data("ToyCNOlist",package="CNORode");
data("ToyModel",package="CNORode");
data("ToyIndices",package="CNORode");

ode_parameters=createLBodeContPars(model,random=TRUE);
%Visualize intial simulation
%simulatedData=plotLBodeFitness(cnolist, model,ode_parameters,indices=indices)
paramsGA = defaultParametersGA()
paramsGA$maxStepSize = 1
paramsGA$popSize = 20
paramsGA$iter = 50
paramsGA$transfer_function = 2

ode_parameters=parEstimationLBode(cnolist,model,ode_parameters=ode_parameters,
    paramsGA=paramsGA)


}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{genetic}
\keyword{algorithm}% __ONLY ONE__ keyword per line
\keyword{CellNOptR}
\keyword{logic}
\keyword{model}
