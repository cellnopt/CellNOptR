\name{createLBodeContPars}
\alias{createLBodeContPars}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
	Create a list with ODE parameter information needed to perform parameter estimation
}
\description{
	Creates a list with the needed the continuous parameters to simulate the model,
	upper and lower bounds for the parameter estimation, parameters names, indices of
	the parameters and others.
}
\usage{
createLBodeContPars(model, LB_n = 1, LB_k = 0.1, LB_tau = 0.01, 
UB_n = 5, UB_k = 0.9, UB_tau = 10, default_n = 3, default_k = 0.5, default_tau = 1,
LB_in = c(), UB_in = c(), opt_n = TRUE, opt_k = TRUE, opt_tau = TRUE, random = FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{model}{
The logic model to be simulated.
}
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
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{

\item{parNames}{An array containing the names of the parameters.}
\item{parValues}{An array containing the values of the parameters, in the same order as the names.}
\item{index_opt_pars}{An array containing the indexes for the parameters to be fitted.}
\item{index_n}{An array containing the indexes of the parameters of type n.}
\item{index_k}{An array containing the indexes of the parameters of type k.}
\item{index_tau}{An array containing the indexes of the parameters of type tau.}
\item{LB}{An array containing the lower bound for each parameter.}
\item{UB}{An array containing the upper bound for each parameter.}
%% ...
}

\author{
David Henriques, Thomas Cokelaer
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

\examples{

data("ToyCNOlist",package="CNORode");
data("ToyModel",package="CNORode");
data("ToyIndices",package="CNORode");
ode_parameters=createLBodeContPars(model, opt_n=FALSE,default_n=2,
random=TRUE,LB_k=0.25,UB_k=0.8,LB_tau=0.01,UB_tau=10);

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
