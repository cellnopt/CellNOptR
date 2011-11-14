\name{get_logic_based_ode_model_simulation}
\alias{get_logic_based_ode_model_simulation}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Simulate the logic-based ODE model
}
\description{
This function simulates a logic-based ODE model and return a list with one matrix for each time point.
The size of each matrix is \( N_{experiments} \times N_{N\_model\_species} \). The input species in the
model are filled with NA values. If the simulation of a particular set of initial conditions fails
the solver will fill the experience row with NA values.
\usage{
get_logic_based_ode_model_simulation(cnolist, model, ode_parameters = NULL, \\
indices = NULL, time = 1, verbose = 0, transfer_function = 3, reltol = 1e-04,\\
 atol = 0.001, maxStepSize = Inf, maxNumSteps = 1e+05, maxErrTestsFails = 50)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{cnolist}{
A list containing the experimental design.
}
  \item{model}{
The logic model to be simulated.
}
  \item{ode_parameters}{
A list with the ODEs parameter information. Obtained with makeParameterList function.
}
  \item{indices}{
Indices to map data in the model. Obtained with indexFinder function from CellNOptR. 
}
  \item{time}{
An integer with the index of the time point to start the simulation. Default is 1.
}
  \item{verbose}{
A logical value. This will trigger a set of comments inside the C code.
}
  \item{transfer_function}{
The type of used transfer. Use 1 for no transfer function, 2 for Hill function and
	3 for normalized Hill function.
}
  \item{reltol}{
The maximum step size allowed to ODE solver.
}
  \item{atol}{
Relative Tolerance for numerical integration.
}
  \item{maxStepSize}{
The maximum number of internal steps between two points being sampled before
	the solver fails.
}
  \item{maxNumSteps}{
Specifies the maximum number of error test failures permitted in attempting one step.
}
  \item{maxErrTestsFails}{
Specifies the maximum number of error test failures permitted in attempting one step.
}
}
\details{
For more details in the configuration of the ODE solver check the CVODES manual. 
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
  \item{1 }{Matrix containing simulation for time-point 1}
  \item{3 }{Matrix containing simulation for time-point 2}
  \item{\( {N_{time\_points}} \) }{Matrix containing simulation for time-point 3}
%% ...
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
    Trans. Math. Softw., 31(3):363–396, September 2005.
    
}
\author{
%%  ~~who you are~~
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as
function
(
		cnolist,				model,					ode_parameters=NULL,	
		indices=NULL,			time=1,					verbose=0, 
		transfer_function=3,	reltol=1e-4,			atol=1e-3,	
		maxStepSize=Inf,		maxNumSteps=100000,		maxErrTestsFails=50
)
{
	adjMat=incidence2Adjacency(model);
	if(is.null(indices))indices <- indexFinder(cnolist,model,verbose=FALSE);
	if(is.null(ode_parameters))ode_parameters=makeParameterList(adjMat,model$namesSpecies);
	sim_function=get_simulation_function(cnolist,model,adjMat,
			indices, ode_parameters1$parValues, time,verbose,
			transfer_function,reltol,atol,maxStepSize,
			maxNumSteps,maxErrTestsFails);
	return(sim_function(cnolist,model,ode_parameters$parValues));
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
