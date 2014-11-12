\name{CNORode}

\alias{CNORode}

\docType{package}

\title{
	Logic based ODE extension for CellNOptR
}
\description{
	This package is used for the simulation and fitting of
	logic based ODE models based on the Odefy approach.
}

\details{
	\tabular{ll}{
		Package: \tab CNORode\cr
		Type: \tab Package\cr
		Version: \tab 1.2.0\cr
		Date: \tab 2012-03-14\cr
		License: \tab GPL-3\cr
		LazyLoad: \tab yes\cr
	}
}

\author{
	David Henriques, Thomas Cokelaer
	Maintainer: David Henriques <dhenriques@ebi.ac.uk>
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

    C. Terfve, T. Cokelaer, A. MacNamara, D. Henriques, E. Goncalves, MK.
    Morris, M. van Iersel, DA Lauffenburger, J. Saez-Rodriguez.
    CellNOptR: a flexible toolkit to train protein signaling networks to data using
    multiple logic formalisms. BMC Systems Biology, 2012, 6:133:
}
\keyword{CNORode}

\seealso{
	\code{\link{CellNOptR}},
	\code{\link{parEstimationLBode}},
	\code{\link{getLBodeModelSim}},
	\code{\link{parEstimationLBode}} 
	\code{\link{plotLBodeFitness}}.
}
