\name{interpretDiscreteGA}
\alias{interpretDiscreteGA}

\title{Interpreter of output of discrete genetic algorithm}

\description{This function takes the integer string output by the discrete
genetic algorithm for training a cFL model and generates the corresponding model
based on the Fuzzy parameters.}

\usage{
interpretDiscreteGA(model,  paramsList, intString, bitString=NULL)
}

\arguments{
  \item{model}{
PKN trained (same model input as to gaDiscrete).
}

  \item{paramsList}{
    List containing parameters (see defaultParametersFuzzy). Only the fuzzy
    parameters are used.
}
  \item{intString}{
    Integer string resulting from gaDiscrete (in bString field of gaDiscreteT1
    output).
}
  \item{bitString}{
    (optional) if you want to cut additional interactions from the model.  Used
    in reduceFuzzy function.
}
}

\details{
After the discreteGA chooses transfer functions from a discrete set of transfer
functions and removing interactions inconsistent with the data, this function
interprets the output and returns an actual model using these transfer functions
as well as a model from which logical redundancy was cut.
}

\value{

\item{model}{The selected initial model based on the provided bitstring.}
\item{simList}{The corresponding data related to Model field}
\item{bitString}{The bitstring corresponding to the Model field}
\item{cutModel}{Same as Model but redundant reactions are also removed.}
\item{cutSimList}{The corresponding data related to cutModel field}
\item{cutBitString}{The corresponding bitstring related to cutModel field.}

}


\author{M.K. Morris}

\seealso{
    \link{gaDiscreteT1}
}

\examples{

    data(ToyModel, package="CellNOptR")
    data(CNOlistToy,package="CellNOptR")
    paramsList = defaultParametersFuzzy()
    \dontrun{

        # preprocessing (see CNORwrapFuzzy or gaDiscreteT1)
        T1opt = gaDiscreteT1(....)  # see CNORwrapFuzzy or gaDiscreteT1 for details
        interpretDiscreteGA(ToyModel, paramsList, T1opt$bString)
}
}

