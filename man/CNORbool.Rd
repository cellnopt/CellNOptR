\name{CNORbool}
\alias{CNORbool}
\title{
CNOR boolean analysis standalone
}
\description{
This function is a wrapper around the whole CNOR analysis, it performs the following steps:

 \enumerate{
  \item Plot the CNOlist;
  \item Call the preprocessing function (Cut the nonc off the model, compress the model and expand the gates);
  \item Optimisation T1 and all TN if providedin MIDAS file. 
  \item Write a HTML report (optional)
 }
}
\usage{
CNORbool(CNOlist, model, paramsList=defaultParameters(),
    compression=TRUE, expansion=TRUE, cutNONC=TRUE, verbose=FALSE)


}

\arguments{
\item{CNOlist}{

}

\item{model}{

}

  \item{paramsList}{
    A list of parameters related to the Genetic Algorithm parameters:
   \enumerate{
     \item sizeFac: default to 1e-04;
     \item NAFac: default to 1;popSize: default to 50;
     \item pMutation: default to 0.5;
     \item maxTime: default to 60;
     \item maxGens: default to 500;
     \item stallGenMax: default to 100;
     \item selPress: default to 1.2;
     \item elitism: default to 5;
     \item relTol: default to 0.1;
     \item verbose: default to FALSE (default to true in the functions used by CNORbool but CNORbool sets them to false by default).
   }
  

  If paramsList is not provided (NA), it is filled internally with the \link{defaultParameters} function.
}
\item{compression}{
compress the model (default TRUE)
}
\item{expansion}{
expand the gates (default TRUE)

}
\item{cutNONC}{
cut the NONC nodes off the model by (default TRUE)
}
\item{verbose}{
FALSE}
}
\details{
If you do not provide a parameters list, you can provide only essential elements,
and all other parameters will be set to their default values. In this case, you
should set paramsList=NA, and provide the following fields: data, model, name, time.
}
\value{
This function does not return anything, it does the analysis, produces all the
plots and puts them in a folder that is in your working directory, and is called "Name".
}
\author{
T. Cokelaer
}
\examples{
#version with paramslist


data(CNOlistToy,package="CellNOptR")
data(ToyModel,package="CellNOptR")

#CNORbool(CNOlist=CNOlistToy, model=ToyModel)

#dontrun
#CNORbool(pathto_MIDASfile, pathto_SIFfile)



}
