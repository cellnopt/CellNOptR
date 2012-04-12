\name{mode2sif}
\alias{mode2sif}

\title{
  Convert model to sif
}

\description{
  This function converts a network form model format to sif format and saves the sif file.
}

\usage{
mode2sif(Model,optimRes=NA,writeSif=FALSE, filename="Model")
}


\arguments{
  \item{Model}{the model, as created by \link{reasSif}}
  \item{optimRes}{a bit string with the reaction of the model to be considered, default considers all reactions}
  \item{writeSif}{if you want to save the sif file set it to TRUE, otherwise FALSE. Default is FALSE}
  \item{filename}{name of the sif file. Default is Model}
  }
  
\details{
This function takes as input, the model and the bit string and saves the corresponding model in sif format.
}

\value{
\item{sifFile}{the corresponding sif}
}

\author{
    F.Eduati
}

\seealso{
\link{sif2graph}
}


