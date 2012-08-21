\name{CNOlist}
\docType{class}

\alias{CNOlist-class}
\alias{CNOlist}
\title{CNOlist objects}
\description{
This function takes as input the filename of a MIDAS file and returns an
instance of CNOlist class.
}


\usage{
CNOlist(MIDASfile, subfield, verbose=FALSE)
}
\arguments{
  \item{MIDASfile}{
filename of a MIDAS file
}
  \item{subfield}{
TRUE or FALSE, specifies if the column headers contain subfields or not i.e. if headers are of the form TR:sthg:sthg or just TR:sthg
}
  \item{verbose}{
(default to FALSE)
}
}
\details{
This function is equivalent to reading the filename with \link{readMIDAS} and
converting with \link{makeCNOlist} and creating a CNOlist class. See the two functions above for details.

This class is not yet used but new functionalities (aftr version 1.3.28) should use it.
}
\value{
An instance of CNOlist class. The methods associated to this class are called cues, inhibitors, stimuli, timepoints, signals.
}
\author{
T. Cokelaer
}
\seealso{
\link{readMIDAS}, \link{makeCNOlist}
}
\examples{
files<-dir(system.file("ToyModel",package="CellNOptR"),full=TRUE)
cnolist  = CNOlist(files[[1]])
cues(cnolist)
}

