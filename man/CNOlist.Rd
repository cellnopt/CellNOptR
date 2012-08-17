\name{CNOlist}
\alias{CNOlist}
\title{
Create an instance of CNOlist class.
}
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
TRUE or FALSE, specifies if the column headers contain subfields or not i.e. if I should look for TR:sthg:sthg or just TR:sthg
}
  \item{verbose}{
(default to FALSE)
}
}
\details{
This function is equivalent to reading the filename with \link{readMIDAS} and
converting with \link{makeCNOlist} and creating a CNOlist class. 

See the two functions above for details.
}
\value{
a CNOlist with fields
\item{cues}{a vector of names of cues}
\item{timuli}{a vector of names of stimuli}
\item{Inhibitors}{a vector of names of inhibitors}
\item{Signals}{a vector of names of signals}
\item{timeSignals}{a vector of times}
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

