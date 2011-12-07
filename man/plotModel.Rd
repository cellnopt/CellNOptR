\name{plotModel}
\alias{plotModel}
\title{
Plot a model
}
\description{
This function can be used to plot a prior model network before any
pre-processing step. However, it can also be completed with information
related to the pre-processing steps (compression. NONC nodes, expansion gates).
It can also be used to plot optimised model given the optimisation bitstring.

This function does not replace writeDot/writeNetwork used by the CNO pipeline.
}
\usage{
plotModel(model, optimRes=NA, indexIntegr=NA, signals=NULL, stimulis=NULL,
inhibitors=NULL, ncno=NULL, compressed=NULL, cnolist=NULL, PDF=FALSE, filename=NULL,SVG=FALSE, PNG=FALSE)

}
\arguments{
  \item{model}{
A model (either the filename or the output of readSif)
}
  \item{optimRes}{a bitstring made of 1 and 0 that must be of same length as the one returned by the GA algorithm (gaBinaryT1).}
  \item{indexIntegr}{a list of weight that must be of same length as the biString returned by the GA algorithm (gaBinaryT1).}
  \item{signals}{a list of nodes belonging to the signals class}
  \item{stimulis}{a list of nodes belonging to the stimulis class}
  \item{inhibitors}{a list of nodes belonging to the inhibitors class}
  \item{ncno}{a list of NCNO nodes.}
  \item{compressed}{a list of compressed nodes}
  \item{cnolist}{the output of makeCNOlist}
  \item{filename}{the filename (without extension) used to write the dot file}
  \item{PDF}{set to TRUE to save the plot in PDF format}
  \item{SVG}{set to TRUE to save the plot in SVG format}
  \item{PNG}{set to TRUE to save the plot in PNG format}
}
\details{
The function plot the model and also save it in a dot file that can be processed
later on. However, you can still save the plot in PNG or PDF or SVG file (one at
a time).

The cnolist variable contains the signals/stimulis/inhibitors so if you provide
a cnolist there is no need to use these arguments. If you decide to use them
they will overwrite the contents of cnolist.

optimRes and indeIntegr are used to weight/colorize the edges.

}
\value{
a graph representation of the model
\item{graph}{a vector of names of cues}
\item{graph}{a vector of names of cues}
\item{graph}{a vector of names of cues}
\item{graph}{a vector of names of cues}
}
%\references{
%J. Saez-Rodriguez, L. G. Alexopoulos, J. Epperlein, R. Samaga, D. A. Lauffenburger, S. Klamt and P. K. Sorger. Discrete logic modeling as a means to link protein signaling networks with functional analysis of
%mammalian signal transduction, Molecular Systems Biology, 5:331, 2009.
%}
\author{
T. Cokelaer
}
\seealso{
readMIDAS
readSif
makeCNOlist
writeNetwork
writeDot
}
\examples{
tmpdir<-tempdir()
setwd(tmpdir)
# Read the model
cpfile<-dir(system.file("ToyModel",package="CellNOptR"),full=TRUE)
file.copy(from=cpfile,to=getwd(),overwrite=TRUE)
model<-readSif(sifFile="ToyPKNMMB.sif")

# the PKN model can be plot easily
res<-plotModel(model)

# optional arguments can be provided to add information to the plot
cpfile<-dir(system.file("ToyModel",package="CellNOptR"),full=TRUE)
file.copy(from=cpfile,to=getwd(),overwrite=TRUE)
data<-readMIDAS(MIDASfile='ToyDataMMB.csv')
cnolist = makeCNOlist(data, subfield=FALSE)
res<-plotModel(model, cnolist=cnolist)
}

