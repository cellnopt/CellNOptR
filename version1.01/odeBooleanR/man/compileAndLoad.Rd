\name{compileAndLoad}
\alias{compileAndLoad}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{compileAndLoad
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
compileAndLoad(cFileName)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{cFileName}{
%%     ~~Describe \code{cFileName} here~~
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
\references{
%% ~put references to the literature/web site here ~
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
function(cFileName)
{  
  actualDir=getwd();
  cFileDir=paste(.find.package("odeBooleanR"),"/models/",sep="")
  setwd(cFileDir);
  try(dyn.unload(paste(cFileName, .Platform$dynlib.ext, sep = "")));
  compileResult=system(paste('R CMD SHLIB ',cFileName,".c",sep=""));
  if(as.logical(compileResult))print("Compilation Failed");
  dyn.load(paste(cFileDir,cFileName, .Platform$dynlib.ext, sep = ""));   
  if(cFileName=="mymod") {library(deSolve);}
  else if(cFileName=="myODE"){library(Rsundials);}
  else {print("A valid File Name was not provided");}
  setwd(actualDir);
  return(compileResult);
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
