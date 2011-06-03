\name{odeSensitivityAnalisys}
\alias{odeSensitivityAnalisys}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{odeSensitivityAnalisys
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
odeSensitivityAnalisys(CNOlist, CNOModel, odePars, times)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{CNOlist}{
%%     ~~Describe \code{CNOlist} here~~
}
  \item{CNOModel}{
%%     ~~Describe \code{CNOModel} here~~
}
  \item{odePars}{
%%     ~~Describe \code{odePars} here~~
}
  \item{times}{
%%     ~~Describe \code{times} here~~
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
function(CNOlist,CNOModel,odePars,times)
{ 
  compileResult=compileAndLoad("myODE");
  
  indexes=indexCNO2ODE(CNOlist,odePars);
  
  TTModel=createOdefyModel(CNOModel$interMat,CNOModel$notMat,CNOModel$namesSpecies);
  
  expressionsInfo=getODEModelExpressions(TTModel)
  
  pderivs=getPartialDerivs(expressionsInfo)
  
  res=computeSensitivities(expressionsInfo,pderivs,CNOlist,indexes,times)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
