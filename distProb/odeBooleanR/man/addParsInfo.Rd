\name{addParsInfo}
\alias{addParsInfo}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
addParsInfo(odePars, initType)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{odePars}{
%%     ~~Describe \code{odePars} here~~
}
  \item{initType}{
%%     ~~Describe \code{initType} here~~
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
function(odePars,initType)
{
  LB=odePars$parValues;
  UB=odePars$parValues;
  LB[odePars$index_k]=0.001;
  LB[odePars$index_tau]=0.001;
  LB[odePars$index_n]=1;
  UB[odePars$index_k]=1;
  UB[odePars$index_tau]=10;
  UB[odePars$index_n]=5;
  
  index_opt=c(); 
  if(length(which(initType=='k'))>0)
  {
        index_opt=c(index_opt,odePars$index_k);
  }
  
  if(length(which(initType=='n'))>0)
  {
        index_opt=c(index_opt,odePars$index_n);
  }
  
  if(length(which(initType=='tau'))>0)
  {
        index_opt=c(index_opt,odePars$index_tau);
  }
  odePars$index_opt=index_opt;
  odePars$LB=LB;
  odePars$UB=UB;
  
  return(odePars);
  
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
