\name{createInSilicoData}
\alias{createInSilicoData}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{createInSilicoData
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
createInSilicoData(paramInfo, y0, times, cFileName)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{paramInfo}{
%%     ~~Describe \code{paramInfo} here~~
}
  \item{y0}{
%%     ~~Describe \code{y0} here~~
}
  \item{times}{
%%     ~~Describe \code{times} here~~
}
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
function(paramInfo,y0,times,cFileName)
{
   if(cFileName=="mymod") 
   {
      out<-lsode(y0, times, func = "derivs", parms=paramInfo$paramValues, jacfunc=NULL, dllname="mymod",initfunc ="initmod",nout = 0);
   }
   else if(cFileName=="myODE")
   {
      out<-cvodes(y0,times,"myODE","rhs",fndata=paramInfo$paramValues,verbose =FALSE);
   }
   else
   {       
      out=NA;
      print("A valid File Name was not provided");
   }
   return(out);
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
