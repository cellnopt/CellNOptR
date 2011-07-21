\name{transitionRules}
\alias{transitionRules}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
transitionRules(intermat, notmat)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{intermat}{
%%     ~~Describe \code{intermat} here~~
}
  \item{notmat}{
%%     ~~Describe \code{notmat} here~~
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
function(intermat,notmat)
{ 
  numSpecies=length(intermat[,1]); 
  result=list();
  for(i in 1:numSpecies)result[i]=list(mat=NULL);
  for(i in 1:numSpecies)
  {
    indexAND=which(intermat[i,]==1);
    vec=c();
    vec=intermat[,indexAND];
    vec[which(vec==1)]=0;
    vec[which(vec==-1)]=1;
    vec[which(notmat[,indexAND]==1)]=-1;
    result[[i]]$mat=as.matrix(vec); 
  }
  return(result);
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
