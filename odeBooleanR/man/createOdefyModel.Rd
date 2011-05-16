\name{createOdefyModel}
\alias{createOdefyModel}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{createOdefyModel
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
createOdefyModel(intermat, notmat, specID)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{intermat}{
%%     ~~Describe \code{intermat} here~~
}
  \item{notmat}{
%%     ~~Describe \code{notmat} here~~
}
  \item{specID}{
%%     ~~Describe \code{specID} here~~
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
function(intermat,notmat,specID)
{
  res=list();
  rules=transitionRules(intermat,notmat);
  count=0;
  for(i in 1:length(specID))
  { 
    tempRules=rules[[i]]$mat; 
     
    numRules=length(tempRules[1,])
    indexInputs=as.matrix(which(tempRules!=0,arr.ind=TRUE));
    indexInputs=unique(indexInputs[,1]);
    numInputs=length(indexInputs);
    nameInputs=specID[indexInputs];
    if(numRules>0)
    { 
      inputTable=createInputTable(numInputs);
      count=count+1
      res[count]=list(inputs=NULL,output=NULL,truthTable=NULL)
      truthTable=logical();
      truthTable[seq(1,2^numInputs)]=FALSE; 
      for(j in 1:numRules)
      {  
        rule=tempRules[,j];
        truth=which(rule!=0)
        rule[which(rule==-1)]=0;
        rule=as.logical(rule);
        stateVecIndex=c();
        
        for(k in 1:length(truth))stateVecIndex[k]=which(indexInputs==truth[k]);
        for(k in 1:2^numInputs)
        {
          if(!any(which((inputTable[k,stateVecIndex]==rule[truth])==FALSE)))
          {
            truthTable[k]=TRUE;
          }
        } 
      }
      res[[count]]$inputs=nameInputs;
      res[[count]]$output=specID[i];
      res[[count]]$truthTable=truthTable;
    }
    }
    return(res)  
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
