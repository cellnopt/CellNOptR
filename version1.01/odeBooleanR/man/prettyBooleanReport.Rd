\name{prettyBooleanReport}
\alias{prettyBooleanReport}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{prettyBooleanReport
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
prettyBooleanReport(boolModel)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{boolModel}{
%%     ~~Describe \code{boolModel} here~~
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
function(boolModel)
{ 
  library(xtable);
  numOutputs=length(boolModel) 
  startString=readChar("latexHeader.txt",nchar=10000);
  repFile="boolReport/report.tex";
  write(startString,file=repFile,sep="",append=FALSE)
  for(i in 1:numOutputs)
  {
    numInputs=length(boolModel[[i]]$inputs);
    ttable=createInputTable(numInputs);
    ttable=cbind(ttable,boolModel[[i]]$truthTable)
    namesTruthTable=c(boolModel[[i]]$inputs,paste(boolModel[[i]]$output,"output",sep="-"));
    colnames(ttable)=namesTruthTable;
    ttable[which(ttable)]=as.integer(1);
    ttable[which(!ttable)]=as.integer(0);
    ttable=xtable(ttable);
    print(ttable, type="latex", file=repFile, append=TRUE,tabular.environment="longtable",floating=FALSE);
    write("\\clearpage",file=repFile,sep="\n",append=TRUE)
  }
  write("\\end{document}",file=repFile,sep="\n",append=TRUE)   
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
