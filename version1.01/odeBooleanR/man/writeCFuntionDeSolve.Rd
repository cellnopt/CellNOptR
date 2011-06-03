\name{writeCFuntionDeSolve}
\alias{writeCFuntionDeSolve}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{writeCFuntionDeSolve
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
writeCFuntionDeSolve(inputs, outputs, strODE, numOutputs, parVecForPrint)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{inputs}{
%%     ~~Describe \code{inputs} here~~
}
  \item{outputs}{
%%     ~~Describe \code{outputs} here~~
}
  \item{strODE}{
%%     ~~Describe \code{strODE} here~~
}
  \item{numOutputs}{
%%     ~~Describe \code{numOutputs} here~~
}
  \item{parVecForPrint}{
%%     ~~Describe \code{parVecForPrint} here~~
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
function(inputs,outputs,strODE,numOutputs,parVecForPrint)
{
   outputsStr=paste(outputs,sep=",",collapse=",");
   outputsStr=paste("double ",outputsStr,";",sep="");
   
   outputAttr=c();
    
   for(i in 1:numOutputs)
   {
      outputAttr[i]=paste("double ",outputs[i],"=y[",i-1,"];",sep="");
   }
   
   for(i in 1:numOutputs)
   {
      strODE[i]=paste("ydot[",i-1,"]=",strODE[i],";",sep="");
   }
       
   parVecForPrintInclude=c();    
   
   for(i in 1:length(parVecForPrint))
   {
     parVecForPrintInclude[i]=paste("#define ",parVecForPrint[i]," parms[",i-1,"]",sep="");
   }
   
   write("#include <R.h>",file="mymod.c",sep="",append=FALSE);
   write("#include <math.h>",file="mymod.c",sep="",append=TRUE);
   tempStr=paste("static double parms[",length(parVecForPrintInclude),"];",sep="");
   write(tempStr,file="mymod.c",sep="",append=TRUE)
   write(parVecForPrintInclude,file="mymod.c",sep="",append=TRUE);
   tempStr=paste("void initmod(void (* odeparms)(int *, double *))\n{\n int N=",length(parVecForPrint),";",sep="")
   write(tempStr,file="mymod.c",sep="",append=TRUE);
   tempStr="odeparms(&N, parms);\n}\n/* Derivatives and 1 output variable */\n void derivs (int *neq, double *t, double *y, double *ydot,\n double *yout, int *ip)\n{";
   write(tempStr,file="mymod.c",sep="",append=TRUE);
   write(outputAttr,file="mymod.c",sep="",append=TRUE);
   write(strODE,file="mymod.c",sep="",append=TRUE);
   write("}",file="mymod.c",sep="",append=TRUE);
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
