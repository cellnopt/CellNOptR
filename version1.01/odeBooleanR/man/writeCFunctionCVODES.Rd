\name{writeCFunctionCVODES}
\alias{writeCFunctionCVODES}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{writeCFunctionCVODES}
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
writeCFunctionCVODES(inputs, outputs, strODE, numOutputs, parVecForPrint)
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
   headerDir=paste(.find.package("odeBooleanR"),"/textHeaders/startStringForFile.txt",sep="");
   startString=readChar(headerDir,nchar=10000);
   outputsStr=paste(outputs,sep=",",collapse=",");
   outputsStr=paste("double ",outputsStr,";",sep="");
   
   declareParVecForPrint=paste(parVecForPrint,sep="",collapse=",");
   declareParVecForPrint=paste("double ",parVecForPrint,";",sep="");
   
   outputAttr=c();

   for(i in 1:numOutputs)
   {
      outputAttr[i]=paste(outputs[i],"=Ith(y,",i,");",sep="");  
   }
   
   for(i in 1:numOutputs)
   {
      strODE[i]=paste("Ith(ydot,",i,")=",strODE[i],";",sep="");        
   }
   
   sizeParVec=length(parVecForPrint);
   count=0;
   
   for(i in 1:length(parVecForPrint))
   {
     parVecForPrint[i]=paste(parVecForPrint[i],"=ptr[",i-1,"];",sep=""); 
   }
   
   cFileDir=paste(.find.package("odeBooleanR"),"/models/myODE.c",sep="");  
   write(startString,file=cFileDir,sep="",append=FALSE)
   write(outputsStr,file=cFileDir,sep="",append=TRUE);
   write(declareParVecForPrint,file=cFileDir,sep="",append=TRUE);
   write(parVecForPrint,file=cFileDir,sep="",append=TRUE);
   write(outputAttr,file=cFileDir,sep="",append=TRUE);
   write(strODE, file=cFileDir, sep="", append = TRUE);  
   write("return(0);\n}",file=cFileDir, sep="", append = TRUE); 
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
