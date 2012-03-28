AddLink <-
function(CueDown, SigUp, Model, Sign){
  
  # Update the vector of the name of reactions

  reaction<-paste(CueDown,"=",SigUp,sep="")
  if (Sign==(-1)){reaction<-paste("!",reaction,sep="")}
    
  if (!any(Model$reacID==reaction)){
    Model$reacID[length(Model$reacID)+1] <- reaction
    
	  tmp<-matrix(data=0,nrow=dim(Model$interMat)[1], ncol=1)
    sou<-match(CueDown,Model$namesSpecies)
	  if (Sign==(-1)){tmp[sou,1]<-1}
	  Model$notMat <- cbind(Model$notMat,tmp)
    
	  tmp<-matrix(data=0,nrow=dim(Model$interMat)[1], ncol=1)
	  sou<-match(CueDown,Model$namesSpecies)
    tar<-match(SigUp,Model$namesSpecies)
    tmp[sou,1]<-(-1)
	  tmp[tar,1]<-1
    Model$interMat <- cbind(Model$interMat,tmp)
  }
      
  colnames(Model$interMat) <- Model$reacID
  colnames(Model$notMat) <- Model$reacID
      
  return(Model)
  
}