transitionRules <-
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

