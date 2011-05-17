createOdefyModel <-
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

