linksRanking <-
function(CNOlist, measErr){
  
  err1 <- measErr[1]
  err2 <- measErr[2]
  
  Nst <- length(CNOlist$namesStimuli)
  Nin <- length(CNOlist$namesInhibitors)
  Nsi <- length(CNOlist$namesSignals)
  
  Lrank <- matrix(data=0, nrow=Nsi*Nst+Nsi*Nin*Nst, ncol=2)
  index <- 1
  for (p in 1:length(CNOlist$namesSignals)){
    #for each measured protein first check the minimum value of k that allows the presence
    #of the direct link stim -> sign and than check the effact of the
    #inhibitor stim -> inh -> sign
    
    #consider as reference the condition with no stimuli and no inhibitors
    ref1_index <- which(apply(CNOlist$valueCues,1,sum) == 0)
    ref1 <- c()
    for (t in 1:length(CNOlist$valueSignals)){
      ref1 <- c(ref1, CNOlist$valueSignals[[t]][ref1_index, p])
    }
    ref1_var <- err1^2+(err2*ref1)^2
    
    ref1 <- ref1[2:length(ref1)] - ref1[1] #subtract the basal because it differs in different experimental conditions
    ref1_var <- ref1_var[2:length(ref1_var)] + ref1_var[1] #add the variance of the basale for error propagation
    ref1_sd <- sqrt(ref1_var) #compute the standard deviation from the variance of the error
    
    for (i in 1:length(CNOlist$namesStimuli)){
      #for each stimulus, verify its effect on the analyzed protein
      #considering as test the condition with that stimulus and no inhibitors
      test1_index <- which(apply(CNOlist$valueStimuli,1,sum) == 1)  #select experiments with only one stimulus..
      test1_index <- intersect(test1_index, which(CNOlist$valueStimuli[,i] == 1)) #..that has to be the one we are interested in..
      test1_index <- intersect(test1_index, which(apply(CNOlist$valueInhibitors,1,sum) == 0)) #..and with no inhibitors
      test1 <- c()
      for (t in 1:length(CNOlist$valueSignals)){
        test1 <- c(test1, CNOlist$valueSignals[[t]][test1_index, p])
      }
      test1_var <- err1^2+(err2*test1)^2
    
      test1 <- test1[2:length(test1)] - test1[1] #subtract the basal because it differs in different experimental conditions
      test1_var <- test1_var[2:length(test1_var)] + test1_var[1] #add the variance of the basale for error propagation
      test1_sd <- sqrt(test1_var) #compute the standard deviation from the variance of the error
      
      check <- test1 - ref1
      check_sd <- sqrt(ref1_sd^2 + test1_sd^2)
      k1 <- max(check/check_sd)
      
      Lrank[index,1] <- paste(CNOlist$namesStimuli[i], '->', CNOlist$namesSignals[p])
      Lrank[index,2] <- k1
      index <- index+1
      
      ref2 <- test1
      ref2_var <- test1_var
      ref2_sd <- test1_sd
      
      for (j in 1:length(CNOlist$namesInhibitors)){
        test2_index <- which(apply(CNOlist$valueStimuli,1,sum) == 1)  #select experiments with only one stimulus..
        test2_index <- intersect(test2_index, which(CNOlist$valueStimuli[,i] == 1)) #..that has to be the one we are interested in..
        test2_index <- intersect(test2_index, which(CNOlist$valueInhibitors[,j] == 1)) #..and with the specific inhibitor I want to test

        test2 <- c()
        for (t in 1:length(CNOlist$valueSignals)){
          test2 <- c(test2, CNOlist$valueSignals[[t]][test2_index, p])
        }
        test2_var <- err1^2+(err2*test2)^2
    
        test2 <- test2[2:length(test2)] - test2[1] #subtract the basal because it differs in different experimental conditions
        test2_var <- test2_var[2:length(test2_var)] + test2_var[1] #add the variance of the basale for error propagation
        test2_sd <- sqrt(test2_var) #compute the standard deviation from the variance of the error
      
        check <- ref2 - test2
        check_sd <- sqrt(ref2_sd^2 + test2_sd^2)
        k2 <- max(check/check_sd)
        
        Lrank[index,1] <- paste(CNOlist$namesStimuli[i], '->', CNOlist$namesInhibitors[j], '->', CNOlist$namesSignals[p])
        Lrank[index,2] <- k2
        index <- index+1
      }
    }    
  }
  
  Lrank <- Lrank[order(as.numeric(Lrank[,2]), decreasing = TRUE),]
  Lrank <- Lrank[1:(min(which(Lrank[,2]<=0))-1),]
  
  write.table(Lrank, file="LinkRanking.txt", sep="\t",row.names=FALSE,col.names=c('link','k min'),quote=FALSE)

  return(Lrank)
  
}