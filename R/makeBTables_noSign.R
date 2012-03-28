#
#  This file is part of the CNO software
#
#  Copyright (c) 2011-2012 - EBI
#
#  File author(s): CNO developers (cno-dev@ebi.ac.uk)
#
#  Distributed under the GPLv2 License.
#  See accompanying file LICENSE.txt or copy at
#      http://www.gnu.org/licenses/gpl-2.0.html
#
#  CNO website: http://www.ebi.ac.uk/saezrodriguez/software.html
#
##############################################################################
# $Id$
makeBTables <-
function(CNOlist, k, measErr){
  #CNOlist should be a list as output from makeCNOlist
  #k is the proportionality constant used to determine the threshold
  #measErr should be a 2 value vector (err1, err2) defining the error model as sd^2=err1^2+(err2*data)^2
  err1 <- measErr[1]
  err2 <- measErr[2]
  
  mat <- list()
  for (p in 1:length(CNOlist$namesSignals)){  
    #for each protein, create a m x n matrix where m is the number of inhibitors and n is the number of stimuli
    mat[[p]] <- matrix(data=0, nrow=length(CNOlist$namesInhibitors), ncol=length(CNOlist$namesStimuli))
    colnames(mat[[p]]) <- CNOlist$namesStimuli
    rownames(mat[[p]]) <- CNOlist$namesInhibitors
    #consider as reference the condition with no stimuli and no inhibitors
    ref1_index <- which(apply(CNOlist$valueCues,1,sum) == 0)
    
    # if the condition with no stimuli and no inhibitors is missing, I assume that the basal is = 0
    if (lenght(ref1_index)==0){ref1=rep(0,length(CNOlist$valueSignals))}
    
    ref1 <- c()
    for (t in 1:length(CNOlist$valueSignals)){
      ref1 <- c(ref1, CNOlist$valueSignals[[t]][ref1_index, p])
	  ref1[is.na(ref1)] <- 0
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
      
      # take care of cases in which single stimulus/single inhibitor conditions are missing
      if (length(test1_index)==0){stop("single-stimulus/single-inhibitor conditions are missing")}
      
      test1 <- c()
      for (t in 1:length(CNOlist$valueSignals)){
        test1 <- c(test1, CNOlist$valueSignals[[t]][test1_index, p])
		    test1[is.na(test1)] <- 0
      }
      test1_var <- err1^2+(err2*test1)^2
    
      test1 <- test1[2:length(test1)] - test1[1] #subtract the basal because it differs in different experimental conditions
      test1_var <- test1_var[2:length(test1_var)] + test1_var[1] #add the variance of the basale for error propagation
      test1_sd <- sqrt(test1_var) #compute the standard deviation from the variance of the error
      
      check <- test1 - ref1
      check_sd <- sqrt(ref1_sd^2 + test1_sd^2)
      
      #the condition with the stimulus and no inhibitor in the new reference to test the effect of the inhibitors
      ref2 <- test1
      ref2_var <- test1_var
      ref2_sd <- test1_sd
      
      if (length(which(check > k*check_sd))>0){
        for (j in 1:length(CNOlist$namesInhibitors)){
          #I'm here only if the stimulus resuled to have effect on the protein, so I put a 1 in each position of the column referred to that stimulus
          mat[[p]][j,i] <- 1
          
          
          test2_index <- which(apply(CNOlist$valueStimuli,1,sum) == 1)  #select experiments with only one stimulus..
          test2_index <- intersect(test2_index, which(CNOlist$valueStimuli[,i] == 1)) #..that has to be the one we are interested in..
          test2_index <- intersect(test2_index, which(CNOlist$valueInhibitors[,j] == 1)) #..and with the specific inhibitor I want to test

          # take care of cases in which single stimulus/single inhibitor conditions are missing
          if (length(test2_index)==0){stop("single-stimulus/single-inhibitor conditions are missing")}
          
          test2 <- c()
          for (t in 1:length(CNOlist$valueSignals)){
            test2 <- c(test2, CNOlist$valueSignals[[t]][test2_index, p])
			      test2[is.na(test2)] <- 0
          }
          test2_var <- err1^2+(err2*test2)^2
    
          test2 <- test2[2:length(test2)] - test2[1] #subtract the basal because it differs in different experimental conditions
          test2_var <- test2_var[2:length(test2_var)] + test2_var[1] #add the variance of the basale for error propagation
          test2_sd <- sqrt(test2_var) #compute the standard deviation from the variance of the error
      
          check <- ref2 - test2
          check_sd <- sqrt(ref2_sd^2 + test2_sd^2)
          if (length(which(check > k*check_sd))>0 && CNOlist$namesSignals[[p]] != CNOlist$namesInhibitors[[j]]){
            mat[[p]][j,i] <- 2
          }     
        }
      }  
    }
  }
  names(mat)<-CNOlist$namesSignals
  return(list(namesSignals=CNOlist$namesSignals, tables=mat))
}
