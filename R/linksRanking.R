#
#  This file is part of the CNO software
#
#  Copyright (c) 2011-2012 - EBI
#
#  File author(s): CNO developers (cno-dev@ebi.ac.uk)
#
#  Distributed under the GPLv3 License.
#  See accompanying file LICENSE.txt or copy at
#      http://www.gnu.org/licenses/gpl-3.0.html
#
#  CNO website: http://www.cellnopt.org
#
##############################################################################
# $Id$
linksRanking <-
function(CNOlist, measErr=c(0.1, 0), savefile=FALSE){
  
  if ((class(CNOlist)=="CNOlist")==FALSE){
  	CNOlist = CellNOptR::CNOlist(CNOlist)
  }	
	
  err1 <- measErr[1]
  err2 <- measErr[2]
  
  Nst <- length(colnames(CNOlist@stimuli))
  Nin <- length(colnames(CNOlist@inhibitors))
  Nsi <- length(colnames(CNOlist@signals[[2]]))
  
  Lrank <- matrix(data=0, nrow=Nsi*Nst+Nsi*Nin*Nst, ncol=2)
  index <- 1
  for (p in 1:length(colnames(CNOlist@signals[[2]]))){
    #for each measured protein first check the minimum value of k that allows the presence
    #of the direct link stim -> sign and than check the effact of the
    #inhibitor stim -> inh -> sign
    
    #consider as reference the condition with no stimuli and no inhibitors
    ref1_index <- which(apply(CNOlist@cues,1,sum) == 0)
    ref1 <- c()
    for (t in 1:length(CNOlist@signals)){
      ref1 <- c(ref1, CNOlist@signals[[t]][ref1_index, p])
    }
    ref1_var <- err1^2+(err2*ref1)^2
    
    ref1 <- ref1[2:length(ref1)] - ref1[1] #subtract the basal because it differs in different experimental conditions
    ref1_var <- ref1_var[2:length(ref1_var)] + ref1_var[1] #add the variance of the basale for error propagation
    ref1_sd <- sqrt(ref1_var) #compute the standard deviation from the variance of the error
    
    for (i in 1:length(colnames(CNOlist@stimuli))){
      #for each stimulus, verify its effect on the analyzed protein
      #considering as test the condition with that stimulus and no inhibitors
      test1_index <- which(apply(CNOlist@stimuli,1,sum) == 1)  #select experiments with only one stimulus..
      test1_index <- intersect(test1_index, which(CNOlist@stimuli[,i] == 1)) #..that has to be the one we are interested in..
      test1_index <- intersect(test1_index, which(apply(CNOlist@inhibitors,1,sum) == 0)) #..and with no inhibitors
      test1 <- c()
      for (t in 1:length(CNOlist@signals)){
        test1 <- c(test1, CNOlist@signals[[t]][test1_index, p])
      }
      test1_var <- err1^2+(err2*test1)^2
    
      test1 <- test1[2:length(test1)] - test1[1] #subtract the basal because it differs in different experimental conditions
      test1_var <- test1_var[2:length(test1_var)] + test1_var[1] #add the variance of the basale for error propagation
      test1_sd <- sqrt(test1_var) #compute the standard deviation from the variance of the error
      
      check <- test1 - ref1
      check_sd <- sqrt(ref1_sd^2 + test1_sd^2)
      k1 <- max(check/check_sd)
      
      Lrank[index,1] <- paste(colnames(CNOlist@stimuli)[i], '->', colnames(CNOlist@signals[[2]])[p])
      Lrank[index,2] <- k1
      index <- index+1
      
      ref2 <- test1
      ref2_var <- test1_var
      ref2_sd <- test1_sd
      
      for (j in 1:length(colnames(CNOlist@inhibitors))){
        test2_index <- which(apply(CNOlist@stimuli,1,sum) == 1)  #select experiments with only one stimulus..
        test2_index <- intersect(test2_index, which(CNOlist@stimuli[,i] == 1)) #..that has to be the one we are interested in..
        test2_index <- intersect(test2_index, which(CNOlist@inhibitors[,j] == 1)) #..and with the specific inhibitor I want to test

        test2 <- c()
        for (t in 1:length(CNOlist@signals)){
          test2 <- c(test2, CNOlist@signals[[t]][test2_index, p])
        }
        test2_var <- err1^2+(err2*test2)^2
    
        test2 <- test2[2:length(test2)] - test2[1] #subtract the basal because it differs in different experimental conditions
        test2_var <- test2_var[2:length(test2_var)] + test2_var[1] #add the variance of the basale for error propagation
        test2_sd <- sqrt(test2_var) #compute the standard deviation from the variance of the error
      
        check <- ref2 - test2
        check_sd <- sqrt(ref2_sd^2 + test2_sd^2)
        k2 <- max(check/check_sd)
        
        Lrank[index,1] <- paste(colnames(CNOlist@stimuli)[i], '->', colnames(CNOlist@inhibitors)[j], '->', colnames(CNOlist@signals[[2]])[p])
        Lrank[index,2] <- k2
        index <- index+1
      }
    }    
  }
  
  Lrank <- Lrank[order(as.numeric(Lrank[,2]), decreasing = TRUE),]
  Lrank <- Lrank[1:(min(which(Lrank[,2]<=0))-1),]
  
  if(savefile==TRUE){
	write.table(Lrank, file="LinkRanking.txt", sep="\t",row.names=FALSE,col.names=c('link','k min'),quote=FALSE)
  }
  return(Lrank)
  
}
