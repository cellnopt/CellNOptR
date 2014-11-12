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
makeBTables <-
function(CNOlist, k=2, measErr=c(0.1, 0), timePoint=NA){
  #CNOlist should be a list as output from makeCNOlist
  #k is the proportionality constant used to determine the threshold
  #measErr should be a 2 value vector (err1, err2) defining the error model as sd^2=err1^2+(err2*data)^2
  
  if ((class(CNOlist)=="CNOlist")==FALSE){
	CNOlist = CellNOptR::CNOlist(CNOlist)
  }	
	
  #if (!is.na(timePoint)){
  #  if (timePoint == "t1"){
	#  CNOlist@signals<-list(t0=CNOlist@signals[[1]], CNOlist@signals[[2]])
  #  }
  #  if (timePoint == "t2"){
	#  CNOlist@signals<-list(t0=CNOlist@signals[[1]], CNOlist@signals[[3]])
  #  }
  #}
  
  if (is.numeric(timePoint)){
    ix <- which(CNOlist@timepoints==timePoint)
    if (length(ix)==0){
      print("You do not have measures for the selected time point, please provide a different time point.")
    }else{
      CNOlist@signals<-list(t0=CNOlist@signals[[1]], CNOlist@signals[[ix]])
    }
  }
  
  err1 <- measErr[1]
  err2 <- measErr[2]
  
  mat <- list()
  NotMatStim <- list()
  NotMatInhib <- list()

  for (p in 1:length(colnames(CNOlist@signals[[2]]))){  
    #for each protein, create a m x n matrix where m is the number of inhibitors and n is the number of stimuli
    mat[[p]] <- matrix(data=0, nrow=length(colnames(CNOlist@inhibitors)), ncol=length(colnames(CNOlist@stimuli)))
    colnames(mat[[p]]) <- colnames(CNOlist@stimuli)
    rownames(mat[[p]]) <- colnames(CNOlist@inhibitors)
    
    NotMatStim[[p]]<-matrix(data=0, nrow=length(colnames(CNOlist@inhibitors)), ncol=length(colnames(CNOlist@stimuli)))
    colnames(NotMatStim[[p]]) <- colnames(CNOlist@stimuli)
    rownames(NotMatStim[[p]]) <- colnames(CNOlist@inhibitors)
    
    NotMatInhib[[p]]<-matrix(data=0, nrow=length(colnames(CNOlist@inhibitors)), ncol=length(colnames(CNOlist@stimuli)))
    colnames(NotMatInhib[[p]]) <- colnames(CNOlist@stimuli)
    rownames(NotMatInhib[[p]]) <- colnames(CNOlist@inhibitors)
    
    #consider as reference the condition with no stimuli and no inhibitors
    ref1_index <- which(apply(CNOlist@cues,1,sum) == 0)
        
    ref1 <- c()
    for (t in 1:length(CNOlist@signals)){
      ref1 <- c(ref1, CNOlist@signals[[t]][ref1_index, p])
      # if the baseline is NA I consider it = 0
	    ref1[is.na(ref1)] <- 0
    }
    
    # if the condition with no stimuli and no inhibitors is missing, I assume that the basal is = 0
    if (length(ref1_index)==0){
      ref1=rep(0,length(CNOlist@signals))
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
      
      # take care of cases in which single stimulus/single inhibitor conditions are missing or NA
      if (length(test1_index)==0 || length(na.omit(test1[-1]))==0){
        
        print(paste(colnames(CNOlist@signals[[2]])[p],": the condition with stimulus", colnames(CNOlist@stimuli)[i] , "is missing or NA"))
        
      }else{
      
        if (is.na(test1[1])){test1[1]<-0}
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
      
        if (length(which(abs(check) > k*check_sd))>0){
          checkStim<-check
          for (j in 1:length(colnames(CNOlist@inhibitors))){
			
			#I'm here only if the stimulus resuled to have effect on the protein, so I put a 1 in each position of the column referred to that stimulus
            mat[[p]][j,i] <- 1
          
            # if test1 < ref1 - k*SD (i.e. abs(check) > k*SD with check < 0) than the regulation has negative sign
            if (any(checkStim<0)){NotMatStim[[p]][j,i] <- 1}
          
            test2_index <- which(apply(CNOlist@stimuli,1,sum) == 1)  #select experiments with only one stimulus..
            test2_index <- intersect(test2_index, which(CNOlist@stimuli[,i] == 1)) #..that has to be the one we are interested in..
            test2_index <- intersect(test2_index, which(CNOlist@inhibitors[,j] == 1)) #..and with the specific inhibitor I want to test

          
            test2 <- c()
            for (t in 1:length(CNOlist@signals)){
              test2 <- c(test2, CNOlist@signals[[t]][test2_index, p])
            }
        
            # take care of cases in which single stimulus/single inhibitor conditions are missing or are all NA
            if (length(test2_index)==0 || length(na.omit(test2[-1]))==0){
            
              print(paste(colnames(CNOlist@signals[[2]])[p],": the condition with stimulus", colnames(CNOlist@stimuli)[i] ,"and inhibitor", colnames(CNOlist@inhibitors)[j]," is missing or NA"))
            
            }else{
            
              if (is.na(test2[1])){test2[1]<-0}
              test2_var <- err1^2+(err2*test2)^2
    
              test2 <- test2[2:length(test2)] - test2[1] #subtract the basal because it differs in different experimental conditions
              test2_var <- test2_var[2:length(test2_var)] + test2_var[1] #add the variance of the basale for error propagation
              test2_sd <- sqrt(test2_var) #compute the standard deviation from the variance of the error
      
              check <- ref2 - test2
              check_sd <- sqrt(ref2_sd^2 + test2_sd^2)
              if (length(which(abs(check) > k*check_sd))>0 && colnames(CNOlist@signals[[2]])[[p]] != colnames(CNOlist@inhibitors)[[j]]){
                mat[[p]][j,i] <- 2
                # if test2 > ref2 + k*SD (i.e. abs(check) > k*SD with check < 0) than the regulation has negative sign
                if (any(check<0)){NotMatInhib[[p]][j,i] <- 1}
              }
            }
          }
        }
      }
    }
  }
  names(mat)<-colnames(CNOlist@signals[[2]])
  names(NotMatStim)<-colnames(CNOlist@signals[[2]])
  names(NotMatInhib)<-colnames(CNOlist@signals[[2]])
  return(list(namesSignals=colnames(CNOlist@signals[[2]]), tables=mat, NotMatStim=NotMatStim, NotMatInhib=NotMatInhib))
}
