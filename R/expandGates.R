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
expandGates<-function(Model, ignoreList=NA){

	# check that Model is a Model list
	if(!is.list(Model)) stop("This function expects as input a Model as output by readSif")
	if(length(Model) == 4) {
		if(all(names(Model) != c("reacID", "namesSpecies","interMat","notMat"))){
			stop("This function expects as input a Model as output by readSif")
			}	
		}
	if(length(Model) == 5) {
		if(all(names(Model) != c("reacID", "namesSpecies","interMat","notMat","speciesCompressed"))) {
			stop("This function expects as input a Model as output by readSif")
		}	
	}	
	
	SplitANDs <- list(initialReac=c("split1","split2"))
	splitR <- 1
	
	# split all the ANDs
	# remove any ANDs 2/3 and save >3 to add later
	# +3 won't get added here again but are part of prior knowledge if contained within PKN
	andToAdd = c()
	remove.and = c()
	reacs2Ignore = c()
	initialReacN <- length(Model$reacID)

  #TODO: move into readSIF ?
  if (initialReacN == 1){
        Model$interMat <- as.matrix(Model$interMat)
   } 

	# which reactions have ignoreList as output?
	if(!is.na(ignoreList[1])) {
		for(s in 1:initialReacN) {
			if(any(Model$interMat[ignoreList,s] == 1)) {
				reacs2Ignore = c(reacs2Ignore, s)
			}
		}
	}
	
	for(r in 1:initialReacN) {
		inNodes <- which(Model$interMat[,r] == -1)
		if(length(inNodes) > 1) {
			if(length(inNodes) > 3) {
				andToAdd = c(andToAdd, r)	
			}
			remove.and = c(remove.and, r)
			
			if(!any(reacs2Ignore == r)) {
				outNode <- which(Model$interMat[,r] == 1)
				newReacs <- matrix(data=0,nrow=dim(Model$interMat)[1],ncol=length(inNodes))
				newReacsNots <- matrix(data=0,nrow=dim(Model$interMat)[1],ncol=length(inNodes))
				newReacs[outNode,] <- 1
				newReacsIDs <- rep("a",length(inNodes))
				for(i in 1:length(inNodes)) {
					newReacs[inNodes[i],i] <- -1
					newReacsNots[inNodes[i],i] <- Model$notMat[inNodes[i],r]
					newReacsIDs[i] <- paste(Model$namesSpecies[inNodes[i]],"=",Model$namesSpecies[outNode],sep="")
					if(Model$notMat[inNodes[i],r] == 1) newReacsIDs[i] <- paste("!",newReacsIDs[i],sep="")
				}
				colnames(newReacs) <- newReacsIDs
				colnames(newReacsNots) <- newReacsIDs
				SplitANDs[[splitR]] <- newReacsIDs
				names(SplitANDs)[splitR] <- Model$reacID[r]
				splitR <- splitR+1
				Model$notMat <- cbind(Model$notMat,newReacsNots)
				Model$interMat <- cbind(Model$interMat,newReacs)
				Model$reacID <- c(Model$reacID,newReacsIDs)
			}
		}
	}
	
	# ***** ADDED *****	
	# remove 'AND' gates that will be made anyway
	# save anything > 3 to add at end
	
	if(length(andToAdd)) {
		toAdd = list()
		toAdd$notMat <- Model$notMat[,andToAdd]
		toAdd$interMat <- Model$interMat[,andToAdd]
		toAdd$reacID <- Model$reacID[andToAdd]
	} else {
		toAdd <- NA	
	}
	
	if(length(remove.and)) {
		Model$notMat <- Model$notMat[,-remove.and]
		Model$interMat <- Model$interMat[,-remove.and]
		Model$reacID <- Model$reacID[-remove.and]
	}
	
	# ***** ADDED *****

	# the list SplitANDs now contains a named element for each AND reac that has been split, 
	# and each element contains a vector with the names of the reactions that result from the split
	# create combinations of ORs
	# the newANDs list will contain an element for each new '&' gate, named by the name of this new and
	# reac, and containing a vector of the names of the reactions from which it was created
	checkExisting<-function(x) {
		all(x == newReac)
	}
	
	newANDs <- list(finalReac=c("or1","or2"))
	ANDsadded <- 1
	total.list = 1:length(Model$namesSpecies)	
	
	for(sp in total.list) {
	#	print(sp)
		inReacsIndex <- which(Model$interMat[sp,] == 1)
		if(length(inReacsIndex) > 1) {
			inReacs <- Model$interMat[,inReacsIndex]
			findInput <- function(x) {
				inNode<-which(x == -1)
			}
			inSp <- apply(inReacs,2,findInput)
			
			# in case there are only 2 inputs, I just need to create a single '&'
			if(length(inSp) == 2) {
				newReac <- rep(0, length(Model$namesSpecies))
				newReac[unlist(inSp)] <- -1
				newReac[sp] <- 1
	
				#check that the reaction doesn't exist
				exist <- apply(Model$interMat,2,checkExisting)
				# if the reaction doesn't exist:
				if(!any(exist)) {
					# 1: add it to interMat
					Model$interMat <- cbind(Model$interMat,newReac)
					# 2: add it to notMat
					newNot <- (Model$notMat[,inReacsIndex[1]] + Model$notMat[,inReacsIndex[2]])
					Model$notMat <- cbind(Model$notMat,newNot)
					# 3: create a new name for it
					newName1 <- ifelse((newNot[inSp[1]]==1),paste("!",Model$namesSpecies[inSp[1]],sep=""),Model$namesSpecies[inSp[1]])
					newName2 <- ifelse((newNot[inSp[2]]==1),paste("!",Model$namesSpecies[inSp[2]],sep=""),Model$namesSpecies[inSp[2]])
					newName <- paste(newName1,"+",newName2,"=",Model$namesSpecies[sp],sep="")
					Model$reacID <- c(Model$reacID,newName)
					colnames(Model$interMat)[dim(Model$interMat)[2]] <- newName
					colnames(Model$notMat)[dim(Model$notMat)[2]] <- newName
					# 4: add this new reac in the newANDs list
					newANDs[[ANDsadded]] <- Model$reacID[inReacsIndex]
					names(newANDs)[[ANDsadded]] <- newName
					ANDsadded <- ANDsadded+1
				}
			} else {
				
				# in case there are 3 inputs, create 3 AND gates with 2 inputs and one AND gate with 3 inputs
				if(length(inSp) == 3) {	
					newReac1 <- rep(0, length(Model$namesSpecies))
					newReac2 <- rep(0, length(Model$namesSpecies))
					newReac3 <- rep(0, length(Model$namesSpecies))
					newReac4 <- rep(0, length(Model$namesSpecies))
					newReac1[sp] <- 1
					newReac2[sp] <- 1
					newReac3[sp] <- 1
					newReac4[sp] <- 1
					newReac1[inSp[1:2]] <- -1
					newReac2[inSp[2:3]] <- -1
					newReac3[inSp[c(1,3)]] <- -1
					newReac4[inSp] <- -1
					newReacs <- cbind(newReac1,newReac2,newReac3,newReac4)
					newNot1 <- (Model$notMat[,inReacsIndex[1]]+Model$notMat[,inReacsIndex[2]])
					newNot2 <- (Model$notMat[,inReacsIndex[2]]+Model$notMat[,inReacsIndex[3]])
					newNot3 <- (Model$notMat[,inReacsIndex[1]]+Model$notMat[,inReacsIndex[3]])
					newNot4 <- (Model$notMat[,inReacsIndex[1]]+Model$notMat[,inReacsIndex[2]]+Model$notMat[,inReacsIndex[3]])
					newNots <- cbind(newNot1,newNot2,newNot3,newNot4)
					
					# go through each one of these reacs and check if they exist
						rexist <- rep(FALSE,4)
						for(i in 1:4) {
							newReac <- newReacs[,i]
							exist <- apply(Model$interMat,2,checkExisting)
							rexist[i] <- any(exist)
						}
					
					# if at least one of the reactions doesn't exist yet
					if(any(rexist == FALSE)) {
						newReacs <- newReacs[,which(rexist == FALSE)]
						# 1:add the new reacs to the interMat
						Model$interMat <- cbind(Model$interMat,newReacs)
						# 2: add to notMat
						newNots<-newNots[,which(rexist == FALSE)]
						Model$notMat<-cbind(Model$notMat,newNots)
						# 3: create the new names
						newNames <- rep("a",4)
						newName1 <- ifelse((newNots[inSp[1],1]==1),paste("!",Model$namesSpecies[inSp[1]],sep=""),Model$namesSpecies[inSp[1]])
						newName2 <- ifelse((newNots[inSp[2],2]==1),paste("!",Model$namesSpecies[inSp[2]],sep=""),Model$namesSpecies[inSp[2]])
						newNames[1] <- paste(newName1,"+",newName2,"=",Model$namesSpecies[sp],sep="")
						newName1 <- ifelse((newNots[inSp[2],2]==1),paste("!",Model$namesSpecies[inSp[2]],sep=""),Model$namesSpecies[inSp[2]])
						newName2 <- ifelse((newNots[inSp[3],3]==1),paste("!",Model$namesSpecies[inSp[3]],sep=""),Model$namesSpecies[inSp[3]])
						newNames[2] <- paste(newName1,"+",newName2,"=",Model$namesSpecies[sp],sep="")
						newName1 <- ifelse((newNots[inSp[1],1]==1),paste("!",Model$namesSpecies[inSp[1]],sep=""),Model$namesSpecies[inSp[1]])
						newName2 <- ifelse((newNots[inSp[3],3]==1),paste("!",Model$namesSpecies[inSp[3]],sep=""),Model$namesSpecies[inSp[3]])
						newNames[3] <- paste(newName1,"+",newName2,"=",Model$namesSpecies[sp],sep="")
						newName1 <- ifelse((newNots[inSp[1],1]==1),paste("!",Model$namesSpecies[inSp[1]],sep=""),Model$namesSpecies[inSp[1]])
						newName2 <- ifelse((newNots[inSp[2],2]==1),paste("!",Model$namesSpecies[inSp[2]],sep=""),Model$namesSpecies[inSp[2]])
						newName3 <- ifelse((newNots[inSp[3],3]==1),paste("!",Model$namesSpecies[inSp[3]],sep=""),Model$namesSpecies[inSp[3]])
						newNames[4]<-paste(newName1,"+",newName2,"+",newName3,"=",Model$namesSpecies[sp],sep="")
						newNames <- newNames[which(rexist == FALSE)]
						colnames(Model$notMat)[grep("newNot",colnames(Model$notMat))] <- newNames
						colnames(Model$interMat)[grep("newReac",colnames(Model$interMat))] <- newNames
						Model$reacID <- c(Model$reacID,newNames)
						
						# 4: add these new reacs in the newANDs list
						newList <- list(c(Model$reacID[inReacsIndex[1]],Model$reacID[inReacsIndex[2]]))
						newList[[2]] <- c(Model$reacID[inReacsIndex[2]],Model$reacID[inReacsIndex[3]])
						newList[[3]] <- c(Model$reacID[inReacsIndex[1]],Model$reacID[inReacsIndex[3]])
						newList[[4]] <- c(Model$reacID[inReacsIndex[1]],Model$reacID[inReacsIndex[2]],Model$reacID[inReacsIndex[3]])
						names(newList) <- newNames
						newList<-newList[!rexist]
						newANDs[ANDsadded:(ANDsadded+sum(!rexist)-1)] <- newList
						names(newANDs)[ANDsadded:(ANDsadded+sum(!rexist)-1)] <- names(newList)
						ANDsadded <- ANDsadded+sum(!rexist)
					}
				} else {
					
					# in case there are 4 inputs or more
					if(length(inSp) >= 4){
		
						# tag the species that enter with a neg (all these reacs should be OR, i.e. they should have 1 input)
						nots <- colSums(Model$notMat[,inReacsIndex])
						nots <- (nots > 0)
						# generate the combinations of 2
						combis <- combn(x=seq(1,length(inSp)),m=2)
						newnot <- matrix(0,nrow=length(Model$namesSpecies),ncol=dim(combis)[2])
						newReacs <- matrix(0,nrow=length(Model$namesSpecies),ncol=dim(combis)[2])
						newReacs[sp,] <- 1
						newNames <- rep("a",dim(combis)[2])
						
						for(i in 1:dim(combis)[2]) {
							newReacs[inSp[combis[,i]],i] <- -1
							newnot[inSp[combis[1,i]],1] <- ifelse(nots[combis[1,i]],1,0)
							newnot[inSp[combis[2,i]],1] <- ifelse(nots[combis[2,i]],1,0)
							newname1 <- ifelse(nots[combis[1,i]],paste("!",Model$namesSpecies[inSp[combis[1,i]]],sep=""),Model$namesSpecies[inSp[combis[1,i]]])
							newname2 <- ifelse(nots[combis[2,i]],paste("!",Model$namesSpecies[inSp[combis[2,i]]],sep=""),Model$namesSpecies[inSp[combis[2,i]]])
							newNames[i] <- paste(newname1,"+",newname2,"=",Model$namesSpecies[sp],sep="")
						}
						
						# generate the combinations of 3
						combis2 <- combn(x=seq(1,length(inSp)),m=3)
						newnot2 <- matrix(0,nrow=length(Model$namesSpecies),ncol=dim(combis2)[2])
						newReacs2 <- matrix(0,nrow=length(Model$namesSpecies),ncol=dim(combis2)[2])
						newReacs2[sp,] <- 1
						newNames2 <- rep("a",dim(combis2)[2])
						for(i in 1:dim(combis2)[2]) {
							newReacs2[inSp[combis2[,i]],i] <- -1
							newnot2[inSp[combis2[1,i]],1] <- ifelse(nots[combis2[1,i]],1,0)
							newnot2[inSp[combis2[2,i]],1] <- ifelse(nots[combis2[2,i]],1,0)
							newnot2[inSp[combis2[3,i]],1] <- ifelse(nots[combis2[3,i]],1,0)
							newname1 <- ifelse(nots[combis2[1,i]],paste("!",Model$namesSpecies[inSp[combis2[1,i]]],sep=""),Model$namesSpecies[inSp[combis2[1,i]]])
							newname2 <- ifelse(nots[combis2[2,i]],paste("!",Model$namesSpecies[inSp[combis2[2,i]]],sep=""),Model$namesSpecies[inSp[combis2[2,i]]])
							newname3 <- ifelse(nots[combis2[3,i]],paste("!",Model$namesSpecies[inSp[combis2[3,i]]],sep=""),Model$namesSpecies[inSp[combis2[3,i]]])
							newNames2[i] <- paste(newname1,"+",newname2,"+",newname3,"=",Model$namesSpecies[sp],sep="")
						}
						
						newnot <- cbind(newnot,newnot2)
						newReacs <- cbind(newReacs,newReacs2)
						newNames <- c(newNames,newNames2)
						colnames(newnot) <- newNames
						colnames(newReacs) <- newNames

						# check which reactions exist
						rexist <- rep(FALSE,length(newNames))
						for(i in 1:length(rexist)) {
							newReac <- newReacs[,i]
							exist <- apply(Model$interMat,2,checkExisting)
							rexist[i] <- any(exist)
						}
						
						if(any(rexist == FALSE)) {
							Model$interMat <- cbind(Model$interMat,newReacs[,which(rexist == FALSE)])
							Model$notMat <- cbind(Model$notMat,newnot[,which(rexist == FALSE)])
							Model$reacID <- c(Model$reacID,newNames[which(rexist == FALSE)])
							newList <- list(c(Model$reacID[inReacsIndex[combis[1,1]]],Model$reacID[inReacsIndex[combis[2,1]]]))
							
							for(i in 2:dim(combis)[2]) {
								newList[[i]] <- c(Model$reacID[inReacsIndex[combis[1,i]]],Model$reacID[inReacsIndex[combis[2,i]]])
							}
							n <- 1	
							for(i in (dim(combis)[2]+1):(dim(combis)[2]+dim(combis2)[2])) {
								newList[[i]] <- c(Model$reacID[inReacsIndex[combis2[1,n]]],Model$reacID[inReacsIndex[combis2[2,n]]],Model$reacID[inReacsIndex[combis2[3,n]]])
								n <- n+1
							}	
							names(newList) <- newNames
							newList <- newList[!rexist]
							newANDs[ANDsadded:(ANDsadded+sum(!rexist)-1)] <- newList
							names(newANDs)[ANDsadded:(ANDsadded+sum(!rexist)-1)] <- names(newList)
							ANDsadded <- ANDsadded+sum(!rexist)
						}
					}	
				}				
			}
		}
	}
	
#	dup.index2 = c(1:length(Model$reacID))
#	dup.index2 = dup.index2[duplicated(Model$reacID)]
#	Model$reacID = Model$reacID[-dup.index2]
#	Model$interMat = Model$interMat[,-dup.index2]
#	Model$notMat = Model$notMat[,-dup.index2]
	
	if(!is.na(toAdd)) {
		Model$notMat = cbind(Model$notMat, toAdd$notMat)
		Model$interMat = cbind(Model$interMat, toAdd$interMat)
		Model$reacID = c(Model$reacID, toAdd$reacID)
	}
	
	ModelExp <- Model
	ModelExp$SplitANDs <- SplitANDs
	ModelExp$newANDs <- newANDs
	return(ModelExp)
}

