# feedbackFinder
# given a model, pick out negative edges that are feedback
# TODO path from output1 is confused by other feedback loops: *** COMEBACK TO ***

#feedbackFinder <- function(Model) {
#	
##findInput <- function(x) {
##	inNode<-which(x == -1)
##}	
#	reac.check = rep(0, length(Model$reacID))
#	names(reac.check) = Model$reacID
#	for (i in 1:length(Model$reacID)) {
#		print(i)
#		# is there a negative edge?
#		if(any(Model$notMat[,i] == 1)) {
#		
#			# what node does the negative edge
#			# originate at
#			output1 = which(Model$notMat[,i]==1)
#			# what is reaction where output1
#			# is a negative input?
#			outReacNeg = which(Model$interMat[output1,]==-1)
#		
#			# what are the output nodes here? 
#			output2 = c()
#			for(j in 1:length(outReacNeg)) {
#				output2 = c(output2, which(Model$interMat[,outReacNeg[j]] == 1))	
#			}
#			
#			# get the outputs of the outputs 
#			outputLoop = output2
#			end.flag = 0
#			
#			while(length(intersect(output1,outputLoop))==0 && end.flag==0) {
#				outputNew = c()
#				for(k in 1:length(outputLoop)) {
#					reac = which(Model$interMat[outputLoop[k],] == -1)
#					for(l in 1:length(reac)) {
#						outputNew = c(outputNew, which(Model$interMat[,reac[l]] == 1))
#					}
#				}
#				
#				outputLoop = outputNew
#				if(length(outputLoop)==0) {
#					end.flag = 1	
#				}
#			}
#			
#			if(length(intersect(output1,outputLoop))) {
#				reac.check[i] = 1
#			}
#		}		
#	}
#	return(reac.check)
#}

# go through each reaction
# pick the output of the neg. edge
# enter while loop
# while output != NOI (node of interest)
# and output node is available
# outputToCheck  = intermat[where input == -1 and output ==1]

# implement Tiernan's algorithm
# takes as input = 'tiernanM'
# may have to run the algorithm from each input?

#Model = tiernanM

feedbackFinder1 <- function(Model) {

V = 1:length(Model$namesSpecies)
Gamma = vector("list",length(V))
names(Gamma) = Model$namesSpecies
for(a in 1:length(V)) {
	pickCols = which(Model$interMat[V[a],]==-1)
	for(b in pickCols) {
		Gamma[[a]] = c(Gamma[[a]], which(Model$interMat[,b]==1))
	}
	if(!is.null(Gamma[[a]])) {
		Gamma[[a]] = unique(Gamma[[a]])
	}
}

H = matrix(0, nrow=length(Model$namesSpecies), ncol=length(Model$namesSpecies))
k = 1
P = 0
P[1] = V[1]
flagT = 1
flagUpdate = c()
j = 1
loops = list()

while (flagT == 1) {

	if(!is.null(Gamma[[P[k]]])) {
		if(length(setdiff(Gamma[[P[k]]],H[P[k],]))) {
			what.diff = setdiff(Gamma[[P[k]]],H[P[k],])
			openEdges = c() 
			for(a in 1:length(what.diff)) {
				openEdges = c(openEdges, which(Gamma[[P[k]]] == what.diff[a]))
			}
		} else {
			openEdges = 1	
		}
	
		if (Gamma[[P[k]]][openEdges[1]] > P[1] && !any(P == Gamma[[P[k]]][openEdges[1]]) && !any(H[P[k],] == Gamma[[P[k]]][openEdges[1]])) {
			flagUpdate = 1	
		} else {
			flagUpdate = 0	
		}
	} else {
		openEdges = 1
		flagUpdate = 0
	}

		
	while (flagUpdate == 1) {
		what.diff = setdiff(Gamma[[P[k]]],H[P[k],])
		openEdges = c() 
		for(a in 1:length(what.diff)) {
			openEdges = c(openEdges, which(Gamma[[P[k]]] == what.diff[a]))
		}
		if (Gamma[[P[k]]][openEdges[1]] > P[1] && !any(P == Gamma[[P[k]]][openEdges[1]]) && !any(H[P[k],] == Gamma[[P[k]]][openEdges[1]])) {
			k = k+1
			P[k] = Gamma[[P[k-1]]][openEdges[1]]
			if(!is.null(Gamma[[P[k]]])) {
				flagUpdate = 1
			} else {
				flagUpdate = 0	
			}	
		} else {
			flagUpdate = 0	
		}	
	}
	
	# check to report circuit here
	if(!is.null(Gamma[[P[k]]])) {
		if(P[1] == Gamma[[P[k]]][openEdges[1]]) {
			loops = c(loops,list(P))	
		}
	}
	
	# check for closure here
	if(k==1) {
		if(P[1] == length(Model$namesSpecies)) {
			# terminate
			flagT = 0	
		} else {
			# go to next node
			P[1] = P[1] + 1
			k = 1
			H = matrix(0, nrow=length(Model$namesSpecies), ncol=length(Model$namesSpecies))
		}
		
	} else {
	#	print(P)
		m = which(H[P[k-1],] == 0)[1]
		H[P[k],] = 0
		H[P[k-1],m] = P[k]
		P[k] = 0			
		k = k-1
	}
}

return(loops)
}



