feedbackWrapper <- function(model) {

	whatLoops = feedbackFinder(model)
	negReacs = c()
	
	if(length(whatLoops) != 0) {
		
		for(b in 1:length(whatLoops)) {
			loop1 = whatLoops[[b]]
			loop1 = loop1[loop1 > 0]
			loop1 = c(loop1, loop1[1])
			for(a in 1:length(loop1)-1) {
				lhs = loop1[a]
				rhs = loop1[a+1]
				LHSreac = which(model$interMat[lhs,] == -1)
				RHSreac = which(model$interMat[rhs,] == 1)
				reac = intersect(LHSreac, RHSreac)
				if(any(model$notMat[,reac] == 1)) {
					negReacs = c(negReacs, reac)	
				}
			}
		}
		negReacs = unique(negReacs)
	#	negEdges = rep(0,length(model$reacID))
	#	negEdges[negReacs] = 1
		return(negReacs)
	} else {
		
		print("No loops found!")
		return(NULL)	
	}
}