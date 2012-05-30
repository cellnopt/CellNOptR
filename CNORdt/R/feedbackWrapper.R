feedbackWrapper <- function(Model) {

	what.loops = feedbackFinder1(Model)
	neg.reacs = c()
	for(b in 1:length(what.loops)) {
		loop1 = what.loops[[b]]
		loop1 = loop1[loop1 > 0]
		loop1 = c(loop1, loop1[1])
		for(a in 1:length(loop1)-1) {
			lhs = loop1[a]
			rhs = loop1[a+1]
			lhs.reac = which(Model$interMat[lhs,] == -1)
			rhs.reac = which(Model$interMat[rhs,] == 1)
			reac = intersect(lhs.reac, rhs.reac)
			if(any(Model$notMat[,reac] == 1)) {
				neg.reacs = c(neg.reacs, reac)	
			}
		}
	}
	neg.reacs = unique(neg.reacs)
	negEdges = rep(0,length(Model$reacID))
	negEdges[neg.reacs] = 1
	return(negEdges)
}