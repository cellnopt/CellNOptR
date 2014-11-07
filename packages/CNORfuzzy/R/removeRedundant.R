#
#  This file is part of the CNO software
#
#  Copyright (c) 2011-2012 - EBI - Massachusetts Institute of Technology
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
# $id: removeRedundant.R 517 2012-02-07 16:57:57Z cokelaer $
removeRedundant = function(model) {
Reacs2Cut = array();
nums = dim(model$interMat)[1] 
for (eachSpec in 1:nums){
    LocalModel <- as.matrix(model$interMat[, model$interMat[eachSpec,]==1])
    if (sum(model$interMat[eachSpec,]==1) > 0) {
    origIndex = which(model$interMat[eachSpec,]==1);
    inputs = matrix(NA,nrow = nums, ncol = nums);
    numInputs = array(1, dim(LocalModel)[2]);
    for (j in 1:dim(LocalModel)[2]) {
		numInputs[j] = length(which(LocalModel[,j] == -1));
        inputs[1:numInputs[j],j] = which(LocalModel[,j] == -1);
    }
    for (eachEdge in 1:dim(LocalModel)[2]) {
        if (numInputs[eachEdge] > 1) {
           if (any(numInputs < numInputs[eachEdge])) {
           	# counterVec looks at inputs with fewer edges
                counterVec = array(0,sum(numInputs < numInputs[eachEdge]));
			    rmVec = array(TRUE,length(numInputs))
			    # only consider removing if it has fewer inputs
			    rmVec[numInputs >= numInputs[eachEdge]] = FALSE
			   	otherEdgeNumInputs = numInputs[rmVec]
                for (eachInput in 1:numInputs[eachEdge]) {
					currEdgeCount = 1
                    for (eachOtherEdge in which(rmVec)) {
						otherEdgesCurr = inputs[1:numInputs[eachOtherEdge],eachOtherEdge]
                    	counterVec[currEdgeCount] = counterVec[currEdgeCount] + any(otherEdgesCurr == inputs[eachInput,eachEdge]);
						currEdgeCount = currEdgeCount + 1
                    }
                }
                if (any(counterVec == otherEdgeNumInputs)) {
                    Reacs2Cut = c(Reacs2Cut,origIndex[eachEdge]);
                }
            }
        }
    }
    }   
}
Reacs2Cut = Reacs2Cut[!is.na(Reacs2Cut)]
return(Reacs2Cut)
}
