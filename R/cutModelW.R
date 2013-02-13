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
# $Id: cutModelW.R 68953 2012-08-30 08:59:47Z t.cokelaer $

cutModelW <- function(model, bString){
    bs = as.logical(bString)
    newmodel <- list()
    newmodel$interMat <- model$interMat[, bs]
    newmodel$notMat <- model$notMat[, bs]
    newmodel$reacID <- model$reacID[bs]
	newmodel$linksWeights <- model$linksWeights[bs]  # added
    newmodel$namesSpecies <- model$namesSpecies

    # could also add the times used in times > T1 if times
    # newmodel$times <- bStringTimes[which(bStringTimes != 0)]

    return(newmodel)
}


