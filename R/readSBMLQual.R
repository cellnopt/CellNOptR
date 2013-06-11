#
#  This file is part of the CNO software
#
#  Copyright (c) 2011-2013 - EMBL-EBI
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

# This function reads a SBMLQual written for storing SIF file. 
# The listOfQualitativeSpecies should contain the species
# From the listOfTransitions, the edgse can be extracted
# TODO: what if there is no transittion for a given specy. Is it possible ? If
# so this is an isolated node.

readSBMLQual <- function(filename){

    library(XML)
    doc = xmlTreeParse(filename)
    r = xmlRoot(doc)

    sif = data.frame(matrix(0, ncol=3))

    species = xmlChildren(xmlChildren(r)[[1]])['listOfQualitativeSpecies']
    nameSpecies = as.vector(unlist(xmlApply(species[[1]], function(x) xmlGetAttr(x, "id"))))

    transitions = xmlChildren(xmlChildren(r)[[1]])['listOfTransitions']

    for (i in seq_along((transitions[[1]]))){
        # somehow we can not loop over the transitions so we use
        # seq_along(length(transitions))
        transition = transitions[[1]][[i]]
        outputs = xmlChildren(transition)$listOfOutputs
        inputs = xmlChildren(transition)$listOfInputs
        signs = as.vector(unlist(sapply(xmlChildren(inputs), function(x) xmlGetAttr(x,"sign"))))
        LHS = as.vector(unlist(sapply(xmlChildren(inputs), function(x) xmlGetAttr(x,"qualitativeSpecies"))))
        RHS = as.vector(unlist(sapply(xmlChildren(outputs), function(x) xmlGetAttr(x,"qualitativeSpecies"))))

        for (j in seq_along(LHS)){

            if (signs[[j]] == "positive"){
                sif = rbind(c(LHS[[j]],1,RHS), sif)
            } else if (signs[[j]] == "negative") {
                sif = rbind(c(LHS[[j]],-1,RHS), sif)
            } else {
                stop("signs must be negative or positive")
            }
        }

    }

    # remove last row (0,0,0) set at the beginning to define the data frame.
   sif = sif[-dim(sif)[[1]],]

   fh = tempfile()
   print(fh)
   write.table(sif,file=fh,
            row.names=FALSE,col.names=FALSE,quote=FALSE,sep="\t")

    sif = readSIF(fh)
    return(sif)

}
