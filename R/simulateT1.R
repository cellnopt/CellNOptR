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
simulateT1<-function(CNOlist, model, bStringT1, simList, indexList){

    Modelcut<-cutModel(model, bStringT1)


    simListcut<-simList
    simListcut$finalCube<-simListcut$finalCube[as.logical(bStringT1),]
    simListcut$ixNeg<-simListcut$ixNeg[as.logical(bStringT1),]
    simListcut$ignoreCube<-simListcut$ignoreCube[as.logical(bStringT1),]
    simListcut$maxIx<-simListcut$maxIx[as.logical(bStringT1)]

    if(is.null(dim(simListcut$finalCube))){
        simListcut$finalCube<-matrix(simListcut$finalCube,ncol=1)
        simListcut$ixNeg<-matrix(simListcut$ixNeg,ncol=1)
        simListcut$ignoreCube<-matrix(simListcut$ignoreCube,ncol=1)
        }

    simRes<-simulatorT1(
        CNOlist=CNOlist,
        model=Modelcut,
        simList=simListcut,
        indexList=indexList)

    return(simRes)
}

