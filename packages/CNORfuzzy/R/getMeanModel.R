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
#  CNO website: http://www.ebi.ac.uk/saezrodriguez/software/cno
#
##############################################################################
getMeanModel<-function(postRefThresh=NULL, allRes, doRefinement=TRUE)
{

    # compute the 
    summary = compileMultiRes(allRes,show=FALSE)
    allFinalMSEs = summary$allFinalMSEs

    # initialize some things
    ReferenceRes = allRes[[1]]
    CNOlist=ReferenceRes$paramsList$data

     if ((class(CNOlist)=="CNOlist")==FALSE){
        CNOlist = CellNOptR::CNOlist(CNOlist)
    }

    indexList<-indexFinder(CNOlist=CNOlist,model=ReferenceRes$processedModel,verbose=ReferenceRes$paramsList$verbose)

    MSE = c()
    common_currIX = max(which(allFinalMSEs[1,]-allFinalMSEs[1,2]<=postRefThresh))
    # for each result, pick the best refined model, simulate it, and save the simulation results
    if(ReferenceRes$paramsList$doRefinement && doRefinement==TRUE){
        for (eachRes in 1:length(allRes)){
            currIX = max(which(allFinalMSEs[eachRes,]-allFinalMSEs[eachRes,2]<=postRefThresh))
            print(currIX)
            if (common_currIX!=currIX){
                print("?skipping a cube because dimensions are not consistent.")
                #next()
            }
            currModel = allRes[[eachRes]]$redRef[[currIX]]$refModel$refinedModel
            currSimList = allRes[[eachRes]]$redRef[[currIX]]$refModel$refinedSimList

            MSE = rbind(MSE, allRes[[eachRes]]$redRef[[currIX]]$MSE)
            if (eachRes == 1){
                kCube = currSimList$kCube

                nR = dim(kCube)[1]
                nC = dim(kCube)[2]

                kCube2 = array(currSimList$kCube, dim=c(nR, nC, 1))
                nCube2 = array(currSimList$nCube, dim=c(nR, nC, 1))
                gCube2 = array(currSimList$gCube, dim=c(nR, nC, 1))

            } else{
                kCube2 = array(kCube2, currSimList$kCube, dim=c(nR, nC, dim(kCube2)[3]+1))
                gCube2 = array(gCube2, currSimList$gCube, dim=c(nR, nC, dim(gCube2)[3]+1))
                nCube2 = array(nCube2, currSimList$nCube, dim=c(nR, nC, dim(nCube2)[3]+1))
            }
          }
        }
        else{
          for (eachRes in 1:length(allRes)){
            currModel = allRes[[eachRes]]$unRef$model
            currSimList = allRes[[eachRes]]$unRef$simList
            MSE = rbind(MSE, allRes[[eachRes]]$unRef$MSE)
            if (eachRes == 1){
                kCube = currSimList$kCube
                nR = dim(kCube)[1]
                nC = dim(kCube)[2]

                kCube2 = array(kCube,  dim=c(nR, nC, 1))
                nCube2 = array(currSimList$nCube,  dim=c(nR, nC, 1))
                gCube2 = array(currSimList$gCube,  dim=c(nR, nC, 1))

            }else{
                kCube2 = array(kCube2, currSimList$kCube, dim=c(nR, nC, dim(kCube2)[3]+1))
                gCube2 = array(gCube2, currSimList$gCube, dim=c(nR, nC, dim(gCube2)[3]+1))
                nCube2 = array(nCube2, currSimList$nCube, dim=c(nR, nC, dim(nCube2)[3]+1))
            }
          }
        }

    model = allRes[[1]]$redRef[[currIX]]$refModel$refinedModel
    names = rownames(allRes[[1]]$redRef[[currIX]]$refModel$refinedSimList$finalCube)
    typeCube = allRes[[1]]$redRef[[currIX]]$refModel$refinedSimList$typeCube

    sens = kCube2
    for (i in 1:dim(kCube2)[1]){
        for (j in 1:dim(kCube2)[2]){
            for (k in 1:dim(kCube2)[3]){
                if (is.na(typeCube[i,j])){
                    sens[i,j,k] = NA
                } else{
                    if (typeCube[i,j]==1){ 
                        sens[i,j,k] = 1 - getEC50(kCube2[i,j,k], nCube2[i,j,k])$solution
                    } 

                    if (typeCube[i,j] == 2) {
                        sens[i,j,k] = 0.5 * gCube2[i,j,k]
                    }
                }
            }
        }
    }   

    # compute mean and std of interest
    meanSensitivity = .mean(sens)
    stdSensitivity = .std(sens)

    # and for debugging
    kCube = .mean(kCube2)
    gCube = .mean(gCube2)
    nCube = .mean(nCube2)

    # set the names
    rownames(kCube) = names
    rownames(nCube) = names
    rownames(gCube) = names
    rownames(stdSensitivity) = names
    rownames(meanSensitivity) = names

return(list(kCube=kCube, gCube=gCube, nCube=nCube2, model=model, names=names,
sens=sens, meanSensitivity=meanSensitivity, stdSensitivity=stdSensitivity,
MSE=MSE))

}



# compute mean and std over all results (on third axis
.mean <- function(data){
    nR = dim(data[,,1])[1]
    nC = dim(data[,,1])[2]
    res = array(NA, nR*nC, dim=c(nR, nC))
    for (i in 1:dim(data)[1]){
        for (j in 1:dim(data)[2]){
            res[i,j] = mean(data[i,j,])
        }
    }
    return(res)
}


.std <- function(data){

    nR = dim(data[,,1])[1]
    nC = dim(data[,,1])[2]
    res = array(NA, nR*nC, dim=c(nR, nC))
    for (i in 1:dim(data)[1]){
        for (j in 1:dim(data)[2]){
            res[i,j] = sqrt(var(data[i,j,]))
        }
    }
    return(res)
}


