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
# $id: compileMultiRes.R 1484 2012-06-14 15:23:41Z schrier $
compileMultiRes = function(allRes,tag=NULL, show=TRUE){

    IxAdj = 0

    N = length(allRes)
    # first create objects to store  number of parameters and MSEs in summary matrices
    for (i in 1:N) {
        Res = allRes[[i]]
        if (!is.null(Res$redRef)){
            if (i==1){
                allFinalMSEs = matrix(NA, N, length(Res$redRef))
                allFinalNumParams = matrix(NA, N,length(Res$redRef))
            }
            for (q in 1:length(Res$redRef)){
                allFinalMSEs[i,q] = Res$redRef[[q]]$MSE
                allFinalNumParams[i,q] = length(Res$redRef[[q]]$finalSet)
            }
        }
        else {
            if (i==1){
               # allFinalMSEs = matrix(NA, N,length(Res$unRef)) This length is the number of fields!!
              allFinalMSEs=matrix(NA,N,length(Res$unRef$MSE))
            }
            for (q in 1:length(Res$unRef$MSE)){
                allFinalMSEs[i,q] = Res$unRef$MSE
            }
        }
    }
    rm(Res)

    # save results in a file if required only
    if (is.null(tag) == FALSE){
        save(allRes,file=paste(tag,'_allRes.RData',sep=''))
        save(allFinalMSEs,file=paste(tag,'_allFinalMSEs.RData',sep=''))
        if (exists('allFinalNumParams')){
            save(allFinalNumParams,file=paste(tag,'_allFinalNumParams.RData',sep=''))
        }
    }

    # now plot what the average MSE and number of parameters would be using different selection thresholds
    if (exists('allFinalNumParams')){

        Thresh2New = c(0.0001,0.0005,0.001,0.002,0.003,0.004,0.005,0.006,0.007,0.008,0.009,0.01,0.013,0.015,0.017,0.02,0.025,0.03,0.05,0.1,0.2,0.3,0.5)

        catExplore = matrix(NA,dim(allFinalMSEs)[1],length(Thresh2New))
        for (i in 1:length(Thresh2New)){
            for (j in 1:dim(allFinalMSEs)[1]){
                currIX = which(allFinalMSEs[j,]-allFinalMSEs[j,2]<=Thresh2New[i])
                catExplore[j,i] = max(currIX)
            }
        }

        AbsNumParams = matrix(NA,dim(allFinalMSEs)[1],length(Thresh2New))
        AbsMSEs = matrix(NA,dim(allFinalMSEs)[1],length(Thresh2New))
        for (i in 1:dim(allFinalMSEs)[1]){
            AbsNumParams[i,] = allFinalNumParams[i,catExplore[i,]]
            AbsMSEs[i,] = allFinalMSEs[i,catExplore[i,]]
        }
        meanMSEs = colMeans(AbsMSEs)
        meanNPs = colMeans(AbsNumParams)

        if (show==TRUE){
            par(mar=c(5,4,4,5)+.1)
            plot(Thresh2New,meanMSEs,type="l",col="red",log="x")
            par(new=TRUE)
            plot(Thresh2New, meanNPs,,type="l",col="blue",xaxt="n",yaxt="n",xlab="",ylab="",log="x")
            axis(4)
            mtext("Mean Number of Parameters",side=4,line=3)
            legend("topleft",col=c("red","blue"),lty=1,legend=c("Mean MSEs","Mean Number of Parameters"))
        }
    }

    if (exists('allFinalNumParams')){
    return(list(allFinalMSEs=allFinalMSEs, allFinalNumParams=allFinalNumParams))}
    else{
      return(list(allFinalMSEs=allFinalMSEs))
    }
}
