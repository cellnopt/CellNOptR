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
# $id: defaultParametersFuzzy.R 675 2012-03-09 11:43:44Z cokelaer $

defaultParametersFuzzy<-function(data=NA, model=NA, nTF=7){


    # data is a CNOlist as returned by makeCNOlist or CNOlist class from
    # CellNOptR
    paramsList<-list()
    paramsList$data<-data
    paramsList$model<-model

    # These are the parameters for the transfer functions. nrow indicates the
    # number of transfer functions.
    #  * column 1 is the parameter G
    #  * column 2 is the parameter N equal to 3 for all
    #  * column 3 is the parameter K
    # The parameters K are chosen so that EC50 is 0.2,0.3,0.4,0.5,0.6,0.7 and 0.5
    Nvalue = 3
    paramsList$type1Funs<-matrix(data = NaN,nrow=nTF,ncol=3)
    paramsList$type1Funs[,1] = 1
    paramsList$type1Funs[,2] = c(rep(Nvalue, nTF-1), 1.01)
    paramsList$type1Funs[,3] = .getk(nTF)
    paramsList$type1Funs[,3] = c(0.2, 0.3, 0.4, 0.55, 0.72,1.03, 68.5098)
 # from matlab code ...
#    paramsList$type1Funs[,3] = c(0.2011, 0.3056, 0.4187, 0.5503, 0.7245,1.03,68.9058)






    paramsList$type2Funs<-matrix(data = NaN,nrow=nTF,ncol=3)
    paramsList$type2Funs[,1]=seq(0.2, 0.8, length.out=nTF)
    paramsList$type2Funs[,2]=1
    paramsList$type2Funs[,3]=1

    #paramsList$type2Funs<-matrix(data = NaN,nrow=6,ncol=3)
    #paramsList$type2Funs[,1]=c(0.2, 0.3, 0.4, 0.5, 0.6, 0.7)
    #paramsList$type2Funs[,2]=1
    #paramsList$type2Funs[,3]=1

    paramsList$redThresh = c(0, 0.0001, 0.0005, 0.001, 0.003, 0.005, 0.01)

    paramsList$doRefinement = TRUE

    # GA parameters
    paramsList$sizeFac<-1e-04
    paramsList$NAFac<-1
    paramsList$popSize<-50
    paramsList$pMutation<-0.5
    paramsList$maxTime<-3*60
    paramsList$maxGens<-500
    paramsList$stallGenMax<-100
    paramsList$selPress<-1.2
    paramsList$elitism<-5
    paramsList$relTol<-0.1
    paramsList$verbose<-FALSE

    # optimisation parameters in getRefinedModel
    paramsList$optimisation = list(algorithm='NLOPT_LN_SBPLX', xtol_abs=0.001,maxEval=10000, maxTime=5*60)

    return(paramsList)
}



.getk <- function(n=20){

   # TC
   # Given N = [3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 1.01]
   # and EC50 = seq(0.2, 0.7, length.out=19) + 0.5
   # we can obtain k using fsolve (from python code see fuzzytools)
   k = c(0.20107819, 0.22960113, 0.258464, 0.28776477, 0.31761979, 0.34816775,
       0.37957511, 0.41204356, 0.44582066, 0.48121482, 0.51861737, 0.55853549,
       0.60164334, 0.64886436, 0.7015105, 0.761532, 0.83200087, 0.91814475,
        1.02988364)
   # last value is harcoded to 68.5098

   interpolated_k = approx(seq(1,19), k, n=n-1)
   results = c(round(interpolated_k$y, 2), 68.5098)
   return(results)

}

# for bookkeeping, may be usefule one day.
.hill_function <- function(k, n=3, EC50=0.5){
    results = (1+k**n)*EC50**n/(EC50**n+k**n) -0.5
    return(matrix(results))

}



