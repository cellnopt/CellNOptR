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
# $Id: test_CNORwrapFuzzy.R 2339 2012-09-06 09:08:02Z cokelaer $
# This is a functional test of the CNORwrapFuzzy function
# author: T.Cokelaer, 2012

test_CNORwrapFuzzy <-function(){

    library(CNORfuzzy)
    # Load some data
    data(CNOlistToy, package="CellNOptR")
    data(ToyModel, package="CellNOptR")
    # Get some default parameters to play with, limiting the duration of the GA
    # algorithm and optimisation step
    paramsList = defaultParametersFuzzy() 
    paramsList$maxTime = 20
    paramsList$verbose = FALSE
    
    nrow = 7
    paramsList$type1Funs = matrix(data = NaN,nrow=nrow,ncol=3)
    paramsList$type1Funs[,1] = 1
    paramsList$type1Funs[,2] = c(3, 3, 3, 3, 3, 3, 1.01)
    paramsList$type1Funs[,3] = c(0.2, 0.3, 0.4, 0.55, 0.72,1.03, 68.5098)
    
    # Default Fuzzy Logic Type2 parameters
    nrow = 7
    paramsList$type2Funs = matrix(data = NaN,nrow=nrow,ncol=3)
    paramsList$type2Funs[,1] = seq(from=0.2, to=0.8, length=nrow)
    paramsList$type2Funs[,2] = 1
    paramsList$type2Funs[,3] = 1
    
    paramsList$redThres = c(0, 0.0001, 0.0005, 0.001, 0.003, 0.005, 0.01)
    

    paramsList$optimisation$algorithm = "NLOPT_LN_SBPLX"
    paramsList$optimisation$xtol_abs = 0.01
    paramsList$optimisation$maxeval = 1000
    paramsList$optimisation$maxtime = 10
    

    results = CNORwrapFuzzy(CNOlist(CNOlistToy), ToyModel, paramsList, verbose=TRUE)
}

