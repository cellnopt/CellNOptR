#
#  This file is part of the CNO software
#
#  Copyright (c) 2011-2013 - EBI
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
# $Id: defaultParametersGA.R 3184 2013-01-21 13:50:31Z cokelaer $
defaultParametersGA <-
 function(){

    params = list()
    params$mutationChance=NA
    params$popSize=200
    params$iters=100
    params$elitism=NA
    params$time = 1
    params$monitor=TRUE
    params$verbose = 0
    params$transfer_function = 3
    params$reltol = 1e-04
    params$atol = 0.001
    params$maxStepSize = Inf 
    params$maxNumSteps = 1e+05
    params$maxErrTestsFails = 50
    params$nan_fac = 1

    return(params)
}
