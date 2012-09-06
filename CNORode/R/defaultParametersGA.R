$Id
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