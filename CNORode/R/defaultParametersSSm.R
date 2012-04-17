defaultParametersSSm <-
function(){

    params = list()
    params$maxeval=Inf
    params$maxtime=100
    params$ndiverse=NULL
    params$dim_refset=NULL
    params$local_solver=NULL
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