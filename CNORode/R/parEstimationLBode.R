parEstimationLBode<-function (cnolist, model, method="ga", 
    ode_parameters = NULL, indices = NULL, paramsGA=NULL, paramsSSM=NULL)
{

    if (method == "essm"){
        if (is.null(paramsSSM)){
            paramsSSM = defaultParametersSSm()
        }
        ode_parameters = parEstimationLBodeSSm(cnolist, model, 
            ode_parameters=ode_parameters, indices=indices,
            maxeval=paramsSSM$maxeval, 
            maxtime=paramsSSM$maxtime,
            ndiverse=paramsSSM$ndiverse, 
            dim_refset=paramsSSM$dim_refset, 
            local_solver=paramsSSM$local_solver, 
            time=paramsSSM$time, 
            verbose=paramsSSM$verbose, 
            transfer_function=paramsSSM$transfer_function, 
            reltol=paramsSSM$reltol, 
            atol=paramsSSM$atol, 
            maxStepSize=paramsSSM$maxStepSize, 
            maxNumSteps=paramsSSM$maxNumSteps, 
            maxErrTestsFails=paramsSSM$maxErrTestsFails, 
            nan_fac=paramsSSM$nan_fac) 


    }
    else if(method=="ga"){
        if (is.null(paramsGA)){
            paramsGA = defaultParametersGA()
        }
        ode_parameters = parEstimationLBodeGA(cnolist, model, 
            ode_parameters=ode_parameters, 
            indices=indices,
            mutationChance=paramsGA$mutationChance, 
            popSize=paramsGA$popSize, 
            iters=paramsGA$iters, 
            elitism=paramsGA$elitism, 
            time=paramsGA$time, 
            monitor=paramsGA$monitor,
            verbose=paramsGA$verbose, 
            transfer_function=paramsGA$transfer_function, 
            reltol=paramsGA$reltol, 
            atol=paramsGA$atol, 
            maxStepSize=paramsGA$maxStepSize, 
            maxNumSteps=paramsGA$maxNumSteps, 
            maxErrTestsFails=paramsGA$maxErrTestsFails, 
            nan_fac=paramsGA$nan_fac) 
    }
    else{
        stop ("method argument must be either 'ga' or 'essm'." )
    }
    return(ode_parameters)
}

defaultParametersGA <- function(){

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

defaultParametersSSm <- function(){

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
