
getEC50 <- function(k, n){

    ConOptions=list(
        algorithm='NLOPT_LN_SBPLX',
        #algorithm="NLOPT_GN_DIRECT",
        xtol_rel=1e-8,
        ftol_rel=1e-8,
        #stopval=0.0000001,
        print_level=0,
        maxeval=10000,
        maxtime=5*60)

    lowBound = 0
    upBound = 10
    guess = 0.27
    res = nloptr(guess,objFun_getEC50,lb=lowBound,ub=upBound,opts=ConOptions,
        k=k, n=n)

    return(list(res=res, solution=res$solution, MSE=res$objective))
}


objFun_getEC50 <- function(x50, k, n){
    res = x50**n * (1. + k**n) - 0.5*(x50**n + k**n) 
    # note that we return the abs value becqause the optimisation looks for the
    # minimum.
    return(abs(res))
}
