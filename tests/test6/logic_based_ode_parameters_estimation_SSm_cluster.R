logic_based_ode_parameters_estimation_SSm_cluster <-function
(
        cnolist,                model,                    ode_parameters=NULL,
        indices=NULL,            maxeval=Inf,            maxtime=100,            
        ndiverse=NULL,            dim_refset=NULL,         local_solver=NULL,      
        time=1,                    verbose=0,                 transfer_function=3,    
        reltol=1e-4,            atol=1e-3,                maxStepSize=Inf,        
        maxNumSteps=100000,        maxErrTestsFails=50,    n_nodes=2,
        n_iter=1
)
{
    
    library(eSSmR)
    library(Rsolnp)
    
    library(CellNOptR)
    library(CNORode)
    library(snow)
    
    cl = makeCluster(n_nodes, type="RCLOUD", pool="R4G");
    res=list();
    tryCatch(
    {
        myLib=.libPaths()[1];
       
        res=clusterEvalQ(cl,.libPaths(c("/net/nas20/rcloud_data/storage/wdir/davidh/RLibs",.libPaths())));
        res=clusterEvalQ(cl,.libPaths());
        res=clusterEvalQ(cl, library(CNORode));
        n_nodes=length(cl);
        print(paste("We could only allocation",n_nodes,"nodes"));
        if(is.null(ode_parameters))
        {
            adjMat=incidence2Adjacency(model);
            ode_parameters=makeParameterList(adjMat,model$namesSpecies,random=TRUE);
        }
        if(is.null(indices))indices <- indexFinder(cnolist,model,verbose=FALSE);
        
        res=clusterApplyLB(cl,1:n_nodes,logic_based_ode_parameters_estimation_SSm,cnolist,model,
            ode_parameters,indices,maxeval,maxtime,ndiverse,dim_refset,local_solver,time,verbose,
            transfer_function,reltol,atol,maxStepSize,maxNumSteps,maxErrTestsFails);

        if(n_iter>0)
        {
            for(i in 1:n_iter)
            {
                for(i in 1:n_nodes)
                {
                    mat=rbind(mat,res[[i]]$Refset_x);
                }
                rand_index=sample(n_nodes*dim_refset);
                mat=mat[rand_index,];
                x_0=list();
                count=0;
                for(j in 1:n_nodes)
                {
                    print(j);
                    count=count+1;
                    x_0[[j]]=mat[count:(count+dim_refset-1),]
                }
                res=clusterApplyLB(cl,1:n_nodes,logic_based_ode_parameters_estimation_SSm,cnolist,model,
                ode_parameters,indices,maxeval,maxtime,ndiverse,dim_refset,local_solver,time,verbose,
                transfer_function,reltol,atol,maxStepSize,maxNumSteps,maxErrTestsFails,x_0); 
            }                     
        }
        }
,error=function(e){print("something happened")},finally={stopCluster(cl);} )
 try(stopCluster(cl)); 
 return(res);
 
}

logic_based_ode_parameters_estimation_SSm <-function
(
        i,                        cnolist,                model,                    ode_parameters=NULL,
        indices=NULL,            maxeval=Inf,            maxtime=100,            
        ndiverse=NULL,            dim_refset=NULL,         local_solver=NULL,      
        time=1,                    verbose=0,                 transfer_function=3,    
        reltol=1e-4,            atol=1e-3,                maxStepSize=Inf,        
        maxNumSteps=100000,        maxErrTestsFails=50,    x_0=NULL
)
{
    
    library(eSSmR)
    library(Rsolnp)
    
    if(is.null(ode_parameters))
    {
        adjMat=incidence2Adjacency(model);
        ode_parameters=makeParameterList(adjMat,model$namesSpecies,random=TRUE);
    }
    if(is.null(indices))indices <- indexFinder(cnolist,model,verbose=FALSE);
    if(!is.null(x_0))
    {
        if(is.list(x_0))x_0=x_0[[i]];
       
    }
    print(x_0)
    
    problem=list();
    problem$f<-get_logic_ode_continuous_objective_function(cnolist,model,ode_parameters,indices,
            time,verbose,transfer_function,reltol,atol,maxStepSize,maxNumSteps,maxErrTestsFails);
    problem$x_L <- ode_parameters$LB[ode_parameters$index_opt_pars];
    problem$x_U <- ode_parameters$UB[ode_parameters$index_opt_pars];
    problem$x_0<- ode_parameters$parValues[ode_parameters$index_opt_pars];
    problem$int_var =0;
    problem$bin_var =0;
    opts=list();
    opts$maxeval=maxeval;
    opts$maxtime=maxtime;
    if(!is.null(local_solver))opts$local_solver=local_solver;
    if(!is.null(ndiverse))opts$ndiverse=ndiverse;      
    if(!is.null(dim_refset))opts$dim_refset=dim_refset; 
    if(!is.null(x_0))problem$x_0=x_0;
    
    return(essR(problem,opts));    
}


