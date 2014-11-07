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
# $Id: createLBodeContPars.R 3184 2013-01-21 13:50:31Z cokelaer $
createLBodeContPars <-function
(
        model,            LB_n=1,            LB_k=0.1,
        LB_tau=1e-2,    UB_n=5,            UB_k=0.9,
        UB_tau=10,        default_n=3,    default_k=0.5,
        default_tau=1,    LB_in=c(),        UB_in=c(),
        opt_n=TRUE,        opt_k=TRUE,        opt_tau=TRUE,
        random=FALSE
)
{
    namesSpecies=model$namesSpecies;
    #convert the graph in incidence matrix to adjacency format
    adjMat=incidence2Adjacency(model);
    #Check what species are dynamic state states and which are inputs to the system
    isState=getStates(adjMat);
    numStates=as.integer(sqrt(length(adjMat)));
    parNames=c();
    parValues=c();
    index_k=c();
    index_n=c();
    index_tau=c();
    LB=c();
    UB=c();
    index_opt_pars=c();

    #Iterate over all species
    count=0;
    for(j in 1:numStates)
    {
        #According to the odefy method, if it a state there must one Tau,
        #and one k and n for each input to that state
        #Important note: Don't change the parameter ordering scheme. These are
        #ordered according to the adjacency matrix. The names are just for guidance
        if(isState[j])
        {
            jCol=adjMat[,j];
            inputs=which(as.logical(jCol));
            numInputs=length(inputs);
            #Add k's and n's for each input
            for(i in 1:numInputs)
            {
                #Add n
                count=count+1;
                parNames[count]=paste(namesSpecies[inputs[i]],"_n_",namesSpecies[j],sep="");
                parValues[count]=default_n;
                index_n=c(index_n,count);
                LB[count]=LB_n;
                UB[count]=UB_n;
                if(opt_n)index_opt_pars=c(index_opt_pars,count);

                #add k
                count=count+1;
                parNames[count]=paste(namesSpecies[inputs[i]],"_k_",namesSpecies[j],sep="");
                parValues[count]=default_k;
                index_k=c(index_k,count);
                LB[count]=LB_k;
                UB[count]=UB_k;
                if(opt_k)index_opt_pars=c(index_opt_pars,count);
            }
            #Add tau
            count=count+1;
            parNames[count]=paste("tau_",namesSpecies[j],sep="");
            parValues[count]=default_tau;
            index_tau=c(index_tau,count);
            LB[count]=LB_tau;
            UB[count]=UB_tau;
            if(opt_tau)index_opt_pars=c(index_opt_pars,count);
        }
    }

    #You pass your predefined vectors of upper al lower bounds
    if(length(LB_in)==length(LB))LB=LB_in;
    if(length(UB_in)==length(UB))UB=UB_in;

    #Create a new uniform random solution if random==TRUE
    if(random)parValues[index_opt_pars]=LB[index_opt_pars]+((UB[index_opt_pars]-LB[index_opt_pars])*runif(length(index_opt_pars)));

    parList=list(parNames=parNames,parValues=parValues,
            index_opt_pars=index_opt_pars,index_n=index_n,
            index_k=index_k,index_tau=index_tau,LB=LB,UB=UB);

    return(parList);
}


