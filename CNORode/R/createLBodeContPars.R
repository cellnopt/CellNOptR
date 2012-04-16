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
#  CNO website: http://www.ebi.ac.uk/saezrodriguez/softare/cno
#
##############################################################################
# $Id$
createLBodeContPars <-function
(		
		model,			LB_n=1,			LB_k=0.1,
		LB_tau=1e-2,	UB_n=5,			UB_k=0.9,
		UB_tau=10,		default_n=3,	default_k=0.5,
		default_tau=1,	LB_in=c(),		UB_in=c(),
		opt_n=TRUE,		opt_k=TRUE,		opt_tau=TRUE,
		random=FALSE
)
{
	namesSpecies=model$namesSpecies;
	adjMat=incidence2Adjacency(model);
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
	
	count=0;
	for(j in 1:numStates)
	{ 
		if(isState[j])
		{
			jCol=adjMat[,j];
			inputs=which(as.logical(jCol));
			numInputs=length(inputs);
			for(i in 1:numInputs)
			{
				count=count+1;
				parNames[count]=paste(namesSpecies[inputs[i]],"_n_",namesSpecies[j],sep="");
				parValues[count]=default_n;
				index_k=c(index_k,count);
				LB[count]=LB_n;
				UB[count]=UB_n;
				if(opt_n)index_opt_pars=c(index_opt_pars,count);
				
				count=count+1;
				parNames[count]=paste(namesSpecies[inputs[i]],"_k_",namesSpecies[j],sep="");
				parValues[count]=default_k;
				index_n=c(index_k,count);
				LB[count]=LB_k;
				UB[count]=UB_k;
				if(opt_k)index_opt_pars=c(index_opt_pars,count);
			}
			count=count+1;
			parNames[count]=paste("tau_",namesSpecies[j],sep="");
			parValues[count]=default_tau;
			index_tau=c(index_tau,count);
			LB[count]=LB_tau;
			UB[count]=UB_tau;
			if(opt_tau)index_opt_pars=c(index_opt_pars,count);
		}
	}
	if(length(LB_in)==length(LB))LB=LB_in;
	if(length(UB_in)==length(UB))UB=UB_in;
	if(random)parValues=LB[index_opt_pars]+((UB[index_opt_pars]-LB[index_opt_pars])*runif(length(index_opt_pars)));
	
	parList=list(parNames=parNames,parValues=parValues,
			index_opt_pars=index_opt_pars,index_n=index_n,
			index_k=index_k,index_tau=index_tau,LB=LB,UB=UB);
	
	return(parList);
}

