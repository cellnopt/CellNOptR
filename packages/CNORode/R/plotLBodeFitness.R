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
# $Id: plotLBodeFitness.R 4446 2014-03-12 15:35:32Z cokelaer $
plotLBodeFitness <-function
(
        cnolist,                    model,                        ode_parameters=NULL,
        indices=NULL,                adjMatrix=NULL,                 time=1,
        verbose=0,                     transfer_function=3,        reltol=1e-4,
        atol=1e-3,                    maxStepSize=Inf,             maxNumSteps=100000,
        maxErrTestsFails=50,        plot_index_signals=NULL,    plot_index_experiments=NULL,
        plot_index_cues=NULL,         colormap="heat",
        plotParams=list(margin=0.1, width=15, height=12,
                  cmap_scale=1, cex=1.6, ymin=NULL), initial_state=0.1
  

)
{

    if (class(cnolist)=="CNOlist"){ cnolist = compatCNOlist(cnolist)}

    if(is.null(plot_index_experiments))plot_index_experiments=1:dim(cnolist$valueCues)[1];
    if(is.null(plot_index_cues))plot_index_cues=1:dim(cnolist$valueCues)[2];
    if(is.null(plot_index_signals))plot_index_signals=1:dim(cnolist$valueSignals[[1]])[2];
    if(is.null(indices))indices=indexFinder(cnolist,model);
    if(is.null(adjMatrix))adjMatrix=incidence2Adjacency(model);
    if(is.null(ode_parameters))ode_parameters=createLBodeContPars(model);
    timeSignals=NULL;

    sim_data=getLBodeDataSim(cnolist,model,
            ode_parameters,indices,timeSignals,time,verbose,
            transfer_function,reltol,atol,maxStepSize,maxNumSteps,
            maxErrTestsFails, initial_state=initial_state);

    times=cnolist$timeSignals;

    sim_data=lapply(sim_data,function(x)x[plot_index_experiments,plot_index_signals]);
    expResults=lapply(cnolist$valueSignals,function(x)x[plot_index_experiments,plot_index_signals]);

    sim_data=lapply(sim_data,function(x)as.matrix(x));
    expResults=lapply(expResults,function(x)as.matrix(x));

    if(dim(cnolist$valueSignals[[1]])[1]!=dim(sim_data[[1]])[1]){
        sim_data=lapply(sim_data,function(x) t(as.matrix(x)));
    }

    namesSignals=cnolist$namesSignals[plot_index_signals];
    namesCues=c(cnolist$namesStimuli,cnolist$namesInhibitors);

    valueCues=cbind(cnolist$valueStimuli,cnolist$valueInhibitors);
    valueCues=as.matrix(valueCues);
    valueCues[which(valueCues>0)]=1;
    valueCues=valueCues[plot_index_experiments,plot_index_cues];
    names(valueCues)=namesCues[plot_index_cues];

    if (colormap=="green"){
        plotOptimResults(simResults=sim_data,expResults=expResults,
        times=times,namesCues=namesCues,namesSignals=namesSignals,
        valueCues=valueCues);
    } else{
        plotOptimResultsPan(sim_data, yInterpol=NULL, xCoords=NULL,
             CNOlist=CNOlist(cnolist), formalism="ode", pdf=FALSE,
             plotParams=plotParams,pdfFileName="", tPt=NULL)
    }

  return(sim_data);
}
