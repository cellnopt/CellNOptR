#
#  This file is part of the CNO software
#
#  Copyright (c) 2011-2013 - EMBL-EBI
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
# $Id: simdata2cnolist.R 4036 2013-09-24 16:09:40Z bernardo $

simdata2cnolist <- function(sim_data, cnolist, model){ 

    if (class(cnolist)=="CNOlist"){cnolist = compatCNOlist(cnolist)}
    adjMatrix=incidence2Adjacency(model);
    states_index=which(as.logical(getStates(adjMatrix)));

    sim_data=lapply(sim_data,function(x) x[,states_index, drop=F]);

    times=cnolist$timeSignals;
    cnolist$valueSignals=sim_data;
 
    cnolist$valueVariances = sim_data # create the data structure
    ## set values to NA
    for (i in seq_along(cnolist$valueVariances)){
        cnolist$valueVariances[[i]][cnolist$valueVariances[[i]]>0] <- NA
    }

    cnolist$namesSignals = model$namesSpecies[states_index];
    cnolist$namesCues = c(cnolist$namesStimuli,cnolist$namesInhibitors);
    cnolist$valueCues = cbind(cnolist$valueStimuli,cnolist$valueInhibitors);
    cnolist$valueCues = as.matrix(cnolist$valueCues);
    cnolist$valueCues[which(cnolist$valueCues>0)]=1;
    colnames(cnolist$valueCues)=cnolist$namesCues;

    return(CNOlist(cnolist))
}
