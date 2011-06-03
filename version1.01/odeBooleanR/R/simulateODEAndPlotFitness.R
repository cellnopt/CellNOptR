simulateODEAndPlotFitness <-
function(CNOlist,CNOModel,odePars)
{
  indexes=indexCNO2ODE(CNOlist,odePars);
  
  simData=simulateData(CNOlist,indexes,odePars)
   
  plotODESimFitness(simData,CNOlist); 
}

