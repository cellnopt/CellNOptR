# authors: TC
# This example provide an example and function to call a C function
# from CNORinterface shared library.


# You can test using :
# R --no-restore --no-save  < example.R
#
# or within a R shell by copy/paste the code that follows


# First read the MIDAS and SIF file to get some data
library(CellNOptR)
library(CNORode)
#m = readMIDAS('initialData.csv')
#cnolist = makeCNOlist(m, subfield=FALSE)

load("CNOlistToyFB.RData")
cnolist=CNOlistToyFB


model = readSIF('ToyModelFeedbackDataGenerator.sif')
#indices <- indexFinder(cnolist, model)

res = preprocessing(cnolist, model)

#indices<- findNONC(model, indices, verbose = TRUE)
#model <- cutNONC(model, indices)
#indices <- indexFinder(cnolist, model)
#model <- compressModel(model, indices)
#indices<- indexFinder(cnolist, model)


ode_parameters = createLBodeContPars(res$model,default_n=3,random=TRUE)

# by default, use GA algorithm with default parameters.
# to overwrite default parameters, use params=defaultParametersGA()
paramsGA = defaultParametersGA()
paramsGA.transfert_functions = 3
ode_parameters = parEstimationLBode(cnolist,res$model,
    ode_parameters=ode_parameters,indices=res$indices, paramsGA=paramsGA)



dev.new();
res = plotLBodeFitness(cnolist,res$model,ode_parameters,res$indices)
 
