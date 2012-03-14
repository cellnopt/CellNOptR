# authors: TC
# This example provide an example and function to call a C function
# from CNORinterface shared library.


# You can test using :
# R --no-restore --no-save  < example.R
#
# or within a R shell by copy/paste the code that follows


# First read the MIDAS and SIF file to get some data


library(CellNOptR)
#install.packages("CNORode_1.0.zip",repos=NULL);
library("CNORode")

s = readSif('model.sif')
m = readMIDAS('initialData.csv')
cnolist = makeCNOlist(m, subfield=FALSE)

indices <- indexFinder(cnolist, s, verbose = TRUE)


ode_parameters = createLBodeContPars(s, default_n=3);

paramsGA = defaultParametersGA()
paramsGA.transfert_functions = 3
ode_parameters = parEstimationLBode(cnolist,s,method="ga",
    ode_parameters=ode_parameters,indices=indices, paramsGA=paramsGA)

