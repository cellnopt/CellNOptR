# authors: TC
# This example provide an example and function to call a C function
# from CNORinterface shared library.


# You can test using :
# R --no-restore --no-save  < example.R
#
# or within a R shell by copy/paste the code that follows


# First read the MIDAS and SIF file to get some data


library(CellNOptR)
install.packages("CNORode_1.0.zip",repos=NULL);
setwd("tests/test4");
#setwd("tests/test1");
library("CNORode")

s = readSif('model.sif')
m = readMIDAS('initialData.csv')
cnolist = makeCNOlist(m, subfield=FALSE)

indices <- indexFinder(cnolist, s, verbose = TRUE)
modelNCNOindices <- findNONC(s, indices, verbose = TRUE)
s <- cutNONC(s, modelNCNOindices);

#results=logic_based_ode_parameters_estimation_SSm(cnolist,s)
logic_based_ode_MINLP_SSm(cnolist,s,ndiverse=10,dim_refset=6)



