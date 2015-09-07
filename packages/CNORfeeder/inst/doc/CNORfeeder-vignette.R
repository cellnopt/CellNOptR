### R code from vignette source 'CNORfeeder-vignette.Rnw'

###################################################
### code chunk number 1: installBio (eval = FALSE)
###################################################
## source("http://bioconductor.org/biocLite.R")
## biocLite(c("RBGL","graph","minet","CellNOptR","igraph","catnet"))


###################################################
### code chunk number 2: installPackage (eval = FALSE)
###################################################
## source("http://bioconductor.org/biocLite.R")
## biocLite("CNORfeeder")


###################################################
### code chunk number 3: installPackage2 (eval = FALSE)
###################################################
## install.packages("path_to_CNORfeeder/CNORfeeder_1.0.0.tar.gz",
##     repos=NULL, type="source")


###################################################
### code chunk number 4: loadLib
###################################################
library(CNORfeeder)


###################################################
### code chunk number 5: getData
###################################################
# load the data already formatted as CNOlist
data(CNOlistDREAM,package="CellNOptR")
# load the model (PKN) already in the CNO format
data(DreamModel,package="CellNOptR")


###################################################
### code chunk number 6: getData
###################################################
BTable <- makeBTables(CNOlist=CNOlistDREAM, k=2, measErr=c(0.1, 0))


###################################################
### code chunk number 7: linkRank
###################################################
Lrank <- linksRanking(CNOlist=CNOlistDREAM, measErr=c(0.1, 0), savefile=FALSE)


###################################################
### code chunk number 8: getData
###################################################
model<-preprocessing(data=CNOlistDREAM, model=DreamModel)


###################################################
### code chunk number 9: integration
###################################################
modelIntegr <- mapBTables2model(BTable=BTable,model=model,allInter=TRUE)


###################################################
### code chunk number 10: integLinks
###################################################
modelIntegr$reacID[modelIntegr$indexIntegr]


###################################################
### code chunk number 11: plotData
###################################################
plotModel(model=modelIntegr, CNOlist=CNOlistDREAM, indexIntegr=modelIntegr$indexIntegr)


###################################################
### code chunk number 12: weight
###################################################
modelIntegrWeight <- weighting(modelIntegr=modelIntegr, PKNmodel=DreamModel,
                               CNOlist=CNOlistDREAM, integrFac=10)


###################################################
### code chunk number 13: weightPPI (eval = FALSE)
###################################################
## data(PPINigraph,package="CNORfeeder")
## data(UniprotIDdream,package="CNORfeeder")
## modelIntegrWeight <- weighting(modelIntegr=modelIntegr, PKNmodel=DreamModel,
##                                CNOlist=CNOlistDREAM, integrFac=10,
##                                UniprotID=UniprotIDdream, PPI=PPINigraph)


###################################################
### code chunk number 14: train
###################################################
initBstring<-rep(1,length(modelIntegrWeight$reacID))
# training to data using genetic algorithm (run longer to obtain better results)
DreamT1opt<-gaBinaryT1W(CNOlist=CNOlistDREAM, model=modelIntegrWeight,
                        initBstring=initBstring, maxGens=2, popSize=5, verbose=FALSE)


###################################################
### code chunk number 15: results
###################################################
# model
plotModel(model=modelIntegrWeight, CNOlist=CNOlistDREAM, bString=DreamT1opt$bString)
# data
cutAndPlotResultsT1(model=modelIntegrWeight, CNOlist=CNOlistDREAM,
                    bString=DreamT1opt$bString)


