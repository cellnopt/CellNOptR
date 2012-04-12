# this is an example of the main steps of the integrated CellNOptR - CNORfeeder pipeline

library(CNORfeeder)

# load the data already formatted as CNOlist
data(CNOlistDREAM,package="CellNOptR")
# load the model (PKN) already in the CNO format
data(DreamModel,package="CellNOptR")
# see CellNOptR documentation to import other data/PKNs)


# load the UniprotID of proteins in the PKN
data(UniprotIDdream,package="CNORfeeder")
# load the curated PIN as igraph 
data(PPINigraph,package="CNORfeeder")

# A. COMPRESSION (see Fig 1A) - CellNOptR
# preprocessing step
res<-preprocessing(Data=CNOlistDREAM, Model=DreamModel)
Model<-res$model  #this is the compressed model


# C. INFERENCE (see Fig 1C) - CNORfeeder
# FEED inference: codified in Boolean Tables
BTable <- makeBTables(CNOlist=CNOlistDREAM, k=2, measErr=c(0.1, 0))


# D. INTEGRATION (see Fig 1D) - CNORfeeder
# integration with the compressed model
modelIntegr <- MapBTables2Model(BTable=BTable,Model=Model,allInter=TRUE)
# see example in ?MapDDN2Model to use other reverse-engineering methods


# E. WEGHTING (see Fig 1E) - CNORfeeder
# proritization of links based on the PIN
# the followig step may take a while:
# if not run, all integrated links will have the same weight in the training step
# resPPIweight <- PPIweight(modelIntegr=modelIntegr,PKNmodel=DreamModel,CNOlist=CNOlistDREAM,UniprotID=UniprotIDdream,PPINigraph=PPINigraph)
# modelIntegr<-resPPIweight$modelIntegr


# B. TRAINING (see Fig 1B) - CellNOptR
DreamFields4Sim<-prep4Sim(modelIntegr)
initBstring<-rep(1,length(modelIntegr$reacID))
# training to data using genetic algorithm (run longer to obtain better results)
DreamT1opt<-gaBinaryT1int(
CNOlist=CNOlistDREAM,
	Model=modelIntegr,
	SimList=DreamFields4Sim,
	indexList=res$indices,
	initBstring=initBstring,
	maxGens=2,
	PopSize=5,
	verbose=FALSE)


# plot the model with selected links in geen (if derived form the PKN) and in purple (if integrated)
plotModel(model=modelIntegr, cnolist=CNOlistDREAM, bString=DreamT1opt$bString, indexIntegr=modelIntegr$indexIntegr)