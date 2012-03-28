library(CellNOptR)


# path to the local copy of the SVN trunk
pathToSVN="/Users/localadmin/CNO_trunk_svn2"


source(paste(pathToSVN,"/CNOR_add_links/R/AddLink.R",sep=""))
source(paste(pathToSVN,"/CNOR_add_links/R/makeBTables.R",sep=""))
source(paste(pathToSVN,"/CNOR_add_links/R/SearchLinkGraph.R",sep=""))
source(paste(pathToSVN,"/CNOR_add_links/R/downCueGraph.R",sep=""))
source(paste(pathToSVN,"/CNOR_add_links/R/MapBTables2Model.R",sep=""))
source(paste(pathToSVN,"/CNOR_add_links/R/upSignalGraph.R",sep=""))
source(paste(pathToSVN,"/CNOR_add_links/R/AddLinkAND.R",sep=""))
source(paste(pathToSVN,"/CNOR_add_links/R/model2sif.R",sep=""))
source(paste(pathToSVN,"/CNOR_add_links/R/sif2graph.R",sep=""))
source(paste(pathToSVN,"/CNOR_add_links/R/AddLinkAND.R",sep=""))
source(paste(pathToSVN,"/CNOR_add_links/R/getFitint.R",sep=""))
source(paste(pathToSVN,"/CNOR_add_links/R/gaBinaryT1int.R",sep=""))



# CNOlist is the original dataset loaded from the csv file
myData<- readMIDAS(MIDASfile = paste(pathToSVN,"/Model_with_data/sampleModels/LiverDREAM/LiverDREAM.csv", sep=""))
CNOlist <- makeCNOlist(dataset = myData, subfield = FALSE)

# generation of the boolean tables codifying the effects of stimuli and inhibitors on each measured protein
BTable <- makeBTables(CNOlist=CNOlist, k=2, measErr=c(0.1, 0))

# model is the original PKN loaded from the sif file
model=readSif(sifFile = paste(pathToSVN,"/Model_with_data/sampleModels/LiverDREAM/LiverDREAM.sif", sep=""))

# compression and expansion of the model
checkSignals(CNOlist,model)
indices<-indexFinder(CNOlist,model,verbose=TRUE)
NCNOindices<-findNONC(model,indices,verbose=TRUE)
NCNOcut<-cutNONC(model,NCNOindices)
indicesNCNOcut<-indexFinder(CNOlist,NCNOcut)
NCNOcutComp<-compressModel(NCNOcut,indicesNCNOcut)
indicesNCNOcutComp<-indexFinder(CNOlist,NCNOcutComp)
NCNOcutCompExp<-expandGates(NCNOcutComp)

Model=NCNOcutCompExp

# integration of inferred links (using Boolean tables) with the Model
#CEN <- MapBTables2Model(BTable=BTable,Model=NULL,allInter=TRUE)
modelIntegr <- MapBTables2Model(BTable=BTable,Model=Model,allInter=TRUE)

# a field $indexIntegr with the index of the integrated links is added to the model,
# the list of links added beyond the PKN can be obtained as
modelIntegr$reacID[modelIntegr$indexIntegr]

