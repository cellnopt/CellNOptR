library(CellNOptR)
dir.create("CNOR_analysis")
setwd("CNOR_analysis")
#If loading the data from the MIDAS file and the model from the sif file, type:
#cpfile<-dir(system.file("ToyModel",package="CellNOptR"),full=TRUE)
#file.copy(from=cpfile,to=getwd(),overwrite=TRUE)
#dataToy<-readMIDAS(MIDASfile='ToyDataMMB.csv')
#CNOlistToy<-makeCNOlist(dataset=dataToy,subfield=FALSE)
#ToyModel<-readSIF(sifFile="ToyPKNMMB.sif")
#If loading data and model directly from the package:
data(CNOlistToy,package="CellNOptR")
data(ToyModel,package="CellNOptR")
CNOlistToy
plotCNOlist(CNOlistToy)
plotCNOlistPDF(CNOlist=CNOlistToy,filename="ToyModelGraph.pdf")
checkSignals(CNOlistToy,ToyModel)
indicesToy<-indexFinder(CNOlistToy,ToyModel,verbose=TRUE)
ToyNCNOindices<-findNONC(ToyModel,indicesToy,verbose=TRUE)
ToyNCNOcut<-cutNONC(ToyModel,ToyNCNOindices)
indicesToyNCNOcut<-indexFinder(CNOlistToy,ToyNCNOcut)
ToyNCNOcutComp<-compressModel(ToyNCNOcut,indicesToyNCNOcut)
indicesToyNCNOcutComp<-indexFinder(CNOlistToy,ToyNCNOcutComp)
ToyNCNOcutCompExp<-expandGates(ToyNCNOcutComp)
resECNOlistToy<-residualError(CNOlistToy)
ToyFields4Sim<-prep4sim(ToyNCNOcutCompExp)
initBstring<-rep(1,length(ToyNCNOcutCompExp$reacID))
ToyT1opt<-gaBinaryT1(CNOlist=CNOlistToy,model=ToyNCNOcutCompExp,simList=ToyFields4Sim,indexList=indicesToyNCNOcutComp,initBstring=initBstring,verbose=TRUE)
cutAndPlotResultsT1(model=ToyNCNOcutCompExp,bString=ToyT1opt$bString,simList=ToyFields4Sim,CNOlist=CNOlistToy,indexList=indicesToyNCNOcutComp,plotPDF=TRUE)
plotFit(OptRes=ToyT1opt)
cutAndPlotResultsT1(model=ToyNCNOcutCompExp,bString=ToyT1opt$bString,simList=ToyFields4Sim,CNOlist=CNOlistToy,indexList=indicesToyNCNOcutComp,plotPDF=TRUE)
pdf("evolFitToyT1.pdf")
plotFit(OptRes=ToyT1opt)
dev.off()
writeScaffold(ModelComprExpanded=ToyNCNOcutCompExp,optimResT1=ToyT1opt,optimResT2=NA,ModelOriginal=ToyModel,CNOlist=CNOlistToy)
writeNetwork(ModelOriginal=ToyModel,ModelComprExpanded=ToyNCNOcutCompExp,optimResT1=ToyT1opt,optimResT2=NA,CNOlist=CNOlistToy)
namesFilesToy<-list(dataPlot="ToyModelGraph.pdf",evolFit1="evolFitToyT1.pdf",evolFit2=NA,SimResults1="ToyNCNOcutCompExpSimResultsT1.pdf",SimResults2=NA,Scaffold="ToyNCNOcutCompExpScaffold.sif",ScaffoldDot="ModelModelComprExpandedScaffold.dot",tscaffold="ToyNCNOcutCompExpTimesScaffold.EA",wscaffold="ToyNCNOcutCompExpweightsScaffold.EA",PKN="ToyModelPKN.sif",PKNdot="ToyModelPKN.dot",wPKN="ToyModelTimesPKN.EA",nPKN="ToyModelnodesPKN.NA")
writeReport(ModelOriginal=ToyModel,ModelOpt=ToyNCNOcutCompExp,optimResT1=ToyT1opt,optimResT2=NA,CNOlist=CNOlistToy,directory="testToy",namesFiles=namesFilesToy,namesdata=list(CNOlist="Toy",model="ToyModel"),resE=resECNOlistToy)
################################################
############The one step version################
library(CellNOptR)
dir.create("CNOR_analysis")
setwd("CNOR_analysis")
#If loading the data from the MIDAS file and the model from the sif file, type:
#cpfile<-dir(system.file("ToyModel",package="CellNOptR"),full=TRUE)
#file.copy(from=cpfile,to=getwd(),overwrite=TRUE)
#dataToy<-readMIDAS(MIDASfile='ToyDataMMB.csv')
#CNOlistToy<-makeCNOlist(dataset=dataToy,subfield=FALSE)
#ToyModel<-readSIF(sifFile="ToyPKNMMB.sif")
#If loading data and model directly from the package:
data(CNOlistToy,package="CellNOptR")
data(ToyModel,package="CellNOptR")
CNORwrap(paramsList=NA,Name="Toy",Namesdata=list(CNOlist="ToyData",model="ToyModel"),data=CNOlistToy,model=ToyModel)
#version 2
pList<-list(data=CNOlistToy,model=ToyModel,sizeFac = 1e-04, NAFac = 1, popSize = 50, pMutation = 0.5, maxTime = 60, maxGens = 500, stallGenMax = 100, selPress = 1.2, elitism = 5, RelTol = 0.1,verbose=TRUE)
CNORwrap(paramsList=pList,Name="Toy",Namesdata=list(CNOlist="ToyData",model="ToyModel"),data=NA,model=NA)
######################################################
############The DREAM data and network################
################################################
############The 2 time points################
library(CellNOptR)
dir.create("CNOR_analysis")
setwd("CNOR_analysis")
######
cpfile<-dir(system.file("ToyModel",package="CellNOptR"),full=TRUE)
file.copy(from=cpfile,to=getwd(),overwrite=TRUE)
dataToy<-readMIDAS(MIDASfile='ToyDataMMB.csv')
CNOlistToy<-makeCNOlist(dataset=dataToy,subfield=FALSE)
CNOlistToy
#Transform data for multiple time points
CNOlistToy2<-CNOlistToy
CNOlistToy2$valueSignals[[3]]<-CNOlistToy2$valueSignals[[2]]
CNOlistToy2$valueSignals[[3]][,6:7]<-0
CNOlistToy2$valueSignals[[2]][which(CNOlistToy2$valueSignals[[2]][,6] > 0),6]<-0.5
CNOlistToy2$valueSignals[[2]][which(CNOlistToy2$valueSignals[[2]][,7] > 0),7]<-0.77118
CNOlistToy2$timeSignals<-c(CNOlistToy2$timeSignals, 100)
#In this model I added a negative fedback between cJun and Jnk (!cJun=Jnk)
#this is the model to use with the data CNOlistToy2
ToyModel2<-readSIF(sifFile="ToyModelMMB2.sif")
#####
data(CNOlistToy2,package="CellNOptR")
data(ToyModel2,package="CellNOptR")
#
plotCNOlist(CNOlistToy2)
plotCNOlistPDF(CNOlist=CNOlistToy2,filename="ToyModelGraphT2.pdf")
checkSignals(CNOlistToy2,ToyModel2)
indicesToy2<-indexFinder(CNOlistToy2,ToyModel2,verbose=TRUE)
ToyNCNOindices2<-findNONC(ToyModel2,indicesToy2,verbose=TRUE)
ToyNCNOcut2<-cutNONC(ToyModel2,ToyNCNOindices2)
indicesToyNCNOcut2<-indexFinder(CNOlistToy2,ToyNCNOcut2)
ToyNCNOcutComp2<-compressModel(ToyNCNOcut2,indicesToyNCNOcut2)
indicesToyNCNOcutComp2<-indexFinder(CNOlistToy2,ToyNCNOcutComp2)
ToyNCNOcutCompExp2<-expandGates(ToyNCNOcutComp2)
#resECNOlistToy2<-residualError(CNOlistToy2)
#ToyFields4Sim2<-prep4sim(ToyNCNOcutCompExp2)
#initBstring2<-rep(1,length(ToyNCNOcutCompExp2$reacID))
#ToyT1opt2<-gaBinaryT1(CNOlist=CNOlistToy2,model=ToyNCNOcutCompExp2,simList=ToyFields4Sim2,indexList=indicesToyNCNOcutComp2,initBstring=initBstring2,MaxTime=18)
#cutAndPlotResultsT1(model=ToyNCNOcutCompExp2,bString=ToyT1opt2$bString,simList=ToyFields4Sim2,CNOlist=CNOlistToy2,indexList=indicesToyNCNOcutComp2,plotPDF=TRUE)
#pdf("evolFitToy2T1.pdf")
#plotFit(OptRes=ToyT1opt2)
#dev.off()
#plotFit(OptRes=ToyT1opt2)
#Optimise T2
#SimToyT12<-simulateT1(CNOlist=CNOlistToy2,model=ToyNCNOcutCompExp2,bStringT1=ToyT1opt2$bString,simList=ToyFields4Sim2,indexList=indicesToyNCNOcutComp2)
#ToyT1opt2T2<-gaBinaryT2(CNOlist=CNOlistToy2,model=ToyNCNOcutCompExp2,simList=ToyFields4Sim2,indexList=indicesToyNCNOcutComp2,bStringT1=ToyT1opt2$bString,SimResT1=SimToyT12,MaxTime=18)
#cutAndPlotResultsT2(model=ToyNCNOcutCompExp2,bStringT1=ToyT1opt2$bString,bStringT2=ToyT1opt2T2$bString,simList=ToyFields4Sim2,CNOlist=CNOlistToy2,indexList=indicesToyNCNOcutComp2,plotPDF=TRUE)
#pdf("evolFitToy2T2.pdf")
#plotFit(OptRes=ToyT1opt2T2)
#dev.off()
#plotFit(OptRes=ToyT1opt2T2)

