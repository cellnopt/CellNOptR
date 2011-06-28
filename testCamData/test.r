#Load the data
library(CellNOptR)
setwd('C:/Users/davidh/Desktop/testCamData')
dataToy<-readMIDAS(MIDASfile='LiverDreamFeedback.csv')
CNOlistToy<-makeCNOlist(dataset=dataToy,subfield=FALSE)
plotCNOlistPDF(CNOlist = CNOlistToy, fileName = "ToyModelGraph.pdf")
#Transform data for multiple time points
#CNOlistToy2<-CNOlistToy
#CNOlistToy2$valueSignals[[3]]<-CNOlistToy2$valueSignals[[2]]
#CNOlistToy2$valueSignals[[3]][,6:7]<-0
#CNOlistToy2$valueSignals[[2]][which(CNOlistToy2$valueSignals[[2]][,6] > 0),6]<-0.5
#CNOlistToy2$valueSignals[[2]][which(CNOlistToy2$valueSignals[[2]][,7] > 0),7]<-0.77118
#CNOlistToy2$timeSignals<-c(CNOlistToy2$timeSignals, 100)