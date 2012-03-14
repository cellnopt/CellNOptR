#
#  This file is part of the CNO software
#
#  Copyright (c) 2011-2012 - EBI
#
#  File author(s): CNO developers (cno-dev@ebi.ac.uk)
#
#  Distributed under the GPLv2 License.
#  See accompanying file LICENSE.txt or copy at
#      http://www.gnu.org/licenses/gpl-2.0.html
#
#  CNO website: http://www.ebi.ac.uk/saezrodriguez/software.html
#
##############################################################################
# $Id: $

plotOptimResultsTimeScale <- function(SimResults=SimResults, yInterpol=yInterpol, xCoords=xCoords, CNOlist=CNOlist, nsplit=1) {

	# check that CNOlist is a CNOlist
	if(!is.list(CNOlist)) {
		stop("This function expects as input a CNOlist as output by makeCNOlist or normaliseCNOlist")
	}
	
	if(all(names(CNOlist) != c(
			"namesCues",
			"namesStimuli",
			"namesInhibitors",
			"namesSignals",
			"timeSignals",
			"valueCues",
			"valueInhibitors",
			"valueStimuli",
			"valueSignals"))) {
	stop("This function expects as input a CNOlist as output by makeCNOlist")
	}
	
	splits <- dim(CNOlist$valueCues)[1]/nsplit	
	splits <- floor(splits)
	CNOlistOriginal <- CNOlist
	SimResultsOriginal <- SimResults
	yInterpolOriginal <- yInterpol
	
	for(i in 1:nsplit) {
		
		if(nsplit > 1) {
			CNOlist <- CNOlistOriginal
			if(i == nsplit) {
				indices <- (indices[length(indices)]+1):dim(CNOlistOriginal$valueCues)[1]
			}
			
			indices <- ((1:splits)+((i-1)*splits))
			CNOlist$valueCues <- CNOlist$valueCues[indices,]
			CNOlist$valueStimuli <- CNOlist$valueStimuli[indices,]
			CNOlist$valueInhibitors <- CNOlist$valueInhibitors[indices,]
		
			for(n in 1:length(CNOlist$valueSignals)) {
				CNOlist$valueSignals[[n]] <- CNOlist$valueSignals[[n]][indices,]
			}
			for(n in 1:length(SimResults)) {
				SimResults[[n]] <- SimResults[[n]][indices,]
			}
			for(n in 1:length(yInterpol)) {
				yInterpol[[n]] <- yInterpol[[n]][indices,]
			}
		}
		
		# set graphical parameters
		par(
			mfrow=c(nr=dim(CNOlist$valueSignals[[1]])[1]+1,nc=dim(CNOlist$valueSignals[[1]])[2]+2),
			cex=1,
			pch=2,
			mar=c(0.5,0.5,0,0),
			oma=c(3,2,2,2), family="Times", mgp=c(3,0.9,0))
		heatCols = heat.colors(1000)
		
		yMax<-max(unlist(lapply(CNOlist$valueSignals,function(x) max(x, na.rm=TRUE))))
		yMin<-min(unlist(lapply(CNOlist$valueSignals,function(x) min(x, na.rm=TRUE))))
		xVal<-CNOlist$timeSignals
		xValS <- xCoords
		xValMax = max(c(xVal, xValS))
		allDiff = matrix(NA, nrow=dim(SimResults)[1], ncol=dim(SimResults)[2])
		for(a in 1:dim(SimResults)[1]) {
			for(b in 1:dim(SimResults)[2]) {
				allDiff[a,b] = sum((SimResults[a,b,]-yInterpol[a,b,])^2)
			}	
		}
		
		diffMax = max(unlist(allDiff))
		
		for(c in 1:dim(CNOlist$valueSignals[[1]])[2]) {
			par(fg="blue",mar=c(0.5,0.5,0.7,0))
			plot(x=xVal, y=rep(-5,length(xVal)), ylim=c(yMin, yMax),xlab=NA,ylab=NA,xaxt="n",yaxt="n")
			text(
				labels=CNOlist$namesSignals[c],
				x = ((xVal[length(xVal)]-xVal[1])/2),
				y = (yMin+((yMax-yMin)/2)),
				cex = 2)
		}
		plot(
			x = xVal, 
			y = rep(-5,length(xVal)), 
			ylim = c(yMin, yMax),
			xlab = NA,ylab=NA,xaxt="n",yaxt="n")
		text(
			labels = "stim",
			x = ((xVal[length(xVal)]-xVal[1])/2),
			y = (yMin+((yMax-yMin)/2)),cex=2.5)
		plot(
			x = xVal, y=rep(-5,length(xVal)), 
			ylim = c(yMin, yMax),
			xlab = NA,ylab=NA,xaxt="n",yaxt="n")
		text(
			labels="inh",
			x=((xVal[length(xVal)]-xVal[1])/2),
			y=(yMin+((yMax-yMin)/2)),cex=2)
		par(fg="black",mar=c(0.5,0.5,0,0))
					
		for(r in 1:dim(CNOlist$valueSignals[[1]])[1]) {
			for(c in 1:dim(CNOlist$valueSignals[[1]])[2]) {
					
				yVal <- lapply(CNOlist$valueSignals, function(x) {x[r,c]})
				yValS <- SimResults[r,c,]
				yValI <- yInterpol[r,c,]
				diff = (1 - (allDiff[r,c] / diffMax)) * 1000
				if(diff<1) {diff=1}
				bgcolor = heatCols[diff]
				flagInter = 0

				plot(x=xVal,y=yVal,ylim=c(yMin, yMax),xlab=NA,ylab=NA,xaxt="n",yaxt="n")
				rect(par("usr")[1], par("usr")[3], par("usr")[2], par("usr")[4], col=bgcolor)
				points(x=xVal,y=yVal,ylim=c(yMin, yMax),xlab=NA,ylab=NA,xaxt="n",yaxt="n",cex=1)
				
				lines(x=xValS, y=yValS, xlim=c(0,xValMax), ylim=c(yMin, yMax), xlab=NA, ylab=NA, xaxt="n", yaxt="n", col="blue", lty=2, lwd=3)
				#points(x=xValS, y=yValS, xlim=c(0,xValMax), ylim=c(yMin, yMax), xlab=NA, ylab=NA, xaxt="n", yaxt="n", col="black")
				
				if(flagInter==1) {
					lines(x=xValS,y=yValI, xlim=c(0,xValMax), ylim=c(yMin, yMax), xlab=NA, ylab=NA, xaxt="n", yaxt="n", col="red", lty=3)
					points(x=xValS ,y=yValI, xlim=c(0,xValMax), ylim=c(yMin, yMax), xlab=NA, ylab=NA, xaxt="n", yaxt="n", col="red")
				}

					
				# add the axis annotations: if we're on the last row, add the x axis
				if(r == dim(CNOlist$valueSignals[[1]])[1]) {
					axis(side=1,at=CNOlist$timeSignals,cex.axis=1)
				}	
					
				# add the axis annotations: if we're on the first column, add the y axis
				if(c == 1)	{
					axis(side=2,at=c(0,0.5,1), labels=c("","0.5",""),las=1, cex.axis=1.5)
				}
			}
			if(r==1) {		
				image(
					t(matrix(1-CNOlist$valueStimuli[r,],nrow=1)),
					col=c("white","white"),xaxt="n",yaxt="n"
				)
			} else {				
			image(
				t(matrix(1-CNOlist$valueStimuli[r,],nrow=1)),
				col=c("black","white"),xaxt="n",yaxt="n")
			}
			if(r == dim(CNOlist$valueSignals[[1]])[1]) {
				axis(
					side=1,
					at=seq(from=0, to=1,length.out=length(CNOlist$namesStimuli)),
					labels=CNOlist$namesStimuli,las=3,cex.axis=1.5)
			}
			
			image(
				t(matrix(CNOlist$valueInhibitors[r,],nrow=1)),
				col=c("white","black"),xaxt="n",yaxt="n")
			
			if(r == dim(CNOlist$valueSignals[[1]])[1]) {
				axis(
					side=1,
					at=seq(from=0, to=1,length.out=length(CNOlist$namesInhibitors)),
					labels=paste(CNOlist$namesInhibitors,"i",sep=""),
					las=3,cex.axis=1.4)
			}
		}
	#	screen(2)
	#	colorlegend(heat.colors(100),LETTERS[1:12], xlim=c(1,2))
	}			
}
		
