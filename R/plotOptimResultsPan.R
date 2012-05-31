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
# $Id: plotOptimResultsPan.R 802 2012-03-22 16:44:12Z cokelaer $

plotOptimResultsPan <- function(SimResults=SimResults, yInterpol=NULL, xCoords=NULL, CNOlist=CNOlist, nsplit=1, formalism=c("ss1","ss2","dt","ode"), pdf=FALSE, pdfFileName="", tPt=NULL) {

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
	formalism <- match.arg(formalism)
	# index valueSignals according to tPt
	valueSignalsI = sapply(c(0,tPt), function(x) which(CNOlist$timeSignals==x))
	
	#####	functions	#####
	
	list2Array = function(x, dim) {
		xUnlist = unlist(x);
		xArray=array(xUnlist,dim=dim)
	}
	
	#####	/functions/	#####
	
	oldPar = par(no.readonly=TRUE)
	if(pdf==TRUE) {
		pdf(file=pdfFileName, width=14.5,height=11)	
	}
	split.screen(c(1,dim(CNOlist$valueSignals[[1]])[2]+3))
	for(a in 1:(dim(CNOlist$valueSignals[[1]])[2]+2)) {
		split.screen(c(dim(CNOlist$valueSignals[[1]])[1]+1,1),a)
	}
		
	# TODO - do i need all these with split.screen?
	par(
		pch=2,
		oma=c(3,2,2,2),
		mgp=c(3,0.9,0),
		family="Times"
	)
	
	heatCols = heat.colors(1000)
	
	# maximum across all data points
	#yMax <- max(unlist(lapply(CNOlist$valueSignals, function(x) max(x,na.rm=TRUE))))
	yMax=1
	# minimum across all data points
	#yMin <- min(unlist(lapply(CNOlist$valueSignals, function(x) min(x,na.rm=TRUE))))
	yMin=0
	# time labels
	xVal <- CNOlist$timeSignals[valueSignalsI]
	if(formalism=="dt") {
		xValS = xCoords	
	} else if (formalism == "ss1") {
		xValS = c(0,tPt[1])
	} else if (formalism == "ss2") {
		xValS = c(0,tPt[1:2])
	} else {
		xValS = xVal
	}
	# latest time point
	xValMax = max(xVal)
		
	# make SimResults array if not already	
	if(!is.array(SimResults)) {
		SimResults = list2Array(SimResults, dim=c(dim(SimResults[[1]]),length(SimResults)))	
	}

	# make valueSignals an array
	valueSignalsArr = list2Array(CNOlist$valueSignals,
	dim=c(dim(CNOlist$valueSignals[[1]]),length(CNOlist$valueSignals)))

	# calculate the MSE
	allDiff = matrix(NA, nrow=dim(SimResults)[1], ncol=dim(SimResults)[2])
	if(formalism != "dt") {
		for(a in 1:dim(SimResults)[1]) {
			for(b in 1:dim(SimResults)[2]) {
				allDiff[a,b] = sum((SimResults[a,b,]-valueSignalsArr[a,b,valueSignalsI])^2)
			}
		}	
	} else {
		for(a in 1:dim(SimResults)[1]) {
			for(b in 1:dim(SimResults)[2]) {
				allDiff[a,b] = sum((SimResults[a,b,]-yInterpol[a,b,])^2)
			}
		}		
	}
	# max difference between sim and exper
	diffMax = max(unlist(!is.na(allDiff)))
	
	# set the count for the split screen window
	count1 = dim(CNOlist$valueSignals[[1]])[2]+4
		
	# plot headers
	for(c in 1:dim(CNOlist$valueSignals[[1]])[2]) {
		screen(count1)
		par(fg="blue",mar=c(0.5,0.5,0.7,0))
		plot(x=xVal, y=rep(-5,length(xVal)), ylim=c(yMin, yMax),
		xlab=NA,ylab=NA,xaxt="n",yaxt="n")
			
		text(
			labels=CNOlist$namesSignals[c],
			x=((xVal[length(xVal)]-xVal[1])/2),
			y=(yMin+((yMax-yMin)/2)),
			cex=1.6
		)
		count1 = count1 + dim(CNOlist$valueSignals[[1]])[1]+1
	}
	
	# stim + inhib	
	screen(count1)
	par(fg="blue",mar=c(0.5,0.5,0.7,0))
	plot(
		x = xVal, 
		y = rep(-5,length(xVal)), 
		ylim = c(yMin, yMax),
		xlab = NA,ylab=NA,xaxt="n",yaxt="n"
	)
	text(
		labels = "Stim",
		x = ((xVal[length(xVal)]-xVal[1])/2),
		y = (yMin+((yMax-yMin)/2)),cex=1.6
	)
	
	count1 = count1 + dim(CNOlist$valueSignals[[1]])[1]+1
	screen(count1)
	par(fg="blue",mar=c(0.5,0.5,0.7,0))
	plot(
		x = xVal, y=rep(-5,length(xVal)), 
		ylim = c(yMin, yMax),
		xlab = NA,ylab=NA,xaxt="n",yaxt="n")
	text(
		labels="Inh",
		x=((xVal[length(xVal)]-xVal[1])/2),
		y=(yMin+((yMax-yMin)/2)),cex=1.6
	)
		
	# new count for plotting results
	countRow = dim(CNOlist$valueSignals[[1]])[2]+4
	
	for(c in 1:dim(CNOlist$valueSignals[[1]])[2]) {
		countRow=countRow+1
		for(r in 1:dim(CNOlist$valueSignals[[1]])[1]) {
				
			screen(countRow)
			par(fg="black",mar=c(0.5,0.5,0,0))
			yVal <- lapply(CNOlist$valueSignals[valueSignalsI], function(x) {x[r,c]})
			yValS <- SimResults[r,c,]
			if(!is.na(allDiff[r,c])) {
				diff = (1 - (allDiff[r,c] / diffMax)) * 1000
			} else {
				diff = 0
			}
			if(diff<1) {diff=1}
			bgcolor = heatCols[diff]
			
			plot(x=xVal,y=yVal,ylim=c(yMin,yMax),xlab=NA,ylab=NA,xaxt="n",yaxt="n",)
			rect(par("usr")[1], par("usr")[3], par("usr")[2], par("usr")[4], col=bgcolor)
			points(x=xVal,y=yVal,ylim=c(yMin, yMax),xlab=NA,ylab=NA,xaxt="n",yaxt="n",pch=2)

			lines(xValS, yValS, xlim=c(0,xValMax), ylim=c(yMin, yMax),
			xlab=NA, ylab=NA, xaxt="n", yaxt="n", col="blue", lty=2, lwd=3)
			points(xValS, yValS, xlim=c(0,xValMax), ylim=c(yMin, yMax),
			xlab=NA, ylab=NA, xaxt="n", yaxt="n", col="blue", pch=16)
			
			# add the axis annotations: if we're on the last row, add the x axis
			if(r == dim(CNOlist$valueSignals[[1]])[1]) {
				axis(side=1,at=CNOlist$timeSignals)
			}	
				
			# add the axis annotations: if we're on the first column, add the y axis
			if(c == 1)	{
				axis(side=2,at=c(0,0.5,1), labels=c("","0.5",""), las=1)
			}
						
			countRow=countRow+1		
		}
	}
		
	sStim = countRow+1
		
	for(c1 in 1:dim(CNOlist$valueSignals[[1]])[1]) {
		screen(sStim)
		par(mar=c(0.5,0.5,0,0))
		if(c1==1) {		
			image(
				t(matrix(1-CNOlist$valueStimuli[c1,],nrow=1)),
				col=c("white","white"),xaxt="n",yaxt="n"
			)
		} else {
			image(
				t(matrix(1-CNOlist$valueStimuli[c1,],nrow=1)),
				col=c("black","white"),xaxt="n",yaxt="n"
			)
		}
		if(c1 == dim(CNOlist$valueSignals[[1]])[1]) {
			axis(
				side=1,
				at=seq(from=0, to=1,length.out=length(CNOlist$namesStimuli)),
				labels=CNOlist$namesStimuli,las=3,cex.axis=1.2
			)	
		}
		sStim = sStim+1
	}
	
	sInhib = sStim+1
	for(i1 in 1:dim(CNOlist$valueSignals[[1]])[1]) {
		screen(sInhib)
		par(mar=c(0.5,0.5,0,0))

		image(
			t(matrix(CNOlist$valueInhibitors[i1,],nrow=1)),
			col=c("white","black"),xaxt="n",yaxt="n"
		)
		if(i1 == dim(CNOlist$valueSignals[[1]])[1]) {
			axis(
				side=1,
				at=seq(from=0, to=1,length.out=length(CNOlist$namesInhibitors)),
				labels=paste(CNOlist$namesInhibitors,"i",sep=""),
				las=3,cex.axis=1.2
			)
		}
		sInhib = sInhib+1			
	}

	screen(dim(CNOlist$valueSignals[[1]])[2]+3)
	splitProp = 1/(dim(CNOlist$valueSignals[[1]])[1]+1)
	sSplit = matrix(c(0,1,(1-splitProp),1,0,1,0,(1-splitProp)),
	ncol=4, byrow=T)
	split.screen(sSplit)
	screen(sInhib)
	par(fg="blue",mar=c(0.5,0.5,0.7,0))
	plot(
		x = xVal, 
		y = rep(-5,length(xVal)), 
		ylim = c(yMin, yMax),
		xlab = NA,ylab=NA,xaxt="n",yaxt="n"
	)
	text(
		labels = "Error",
		x = ((xVal[length(xVal)]-xVal[1])/2),
		y = (yMin+((yMax-yMin)/2)),cex=1.6
	)
	
	screen(sInhib+1)	
	colbar = heat.colors(100)
	labels = c(1,0.5,0)
	len <- length(colbar)
	rhs <- 0.6
	rhs2 <- rhs + rhs/10
	at <- c(0, 0.5, 1)
	
	par(mai=c(0,0.2,0,0))
	plot.new()
	yyy <- seq(0,1,length=len+1)
	rect(0, yyy[1:len], rep(rhs, len), yyy[-1],
	col = colbar, border = colbar)
	rect(0, 0, rhs, 1, border="black")
	segments(rhs, at, rhs2, at)
	text(x=rhs2, y=at, labels=labels, pos=4, offset=0.2)
	if(pdf==TRUE) {
		dev.off()
	}
	close.screen(all.screens=TRUE)
	par(oldPar)
}			



	########## TO ADD ##########


#	splits <- dim(CNOlist$valueCues)[1]/nsplit	
#	splits <- floor(splits)
#	CNOlistOriginal <- CNOlist
#	SimResultsOriginal <- SimResults
	
# 	for(i in 1:nsplit) {
# 		if(nsplit > 1) {
# 			CNOlist <- CNOlistOriginal
# 		if(i == nsplit) {
#			indices <- (indices[length(indices)]+1):dim(CNOlistOriginal$valueCues)[1]
# 		}

# 		indices <- ((1:splits)+((i-1)*splits))
# 		CNOlist$valueCues <- CNOlist$valueCues[indices,]
# 		CNOlist$valueStimuli <- CNOlist$valueStimuli[indices,]
# 		CNOlist$valueInhibitors <- CNOlist$valueInhibitors[indices,]

# 		for(n in 1:length(CNOlist$valueSignals)) {
# 			CNOlist$valueSignals[[n]] <- CNOlist$valueSignals[[n]][indices,]
# 		}
#		for(n in 1:length(SimResults)) {
#			SimResults[[n]] <- SimResults[[n]][indices,]
# 		}
# 	}

