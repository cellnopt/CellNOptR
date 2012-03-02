plotOptimResultsTimeScale <- function(SimResults=SimResults, yInterpol=yInterpol, xCoords=xCoords, expResults=expResults, times=times, namesCues=namesCues, namesSignals=namesSignals, valueCues=valueCues) {

	# set graphical parameters	
	par(mfrow=c(nr=dim(SimResults)[1]+1, nc=dim(SimResults)[2]+1), cex=0.5, pch=20, mar=c(0.5,0.5,0,0), oma=c(3,2,2,2))
	yMax <- max(max(unlist(lapply(expResults, function(x) max(x, na.rm=TRUE)))), 1)
	yMin <- min(min(unlist(lapply(expResults, function(x) min(x, na.rm=TRUE)))), 0)
	xVal <- times
	xValS <- xCoords
	xValMax = max(c(xVal, xValS))
	allDiff = matrix(NA, nrow=dim(SimResults)[1], ncol=dim(SimResults)[2])
	for(a in 1:dim(SimResults)[1]) {
		for(b in 1:dim(SimResults)[2]) {
			allDiff[a,b] = sum((SimResults[a,b,]-yInterpol[a,b,])^2)
		}	
	}
	diffMax = max(unlist(allDiff))

	# write the boxes on top with the labels of the signals
	for(c in 1:dim(expResults[[1]])[2]){
		par(fg="blue", mar=c(0.5,0.5,0.5,0))
		plot(x=xVal, y=rep(-5, length(xVal)), ylim=c(yMin, yMax), xlab=NA, ylab=NA, xaxt="n", yaxt="n")
		text(labels=as.character(namesSignals[c]), x=((xVal[length(xVal)]-xVal[1])/2), y=(yMin+((yMax-yMin)/2)), cex=2)
	}
	
	plot(x=xVal, y=rep(-5,length(xVal)), ylim=c(yMin, yMax), xlab=NA, ylab=NA, xaxt="n", yaxt="n")
	text(labels="Cues", x=((xVal[length(xVal)]-xVal[1])/2), y=(yMin+((yMax-yMin)/2)), cex=2)
	par(fg="black", mar=c(0.5,0.5,0,0))	

	# go through each elements of the results matrix (i.e. one plot for each signal (column) for each condition (row)	
	for(r in 1:dim(expResults[[1]])[1]) {
		for(c in 1:dim(expResults[[1]])[2]) {

			# determine the simulated and experimental values 
			yVal <-  as.numeric(lapply(expResults, function(x) {x[r,c]}))
			yValS <- SimResults[r,c,]
			yValI <- yInterpol[r,c,]
			
			diff = 1 - (allDiff[r,c] / diffMax)
			bgcolor = rgb(diff,0,0,alpha=0.5)
	
			# plot
			# flag for showing interpolated data
			flagInter = 0
			if(r == dim(expResults[[1]])[1] && c == 1){
				plot(x=xVal, y=yVal, xlim=c(0,xValMax), ylim=c(yMin, yMax), xlab=NA, ylab=NA, xaxt="n", yaxt="n")
				rect(par("usr")[1], par("usr")[3], par("usr")[2], par("usr")[4], col=bgcolor)
				points(x=xVal, y=yVal, xlim=c(0,xValMax), ylim=c(yMin, yMax), xlab=NA, ylab=NA, xaxt="n", yaxt="n")
				lines(x=xVal, y=yVal, xlim=c(0,xValMax), ylim=c(yMin, yMax), xlab=NA, ylab=NA, xaxt="n", yaxt="n")
				lines(x=xValS,y=yValS, xlim=c(0,xValMax), ylim=c(yMin, yMax), xlab=NA, ylab=NA, xaxt="n", yaxt="n", col="blue", lty=2)
				points(x=xValS ,y=yValS, xlim=c(0,xValMax), ylim=c(yMin, yMax), xlab=NA, ylab=NA, xaxt="n", yaxt="n", col="blue")
				if(flagInter==1) {
					lines(x=xValS,y=yValI, xlim=c(0,xValMax), ylim=c(yMin, yMax), xlab=NA, ylab=NA, xaxt="n", yaxt="n", col="red", lty=3)
					points(x=xValS ,y=yValI, xlim=c(0,xValMax), ylim=c(yMin, yMax), xlab=NA, ylab=NA, xaxt="n", yaxt="n", col="red")
				}
				axis(side=1,at=xVal)
				axis(side=2,at=c(-0.5,0,0.5))
			}	
	
			if(r == dim(expResults[[1]])[1] && c != 1){
				plot(x=xVal, y=yVal, xlim=c(0,xValMax), ylim=c(yMin, yMax), xlab=NA, ylab=NA, xaxt="n", yaxt="n")
				rect(par("usr")[1], par("usr")[3], par("usr")[2], par("usr")[4], col=bgcolor)
				points(x=xVal, y=yVal, xlim=c(0,xValMax), ylim=c(yMin, yMax), xlab=NA, ylab=NA, xaxt="n", yaxt="n")
				lines(x=xVal, y=yVal, xlim=c(0,xValMax), ylim=c(yMin, yMax), xlab=NA, ylab=NA, xaxt="n", yaxt="n")
				lines(x=xValS,y=yValS, xlim=c(0,xValMax), ylim=c(yMin, yMax), xlab=NA, ylab=NA, xaxt="n", yaxt="n", col="blue", lty=2)
				points(x=xValS ,y=yValS, xlim=c(0,xValMax), ylim=c(yMin, yMax), xlab=NA, ylab=NA, xaxt="n", yaxt="n", col="blue")
				if(flagInter==1) {		
					lines(x=xValS,y=yValI, xlim=c(0,xValMax), ylim=c(yMin, yMax), xlab=NA, ylab=NA, xaxt="n", yaxt="n", col="red", lty=3)
					points(x=xValS ,y=yValI, xlim=c(0,xValMax), ylim=c(yMin, yMax), xlab=NA, ylab=NA, xaxt="n", yaxt="n", col="red")
				}
				axis(side=1, at=xVal)
			}
				
			if(c == 1 && r == 1){
				plot(x=xVal, y=yVal, xlim=c(0,xValMax), ylim=c(yMin, yMax), xlab=NA, ylab=NA, xaxt="n", yaxt="n")
				rect(par("usr")[1], par("usr")[3], par("usr")[2], par("usr")[4], col=bgcolor)
				points(x=xVal, y=yVal, xlim=c(0,xValMax), ylim=c(yMin, yMax), xlab=NA, ylab=NA, xaxt="n", yaxt="n")
				lines(x=xVal, y=yVal, xlim=c(0,xValMax), ylim=c(yMin, yMax), xlab=NA, ylab=NA, xaxt="n", yaxt="n")
				lines(x=xValS,y=yValS, xlim=c(0,xValMax), ylim=c(yMin, yMax), xlab=NA, ylab=NA, xaxt="n", yaxt="n", col="blue", lty=2)
				points(x=xValS ,y=yValS, xlim=c(0,xValMax), ylim=c(yMin, yMax), xlab=NA, ylab=NA, xaxt="n", yaxt="n", col="blue")
				if(flagInter==1) {		
					lines(x=xValS,y=yValI, xlim=c(0,xValMax), ylim=c(yMin, yMax), xlab=NA, ylab=NA, xaxt="n", yaxt="n", col="red", lty=3)
					points(x=xValS ,y=yValI, xlim=c(0,xValMax), ylim=c(yMin, yMax), xlab=NA, ylab=NA, xaxt="n", yaxt="n", col="red")
				}
				axis(side=2, at=c(-0.5,0,0.5))
			}
	
			if(c == 1 && r != 1 && r != dim(expResults[[1]])[1]){
				plot(x=xVal, y=yVal, xlim=c(0,xValMax), ylim=c(yMin, yMax), xlab=NA, ylab=NA, xaxt="n", yaxt="n")
				rect(par("usr")[1], par("usr")[3], par("usr")[2], par("usr")[4], col=bgcolor)
				points(x=xVal, y=yVal, xlim=c(0,xValMax), ylim=c(yMin, yMax), xlab=NA, ylab=NA, xaxt="n", yaxt="n")
				lines(x=xVal, y=yVal, xlim=c(0,xValMax), ylim=c(yMin, yMax), xlab=NA, ylab=NA, xaxt="n", yaxt="n")
				lines(x=xValS,y=yValS, xlim=c(0,xValMax), ylim=c(yMin, yMax), xlab=NA, ylab=NA, xaxt="n", yaxt="n", col="blue", lty=2)
				points(x=xValS ,y=yValS, xlim=c(0,xValMax), ylim=c(yMin, yMax), xlab=NA, ylab=NA, xaxt="n", yaxt="n", col="blue")
				if(flagInter==1) {		
					lines(x=xValS,y=yValI, xlim=c(0,xValMax), ylim=c(yMin, yMax), xlab=NA, ylab=NA, xaxt="n", yaxt="n", col="red", lty=3)
					points(x=xValS ,y=yValI, xlim=c(0,xValMax), ylim=c(yMin, yMax), xlab=NA, ylab=NA, xaxt="n", yaxt="n", col="red")
				}
				axis(side=2, at=c(-0.5,0,0.5))
			}
				
			if(r == 1 && c!= 1){
				plot(x=xVal, y=yVal, xlim=c(0,xValMax), ylim=c(yMin, yMax), xlab=NA, ylab=NA, xaxt="n", yaxt="n")
				rect(par("usr")[1], par("usr")[3], par("usr")[2], par("usr")[4], col=bgcolor)
				points(x=xVal, y=yVal, xlim=c(0,xValMax), ylim=c(yMin, yMax), xlab=NA, ylab=NA, xaxt="n", yaxt="n")
				lines(x=xVal, y=yVal, xlim=c(0,xValMax), ylim=c(yMin, yMax), xlab=NA, ylab=NA, xaxt="n", yaxt="n")
				lines(x=xValS,y=yValS, xlim=c(0,xValMax), ylim=c(yMin, yMax), xlab=NA, ylab=NA, xaxt="n", yaxt="n", col="blue", lty=2)
				points(x=xValS ,y=yValS, xlim=c(0,xValMax), ylim=c(yMin, yMax), xlab=NA, ylab=NA, xaxt="n", yaxt="n", col="blue")
				if(flagInter==1) {		
					lines(x=xValS,y=yValI, xlim=c(0,xValMax), ylim=c(yMin, yMax), xlab=NA, ylab=NA, xaxt="n", yaxt="n", col="red", lty=3)
					points(x=xValS ,y=yValI, xlim=c(0,xValMax), ylim=c(yMin, yMax), xlab=NA, ylab=NA, xaxt="n", yaxt="n", col="red")
				}
			}
			
			if(c != 1 && r != 1 && r!=dim(expResults[[1]])[1] && c != dim(expResults[[1]])[1]+1){
				plot(x=xVal, y=yVal, xlim=c(0,xValMax), ylim=c(yMin, yMax), xlab=NA, ylab=NA, xaxt="n", yaxt="n")
				rect(par("usr")[1], par("usr")[3], par("usr")[2], par("usr")[4], col=bgcolor)
				points(x=xVal, y=yVal, xlim=c(0,xValMax), ylim=c(yMin, yMax), xlab=NA, ylab=NA, xaxt="n", yaxt="n")
				lines(x=xVal, y=yVal, xlim=c(0,xValMax), ylim=c(yMin, yMax), xlab=NA, ylab=NA, xaxt="n", yaxt="n")
				lines(x=xValS,y=yValS, xlim=c(0,xValMax), ylim=c(yMin, yMax), xlab=NA, ylab=NA, xaxt="n", yaxt="n", col="blue", lty=2)
				points(x=xValS ,y=yValS, xlim=c(0,xValMax), ylim=c(yMin, yMax), xlab=NA, ylab=NA, xaxt="n", yaxt="n", col="blue")
				if(flagInter==1) {		
					lines(x=xValS,y=yValI, xlim=c(0,xValMax), ylim=c(yMin, yMax), xlab=NA, ylab=NA, xaxt="n", yaxt="n", col="red", lty=3)
					points(x=xValS ,y=yValI, xlim=c(0,xValMax), ylim=c(yMin, yMax), xlab=NA, ylab=NA, xaxt="n", yaxt="n", col="red")
				}
			}
		}
	
		# plot the image plots that tell presence/absence of cues
		image(t(matrix(as.numeric(valueCues[r,]), nrow=1)), col=c("white","black"), xaxt="n", yaxt="n")
		if(r == dim(expResults[[1]])[1]) {
			axis(side=1, at=seq(from=0, to=1, length.out=length(namesCues)), labels=as.character(namesCues), las=3)
		}
	}	
}
