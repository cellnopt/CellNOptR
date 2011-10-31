plotCNOlistPDF <-
function(CNOlist,filename){
	pdf(file=filename,width=14,height=7)
	plotCNOlist(CNOlist)
	dev.off()
	}

