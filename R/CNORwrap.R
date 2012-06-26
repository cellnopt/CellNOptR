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
# $Id$
#This function is a wrapper around the whole CNOR analysis, it performs the following steps:
#1.Plot the CNOlist
#2.Checks data to model compatibility
#3.Find the indices, in the model, of the species that are inh/stim/sign
#4.Find the indices of the non-osb/non-contr
#5.Cut the nonc off the model
#6.Recompute the indices
#7.Compress the model
#8.Recompute the indices
#9.Expand the gates
#10.Compute the residual error
#11.Prepare for simulation
#12.Optimisation t1
#13.Plot simulated and experimental results
#14.Plot the evolution of fit
#15.Optimise t2 (not implemented in this version)
#16.Write the scaffold and PKN
#17.Write the report

CNORwrap<-function(paramsList=NA, data=NA, model=NA, name, namesData=NA, time=1)
{

    # aliases
    Name = name
    NamesData = namesData
    Time = time

    # if paramsList empty, we will fill it with
    # default values and the Data and model provided
	if(is.na(paramsList[1])==TRUE){
        #print("paramList not provided")
        # Data must be provided
        # is.na raise warning, so let us use is.list
        if (is.na(data[1])==TRUE){
            stop("if paramsList not provided, Data must be provided")
        }
        # model must be provided
        if (is.na(model[1])==TRUE){
            stop("if paramsList not provided, model must be provided")
        }
        paramsList = defaultParameters(data, model)
    }

    # if paramsList is provided, and model and Data are NA, then we must find
    # Data and model in paramsList
    else if (is.na(paramsList[1])==FALSE){
        #print("paramsList  provided")

        if (is.na(data[1])==TRUE && is.na(paramsList$data[1])==TRUE){
            stop("Data must be provided either in paramsList or Data argument")
        }
        if (is.na(model[1])==TRUE && is.na(paramsList$model[1])==TRUE){
            stop("model must be provided either in paramsList or model argument")
        }

        # Data must be provided
        if (is.na(data[1])==FALSE){
            #print("Data provided")
            if (is.na(paramsList$data[1])==FALSE){
                warning("overwritting paramsList$data with provided Data")
            }
            paramsList$data = data
        }
        # model must be provided
        if (is.na(model[1])==FALSE){
            #print("model provided")
            if (is.na(paramsList$model[1])==FALSE){
                warning("overwritting paramsList$Model with provided Model")
            }
            paramsList$model = model
        }

    }

    if (is.list(NamesData)==FALSE){
        if (is.na(NamesData)==TRUE){
            NamesData <- list(
                CNOlist=paste(Name, "Data",sep=""),
               model=paste(Name,"Model", sep=""))
        }
        else{
            stop("NamesData must be a list or kept to the default values (NA)")
        }
    }

    return(paramsList)

    #1.Plot the CNOlist
	plotCNOlist(paramsList$data)
	plotCNOlistPDF(
		CNOlist=paramsList$data,
		filename=paste(Name,"DataPlot.pdf",sep="")
		)

    #2. Checks data to model compatibility
	checkSignals(CNOlist=paramsList$data,model=paramsList$model)

    #3. Find the indices, in the model, of the species that are inh/stim/sign
	Indices<-indexFinder(
		CNOlist=paramsList$data,
		model=paramsList$model,
		verbose=paramsList$verbose)

    #4. Find the indices of the non-osb/non-contr
	NCNOindices<-findNONC(
		model=paramsList$model,
		indexes=Indices,
		verbose=paramsList$verbose)

    #5.Cut the nonc off the model
	NCNOcut<-cutNONC(model=paramsList$model, NONCindexes=NCNOindices)

    #6.Recompute the indices
	IndicesNCNOcut<-indexFinder(CNOlist=paramsList$data,model=NCNOcut)

    #7.Compress the model
	NCNOcutComp<-compressModel(model=NCNOcut,indexes=IndicesNCNOcut)

    #8.Recompute the indices
	IndicesNCNOcutComp<-indexFinder(CNOlist=paramsList$data,model=NCNOcutComp)

    #9.Expand the gates
	NCNOcutCompExp<-expandGates(model=NCNOcutComp)

    #10.Compute the residual error
	resE<-residualError(CNOlist=paramsList$data)

    #11.Prepare for simulation
	fields4Sim<-prep4sim(model=NCNOcutCompExp)

    #12.Optimisation t1
	initBstring<-rep(1,length(NCNOcutCompExp$reacID))
	T1opt<-gaBinaryT1(CNOlist=paramsList$data,
		model=NCNOcutCompExp,
		simList=fields4Sim,
		indexList=IndicesNCNOcutComp,
		initBstring=initBstring,
		sizeFac=paramsList$sizeFac,
		NAFac=paramsList$NAFac,
		popSize=paramsList$popSize,
		pMutation=paramsList$pMutation,
		maxTime=paramsList$maxTime,
		maxGens=paramsList$maxGens,
		stallGenMax=paramsList$stallGenMax,
		selPress=paramsList$selPress,
		elitism=paramsList$elitism,
		relTol=paramsList$relTol,
		verbose=paramsList$verbose)

#13.Plot simulated and experimental results
	cutAndPlotResultsT1(
		model=NCNOcutCompExp,
		bString=T1opt$bString,
		simList=fields4Sim,
		CNOlist=paramsList$data,
		indexList=IndicesNCNOcutComp,
		plotPDF=TRUE)

#14.Plot the evolution of fit
	pdf(paste(Name,"evolFitT1.pdf",sep=""))
	plotFit(optRes=T1opt)
	dev.off()
	plotFit(optRes=T1opt)

#15.Optimise t2
	if(Time==2){

		SimT1<-simulateT1(
			CNOlist=paramsList$data,
			model=NCNOcutCompExp,
			bStringT1=T1opt$bString,
			simList=fields4Sim,
			indexList=IndicesNCNOcutComp)

		T2opt<-gaBinaryT2(
			CNOlist=paramsList$data,
			model=NCNOcutCompExp,
			simList=fields4Sim,
			indexList=IndicesNCNOcutComp,
			bStringT1=T1opt$bString,
			simResT1=SimT1,
			sizeFac=paramsList$sizeFac,
			NAFac=paramsList$NAFac,
			popSize=paramsList$popSize,
			pMutation=paramsList$pMutation,
			maxTime=paramsList$maxTime,
			maxGens=paramsList$maxGens,
			stallGenMax=paramsList$stallGenMax,
			selPress=paramsList$selPress,
			elitism=paramsList$elitism,
			relTol=paramsList$relTol,
			verbose=paramsList$verbose)

		cutAndPlotResultsT2(model=NCNOcutCompExp,bStringT1=T1opt$bString,
            bStringT2=T2opt$bString,simList=fields4Sim,CNOlist=paramsList$data,
            indexList=IndicesNCNOcutComp,plotPDF=TRUE)

		pdf(paste(Name,"evolFitT2.pdf",sep=""))
		plotFit(optRes=T2opt)
		dev.off()
		plotFit(optRes=T2opt)

		}else{

			T2opt<-NA

			}
#16.Write the scaffold and PKN
#and
#17.Write the report
	if(Time==2){


		writeScaffold(
			modelComprExpanded=NCNOcutCompExp,
			optimResT1=T1opt,
			optimResT2=T2opt,
			modelOriginal=paramsList$model,
			CNOlist=paramsList$data)
		writeNetwork(
			modelOriginal=paramsList$model,
			modelComprExpanded=NCNOcutCompExp,
			optimResT1=T1opt,
			optimResT2=T2opt,
			CNOlist=paramsList$data)
		namesfiles<-list(
			dataPlot=paste(Name,"DataPlot.pdf",sep=""),
			evolFitT1=paste(Name,"evolFitT1.pdf",sep=""),
			evolFitT2=paste(Name,"evolFitT2.pdf",sep=""),
			simResults2="NCNOcutCompExpSimResultsT1T2.pdf",
			simResults1="NCNOcutCompExpSimResultsT1.pdf",
			scaffold="Scaffold.sif",
			scaffoldDot="Scaffold.dot",
			tscaffold="TimesScaffold.EA",
			wscaffold="weightsScaffold.EA",
			PKN="PKN.sif",
			PKNdot="PKN.dot",
			wPKN="TimesPKN.EA",
			nPKN="nodesPKN.NA")
		writeReport(
			modelOriginal=paramsList$model,
			modelOpt=NCNOcutCompExp,
			optimResT1=T1opt,
			optimResT2=T2opt,
			CNOlist=paramsList$data,
			directory=Name,
			namesFiles=namesfiles,
			namesData=NamesData,
			resE=resE)

		}else{

			writeScaffold(
				modelComprExpanded=NCNOcutCompExp,
				optimResT1=T1opt,
				optimResT2=NA,
				modelOriginal=paramsList$model,
				CNOlist=paramsList$data)
			writeNetwork(
				modelOriginal=paramsList$model,
				modelComprExpanded=NCNOcutCompExp,
				optimResT1=T1opt,
				optimResT2=NA,
				CNOlist=paramsList$data)
			namesfiles<-list(
				dataPlot=paste(Name,"DataPlot.pdf",sep=""),
				evolFitT1=paste(Name,"evolFitT1.pdf",sep=""),
				evolFitT2=NA,simResults2=NA,
				simResults1="NCNOcutCompExpSimResultsT1.pdf",
				scaffold="Scaffold.sif",
				scaffoldDot="Scaffold.dot",
				tscaffold="TimesScaffold.EA",
				wscaffold="weightsScaffold.EA",
				PKN="PKN.sif",
				PKNdot="PKN.dot",
				wPKN="TimesPKN.EA",
				nPKN="nodesPKN.NA")
			writeReport(
				modelOriginal=paramsList$model,
				modelOpt=NCNOcutCompExp,
				optimResT1=T1opt,
				optimResT2=NA,
				CNOlist=paramsList$data,
				directory=Name,
				namesFiles=namesfiles,
				namesData=NamesData,
				resE=resE)

			}
}
