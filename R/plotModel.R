plotModel <- function(model, optimRes=NA, indexIntegr=NA, signals=NULL, stimulis=NULL, inhibitors=NULL,
    ncno=NULL, compressed=NULL, cnolist=NULL, PDF=FALSE, filename=NULL, SVG=FALSE, PNG=FALSE){
#  "$Id$"

# Purpose: Plot a model at various stages (e.g., before running the simulation.)
# Argument: this function takes as an input either an object as returned by readSif
#  function or the filename of a sif file.
# Optional arguments: In order to color the node, you must provide optional
#  arguments. They can be a CNOList, or list of signals, stimulis, inhibitors. 
#  compressed arguments can also be provided. Finally, optional arguments allows
#  to save the plot into either a PDF or SVG format.
# authors: tc, 2011
# example: 
#  > plotModel2("ToyModelPKN.sif")
#  > s = readSif("ToyModelPKN.sif"); plotModel2(s)
  
# g<-plotModel(model="Model.sif",optimRes=optimRes, indexIntegr=indexIntegr, cnolist=CNOlist)
   

# Example2:
# filename = "Model_with_data/sampleModels/Toy_PKN_MMB/ToyPKNMMB.sif"
# g = plotModel(model, cnolist=cnolist, filename="myPKN4.dot", optimRes=T1opt)

    library(Rgraphviz)
    library(RBGL)
    
    
    
    if (typeof(model) == "character"){
        print("Reading the SIF file and parsing")
        raw = read.table(model)  # read the PKN data
        # build the unique vertices from the column 1 and 3 of the SIF file
        vertices = unique(c(as.character(raw$V1), as.character(raw$V3)))
        # some aliases
        v1 = raw$V1
        v2 = raw$V3
        edges = raw$V2
        BStimes<-rep(1,length(edges))
        Integr<-rep(0,length(edges))
    }
    else if (typeof(model)=="list" && any("namesSpecies" == names(model))){ 
        # is it the output of readSif ?
        # build the unique vertices from the nameSpecies
        vertices = model$namesSpecies
        mysplit = function(x){strsplit(x, "=")}
        reacs = mysplit(model$reacID) # separate the reactions into left and right parts
        tmp <- unlist(mysplit(model$reacID))
        reacs = t(matrix(unlist(mysplit(model$reacID)), ncol=length(tmp)/2)) # reordering
        
        if (is.na(optimRes[1])){
            optimBStimes<-rep(1,dim(reacs)[1])
        }else{
            optimBStimes<-optimRes$bString
        }
    
        optIntegr<-rep(0,length(optimBStimes))
        if (!is.na(indexIntegr[1])){
            optIntegr[indexIntegr]<-1
        }
    
        BStimes<-vector()
        Integr<-vector()
        
        
        CountReac<-1
        CountAnds<-1
        mysplit2 = function(x){strsplit(x, "+", fixed=TRUE)}
        v1<-vector()
        v2<-vector()
        edges<-vector()
        for (i in 1:dim(reacs)[1]){
          inputs<-unlist(strsplit(reacs[i,1],"+", fixed=TRUE))
          if (length(inputs)==1){
            v1[CountReac] = reacs[i,1]
            edges[CountReac] = 1
            v2[CountReac] = reacs[i,2]
            if (length(grep("!", v1))){
                v1[CountReac] = sub("!", "", v1[CountReac])
                edges[CountReac] = -1
            }
            BStimes[CountReac]<-optimBStimes[i]
            Integr[CountReac]<-optIntegr[i]
            CountReac<-CountReac+1
          }else{
            for (j in 1:length(inputs)){
              v1[CountReac]<-inputs[j]
              edges[CountReac] = 1
              v2[CountReac]<-paste("and",CountAnds,sep="")
              if (length(grep("!", v1[CountReac]))){
                v1[CountReac] = sub("!", "", v1[CountReac])
                edges[CountReac] = -1
              }
              BStimes[CountReac]<-optimBStimes[i]
              Integr[CountReac]<-optIntegr[i]
              CountReac<-CountReac+1
            }
            v1[CountReac]<-paste("and",CountAnds,sep="")
            edges[CountReac] = 1
            v2[CountReac] = reacs[i,2]
            BStimes[CountReac]<-optimBStimes[i]
            Integr[CountReac]<-optIntegr[i]
            CountReac<-CountReac+1
            vertices<-c(vertices,paste("and",CountAnds,sep=""))
            CountAnds<-CountAnds+1

          }
        }
        

    }
    
    

    if (is.null(cnolist) == FALSE){ # if a cnolist is provided, fill
                                    # signals/stimulis/inhitors
        stimulis <- cnolist$namesStimuli
        signals <- cnolist$namesSignals
        inhibitors <- cnolist$namesInhibitors
        inhibitors <- cnolist$namesInhibitors
    }


    # build the edges. IGraph does not use names for the vertices but ids 
    # that starts at zero. Let us build a data.frame to store the correspondence
    # between the ids and names.
    l = length(vertices) - 1

    # build the graph now
    g <- new("graphNEL", nodes=vertices, edgemode="directed")
    weights = rep(1, l)
    for (i in 1:length(v1)){
        g <- addEdge(as.character(v1[i]), as.character(v2[i]), g, weights[i])
    }

    # ----------------------------------------------- Build the node attributes list
    nodeAttrs <- list()
    fillcolor <- list()
    color <- list()
    style <- list()  # the style of what is inside the node. For the lty (line
                     # style, this should be different
    lty <- list()
    height <-list()
    label <- list()
    width <- list()
    fixedsize <- list()
    # default. Must be filled in case no signals/stimulis/cnolist are provide
    # default. Must be filled in case no signals/stimulis/cnolist are providedd
    for (s in vertices){
        color[s] <- "black"
        fillcolor[s] <- "white"
        style[s] <- "filled, bold"
        lty[s] <- "solid"
        label[s] <- s
        height[s] <- 1
        width[s] <-1
        fixedsize[s] <- FALSE
    }
    #lty <- list()
    # color first (fill and contour)
    for (s in stimulis){ 
        fillcolor[s] <- "olivedrab3"; 
        color[s] <- "olivedrab3";
        style[s] <- "filled"
        lty[s] <- "solid"
    }
    for (s in signals){ #must be before the inhibitors to allow bicolors
        fillcolor[s] <- "SkyBlue2"; 
        color[s] <-"SkyBlue2";
    #    style[s] <- "filled"
        lty[s] <- "solid"
    }
    for (s in inhibitors){
        if (length(grep(s, signals))>=1){
            fillcolor[s] <- "SkyBlue2"
            style[s] <- "filled,bold,diagonals"
            lty[s] <- "solid"
            color[s] <-"orangered"
        }
        else{
            fillcolor[s] <- "orangered"; 
            color[s] <-"orangered";
     #       style[s] <- "filled"
        lty[s] <- "solid"
        }
    }
    for (s in compressed){ 
        fillcolor[s] <- "white"; 
        color[s] <- "black"; 
        style[s] <- "dashed,bold"; 
        lty[s]='dashed';
    }
    for (s in ncno){ 
        fillcolor[s] <- "white"; 
        color[s] <- "grey"; 
        fillcolor[s] = "grey"
#        style[s] <- "filled"; 
        lty[s]='dashed';
    }
    for (s in vertices){
        if (length(grep("and", s))>=1){
            color[s] = "black"
            fillcolor[s] = "black"
            width[s]=0.10
            height[s]=0.10
            fixedsize[s]=TRUE
 #           style[s]="filled"
            label[s] = ""
            lty[s]= 'solid'
        }
    }
    nodeAttrs <- list(
        fillcolor=fillcolor,
        color=color,
        label=label,
        width=width,
        height=height,
        style=style,lty=lty,
        fixedsize=fixedsize)


    edgeAttrs <- list()
    arrowhead <- list()
    edgecolor <- list()
    
    ######## MODIFIED: , selected links in black other in  ######
    # negative links not colored
    # unselected in grey
    # selected links in black
    # integrated selected links in blue
    
    for (i in 1:length(edges)){
       edgename = paste(v1[i], "~", v2[i], sep="")
       edgecolor[edgename] <- "grey90"
       if (edges[i] == 1){
           arrowhead[edgename] <- "normal"
           #edgecolor[edgename] <- "black"
       }
       else if (edges[i] == -1){
           arrowhead[edgename] <- "tee"
           #edgecolor[edgename] <- "red"
       }
       if (BStimes[i] == 1){
           edgecolor[edgename] <- "black"
       }

    }


    indexI<-intersect(which(Integr==1), which(BStimes==1))
    edgecolor[indexI]<-"red"
    
    
    edgeAttrs <- list(color=edgecolor,arrowhead=arrowhead)
    # --------------------------- the ranks computation for the layout

    if (is.null(signals) || is.null(stimulis)){
       clusters = NULL
    }
    else{
        clusters = create_layout(g, signals, stimulis)
   } 

    # ------------------------------ general attributes
    # for some obscure reasons, must set style="filled" right her even though it
    # is then overwritten by nodesAttrs$style later on otherwise the output.dot
    # does not contain any style option
    attrs <-
        list(node=list(shape="ellipse",fontsize=22,fontname="Helvetica",style="filled,bold"),
             edge=list(style="solid", penwidth="2.0", weight="1.0"),
             graph=list(splines=TRUE,size="18.5,11",bgcolor="white",ratio="fill"))
    #other options
    #in graph: pad="0.5,5"))

   # --------------------------------------- some output 
   if (PDF==TRUE){ 
        if (is.null(filename)){
            pdf("pkn.pdf")
        }
        else {
            pdf(paste(filename, ".pdf", sep=""))
        }
    }
    if (PNG==TRUE){ 
        if (is.null(filename)){
            png("pkn.png")
        }
        else {
            png(paste(filename, ".png", sep=""))
        }
    }
    if (SVG==TRUE){
        if (is.null(filename)){
            svg("pkn.svg", pointsize=22, width=10, height=10)
        }
        else {
            svg(paste(filename, ".svg", sep=""),pointsize=22, width=10, height=10)
        }
    }

    if (is.null(filename)){
        output_dot = "tst.dot"
    }
    else{
        output_dot = paste(filename, ".dot", sep="")
    }

    # ------------------ plotting
    recipEdges="distinct"
    if (is.null(clusters)){
        plot(g,"dot",attrs=attrs,nodeAttrs=nodeAttrs,edgeAttrs=edgeAttrs,recipEdges=recipEdges)
        toDot(g, output_dot, nodeAttrs=nodeAttrs,edgeAttrs=edgeAttrs,attrs=attrs, recipEdges=recipEdges)
    }
    else{
        copyg <- g
print(nodeAttrs$lty)
#        plot(g,"dot",attrs=attrs,nodeAttrs=nodeAttrs,edgeAttrs=edgeAttrs,subGList=clusters,recipEdges=recipEdges)

        nodeRenderInfo(g) <- list(
            fill=nodeAttrs$fillcolor, 
            col=nodeAttrs$color,
            style=nodeAttrs$style,
#            lty=nodeAttrs$lty,
            #iagonals,bold",
#            lWidth=nodeAttrs$width,
            lwd=2,
#            height=2,
            textLabel=nodeAttrs$label,
#            textCol="red",
            fontsize=22,
            shape="ellipse",
            cex=0.5
            )
#othetr possible arguments:    fontsize, cex(expansion factor)

#arrowhead: Object of class "character" The style of arrowh#ead for this edge.
#arrowsize: Object of class "character" A scale factor for the length of the
#arrow heads & tails
#arrowtail: Object of class "character" The style of arrowtail for this edge.
#color: Object of class "character" The edge color.
#dir: Object of class "character" The edge direction.
#ep: Object of class "xyPoint" The end point of the edge.
#head: Object of class "character" The head node for this edge.
#lty: Object of class "character" The edge line type. e.g. (dashed)
#lwd: Object of class "numeric" The edge line width.
#sp: Object of class "xyPoint" The starting point of the edge.
#splines: Object of class "list" A list of BezierCurve objects
#tail: Object of class "character" The tail node for this edge.
#txtLabel: Object of class "character" The edge text label.

        edgeRenderInfo(g) <- list(
            color=edgeAttrs$color,
#           arrowhead=edgeAttrs$arrowhead,
#            arrowsize=4,
#            head=v2,
#            tail=v1,
#            label="test",
##            dir=-1,
            lwd=2
            #lty="dashed"
        )

#        graphRenderInfo(g) <-  list(
#                    node=list(
#                        shape="circle",
#                        fontsize=22,
#                        fontname="Helvetica"),
#                    edge=list(style="dashed"),
#                    graph=list(
#                        splines=FALSE,
#                        #size="18.5,11",
#                        bgcolor="white"#, ratio="fill"))
#                ))

        x <- layoutGraph(g, layoutType="dot",  subGList=clusters,recipEdges=recipEdges)
        renderGraph(x)
#        title("plotModel draft version. Use the dot file for optimal
#rendering",outer=FALSE)

        toDot(copyg,output_dot,nodeAttrs=nodeAttrs,edgeAttrs=edgeAttrs,subGList=clusters,attrs=attrs,recipEdges=recipEdges)
    }
    #title(main="Prior Knowledge Network",  xlab="created by plotModel")
    if (PDF==TRUE || SVG==TRUE || PNG==TRUE) { dev.off()}

    return(list(graph=g, attrs=attrs, nodeAttrs=nodeAttrs, edgeAttrs=edgeAttrs, clusters=clusters))
}


# if a cnolist, or at least signals/stimuli, then we can create ranking for the
# layout 
create_layout <- function(g, signals, stimulis){

    # this algo will tell us the distance between vertices
    # distMatrix columns contains the distance from a vertex to the others
    distMatrix <- floyd.warshall.all.pairs.sp(g)
    distMatrix[is.infinite(distMatrix) == TRUE] <- -Inf # convert Inf to be able
                                                        # to use the max
    # we will need to know the sinks
    sinks  <- signals[degree(g, signals)$outDegree == 0]
    # compute max rank for each column
    ranks <- apply(distMatrix, 2, max)
    mrank = max(ranks, na.rm=TRUE)-1  # -1 because we already know the sinks

    clusters <- list()
    if (mrank >= 1){ # otherwise, nothing to do. 
        # for each different rank select the names of the column to create a
        # cluster
        for (rank in 1:mrank){  # starts at 1 although ranks starts at 0.
                                # However, 0 corresponds to the source that are know already
            # nodes found for a particular rank may contain a sink, that should
            # be removed
            nodes2keep = NULL
            nodes <- names(ranks[which(ranks==rank)])
            for (n in nodes){
                if (any(n==sinks) == FALSE){ nodes2keep <- c(nodes2keep, n)}
            }

            # may fail sometimes
            tryCatch({
                thisCluster <- subGraph(nodes2keep, g)
                thisGraph <-list(graph=thisCluster,cluster=FALSE,attrs=c(rank="same"))
                clusters[[length(clusters)+1]] =  thisGraph},
             error=function(e){})

        }
    }
    # first the stimulis
    tryCatch(
        {
            clusterSource <- subGraph(stimulis, g)
            clusterSource<-list(graph=clusterSource,cluster=FALSE,attrs=c(rank="source"))
        },
        error=function(e){}
    )

    # then the signals keeping only those with outDegree==0
    tryCatch(
        {
            clusterSink <- subGraph(signals[degree(g, signals)$outDegree == 0], g)
            clusterSink <- list(graph=clusterSink, cluster=FALSE,
            attrs=c(rank="sink"))
        }, error=function(e){}
    )

    tryCatch(
        {clusters[[length(clusters)+1]] = clusterSource},
         error=function(e){}
    )
    tryCatch(
        {clusters[[length(clusters)+1]] = clusterSink}, 
        error=function(e){}
    )

    return(clusters)
}
