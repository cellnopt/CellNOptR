downCue <-
function(AllCues,cue,Model,CueDown=vector()){
  # creates an array of nodes downstream the given cue (but upstream the previous cue)
  
  # 1. find reactions in which the investigated cue is the source node
  ix <- which(Model$interMat[rownames(Model$interMat)==cue] == -1)
  # 2. check if the investigated node is already in the array of downstream nodes
  #    if it is not, it adds it to the array
  if (sum(CueDown == cue) == 0){
    CueDown[length(CueDown)+1] <- cue
  }
  if (length(ix) > 0) {
    # 3. finds which species are target nodes in those reactions (see point 1.)
    index <- unique(as.matrix(which(Model$interMat[,ix]==1, arr.ind=TRUE))[,1])
    # 4. and keeps seraching/adding downstream nodes...
    for (i in 1:length(index)){
      #...until the next cue is reached
      if (sum(rownames(Model$interMat)[index[i]] == AllCues) == 0){
        cue <- rownames(Model$interMat)[index[i]]
        CueDown <- downCue(AllCues = AllCues, cue = cue, Model = Model, CueDown = CueDown)
      }
    }
  }
  return(CueDown)
}