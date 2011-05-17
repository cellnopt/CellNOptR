bin2dec <-
function(binaryvector){sum(2^(which(rev(binaryvector)==TRUE)-1))}

