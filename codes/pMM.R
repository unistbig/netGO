pMM = function(GeneIdx){
  #alpha = 0
  alpha = 1
  A = length(GeneIdx)
  sapply(1:length(C2GS), 
    function(i){
      B = C2GSIDX[[i]]
      U = setdiff(GeneIdx, B) # most killing part, 120 ms
      if(length(U)==0){
        PPI = 0
      }
      else{
        PPI = sum(GenesetV[U, i]) / length(U) / length(B) # 20 ms
      }
      Ovl = (A-length(U))/A
      
      return(1-(Ovl+PPI*alpha)) 
    })
}

