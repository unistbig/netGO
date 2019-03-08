#' @export
#'
#'


# Assume normalized with 0 ~ 1 scale
getPPI = function(PPIMatrix){

  # ppi distribution
  PPISum = sapply(1:nrow(PPIMatrix),function(i){sum(PPIMatrix[i,],na.rm = T)})

}
