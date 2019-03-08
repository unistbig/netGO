#' @export
GetPPIIndex = function(genes, PPI){ unlist(sapply(genes, function(i){ which(i==rownames(PPI)) })) }