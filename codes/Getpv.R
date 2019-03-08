library(doParallel)
library(foreach)
library(parallel)

GetStringIndex= function(gene){ which(gene == rownames(string)) }

Getpv = function(genes){	
  od = pMM(sapply(genes, GetStringIndex))
  SampleGeneLevel = GetSampleGeneLevel(genes)
  
  sim = function(){ as.numeric(pMM(Resample(SampleGeneLevel)) <= od) }
  
  numCores = detectCores()
  cl = makeCluster(numCores-1)
  registerDoParallel(cl)
  pv = ( foreach(i=1:10000, .combine = '+',.inorder = FALSE) %dopar% { sim() } + 1 ) / 10001
  stopCluster(cl)
  return(pv)
}
