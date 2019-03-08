Getpvh = function(genes){
  A = length(unique(unlist(C2GS)))
  pv = rep(0,length(C2GS))
  
  for(i in 1:length(C2GS)){
    q = length(intersect(C2GS[[i]],genes))
    m = length(C2GS[[i]])
    n = A-m
    k = length(genes)
    pv[i] = 1-phyper(q-1,m,n,k)
  }
  return(pv)
}