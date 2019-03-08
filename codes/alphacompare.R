# alpha comparing

pMM = function(GeneIdx){
  alpha = 0
  A = length(GeneIdx)
  sapply(1:length(C2GS),
         function(i){
           B = C2GSIDX[[i]]
           U = setdiff(GeneIdx, B) # most killing part, 120 ms
           if(length(U)==0){ PPI = 0 }
           else{ PPI = sum(GenesetV[U, i]) / length(U) / length(B) # 20 ms
           }
           Ovl = (A-length(U))/A
           return(1-(Ovl+PPI*alpha))
         })
}

od = pMM(unlist(sapply(genes, GetStringIndex)))
SampleGeneLevel = GetSampleGeneLevel(genes)
cl = makeCluster(numCores-1)
registerDoParallel(cl)
pv0 = ( foreach(i=1:10000, .combine = '+',.inorder = FALSE) %dopar% { sim() } + 1 ) / 10001
stopCluster(cl)
names(pv0) = names(C2GS)

od = pMM(unlist(sapply(genes, GetStringIndex)))
SampleGeneLevel = GetSampleGeneLevel(genes)
cl = makeCluster(numCores-1)
registerDoParallel(cl)
pv10 = ( foreach(i=1:10000, .combine = '+',.inorder = FALSE) %dopar% { sim() } + 1 ) / 10001
stopCluster(cl)
names(pv10) = names(C2GS)

pvh = Getpvh(genes)
names(pvh) = names(C2GS)

k = cbind(
  rank(pv0, ties.method = 'first'),
  rank(pv10, ties.method = 'first'),
  rank(pvh, ties.method = 'first'),
  sapply(names(C2GS), function(i){i %in% ans}),
  k2 = sapply(1:length(C2GS), function(i){
    GeneIdx = unlist(sapply(genes, GetStringIndex))
    B = C2GSIDX[[i]]
    U = setdiff(GeneIdx, B) # most killing part, 120 ms

    if(length(U)==0 || length(B)==0){
      PPI = 0
    }
    else{
      PPI = sum(GenesetV[U, i]) / length(U) / length(B) # 20 ms
    }
    return(PPI*1000)
  })
  ,
  round(sapply(1:length(C2GS), function(i){length(intersect(genes,C2GS[[i]]))})/sizes*100,4),
  round(sapply(1:length(C2GS),
               function(i){length(intersect(C2GSIDX[[i]], unlist(sapply(genes, GetStringIndex))))})/sizes*100,4),
  sizes
)
colnames(k) = c(
  "R1", "R2","R3", "R4", "R5", "TRUE","PPIAVE","OVLRATIO",
                "OVLRATIO2", "SIZE")
datatable(
  k[which(k[,1]<=100 | k[,2]<=100 | k[,3]<=100 | k[,4]<=100 | k[,5]<=100) ,]
)

ans = mycANS
genes = myc[1:20]

ans = p53ANS
genes = p53b[1:20]

ans = brcaANS
genes = brca[1:20]
