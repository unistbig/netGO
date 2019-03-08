
GetPPILV = function(gene, n){ # 4 Group
  ppi = n[gene]
  v = summary(n)
  if(is.na(ppi)){return(0)}
  if(ppi<=v[2]){return(1)} # 1Q
  if(ppi<=v[3]){return(2)} # 2Q
  if(ppi<=v[5]){return(3)} # 3Q
  return(4)     
}


GeneLevel = sapply(rownames(string), GetPPILV, n = StringNetworkSum)
GeneLevelIndex = GeneLevel
names(GeneLevelIndex) = 1:length(GeneLevel)

save(GeneLevel,file='GeneLevel2.RData')
save(GeneLevelIndex,file='GeneLevelIndex2.RData')

save(GeneLevel,file='GeneLevel_hnet.RData')
save(GeneLevelIndex,file='GeneLevelIndex_hnet.RData')

save(GeneLevel,file='GeneLevel_hippie.RData')
save(GeneLevelIndex,file='GeneLevelIndex_hippie.RData')

save(GeneLevel,file='GeneLevel_string7.RData')
save(GeneLevelIndex,file='GeneLevelIndex_string7.RData')

save(GeneLevel,file='GeneLevel_string8.RData')
save(GeneLevelIndex,file='GeneLevelIndex_string8.RData')
