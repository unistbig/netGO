# 
GetStringIndex= function(gene){ which(gene == rownames(string)) }
GethnetIndex = function(gene){ which(gene == rownames(hnet)) }
GethippieIndex = function(gene){ which(gene == rownames(hippie)) }
GetV = function(geneA, geneB){ return(sum(string[geneA, geneB])) }
GetV = function(geneA, geneB){ return(sum(hnet[geneA, geneB])) }
GetV = function(geneA, geneB){ return(sum(hippie[geneA, geneB])) }
BuildGenesetV = function(){
  GenesetV = matrix(,nrow(string), length(C2GS))
  
  for(i in 1:nrow(string)){
    if(i%%10==0){print(i)} # CHECK
    
    SubGeneset = rep(0,length(C2GS))
    ThisGene = i # GetStringIndex(rownames(string)[i])
    
    for(j in 1:length(C2GS)){
      ThisGeneset = C2GSIDX[[j]]
      GenesetPart = setdiff(ThisGeneset, ThisGene)
      v = GetV(GenesetPart, ThisGene)
      if(v!=0){SubGeneset[j] = v}
    }
    GenesetV[i,] = SubGeneset
  } 
  
  return(GenesetV)
}


GenesetV = BuildGenesetV() # 13948 * 3795
rownames(GenesetV) = rownames(string)
colnames(GenesetV) = names(C2GS)
save(GenesetV, file='GenesetV_hnet.RData')
save(GenesetV, file='GenesetV_hippie.RData')
save(GenesetV, file='GenesetV_string8.RData')

