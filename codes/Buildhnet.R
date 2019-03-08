# build humannet
tab = read.delim("HumanNet.v1.join.txt", header = F) 
dim(tab)
library(biomaRt)
library(org.Hs.eg.db)
library(annotate)
a = getSYMBOL(as.character(tab[,1]),data = 'org.Hs.eg')
b = getSYMBOL(as.character(tab[,2]),data = 'org.Hs.eg')
genes = sort(union(a,b))
hnet = matrix(0,length(genes), length(genes))
rownames(hnet) = colnames(hnet) = genes
for(i in 1:nrow(tab)){
  geneA = getSYMBOL(as.character(tab[i,1]), data = 'org.Hs.eg')
  geneB = getSYMBOL(as.character(tab[i,2]), data = 'org.Hs.eg')
  hnet[geneA,geneB] = hnet[geneB,geneA] = tab[i,3]
  
}
