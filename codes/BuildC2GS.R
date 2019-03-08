# Build C2GS.R

library(GSA)
tmp = GSA.read.gmt('c2.all.v6.2.symbols.gmt') # Downloaded.
C2GS = tmp$genesets
names(C2GS) = tmp$geneset.names

# Apply gene size threshold
C2GS = C2GS[which(sapply(C2GS,length) >= 15 & sapply(C2GS,length) <= 500 )] # 3795 Geneset;

save(C2GS,file="C2GS.RData")

C2GENES = sort(unique(unlist(C2GS))) # 20019 Genes
C2AnnotationSum = as.numeric(table(unlist(C2GS)))
names(C2AnnotationSum)=C2GENES

save(C2AnnotationSum,file='C2AnnotationSum.RData')

C2GSIDX = list()
for(i in 1:length(C2GS)){
  C2GSIDX[[i]] = unlist(sapply(C2GS[[i]], function(k){
    which(k == rownames(string))
  }), use.names = FALSE)
}

save(C2GSIDX,file='C2GSIDX.RData')
