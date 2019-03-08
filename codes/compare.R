l1 = which(rank(pv, ties.method = 'first')<=100)
l2 = which(rank(pvh, ties.method = 'first')<=100)

myf = function(pv){  sapply(1:100, function(i){which(rank(pv,ties.method = 'first')==i)}) }

names(C2GS)[myf(pv)]
names(C2GS)[myf(pvh)]

myf2 = function(i){ mean(StringNetworkSum[C2GS[[i]]], na.rm = T) }         
myf3 = function(i){ length(C2GS[[i]]) }
sums = sapply(1:length(C2GS), function(i){sum(StringNetworkSum[C2GS[[i]]] , na.rm = T)})
means = sapply(1:length(C2GS), function(i){mean(StringNetworkSum[C2GS[[i]]] , na.rm = T)})
medians = sapply(1:length(C2GS), function(i){median(StringNetworkSum[C2GS[[i]]] , na.rm = T)})
b1 = sapply(1:100, function(i){myf2(myf(pv)[i])})
b2 = sapply(1:100, function(i){myf2(myf(pvh)[i])})
sizes = sapply(1:length(C2GS), myf3)
b3 = sapply(1:100, function(i){myf3(myf(pv)[i])})
b4 = sapply(1:100, function(i){myf3(myf(pvh)[i])})
cbind(l1,l2,sapply(l1,myf2), sapply(l2,myf2) )
union(l1,l2)
s = union(l1,l2)
sapply(s, myf2)
sapply(s, function(i){which(i==rank(pv, ties.method = 'first'))})
sapply(s, function(i){which(i==rank(pvh, ties.method = 'first'))})
k = cbind(round(pv,4),round(pvh,4), 
          round(p.adjust(pv,'fdr'),4),round(p.adjust(pvh,'fdr'),4),
          r1 = rank(p.adjust(pv,'fdr'), ties.method = 'first'), 
          r2 = rank(p.adjust(pvh,'fdr'), ties.method = 'first'), 
          sapply(names(C2GS), function(i){i %in% ans}),
          k2 = sapply(1:length(C2GS), function(i){
            GeneIdx = unlist(sapply(genes, GethippieIndex))
            B = C2GSIDX[[i]]
            U = setdiff(GeneIdx, B) # most killing part, 120 ms
            
            if(length(U)==0 || length(B)==0){ PPI = 0 }
            else{ PPI = sum(GenesetV[U, i]) / length(U) / length(B) } # 20 ms 
            return(PPI*1000)
          })
          ,
          sapply(1:length(C2GS), function(i){length(intersect(genes,C2GS[[i]]))}),
          sizes
          )
colnames(k ) = c("netGO", "Hyper","netGO", "Hyper",  "R1", "R2","TRUE","PPIAVE","OVL", "SIZE")
datatable( k[which(k[,5]<=100 | k[,6]<=100) ,] )
datatable( k[which(k[,7]==1) ,] )

plot(RankLine(pv),type='o', col = '#fd79a8',lwd = 2,
     main = 'TRUE in Rank 100 (Dot with Qvalue < 0.25)' , 
     ylab = 'True Positive Genesets', xlab = 'Rank')
lines(RankLine(pvh)[1:8], col = '#74b9ff', lwd = 2, type='o')
lines(RankLine(pvh), col = '#74b9ff', lwd = 2)
abline(a=0,b=length(ans)/3795, lwd = 2, col = '#2ecc71')
legend(0,10,col = c('#fd79a8', '#74b9ff', '#2ecc71'), lwd = 2, c('netGO', 'Hyper', 'Random'))
