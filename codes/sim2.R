load("hnet.RData")
load("hnetSum.RData")
load("C2GSIDX-hnet.RData")
load("GeneLevel_hnet.RData")
load("GeneLevelIndex_hnet.RData")
load("GenesetV_hnet.RData")

load("string.RData")
load("StringNetworkSum.RData")
load("C2GSIDX.RData")
load("GeneLevel2.RData")
load("GeneLevelIndex2.RData")
load("GenesetV2.RData")

load("string5.RData")
load("String5NetworkSum.RData")
load("C2GSIDX_string5.RData")
load("GeneLevel_string5.RData.RData")
load("GeneLevelIndex_string5.RData")
load("GenesetV_.RData")

load("string7.RData")
load("String7NetworkSum.RData")
load("C2GSIDX_string7.RData")
load("GeneLevel_string5.RData.RData")
load("GeneLevelIndex_string5.RData")
load("GenesetV_string7.RData")

#od = pMM(unlist(sapply(genes, GethippieIndex)))
#od = pMM(unlist(sapply(genes, GethnetIndex)))
od = pMM(unlist(sapply(genes, GetStringIndex)))

SampleGeneLevel = GetSampleGeneLevel(genes)
numCores = detectCores()
cl = makeCluster(numCores-1)
registerDoParallel(cl)
pv = ( foreach(i=1:10000, .combine = '+',.inorder = FALSE) %dopar% { sim() } + 1 ) / 10001
stopCluster(cl)
pvh = Getpvh(genes)
names(pv) = names(C2GS)
names(pvh) = names(C2GS)
RankLine(pv)
RankLine(pvh)
RankLine(pv, type='p')
RankLine(pvh, type='p')


#lines(RankLine(pvString)[1:min(c(length(which(p.adjust(pv,'fdr')<0.25)),100))], type='o',col = '#fd79a8',lwd = 2)

# NETGOS, red, pink, purple
plot(RankLine(pvString), 
     type='l', col = '#e74c3c',lwd = 5,
     main = paste0('TRUE in Rank 100 (Dot with Qvalue < 0.25), ',length(genes) ,' genes' ), 
     ylab = 'True Positive Genesets', xlab = 'Rank', ylim  = c(0,12))
lines(RankLine(pvHnet),col = '#FDA7DF', lwd = 5)
lines(RankLine(pvHippie),col = '#D980FA', lwd = 5)
# Binoxes - sky, blue
lines(RankLine(binoxHnet),col = '#74b9ff', lwd = 5)
lines(RankLine(binoxString),col = '#0652DD', lwd = 5)
# Neas - light green, deep green, 
lines(RankLine(neahnet),col = '#C4E538', lwd = 5) 
lines(RankLine(neahippie),col = '#A3CB38', lwd = 5)
lines(RankLine(neastring),col = '#009432', lwd = 5)

# Hyper - yellow
lines(RankLine(pvh), col = '#FFC312', lwd = 5)
# TRUE - black
abline(a=0,b=length(ans)/3795, lwd = 5, col = '#2f3542')

#lines(RankLine(pvh)[1:min(c(length(which(p.adjust(pvh,'fdr')<0.25)),100))], col = '#74b9ff', lwd = 2, type='o')

legend(0,12,col = 
         c('#e74c3c', '#FDA7DF', '#D980FA', '#74b9ff','#0652DD','#C4E538','#A3CB38','#009432','#FFC312','#2f3542'), 
       lwd = 5, 
       c('netGOString','netGOHnet','netGOHippie','binoxHnet', 'binoxString','neaHnet', 'neaHippie','neaString', 'Hyper', 'Random')
       )
     
genes = p53b[1:15] # 12/8, 12/10
genes = p53b[1:20] # 11/4, 11/9
genes = p53b[1:30] # 6/4, 6/10
genes = p53b[1:50] # 6/5, 6/9

genes = myc[1:5] # 5/0,5/4
genes = myc[1:10] # 5/1, 4/5
genes = myc[1:15] # 7/4, 6/9
genes = myc[1:20] # 5/5, 5/9
genes = myc[1:30] # 3/4,3/6
genes = myc[1:50] # 6/7, 6/9

genes = brca[1:5] 
genes = brca[1:10]
genes = brca[1:15] 
genes = brca[1:20]
genes = brca[1:30]
genes = brca[1:50]

