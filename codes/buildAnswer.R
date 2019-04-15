tab = read.delim('GSE120568.txt', stringsAsFactors = F)
tab = tab[which(abs(tab[,'logFC']) >=1), ]
tab = tab[,-8]
tab = read.delim('GSE120568.txt', stringsAsFactors = F)
g = unique(tab[,4])
g = intersect(g,rownames(hnet))
g = g[1:100]
myc = g
write(paste(g, 1:100, sep = '\t'),ncolumns = 1,file='genes.rnk')

genes = g[1:5]
genes = g[1:10]
genes = g[1:30]
genes = g[1:50]

mycANS = names(C2GS)[union(grep("_MYC", names(C2GS)), grep("MYC_", names(C2GS)))]
mycANS = mycANS[-4]
mycANS = mycANS[-c(6:11)]

# MYCOSIS, MYCOBACTERIUM removed

ans = mycANS

teloANS = union(names(C2GS)[grep("_TEL", names(C2GS))], names(C2GS)[grep("_TERT", names(C2GS))])

# http://science.sciencemag.org/content/284/5419/1431.2
# p53, FAS, Apoptosis
p53ANS = unique(c(names(C2GS)[grep("P53", names(C2GS))], names(C2GS)[grep("FAS", names(C2GS))], names(C2GS)[grep("APOPT", names(C2GS))]))
p53ANS = names(C2GS)[grep("P53", names(C2GS))] # p53 only

# GSE120559 MYC , MYC GENESETS;
genes = g[1:5] # 10/1
genes = g[1:10] # 1/0 ???
genes = g[1:15] # 4/1
genes = g[1:30] # 8/3
genes = g[1:50] # 6/6
genes = g[1:100] # 11/10

# GSE120568 - MYC GENESETS;
genes = g[1:5] # 11/5
genes = g[1:10] # 18/4 ???
genes = g[1:15] # 17/4
genes = g[1:30] # 15/3
genes = g[1:50] # 10/2
genes = g[1:100] # 7/4

# GSE35896 Colon cancer p53+-, rank plot with QV
genes = p53b[1:5] # 11/9
genes = p53b[1:15] # 17/13
genes = p53b[1:20] # 10/7
genes = p53b[1:30] # 8/6
genes = p53b[1:50] # 9/5

# GSE3744, breast cancer
brcaANS = union(names(C2GS)[grep("BRCA", names(C2GS))], names(C2GS)[grep("BREAST", names(C2GS))])
ans = brcaANS
genes = brca[1:5] # 8/1
genes = brca[1:10]
genes = brca[1:15] # 11/4
genes = brca[1:20]
genes = brca[1:30]
