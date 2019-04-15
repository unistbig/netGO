# BUILD GENES

# brca 20 genes.
filename = 'GSE3744BRCA.txt'
tab = read.delim(filename, header = T, stringsAsFactors = F)
genes = unique(tab[,2])
genes = genes[which(genes!='')]
genes[1] = 'FIGF'
genes[8] = 'TNXB'
brca = genes[1:20]

# p53 20 genes.
filename = 'GSE35896P53.txt'
tab = read.delim(filename, header = T, stringsAsFactors = F)
genes = unique(tab[,2])
genes = genes[which(genes!='')]
p53 = genes[1:20]

# myc 20 genes. ( THIS )
filename = 'GSE120568MYC.txt'
tab = read.delim(filename, header = T, stringsAsFactors = F)
genes = unique(tab[,2])
genes = genes[which(genes!='')]
myc = genes[1:20]

# binox 20 genes

# GENE1\tBRCA
# GENE2\tBRCA

filename = 'GSE98699BRCA1.txt'
tab = read.delim(filename, header = T, stringsAsFactors = F)
genes = unique(tab[,2])
genes = genes[which(genes!='')]
brca = genes[1:20]
save(brca,file='brca1.RData')

