

numCores = detectCores()
od = pMM(unlist(sapply(genes, GethnetIndex)))
SampleGeneLevel = GetSampleGeneLevel(genes)
cl = makeCluster(numCores-1)
registerDoParallel(cl)
pv = ( foreach(i=1:10000, .combine = '+',.inorder = FALSE) %dopar% { sim() } + 1 ) / 10001
stopCluster(cl)
pvh = Getpvh(genes)
names(pv) = names(C2GS)
names(pvh) = names(C2GS)

# hnet p53
ans = p53ANS
genes = p53b[1:5] # 2/4
genes = p53b[1:10] # 9/8, 10/10+
genes = p53b[1:15] # 10/8, 10/10+
genes = p53b[1:20] # 9/4, 9/9
genes = p53b[1:30] # 7/4, 7/10
genes = p53b[1:50] # 6/5, 6/9

# hnet breast
ans = brcaANS
genes = brca[1:5] # 6/1 , 7+/7
genes = brca[1:10] # 10/3, 10/9
genes = brca[1:15] # 5/4,  9/7
genes = brca[1:20] # 7/6, 11/10
genes = brca[1:30] # 12/10,12/10
genes = brca[1:50] # 14/13 , 14/13

# hnet myc
tab = read.delim('GSE120568.txt', stringsAsFactors = F)
g = unique(tab[,4])
g = intersect(g,rownames(hnet))
g = g[1:100]
myc = g

ans = mycANS # 48 answers.
genes = myc[1:5] # 5/0,5/4
genes = myc[1:10] # 8/1, 8/5
genes = myc[1:15] # 9/4, 10/9 #########################
genes = myc[1:20] # 10/5, 10/9 
genes = myc[1:30] # 7/4, 7/6
genes = myc[1:50] # 11/7, 11/9

# hippie myc

genes = myc[1:5] # 5/0,5/4
genes = myc[1:10] # 5/1, 4/5
genes = myc[1:15] # 7/4, 6/9
genes = myc[1:20] # 5/5, 5/9
genes = myc[1:30] # 3/4,3/6
genes = myc[1:50] # 6/7, 6/9
