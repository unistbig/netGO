# handle nea, binox results

load("neaBRCA.Rdata")
load("neaBRCAH.Rdata")
load("neaMYC.Rdata")
load("neaMYCH.Rdata")
load("neaP53.Rdata")
load("neaP53H.Rdata")

GetNeaPv = function(NeaRes){ 2*pnorm(-abs(NeaRes$zscore)) }
nBRCA = GetNeaPv(neaBRCA)
nBRCAH = GetNeaPv(neaBRCAH)
nBRCA1 = GetNeaPv(neaBRCA1)
nBRCA1H = GetNeaPv(neaBRCA1H)

nMYC = GetNeaPv(neaMYC)
nMYCH = GetNeaPv(neaMYCH)

nP53 = GetNeaPv(neaP53)
nP53H = GetNeaPv(neaP53H)
save(BRCA, file='nBRCA.RData')
save(BRCAH, file='nBRCAH.RData')

save(P53, file='nP53.RData')
save(P53H, file='nP53H.RData')

save(MYC, file='nMYC.RData')
save(MYCH, file='nMYCH.RData')

ReadBinoxRes = function(filename){
  tab = readLines(filename)
  filename = paste0('2',filename)
  write(paste(strsplit(tab[1],'\t')[[1]][3:7],collapse = '\t'),file=filename)

  for(i in seq(from = 1, to = (length(tab)-1)/3, by = 1 )){
    write( paste( tab[i*3],strsplit(tab[i*3+1],'\t')[[1]][c(2,4)], collapse = ''), file=filename,append = T)
  }
  tab = read.table(filename, header = TRUE,stringsAsFactors = FALSE)
  colnames(tab) = c("GENESET","Pvalue","RelationType")

  # relationType '-' -> Reverse pvalue
  tab[which(tab[,4]=='-'),2] = 1-tab[which(tab[,4]=='-'),2]

  pv = tab[,2]
  names(pv) = tab[,1]
  return(pv)
}

bBRCA = ReadBinoxRes('BinoxBRCARes.txt')
bBRCAH = ReadBinoxRes('BinoxBRCAResH.txt')
bP53 = ReadBinoxRes('BinoxP53Res.txt')
bP53H = ReadBinoxRes('BinoxP53ResH.txt')
bMYC = ReadBinoxRes('BinoxMYCRes.txt')
bMYCH = ReadBinoxRes('BinoxMYCResH.txt')
bBRCA1 = ReadBinoxRes('BinoxBRCA1Res.txt')
bBRCA1H = ReadBinoxRes('BinoxBRCA1ResH.txt')

save(bBRCA, file='bBRCA.RData')
save(bBRCAH, file='bBRCAH.RData')

save(bBRCA1, file='bBRCA1.RData')
save(bBRCA1H, file='bBRCAH1.RData')

save(bMYC, file='bMYC.RData')
save(bMYCH, file='bMYCH.RData')
save(bP53, file='bP53.RData')
save(bP53H, file='bP53H.RData')
