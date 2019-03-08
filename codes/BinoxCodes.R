# Example Usage
# WriteBinoxGene(kras[1:50],'KRAS')
# result -> BinoxKRAS.tsv

#' @export
WriteBinoxGene = function(genes, keyword){
	filename = paste0('Binox',keyword,'.tsv')
	write(paste(genes,keyword, sep='\t'),file=filename)
}

# write binox network
# Example Usage
# load("string.RData")
# PPI = string;
# filename = 'binoxString.tsv'
# buildBinoxNetwork(filename, PPI)

#' @export
buildBinoxNetwork = function(filename, PPI, Cutoff = 0){
  rn = rownames(PPI)
  for(i in 1:nrow(PPI)){
    for(j in (i+1):ncol(PPI)){
      if(PPI[i,j]>Cutoff){
        write( paste(rn[i], rn[j], PPI[i,j], sep='\t'), file = filename, append = TRUE )
      }
    }
  }
}

# Example Usage
# Binoxpv= ReadBinoxRes('BinoxResult_p53.txt')
#
#' @export
ReadBinoxRes = function(filename){
  tab = readLines(filename)
  filename = paste0('2',filename)
  write(paste(strsplit(tab[1],'\t')[[1]][3:7],collapse = '\t'),file=filename)

  for(i in 1:(length(tab)-1)/3){
    write( paste( strsplit(tab[i*3],'\t')[[1]][2], tab[i*3+1], collapse = ''), file=filename,append = T)
  }
  tab = read.table(filename, header = TRUE,stringsAsFactors = FALSE)
  colnames(tab) = c("GENESET","Pvalue","FDR","RelationType",'PFC')

  # relationType '-' -> Reverse pvalue
  tab[which(tab[,4]=='-'),2] = 1-tab[which(tab[,4]=='-'),2]

  pv = tab[,2]
  names(pv) = tab[,1]
  return(pv)
}

#' @export 

# transpose c2gs to binox form. 
filename = 'BinoxPathwaySymbol_15_500.tsv'
for(i in 1:length(C2GS)){
  thisgs = C2GS[[i]]
  for(j in 1:length(thisgs)){
    write(paste(thisgs[j], names(C2GS)[i], sep = '\t'), file =filename, append = T )
  }
}

filename = 'brca20genes.tsv'

write(paste(genes, 'BRCA', sep = '\t'), file= filename)
r