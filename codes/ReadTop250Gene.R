


# Assume column 3 : GeneSymbol
ReadTop250Gene = function(filename){
	tab = read.table(filename, header = T, stringsAsFactors = FALSE)
	tab = tab[which(tab[,3]!=''),]  # Remove Empty GeneSymbol
	tab = tab[!duplicated(tab[,3]),] # Remove Duplicated
	lists = grep('///',tab[,3])
	
### HANDLE THIS ONE
### LOAD ROWNAMES ONLY

	# Set GeneSymbol with intersect ; string
	load('string.RData')
	for(i in 1:length(lists)){
		txt = strsplit(tab[lists[i],3],'///')[[1]]
		txt = intersect(rownames(string), txt)
		if(length(txt)==0){tab[lists[i],3]=''}
		else{ tab[lists[i],3] = txt[1] }
	}

	tab = tab[which(tab[,3]!=''),3]
	
	return(tab[1:250])
}


# Example
# Genes : ReadTop250Gene(filename)
