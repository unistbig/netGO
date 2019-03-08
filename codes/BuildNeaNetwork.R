# Build and Return that object
# Example 
# network = BuildNeaNetwork('NeaNetworkString.txt')

BuildNeaNetwork = function(filename, Cutoff = .75){	
	rn = rownames(PPI)
	for(i in 1:nrow(PPI)){
		for(j in (i+1):ncol(PPI)){ 
			if(PPI[i,j]>Cutoff){ 
				write( paste(rn[i], rn[j], sep=' '), file = filename, append = TRUE ) 
			} 
		}
	}
	network = readLines(filename)
	network
}
