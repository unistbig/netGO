RankLine = function(pv, type = 'q'){
	if(type=='q'){qv = p.adjust(pv,'fdr')}
	else{qv = pv}
	if(type=='n'){
		res = rep(0,100)
		for(i in 1:100){
			targets = names(C2GS)[which(rank(qv, ties.method = 'first')<=i)] 
			res[i] = length(intersect(ans,targets))
		}
		res
	}
  
	else{
		res = rep(0,100)
		for(i in 1:100){
			targets = names(pv)[which(rank(qv, ties.method = 'first')<=i)]
			res[i] = length(intersect(ans,targets))
		}
		res
	}  
}


# each pv, pvh , pvnea, pvbinox
RankPlot = function(main, pv1, pv2, pv3, pv4, y = 12){
	plot(RankLine(pv1), type='l',ylim = c(0,y), lwd = 2, ylab = 'count', xlab = 'rank', col = '#7ed6df', main = main)
	lines(RankLine(pv2),col='#eb4d4b', lwd = 2) 
	lines(RankLine(pv3, type='b'), col='#6ab04c', lwd = 2)
	lines(RankLine(pv4),col='#be2edd', lwd = 2)
	legend(0,y,c('netGO', 'Hyper', 'BinoX', 'NEA'), lwd = 2, col = c("#7ed6df",'#eb4d4b','#6ab04c','#be2edd'))
}

# Example Usage
# Rankplot("KRAS 5 Genes", pv, pvh, pvnea, pvbinox, y = 12)

