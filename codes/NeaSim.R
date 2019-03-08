# nea Usage , load nea library;
NeaRes = nea(ags = genes, fgs = C2GS, network = network)

GetNeaPv = function(NeaRes){ 2*pnorm(-abs(NeaRes$zscore)) }

NeaSim = function(ags, fgs, network, filename){
	NeaRes = nea(ags, fgs, network)	
	save(NeaRes, file=filename) # Save for Later Usage
	GetNeaPv(NeaRes)
}

# NOTICE fgs for NeaSim should be 15 <= size <= 500

GetNeaPv = function(NeaRes){ 2*pnorm(-abs(NeaRes$zscore)) }
# brca 20 genes

# string nea

load('neastring.RData')
NeaResstring = nea(ags =ags, fgs = C2GS, network = network)
save(NeaResstring, file='nearesstring.RData')
