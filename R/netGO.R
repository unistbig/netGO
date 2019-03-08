#'  @import doParallel
#'  @import foreach
#'  @import parallel

getHyperPvalue = function(genes, genesets){
  A = length(unique(unlist(genesets)))
  pv = rep(0,length(genesets))

  for(i in 1:length(genesets)){
    q = length(intersect(genesets[[i]],genes))
    m = length(genesets[[i]])
    n = A-m
    k = length(genes)
    pv[i] = 1-phyper(q-1,m,n,k)
  }
  names(pv) = names(genesets)
  return(pv)
}

IndexGenes = function(genes, rn){ unlist(sapply(genes, function(i){ which(i==rn) }, USE.NAMES = F)) }

GetSamplePPILV = function(genes, PPILV){ table(PPILV[genes]) }

Resample = function(SamplePPILV, PPILV){
  res = c()
  if(!is.na(SamplePPILV['0'])){ res = c(res, sample(which(PPILV == '0'), size = SamplePPILV['0'])) } # NA
  if(!is.na(SamplePPILV['1'])){ res = c(res, sample(which(PPILV == '1'), size = SamplePPILV['1'])) } # 1Q
  if(!is.na(SamplePPILV['2'])){ res = c(res, sample(which(PPILV == '2'), size = SamplePPILV['2'])) } # 2Q
  if(!is.na(SamplePPILV['3'])){ res = c(res, sample(which(PPILV == '3'), size = SamplePPILV['3'])) } # 3Q
  if(!is.na(SamplePPILV['4'])){ res = c(res, sample(which(PPILV == '4'), size = SamplePPILV['4'])) } # 4Q
  return(names(res))
}

getPPI = function(PPI){
  PPISum = sapply(1:nrow(PPI),function(i){sum(PPI[i,],na.rm = T)})

  rn = rownames(PPI)

  names(PPISum) = rn
  res = rep(0,length(rn))
  v = summary(PPISum)
  for(i in 1:length(rn)){
    if(is.na(PPISum[rn[i]])){res[i] = 0}
    else if(PPISum[rn[i]]<=v[2]){res[i] = 1} # 1Q
    else if(PPISum[rn[i]]<=v[3]){res[i] = 2} # 2Q
    else if(PPISum[rn[i]]<=v[5]){res[i] = 3} # 3Q
    else {res[i] = 4}
  }
  PPILV = res
  names(PPILV) = rn
  return(PPILV)
}

BuildGenesetsI = function(rn, genesets){
  genesetsI = list()
  for(i in 1:length(genesets)){ genesetsI[[i]] = unlist(sapply(genesets[[i]], function(k){which(k == rn)}), use.names = FALSE)}
  names(genesetsI) = names(genesets)
  genesetsI
}

pMM = function(genes, genesI,  genesets, genesetsI, genesetV){
  ovl = sapply(1:length(genesets), function(i){length(intersect(genesets[[i]], genes))/length(genes)})
  alpha = 1
  net = sapply(1:length(genesets), function(i){
    B = genesetsI[[i]]
    U = setdiff(genesI, B)
    if(length(U)==0){ v = 0 }
    else{ v = sum(genesetV[U, i]) / length(U) / length(B) }
    v
  })
  return(1-(ovl+net*alpha))
}

getPvalue = function(genes, genesets, PPI, genesetV){
  require(foreach)
  genesI = IndexGenes(genes,rownames(PPI))
  genesetsI = BuildGenesetsI(rownames(PPI), genesets)
  PPILV = getPPI(PPI)

  od = pMM(genes, genesI, genesets, genesetsI, genesetV)
  SamplePPILV = GetSamplePPILV(genes, PPILV)


  sim = function(){
    sampled = Resample(SamplePPILV, PPILV)
    sampledI = IndexGenes(sampled,rownames(PPI))
    as.numeric(pMM(sampled, sampledI, genesets, genesetsI, genesetV) <= od)
  }


  numCores = parallel::detectCores()
  cl = parallel::makeCluster(numCores-1)
  on.exit(parallel::stopCluster(cl))
  doParallel::registerDoParallel(cl)
  pv = foreach(i=1:10000,.combine = '+',.inorder = FALSE) %dopar% { sim() }
  pv = (pv+1)/10001

  names(pv) = names(genesets)
  return(pv)

}


#' @export
BuildGenesetV = function(PPI, genesets){
  GenesetV = matrix(0,nrow(PPI), length(genesets))
  GetV = function(geneA, geneB){return(sum(PPI[geneA,geneB]))}
  for(i in 1:nrow(PPI)){
    SubGeneset = rep(0,length(genesets))
    ThisGene = i
    for(j in 1:length(genesets)){
      ThisGeneset = genesets[[j]]
      GenesetPart = setdiff(ThisGeneset, ThisGene)
      v = GetV(ThisGene, GenesetPart)
      if(v!=0){SubGeneset[j] = v}
    }
    GenesetV[i,] = SubGeneset
  }
  rownames(GenesetV) = rownames(PPI)
  colnames(GenesetV) = names(genesets)
  return(GenesetV)
}

#' @export
netGO = function(genes, genesets, PPI, genesetV){
  pvh = getHyperPvalue(genes, genesets)
  pv = getPvalue(genes, genesets, PPI, genesetV)
  return(list(pv = pv, pvh = pvh))
}

