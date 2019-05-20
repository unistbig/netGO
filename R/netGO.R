#' @import doParallel
#' @import foreach
#' @import parallel
#'

#' @export
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

#' @export
IndexGenes = function(genes, rn){ unlist(sapply(genes, function(i){ which(i==rn) }, USE.NAMES = F)) }

#' @export
GetSamplePPILV = function(genes, PPILV){ table(PPILV[genes]) }

#' @export
Resample = function(SamplePPILV, PPILV){
  res = c()
  if(!is.na(SamplePPILV['0'])){ res = c(res, sample(which(PPILV == '0'), size = SamplePPILV['0'])) } # NA
  if(!is.na(SamplePPILV['1'])){ res = c(res, sample(which(PPILV == '1'), size = SamplePPILV['1'])) } # 1Q
  if(!is.na(SamplePPILV['2'])){ res = c(res, sample(which(PPILV == '2'), size = SamplePPILV['2'])) } # 2Q
  if(!is.na(SamplePPILV['3'])){ res = c(res, sample(which(PPILV == '3'), size = SamplePPILV['3'])) } # 3Q
  if(!is.na(SamplePPILV['4'])){ res = c(res, sample(which(PPILV == '4'), size = SamplePPILV['4'])) } # 4Q
  # 10 Group added
  if(!is.na(SamplePPILV['5'])){ res = c(res, sample(which(PPILV == '5'), size = SamplePPILV['5'])) } # 4Q
  if(!is.na(SamplePPILV['6'])){ res = c(res, sample(which(PPILV == '6'), size = SamplePPILV['6'])) } # 4Q
  if(!is.na(SamplePPILV['7'])){ res = c(res, sample(which(PPILV == '7'), size = SamplePPILV['7'])) } # 4Q
  if(!is.na(SamplePPILV['8'])){ res = c(res, sample(which(PPILV == '8'), size = SamplePPILV['8'])) } # 4Q
  if(!is.na(SamplePPILV['9'])){ res = c(res, sample(which(PPILV == '9'), size = SamplePPILV['9'])) } # 4Q
  if(!is.na(SamplePPILV['10'])){ res = c(res, sample(which(PPILV == '10'), size = SamplePPILV['10'])) } # 4Q
  return(names(res))
}

#' @export
getPPI = function(PPI){
  PPISum = sapply(1:nrow(PPI),function(i){sum(PPI[i,],na.rm = T)})

  rn = rownames(PPI)

  names(PPISum) = rn
  res = rep(0,length(rn))
  # 4 Group
  #v = summary(PPISum)
  #for(i in 1:length(rn)){
   #if(is.na(PPISum[rn[i]])){res[i] = 0}
    #else if(PPISum[rn[i]]<=v[2]){res[i] = 1} # 1Q
    #else if(PPISum[rn[i]]<=v[3]){res[i] = 2} # 2Q
    #else if(PPISum[rn[i]]<=v[5]){res[i] = 3} # 3Q
    #else {res[i] = 4}
  #}
  # 10 Group
  v = unname(quantile(PPISum, probs = seq(0,1,1/9)))
  for(i in 1:length(rn)){
    a = unname(PPISum[rn[i]])
    if(is.na(a)){res[i] = 0;}
    else if(a <= v[2]){res[i] = 1}
    else if(a <= v[3]){res[i] = 2}
    else if(a <= v[4]){res[i] = 3}
    else if(a <= v[5]){res[i] = 4}
    else if(a <= v[6]){res[i] = 5}
    else if(a <= v[7]){res[i] = 6}
    else if(a <= v[8]){res[i] = 7}
    else if(a <= v[9]){res[i] = 8}
    else if(a <= v[10]){res[i] = 9}
  }

  PPILV = res
  names(PPILV) = rn
  return(PPILV)
}

#' @export
BuildGenesetsI = function(rn, genesets){
  genesetsI = list()
  for(i in 1:length(genesets)){ genesetsI[[i]] = unlist(sapply(genesets[[i]], function(k){which(k == rn)}), use.names = FALSE)}
  names(genesetsI) = names(genesets)
  genesetsI
}

#' @export
pMM = function(genes, genesI,  genesets, genesetsI, genesetV){
  ovl = sapply(1:length(genesets), function(i){length(intersect(genesets[[i]], genes))/length(genes)})
  alpha = 1
  net = sapply(1:length(genesets), function(i){
    B = genesetsI[[i]]
    U = setdiff(genesI, B)
    if(length(U)==0){ v = 0 }
    else{
      v2 = genesetV[U, i]
      #v = sum(genesetV[U, i]) / length(U) / length(B)
      v2[which(v2>1.5)] = 1.5
      v = sum(v2)
      v = v / length(U) / length(B)
      return(v)
    }
    v
  })
  return(1-(ovl+net*alpha))
}
#' @export
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
  # added genesets indexing ; too slow as symbol
  for(i in 1:length(genesets)){
    v = intersect(rownames(PPI), genesets[[i]])
    v = sapply(1:length(v), function(j){which(rownames(PPI)==v[[j]])})
    genesets[[i]]=v
  }

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

#' @export
netGOVis = function(obj, brca, genesets, PPI, R = 20, Q = NULL){
  # arg1 -> inst/FOLDERNAME ( GScluster )
  # arg2 -> ...?

  # code to run netGOVis
  appDir = system.file("netGO", package = 'netGO')
  if(appDir ==''){
    stop(
      "Could not find shinyCyJS Directory, Try re-installing 'shinyCyJS'.",
      call. = FALSE
    )
  }

  .GlobalEnv$.obj = obj
  .GlobalEnv$.brca = brca
  .GlobalEnv$.PPI = PPI
  .GlobalEnv$.genesets = genesets
  .GlobalEnv$.Q = Q
  .GlobalEnv$.R = R

  on.exit(rm(list=c('.obj', '.brca', '.PPI','.genesets','.R','.Q'),
             envir=.GlobalEnv))

  shiny::runApp(
    appDir,
    launch.browser = TRUE,
    display.mode ='normal'
  )
}
