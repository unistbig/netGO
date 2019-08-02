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
IndexGenes = function(genes, rn){ sort(unlist(sapply(genes, function(i){ which(i==rn) }, USE.NAMES = F))) }

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
BuildGenesetsI = function(rn, genesets){ lapply(genesets, IndexGenes, rn = rn) }

DIF = function(A,B){
  V = setdiff(A,B);
  if(L(V)==0){return(0)}
  return(V)
}
INT = function(A,B){intersect(A,B)}
L = function(A){length(A)}
UNI = function(A,B){union(A,B)}

#' @export
pMM = function(genes, genesI, genesets, genesetsI, LGS, PPI ){
  ALPHA = 1
  P = function(A,B){ sum(PPI[A,B]) }

  OVL = sapply(1:L(genesets), function(i){L(INT(genesets[[i]], genes))}) / L(genes)

  NET = sapply(1:L(genesets), function(i){
    g = L(genes)
    A = genesI
    B = genesetsI[[i]]
    D = DIF(A,B)

    w = LGS[i] / ( g + LGS[i] - g * OVL[i] )
    v = ( P(D, B) + (w-1) * P(D, INT(A,B)) ) / ( LGS[i] + (1-w) * g * OVL[i] )

    v
  })

  return(1-(OVL+ALPHA*NET))
}

#' @export
pMM2 = function(genes, genesets, genesI, genesetV, RS, alpha){
  OVL = sapply(1:L(genesets), function(i){L(INT(genes, genesets[[i]]))})
  NET = sapply(1:L(genesets), function(i){
    # sum(genesetV[genesI,i]) / sqrt(  (sum(RS[genesI]) * length(genesets[[i]]) ) )
    sum(genesetV[genesI,i]) / (sum(RS[genesI])^(alpha) / (length(genesets[[i]]))^(1-alpha))
    })
  1-(OVL+NET)/L(genes)
}

#' @export
getPvalue = function(genes, genesets, PPI, genesetV, alpha, REP = 10000){

  LGS = sapply(1:L(genesets), function(i){length(genesets[[i]])})

  RS = sapply(1:nrow(PPI), function(i){sum(PPI[,i])})
  names(RS) = rownames(PPI)
  sim = function(){
    sampled = Resample(SamplePPILV, PPILV)
    sampledI = IndexGenes(sampled,rownames(PPI))
    as.numeric(pMM(sampled, sampledI, genesets, genesetsI, LGS, PPI) <= od)
  }

  sim2 = function(){
    sampled = Resample(SamplePPILV, PPILV)
    sampledI = IndexGenes(sampled,rn)
    as.numeric(pMM2(sampled, genesets, sampledI, genesetV, RS, alpha) <= od)
  }
  genesI = IndexGenes(genes,rownames(PPI)) # 0

  #genesetsI = BuildGenesetsI(rownames(PPI), genesets) # 30 second

  PPILV = getPPI(PPI) # 13 second

  # od = pMM(genes, genesI, genesets, genesetsI, LGS, PPI) # 0.1 second
  od = pMM2(genes, genesets, genesI, genesetV, RS, alpha)

  SamplePPILV = GetSamplePPILV(genes, PPILV) # 0 second

  numCores = parallel::detectCores()
  cl = parallel::makeCluster(numCores-1)
  on.exit(parallel::stopCluster(cl))
  doParallel::registerDoParallel(cl)
  rn = rownames(PPI)
  # ~ 1 second
  #pv = foreach( i=1:REP, .inorder = FALSE, .combine = '+', .verbose = TRUE ) %dopar% { sim() }
  pv = foreach(i = 1:REP, .inorder = FALSE, .combine = '+', .noexport = 'PPI') %dopar% {sim2()}
  #pv = Reduce(`+`,pv)

  pv = (pv+1)/(REP+1)
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
netGO = function(genes, genesets, PPI, genesetV, alpha = 0.5){

  require(foreach)
  require(parallel)
  require(doParallel)

  pvh = getHyperPvalue(genes, genesets)
  pv = getPvalue(genes, genesets, PPI, genesetV, alpha)
  return(list(pv = pv, pvh = pvh))
}

#' @export
netGOVis = function(obj, genes, genesets, PPI, R = 50, Q = NULL){
  suppressPackageStartupMessages('')
  # arg1 -> inst/FOLDERNAME ( GScluster )
  # arg2 -> ...?
  require(shinyCyJS)
  # code to run netGOVis
  appDir = system.file("netGO", package = 'netGO')
  if(appDir ==''){
    stop(
      "Could not find shinyCyJS Directory, Try re-installing 'shinyCyJS'.",
      call. = FALSE
    )
  }

  .GlobalEnv$.obj = obj
  .GlobalEnv$.genes = genes
  .GlobalEnv$.PPI = PPI
  .GlobalEnv$.genesets = genesets
  .GlobalEnv$.Q = Q
  .GlobalEnv$.R = R

  on.exit(rm(list=c('.obj', '.genes', '.PPI','.genesets','.R','.Q'),
             envir=.GlobalEnv))

  shiny::runApp(
    appDir,
    launch.browser = TRUE,
    display.mode ='normal'
  )
}
