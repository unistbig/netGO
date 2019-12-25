library(foreach)
library(parallel)
library(doParallel)
library(shinyCyJS)

getHyperPvalue <- function(genes, genesets) {
  A <- length(unique(unlist(genesets)))
  pv <- rep(0, length(genesets))

  for (i in 1:length(genesets)) {
    q <- length(intersect(genesets[[i]], genes))
    m <- length(genesets[[i]])
    n <- A - m
    k <- length(genes)
    pv[i] <- 1 - phyper(q - 1, m, n, k)
  }
  names(pv) <- names(genesets)
  return(pv)
}

IndexGenes <- function(genes, rn) {
  sort(unlist(sapply(genes, function(i) {
    which(i == rn)
  }, USE.NAMES = F)))
}

GetSamplenetworkLV <- function(genes, networkLV) {
  table(networkLV[genes])
}

Resample <- function(SamplenetworkLV, networkLV, category) {
  uns <- unique(names(SamplenetworkLV))
  res <- c()
  for (i in 1:length(uns)) {
    res <- c(res, sample(which(networkLV == uns[i]), size = SamplenetworkLV[uns[i]]))
  }
  return(names(res))
}

getnetwork <- function(network, category) {
  networkSum <- sapply(1:nrow(network), function(i) {
    sum(network[i, ], na.rm = T)
  })

  rn <- rownames(network)

  names(networkSum) <- rn
  res <- rep(0, length(rn))
  category <- category - 1
  # 10 Group
  v <- unname(quantile(networkSum, probs = seq(0, 1, 1 / category)))
  for (i in 1:length(rn)) {
    a <- unname(networkSum[rn[i]])
    if (is.na(a)) {
      res[i] <- 0
    } else {
      for (j in (category + 1):2) {
        if (a <= v[j]) {
          res[i] <- j - 1
        }
      }
    }
  }

  networkLV <- res
  names(networkLV) <- rn
  return(networkLV)
}

BuildGenesetsI <- function(rn, genesets) {
  lapply(genesets, IndexGenes, rn = rn)
}

DIF <- function(A, B) {
  V <- setdiff(A, B)
  if (L(V) == 0) {
    return(0)
  }
  return(V)
}
INT <- function(A, B) {
  intersect(A, B)
}
L <- function(A) {
  length(A)
}
UNI <- function(A, B) {
  union(A, B)
}

pMM <- function(genes, genesets, genesI, genesetV, RS, beta) { # for netGO (network only)
  NET <- sapply(1:L(genesets), function(i) {
    sum(genesetV[genesI, i]) / (sum(RS[genesI])^(1 - beta) * (length(genesets[[i]]))^(beta))
  })
  1 - NET / L(genes) # netGO
}

pMM2 <- function(genes, genesets, genesI, genesetV, RS, alpha, beta) { # for netGO+ (overlap + network)
  OVL <- sapply(1:L(genesets), function(i) {
    L(INT(genes, genesets[[i]]))
  })

  NET <- sapply(1:L(genesets), function(i) {
    sum(genesetV[genesI, i]) / (sum(RS[genesI])^(1 - beta) * (length(genesets[[i]]))^(beta))
  })
  1 - (OVL + NET * alpha) / L(genes) # netGO Plus
}

getValues <- function(genes, genesets, genesI, genesetV, RS, alpha, beta) {
  OVL <- sapply(1:L(genesets), function(i) {
    L(INT(genes, genesets[[i]]))
  })
  NET <- sapply(1:L(genesets), function(i) {
    sum(genesetV[genesI, i]) / (sum(RS[genesI])^(1 - beta) * (length(genesets[[i]]))^(beta))
  })
  OVL <- OVL / L(genes)
  NET <- alpha * NET / L(genes)
  list(OVL = OVL, NET = NET)
}

getPvalue <- function(genes, genesets, network, genesetV, alpha, beta, nperm, category, plus, verbose) {
  numCores <- parallel::detectCores()
  cl <- parallel::makeCluster(numCores - 1)
  on.exit(parallel::stopCluster(cl))
  additional <- FALSE
  LGS <- sapply(1:L(genesets), function(i) {
    length(genesets[[i]])
  })
  RS <- sapply(1:nrow(network), function(i) {
    sum(network[, i])
  })
  names(RS) <- rownames(network)
  if(verbose){
    cat("Indexing genes\n")
  }
  genesI <- IndexGenes(genes, rownames(network))

  if(verbose){
    cat("Build category\n")
  }
  networkLV <- getnetwork(network, category)
  SamplenetworkLV <- GetSamplenetworkLV(genes, networkLV)

  if(verbose){
    cat("Parallel functions load\n")
  }

  doSNOW::registerDoSNOW(cl)
  pb <- txtProgressBar(min = 0, max = nperm, width = 20)
  progress <- function(n) setTxtProgressBar(pb, n)
  opts <- list(progress = progress)
  rn <- rownames(network)
  if (plus) {
    if(verbose){
      cat("netGO skipped\n")
    }
  } else {
    # sim, od, pMM, pv -> netGO
    sim <- function() {
      sampled <- Resample(SamplenetworkLV, networkLV, category)
      sampledI <- IndexGenes(sampled, rn)
      pMM(sampled, genesets, sampledI, genesetV, RS, beta) <= od
    }
    od <- pMM(genes, genesets, genesI, genesetV, RS, beta)
    cnt <- 0
    if(verbose){
      cat("netGO Calculation start\n")
    }
    cat("\nProgress - each = means 5%\n")
    if (nperm > 50000) {
      if(verbose){
        cat("nperm > 50000, this will take time\n")
      }

      pv <- foreach::foreach(
        i = 1:(nperm / 10), .inorder = FALSE,
        .combine = "+", .noexport = "network",
        .options.snow = opts
      ) %dopar% {
        cnt <- cnt + 1
        setTxtProgressBar(pb, cnt)
        return(sim())
      }
      for (j in 1:9) {
        pv <- pv + foreach::foreach(
          i = 1:(nperm / 10), .inorder = FALSE,
          .combine = "+", .noexport = "network",
          .options.snow = opts
        ) %dopar% {
          cnt <- cnt + 1
          setTxtProgressBar(pb, cnt)
          return(sim())
        }
      }
    }
    else {
      pv <- foreach::foreach(
        i = 1:nperm, .inorder = FALSE,
        .combine = "+", .noexport = "network",
        .options.snow = opts
      ) %dopar% {
        cnt <- cnt + 1
        setTxtProgressBar(pb, cnt)
        return(sim())
      }
    }

    names(pv) <- names(genesets)
    nperms <- rep(nperm, length(genesets))
    names(nperms) <- names(pv)
    pv <- (pv + 1) / (nperms + 1)
    if(verbose){
      cat("\n")
      cat("netGO Calculation finished\n")
    }
  }

  # sim2, od2, pMM2, pv2 -> netGO+
  sim2 <- function() {
    sampled <- Resample(SamplenetworkLV, networkLV, category)
    sampledI <- IndexGenes(sampled, rn)
    pMM2(sampled, genesets, sampledI, genesetV, RS, alpha, beta) <= od2
  }
  od2 <- pMM2(genes, genesets, genesI, genesetV, RS, alpha, beta)
  cnt <- 0
  if(verbose){
    cat("netGO+ Calculation start\n")
  }
  cat("\nProgress - each = means 5%\n")
  if (nperm > 50000) {
    if(verbose){cat("nperm > 50000, this will take time\n")}
    pv <- foreach::foreach(
      i = 1:(nperm / 10), .inorder = FALSE,
      .combine = "+", .noexport = "network",
      .options.snow = opts
    ) %dopar% {
      cnt <- cnt + 1
      setTxtProgressBar(pb, cnt)
      return(sim2())
    }
    for (i in 1:9) {
      pv <- pv + foreach::foreach(
        i = 1:(nperm / 10), .inorder = FALSE,
        .combine = "+", .noexport = "network",
        .options.snow = opts
      ) %dopar% {
        cnt <- cnt + 1
        setTxtProgressBar(pb, cnt)
        return(sim2())
      }
    }
  }
  else {
    pv2 <- foreach::foreach(
      i = 1:nperm, .inorder = FALSE,
      .combine = "+", .noexport = "network",
      .options.snow = opts
    ) %dopar% {
      cnt <- cnt + 1
      setTxtProgressBar(pb, cnt)
      return(sim2())
    }
  }
  close(pb)

  names(pv2) <- names(genesets)
  nperms <- rep(nperm, length(genesets))
  names(nperms) <- names(pv2)

  pv2 <- (pv2 + 1) / (nperms + 1)
  if(verbose){
    cat("netGO+ Calculation finished\n")
  }

  # not use
  if (additional) {
    if(verbose){
      cat("Additional calculation start")
    }

    R <- rank(pv, ties.method = "min")
    higher <- c()
    for (i in 1:80) {
      A <- names(which(R == i))
      if (length(A) > 1) {
        higher <- c(higher, A)
      }
      if (length(higher) > 80) {
        break
      }
    }
    sim3 <- function() {
      sampled <- Resample(SamplenetworkLV, networkLV, category)
      sampledI <- IndexGenes(sampled, rn)
      as.numeric(pMM2(
        sampled, genesets[higher], sampledI,
        genesetV, RS, alpha, beta
      ) <= od2)
    }
    v <- 0
    if (length(higher) > 80) {
      cat("Additional permutation")
      while (length(higher) > 80 & v < 10) {
        v <- v + 1
        pv <- pv * (nperms + 1) - 1
        cat(".")
        od2 <- pMM2(
          genes, genesets[higher], genesI,
          genesetV, RS, alpha, beta
        )
        additional_pv <- foreach(
          i = 1:nperm, .inorder = FALSE,
          .combine = "+", .noexport = "network"
        ) %dopar% {
          sim3()
        }
        nperms[higher] <- nperms[higher] + nperm
        pv[higher] <- pv[higher] + additional_pv
        pv <- (pv + 1) / (nperms + 1)
        R <- rank(pv, ties.method = "min")
        higher <- c()
        for (i in 1:80) {
          A <- which(R == i)
          if (length(A) > 1) {
            higher <- c(higher, A)
          }
          if (length(higher) > 80) {
            break
          }
        }
      }
    }
  }

  values <- getValues(genes, genesets, genesI, genesetV, RS, alpha, beta)

  if (plus) {
    return(list(pv2 = pv2, values = values))
  }

  return(list(pv = pv, pv2 = pv2, values = values))
}

#' @export
BuildGenesetV <- function(network, genesets) {
  GenesetV <- matrix(0, nrow(network), length(genesets))
  GetV <- function(geneA, geneB) {
    return(sum(network[geneA, geneB]))
  }
  # added genesets indexing ; too slow as symbol
  for (i in 1:length(genesets)) {
    v <- intersect(rownames(network), genesets[[i]])
    v <- sapply(1:length(v), function(j) {
      which(rownames(network) == v[[j]])
    })
    genesets[[i]] <- v
  }

  for (i in 1:nrow(network)) {
    SubGeneset <- rep(0, length(genesets))
    ThisGene <- i
    for (j in 1:length(genesets)) {
      ThisGeneset <- genesets[[j]]
      GenesetPart <- setdiff(ThisGeneset, ThisGene)
      v <- GetV(ThisGene, GenesetPart)
      if (v != 0) {
        SubGeneset[j] <- v
      }
    }
    GenesetV[i, ] <- SubGeneset
  }
  rownames(GenesetV) <- rownames(network)
  colnames(GenesetV) <- names(genesets)
  return(GenesetV)
}

#' @export
netGO <- function(genes, genesets, network, genesetV,
                  alpha = 20, beta = 0.5, nperm = 10000, category = NULL,
                  pvalue = TRUE, plus = TRUE, verbose = TRUE) {

  pvh <- getHyperPvalue(genes, genesets)
  if(verbose){
   cat('Fisher Pvalue Calculation finished\n')
  }

  if(is.null(category)) {
    category <- 2000
    if (nrow(genesetV) < 10000) {
      category <- 1000
    }
    if(verbose){ cat(paste0(category,' genes in each category\n') )}
    category <- ceiling(nrow(genesetV) / category)
  }
  else {
    category <- ceiling(nrow(genesetV) / category)
  }
  obj <- getPvalue(genes, genesets, network, genesetV, alpha, beta, nperm, category, plus, verbose)

  if (plus) { # + only
    qv2 <- p.adjust(obj$pv2, "fdr") # netGO+
    qvh <- p.adjust(pvh, "fdr")
    ret <- data.frame(qv2, qvh)
    ret <- cbind(rownames(ret), ret)
    colnames(ret) <- c("gene-set", "netGO+Q", "FisherQ")
    if (pvalue) {
      tmp <- colnames(ret)
      ret <- cbind(ret, obj$pv2, pvh)
      colnames(ret) <- c(tmp, "netGO+P", "FisherP")
    }
  }

  if (!plus) {
    qv <- p.adjust(obj$pv, "fdr") # netGO
    qv2 <- p.adjust(obj$pv2, "fdr") # netGO+
    qvh <- p.adjust(pvh, "fdr")
    ret <- data.frame(qv, qv2, qvh)
    ret <- cbind(rownames(ret), ret)
    colnames(ret) <- c("gene-set", "netGOQ", "netGO+Q", "FisherQ")
    if (pvalue) {
      tmp <- colnames(ret)
      ret <- cbind(ret, obj$pv, obj$pv2, pvh)
      colnames(ret) <- c(tmp, "netGOP", "netGO+P", "FisherP")
    }
  }

  values <- obj$values
  tmp <- colnames(ret)
  ret <- cbind(ret, values$OVL, values$NET)
  colnames(ret) <- c(tmp, "OverlapScore", "NetworkScore")
  rownames(ret) <- NULL
  ret <- ret[order(ret$`netGO+Q`), ] # set order as netGO+ Qvalue
  return(ret)
}

#' @export
netGOVis <- function(obj, genes, genesets, network, R = 50, Q = 0.25) {
  suppressPackageStartupMessages("")

  appDir <- system.file("netGO", package = "netGO")
  if (appDir == "") {
    stop(
      "Could not find shinyCyJS Directory, Try re-installing 'shinyCyJS'.",
      call. = FALSE
    )
  }

  .GlobalEnv$.obj <- obj
  .GlobalEnv$.genes <- genes
  .GlobalEnv$.network <- network
  .GlobalEnv$.genesets <- genesets
  .GlobalEnv$.Q <- Q
  .GlobalEnv$.R <- R

  on.exit(rm(list = c(".obj", ".genes", ".network", ".genesets", ".R", ".Q"), envir = .GlobalEnv))

  shiny::runApp(appDir, launch.browser = TRUE, display.mode = "normal")
}

#' @export
DownloadExampleData <- function() {
  filelist <- c(
    "networkString.RData", "brca.RData", "brcaresult.RData", "c2gs.RData",
    "genesetVString1.RData", "genesetVString2.RData"
  )
  filelist <- setdiff(filelist, dir())

  urls <- paste0("https://github.com/unistbig/netGO-Data/raw/master/Human/", filelist)
  if (length(filelist)) {
    for (i in 1:length(urls)) {
      download.file(urls[i], filelist[i])
    }
  }

  filelist <- c(
    "networkString.RData", "brca.RData", "brcaresult.RData", "c2gs.RData",
    "genesetVString1.RData", "genesetVString2.RData"
  )

  for (i in 1:length(filelist)) {
    cat(paste0("Loading ", filelist[i], "\n"))
    load(filelist[i], envir = .GlobalEnv)
  }

  genesetV <- rbind(genesetV1, genesetV2)
  rm(genesetV1, genesetV2, envir = .GlobalEnv)
  assign("genesetV", genesetV, envir = .GlobalEnv)
}

#' @export
exportTable <- function(obj, type = "", R = 50, Q = 0.25) {

  sigIdx <- function(obj, R, Q) {
    pv = obj$`netGO+P`
    pvh = obj$FisherP
    qv <- obj$`netGO+Q`
    qvh <- obj$FisherQ

    names(pv) <- names(pvh) <- obj$`gene-set`
    names(qv) <- names(qvh) <- obj$`gene-set`

    if (!is.null(Q)) {
      idx <- which(qv <= Q | qvh <= Q)
    }
    if (!is.null(R)) {
      idx2 <- which(rank(pv, ties.method = "first") <= R |
                      rank(pvh, ties.method = "first") <= R)
      if(!exists('idx')){return(idx2)}

      idx = union(names(idx),names(idx2))
      idx = sapply(idx, function(i){which(obj$`gene-set`==i)})
    }
    return(idx)
  }
  si <- sigIdx(obj, R, Q)


  myTab <- cbind(names(si), round(cbind(obj$`netGO+Q`, obj$FisherQ)[si, ], 4))
  myTab <- data.frame(myTab, stringsAsFactors = FALSE)
  colnames(myTab) <- c("Gene-set name", "netGO+ q-value", "Fisher q-value")
  if(type=='D'){
    colnames(myTab) <- c("Gene-set name", "netGO+<br>q-value", "Fisher<br>q-value")
  }


  if(!is.null(obj$netGOQ)){ # netGO and netGO+
    myTab <- cbind(names(si), round(cbind(obj$`netGO+Q`, obj$`netGOQ`, obj$FisherQ)[si, ], 4))
    myTab <- data.frame(myTab, stringsAsFactors = FALSE)
    colnames(myTab) <- c("Gene-set name", "netGO+ q-value", "netGO q-value", "Fisher q-value")
    if(type=='D'){
      colnames(myTab) <- c("Gene-set name", "netGO+<br>q-value", "netGO<br>q-value", "Fisher<br>q-value")
    }

  }

  for(i in 2:ncol(myTab)){
    myTab[,i] = as.numeric(myTab[,i])
  }

  rownames(myTab) <- myTab[, 1]
  myTab <- myTab[order(myTab[, 2]), ] # sort by netGO+Q

  D <- myTab

  if (type == "D") {

    D <- datatable(
      myTab,
      rownames = FALSE,
      extensions = c("Scroller", "Buttons"),
      options = list(
        processing = TRUE,
        order = list(list(1, "asc")),
        deferRender = TRUE,
        scrollY = "34vh",
        scroller = TRUE,
        scrollX = TRUE,
        dom = "ltipr",
        autoWidth = TRUE,
        columnDefs = list(
          list(width = "60px", targets = c(1:(ncol(myTab)-1))),
          list(width = "100%", targets = 0)
        )
      ),
      selection = "single",
      escape = FALSE
    )
  }

  return(D)
}

#' @export
exportGraph <- function(genes, geneset, network) {
  res <- list()
  isobj <- getIntersectPart(genes, geneset, network)

  if (length(setdiff(genes, geneset))) {
    res[[length(res) + 1]] <-
      buildElems(
        data.frame(
          id = setdiff(genes, geneset),
          bgColor = "#FFFFFF",
          borderColor = "#48DBFB",
          borderWidth = 2,
          fontSize = 10,
          width = 60, height = 20, opacity = 1, stringsAsFactors = FALSE
        ), "Node"
      )
  }

  if (length(intersect(genes, geneset))) {
    res[[length(res) + 1]] <-
      buildElems(
        data.frame(
          id = intersect(genes, geneset),
          bgColor = "#FFFFFF",
          borderColor = "#03CB5D",
          borderWidth = 2,
          fontSize = 10,
          width = 60, height = 20, opacity = 1, stringsAsFactors = FALSE
        ), "Node"
      )
  }

  if (nrow(isobj$res)) {
    res[[length(res) + 1]] <- buildElems(data.frame(
      id = isobj$nodes,
      bgColor = "#FFFFFF",
      borderColor = "#FCCE00",
      borderWidth = 2,
      fontSize = 10,
      width = 60, height = 20, opacity = 1, stringsAsFactors = FALSE
    ), "Node")

    res[[length(res) + 1]] <- buildElems(isobj$res, "Edge")
  }
  return(c(res[[1]], res[[2]], res[[3]], res[[4]]))
}

getIntersectPart <- function(gene, geneset, network) {
  elements <- list()
  gs <- intersect(geneset, rownames(network))
  g <- intersect(gene, rownames(network))

  edges <- list()

  res <- data.frame(source = "", target = "", width = "", lineColor = "", opacity = "", stringsAsFactors = FALSE)
  for (i in 1:length(g)) {
    E <- network[g[i], names(which(network[g[i], gs] > 0))]
    if (length(E)) {
      if (length(E) == 1) {
        n <- names(which(network[g[i], gs] > 0))
        res <- rbind(
          res,
          data.frame(
            source = g[i],
            target = n, width = (E + .5) * 3, lineColor = "#4B4B4B", opacity = 1
          )
        )
        next
      }
      res <- rbind(
        res,
        data.frame(source = g[i], target = names(E), width = (unname(E) + .5) * 3, lineColor = "#4B4B4B", opacity = 1)
      )
    }
  }
  res <- res[-1, ]
  nodes <- union(res[, 1], res[, 2])
  return(list(res = res, nodes = nodes))
}

#' @export
exportGraphTxt <- function(gene, geneset, network) {
  rn <- rownames(network)
  g <- intersect(gene, rn)
  gs <- intersect(geneset, rn)

  res <- data.frame(geneA = "", geneB = "", strength = "", stringsAsFactors = FALSE)
  for (i in 1:length(g)) {
    E <- network[g[i], names(which(network[g[i], gs] > 0))]
    if (length(E)) {
      if (length(E) == 1) {
        n <- names(which(network[g[i], gs] > 0))
        res <- rbind(res, data.frame(geneA = g[i], geneB = n, strength = E, stringsAsFactors = FALSE))
        next
      }
      res <- rbind(res, data.frame(geneA = g[i], geneB = names(E), strength = unname(E), stringsAsFactors = FALSE))
    }
  }

  res <- res[-1, ]

  type <- rep("Inter", nrow(res))
  for (i in 1:nrow(res)) {
    if (res[i, 1] %in% g && res[i, 2] %in% g) {
      type[i] <- "Inner"
    }
  }
  res <- cbind(res, type)
  return(res)
}
