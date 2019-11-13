library(shinyCyJS)

#' @export
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

#' @export
IndexGenes <- function(genes, rn) {
  sort(unlist(sapply(genes, function(i) {
    which(i == rn)
  }, USE.NAMES = F)))
}

#' @export
GetSamplenetworkLV <- function(genes, networkLV) {
  table(networkLV[genes])
}

#' @export
Resample <- function(SamplenetworkLV, networkLV) {
  res <- c()
  if (!is.na(SamplenetworkLV["0"])) {
    res <- c(res, sample(which(networkLV == "0"), size = SamplenetworkLV["0"]))
  } # NA
  if (!is.na(SamplenetworkLV["1"])) {
    res <- c(res, sample(which(networkLV == "1"), size = SamplenetworkLV["1"]))
  } # 1Q
  if (!is.na(SamplenetworkLV["2"])) {
    res <- c(res, sample(which(networkLV == "2"), size = SamplenetworkLV["2"]))
  } # 2Q
  if (!is.na(SamplenetworkLV["3"])) {
    res <- c(res, sample(which(networkLV == "3"), size = SamplenetworkLV["3"]))
  } # 3Q
  if (!is.na(SamplenetworkLV["4"])) {
    res <- c(res, sample(which(networkLV == "4"), size = SamplenetworkLV["4"]))
  } # 4Q
  # 10 Group added
  if (!is.na(SamplenetworkLV["5"])) {
    res <- c(res, sample(which(networkLV == "5"), size = SamplenetworkLV["5"]))
  } # 4Q
  if (!is.na(SamplenetworkLV["6"])) {
    res <- c(res, sample(which(networkLV == "6"), size = SamplenetworkLV["6"]))
  } # 4Q
  if (!is.na(SamplenetworkLV["7"])) {
    res <- c(res, sample(which(networkLV == "7"), size = SamplenetworkLV["7"]))
  } # 4Q
  if (!is.na(SamplenetworkLV["8"])) {
    res <- c(res, sample(which(networkLV == "8"), size = SamplenetworkLV["8"]))
  } # 4Q
  if (!is.na(SamplenetworkLV["9"])) {
    res <- c(res, sample(which(networkLV == "9"), size = SamplenetworkLV["9"]))
  } # 4Q
  if (!is.na(SamplenetworkLV["10"])) {
    res <- c(res, sample(which(networkLV == "10"), size = SamplenetworkLV["10"]))
  } # 4Q
  return(names(res))
}

#' @export
getnetwork <- function(network) {
  networkSum <- sapply(1:nrow(network), function(i) {
    sum(network[i, ], na.rm = T)
  })

  rn <- rownames(network)

  names(networkSum) <- rn
  res <- rep(0, length(rn))

  # 10 Group
  v <- unname(quantile(networkSum, probs = seq(0, 1, 1 / 9)))
  for (i in 1:length(rn)) {
    a <- unname(networkSum[rn[i]])
    if (is.na(a)) {
      res[i] <- 0
    }
    else if (a <= v[2]) {
      res[i] <- 1
    }
    else if (a <= v[3]) {
      res[i] <- 2
    }
    else if (a <= v[4]) {
      res[i] <- 3
    }
    else if (a <= v[5]) {
      res[i] <- 4
    }
    else if (a <= v[6]) {
      res[i] <- 5
    }
    else if (a <= v[7]) {
      res[i] <- 6
    }
    else if (a <= v[8]) {
      res[i] <- 7
    }
    else if (a <= v[9]) {
      res[i] <- 8
    }
    else if (a <= v[10]) {
      res[i] <- 9
    }
  }

  networkLV <- res
  names(networkLV) <- rn
  return(networkLV)
}

#' @export
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

#' @export
pMM <- function(genes, genesI, genesets, genesetsI, LGS, network) {
  ALPHA <- 1
  P <- function(A, B) {
    sum(network[A, B])
  }

  OVL <- sapply(1:L(genesets), function(i) {
    L(INT(genesets[[i]], genes))
  }) / L(genes)

  NET <- sapply(1:L(genesets), function(i) {
    g <- L(genes)
    A <- genesI
    B <- genesetsI[[i]]
    D <- DIF(A, B)

    w <- LGS[i] / (g + LGS[i] - g * OVL[i])
    v <- (P(D, B) + (w - 1) * P(D, INT(A, B))) / (LGS[i] + (1 - w) * g * OVL[i])

    v
  })

  return(1 - (OVL + ALPHA * NET))
}

#' @export
pMM2 <- function(genes, genesets, genesI, genesetV, RS, alpha, beta) {
  OVL <- sapply(1:L(genesets), function(i) {
    L(INT(genes, genesets[[i]]))
  })
  NET <- sapply(1:L(genesets), function(i) {
    sum(genesetV[genesI, i]) / (sum(RS[genesI])^(1 - beta) * (length(genesets[[i]]))^(beta))
  }) * alpha
  1 - (OVL + NET) / L(genes)
}

#' @export
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

#' @export
getPvalue <- function(genes, genesets, network, genesetV, alpha, beta, nperm) {
  library(parallel)

  numCores <- parallel::detectCores()
  cl <- parallel::makeCluster(numCores - 1)

  on.exit(parallel::stopCluster(cl))


  library(foreach)
  library(doSNOW)
  library(doParallel)

  additional <- FALSE
  LGS <- sapply(1:L(genesets), function(i) {
    length(genesets[[i]])
  })

  RS <- sapply(1:nrow(network), function(i) {
    sum(network[, i])
  })
  names(RS) <- rownames(network)
  sim <- function() {
    sampled <- Resample(SamplenetworkLV, networkLV)
    sampledI <- IndexGenes(sampled, rownames(network))
    as.numeric(pMM(sampled, sampledI, genesets, genesetsI, LGS, network) <= od)
  }

  sim2 <- function() {
    sampled <- Resample(SamplenetworkLV, networkLV)
    sampledI <- IndexGenes(sampled, rn)
    as.numeric(pMM2(sampled, genesets, sampledI, genesetV, RS, alpha, beta) <= od)
  }

  genesI <- IndexGenes(genes, rownames(network))

  networkLV <- getnetwork(network)

  od <- pMM2(genes, genesets, genesI, genesetV, RS, alpha, beta)

  SamplenetworkLV <- GetSamplenetworkLV(genes, networkLV)

  print("Parallel function loads")
  doSNOW::registerDoSNOW(cl)
  cnt = 0
  pb <- txtProgressBar(min = 0, max = nperm, width = 100)
  progress <- function(n) setTxtProgressBar(pb,n)
  opts <- list(progress = progress)

  rn <- rownames(network)
  print('Calculation start')
  print('Progress - each = mean 1%')

  if (nperm > 50000) {
    print("nperm > 50000")
    pv <- foreach::foreach(i = 1:(nperm / 10), .inorder = FALSE, .combine = "+", .noexport = "network", .options.snow = opts) %dopar% {
      cnt <- cnt + 1
      setTxtProgressBar(pb, cnt)
      return(sim2())
    }
    for (i in 1:9) {
      pv <- pv + foreach::foreach(i = 1:(nperm / 10), .inorder = FALSE, .combine = "+", .noexport = "network", .options.snow = opts) %dopar% {
        cnt <- cnt + 1
        setTxtProgressBar(pb, cnt)
        return(sim2())
      }
    }
  }
  else {
    pv <- foreach::foreach(i = 1:nperm, .inorder = FALSE, .combine = "+", .noexport = "network", .options.snow = opts) %dopar% {

      cnt <- cnt + 1
      setTxtProgressBar(pb, cnt)
      return(sim2())
    }
  }
  close(pb)
  print('Calculation finished')

  names(pv) <- names(genesets)

  nperms <- rep(nperm, length(genesets))
  names(nperms) <- names(pv)

  pv <- (pv + 1) / (nperms + 1)

  # more permutation with higher ranks
  if (additional) {
    print('Additional calculation start')
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
      sampled <- Resample(SamplenetworkLV, networkLV)
      sampledI <- IndexGenes(sampled, rn)
      as.numeric(pMM2(sampled, genesets[higher], sampledI, genesetV, RS, alpha, beta) <= od2)
    }
    v <- 0

    if (length(higher) > 80) {
      print("Additional permutation")
      while (length(higher) > 80 & v < 10) {
        v <- v + 1
        pv <- pv * (nperms + 1) - 1
        cat(".")
        od2 <- pMM2(genes, genesets[higher], genesI, genesetV, RS, alpha, beta)

        additional_pv <- foreach(i = 1:nperm, .inorder = FALSE, .combine = "+", .noexport = "network") %dopar% {
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
  return(list(pv = pv, values = values))
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
netGO <- function(genes, genesets, network, genesetV, alpha = 20, beta = 0.5, nperm = 10000) {
  pvh <- getHyperPvalue(genes, genesets)
  obj <- getPvalue(genes, genesets, network, genesetV, alpha, beta, nperm)
  pv <- obj$pv
  values <- obj$values

  ret <- data.frame(pv, pvh, p.adjust(pv, "fdr"), p.adjust(pvh, "fdr"))
  ret <- cbind(rownames(ret), ret)
  ret <- cbind(ret, values$OVL, values$NET)
  colnames(ret) <- c("gene-set", "netGOP", "FisherP", "netGOFDR", "FisherFDR", "overlap_score", "network_score")
  rownames(ret) <- NULL
  ret <- ret[order(ret$netGOP), ]
  return(ret)
}

#' @export
netGOVis <- function(obj, genes, genesets, network, R = 50, Q = NULL) {
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
  NowDir <- getwd()
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
    "networkString.RData", "brca.RData", "brcaresult.RData", "c2gs.RData","genesetVString1.RData", "genesetVString2.RData")
  for (i in 1:length(filelist)) {
    print(paste0("Loading ", filelist[i]))
    load(filelist[i], envir = .GlobalEnv)
  }
  genesetV <- rbind(genesetV1, genesetV2)
  rm(genesetV1,genesetV2, envir = .GlobalEnv)
  assign("genesetV", genesetV, envir = .GlobalEnv)
}
