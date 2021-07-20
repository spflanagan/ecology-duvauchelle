
# =========================================
### Calculate turnover (beta-1)/(N-1). Requires a measure of
### multiplicative beta diversity and sample size N.
# =========================================

turnover <- function(beta, N) {
  turn <- (beta - 1)/(N - 1)
  return(turn)
}


# =========================================
### Calculate 'local' (Sorensen-type) dissimilarity for
### any community size and any order q. Requires the order
### (q), multiplicative beta diversity estimate, and
### sample size (N)
# =========================================

local <- function(q, beta, N = 2) {
  if (q != 1) {
    CqN <- 1 - ((1/beta)^(q - 1) - (1/N)^(q - 1))/
      (1 - (1/N)^(q - 1))
  } else {
    CqN <- log(beta)/log(N)
  }
  return(CqN)
}


# =========================================
### Calculate 'regional' (Jaccard-type) dissimilarity for
### any community size and any order q. Requires the order
### (q), multiplicative beta diversity estimate, and
### sample size (N)
# =========================================

regional <- function(q, beta, N = 2) {
  if (q != 1) {
    UqN <- 1 - ((1/beta)^(1 - q) - (1/N)^(1 - q))/
      (1 - (1/N)^(1 - q))
  } else {
    UqN <- log(beta)/log(N)
  }
  return(UqN)
}


# =========================================
### Given a diversity order and community matrix, this
### function calculates alpha and gamma diversities
# =========================================

divPart <- function(q = 1, abundances) {
  pz <- abundances/sum(abundances)
  pg <- colSums(abundances)/sum(abundances)
  
  if (q == 1) {
    pq <- ifelse(pz > 0, pz * log(pz), 0) 
    alpha <- exp(-sum(colSums(pq)) - log(nrow(pq)))
    
    gamma <- exp(-1 * sum(pg[pg > 0] * log(pg[pg > 0])))
  } else {
    pq <- ifelse(pz > 0, pz^q, 0)
    alpha <- sum(colSums(pq))^(1/(1 - q)) * (1/nrow(pq))
    
    pGam <- ifelse(pg > 0, pg^q, 0)
    gamma <- sum(pGam)^(1/(1 - q))
  } 
  
  return(list(gamma = gamma,alpha = alpha))
}


# =========================================
### This function calculates pairwise equivalents of
### multiplicative beta diversity and local and regional
### dissimilarity. Requires abundances or incidences for a
### pair of sites and a sequence of diversity orders (q)
# =========================================

pairwiseFunc <- function(pairs, q.seq = q.seq) {
  divProf <- sapply(q.seq, divPart, abundances = pairs)
  
  beta <- unlist(divProf["gamma", ])/unlist(divProf["alpha",])
  loc <- mapply(local, q = q.seq, beta = beta, N = 2)
  reg <- mapply(regional, q = q.seq, beta = beta, N = 2)
  
  return(data.frame(beta = beta, local = loc, region = reg))
}

# =========================================
### noHier calculates pairwise and N-community
### diversities for a given sequence of q. It requires an
### nObs x nSpp matrix of abundances and a sequence of
### diveristy orders. It outputs alpha, multiplicative
### beta, gamma, and if pairwise = TRUE, mean pairwise beta 
### (as turnover), and mean local and regional 
### dissimilarities. The output is a dataframe with ascending 
### diversity order as rows and diversity metric as columns.
# =========================================

noHier <- function(comm, q.seq, pairwise = pairwise) {
  Diversity.seq <- sapply(q.seq, divPart, abundances = comm)
  gamma <- unlist(Diversity.seq["gamma", ])
  alpha <- unlist(Diversity.seq["alpha", ])
  mBeta <- gamma/alpha
  
  if(pairwise == FALSE) {
    commDiv <- data.frame(alpha = alpha, mBeta = mBeta, 
                          gamma = gamma)
  }
  
  if(pairwise == TRUE) {
    if (nrow(comm) == 2) {
      pbeta <- mBeta - 1
      
    } else {
      target <- combn(1:nrow(comm), 2)
      pairs <- vector("list", length = ncol(target))
      
      for (i in 1:ncol(target)) {
        pairs[[i]] <- rbind(comm[target[1, i], ], 
                            comm[target[2, i], ])
      }
      
      betaPairs <- sapply(pairs, pairwiseFunc, q.seq = q.seq)
      pbeta <- do.call("cbind", betaPairs["beta",]) - 1
      
    }
    commDiv <- data.frame(pBeta = t(pbeta))
    nms <- combn(rownames(comm),2)
    rownames(commDiv) <- apply(nms,2,paste, collapse="_")
  }
  return(commDiv)
}


# =========================================
### patchDiversity calculates N-community and pairwise
### diversities (if pairwise = true) for a given sequence
### of q within each group of the hierarchy It requires an
### nObs x nSpp matrix of abundances, a vector indicating
### which rows correspond to which group, and a sequence
### of diveristy orders. For each patch, it outputs alpha,
### multiplicative beta, gamma, and if pairwise = true,
### mean pairwise beta (as turnover) as a dataframe with ascending
### diversity order as rows and diversity metric as
### columns.
# =========================================

patchDiversity <- function(comm, hier = hier, q.seq = q.seq, pairwise = pairwise) {
  patchNum <- unique(as.numeric(hier))
  commDat <- vector("list", max(patchNum))
  for (p in 1:max(patchNum)) {
    patchComm <- comm[as.numeric(hier) == patchNum[p], ]
    commDat[[p]] <- patchComm
  }
  
  patchList <- lapply(commDat, FUN = noHier, q.seq = q.seq,
                      pairwise = pairwise)
  patchDiv <- do.call(rbind, patchList)
  
  if(pairwise == FALSE) {
    names(patchDiv) <- paste(c("alphaP", "mBetaP", "gammaP"),
                             rep(patchNum, each = 3), sep = "")
  } else {
    colnames(patchDiv) <- paste("q", q.seq, sep="")
  }
  
  return(patchDiv)
}


# =========================================
### RegionDiversity calculates pairwise and N-community
### diversities for a given sequence of q for the highest
### level of the hierarchy. Measures of pairwise
### beta/dissimilarity are among-group. It requires an
### nObs x nSpp matrix of abundances, a vector indicating
### which rows correspond to which group, and a sequence
### of diveristy orders. For each patch, it outputs alpha,
### multiplicative beta, gamma, and if pairwise = TRUE,
### mean pairwise beta (as turnover) and mean local and
### regional dissimilarities. The output returns a
### dataframe with ascending diversity order as rows and
### diversity metric as columns.
# =========================================

regionDiversity <- function(comm, hier = hier, q.seq = q.seq,
                            pairwise = pairwise) {
  Diversity.seq <- sapply(q.seq, divPart, abundances = comm)
  gammaR <- unlist(Diversity.seq["gamma", ])
  alphaR <- unlist(Diversity.seq["alpha", ])
  betaR <- gammaR / alphaR
  
  if (pairwise == FALSE) {
    regionDiv <- data.frame(alphaR = alphaR, mBetaR = betaR, 
                            gammaR = gammaR)
  } else {
    #------Regional Beta------------
    nObs <- nrow(comm)
    a <- rep(hier, nObs)
    b <- rep(unique(hier), times = as.vector(table(hier)) * nObs)
    combs <- expand.grid(1:nObs, 1:nObs)
    target <- subset(combs, Var1 < Var2 & a != b)
    pairs <- vector("list", length = nrow(target))
    for (i in 1:nrow(target)) {
      pairs[[i]] <- rbind(comm[target$Var1[i], ], 
                          comm[target$Var2[i], ])
    }
    betaPairs <- sapply(pairs, pairwiseFunc, q.seq = q.seq)
    
    pbeta <- rowMeans(do.call("cbind", betaPairs["beta", ])) - 1
    loc <- rowMeans(do.call("cbind", betaPairs["local", ]))
    reg <- rowMeans(do.call("cbind", betaPairs["region", ]))
    
    regionDiv <-
      data.frame(alphaR = alphaR, mBetaR = betaR, 
                 gammaR = gammaR, pBetaR = pbeta, lDiffR = loc, 
                 rDiffR = reg)
  }
  
  return(regionDiv)
}

noHierGeo <- function(geo) {
  if (nrow(geo) == 2) {
    geoDist <- data.frame(distance = distVincentyEllipsoid(geo))
    rownames(geoDist) <- paste(rownames(geo), collapse="_")
    
  } else {
    target <- combn(1:nrow(geo), 2)
    pairs <- vector("list", length = ncol(target))
    
    for (i in 1:ncol(target)) {
      pairs[[i]] <- rbind(geo[target[1, i], ], 
                          geo[target[2, i], ])
    }
    
    geoPairs <- sapply(pairs, distVincentyEllipsoid)
    geoDist <- data.frame(distance = geoPairs)    
    nms <- combn(rownames(geo),2)
    rownames(geoDist) <- apply(nms,2,paste, collapse="_")
  }
  
  return(geoDist)
}


# =========================================
### patchDistances calcualates pairwise geographic distances among sites. It requires longitude and latitude, as well as a vector of groupings called hier. 
# =========================================

patchDistances <- function(geo=geo, hier = hier) {
  patchNum <- unique(as.numeric(hier))
  geoData <- vector("list", max(patchNum))
  for (p in 1:max(patchNum)) {
    patchLocation <- geo[as.numeric(hier) == patchNum[p], ]
    geoData[[p]] <- patchLocation
  }
  
  patchList <- lapply(geoData, FUN = noHierGeo)
  patchDist <- do.call(rbind, patchList)
  
  
  
  return(patchDist)
}

############################
#This function calcualates pairwise turnover for a series of values of q, and calculates geographic distances. It requires a community matrix, a set of geographic coordinates, a vector of q's desired, and a vector of land uses.
############################
microDivNFH <- function(comm=comm, locations=locations, q.seq=q.seq, pairwise=TRUE, landUse = landuse) {
  # get community set up
  comm <- comm[grep("Hub|Near|Far", rownames(comm)), ]
  comm <- comm[order(rownames(comm)),]
  
  # set up plotIDs for patch partitioning wankery
  plots <- gsub("(.*)\\-.*-.*", "\\1" , rownames(comm))
  
  # get geographic locations set up
  geo <- locations[match(plotName(comm), locations$sID), c("sID", "fLong", "fLat")]
  geoDat <- geo[, c("fLong", "fLat")]
  rownames(geoDat) <- geo$sID
  
  # calculate within-patch geographical distances
  geoDistances <- patchDistances(geo=geoDat, hier=as.factor(plots))
  
  # get plot names into format to match up with other data 
  # if necessary
  matchPlots <- gsub("(.*\\w)_.*\\w$", "\\1" , rownames(geoDistances))
  
  # like landuse
  land <- landuse$Landuse[match(matchPlots, landuse$PlotName)]
  
  
  # calculate within-patch pairwise turnover
  div <- patchDiversity(comm = comm, hier = as.factor(plots), q.seq=q.seq, pairwise=pairwise)
  
  # smoosh them all together into one dataframe  
  out <- cbind(div, geoDistances, land, plots=gsub("(.*)\\w", "\\1" , matchPlots))
  
  return(out)
}




# =========================================
### This function is a wrapper to calculate the
### diversities for communities with or without
### hierarchical structure. It requires X which can be
### either a list of a community matrix and a hierarchical
### structure, or just a community matrix. The diveristies
### of the hiersearchy can be specified (i.g., just patch or
### just region or both). It requires a sequence of
### diversity orders (q), and can calculate diversities on
### absolute abundances or relative abundances. If pairwise 
### is true, divHier will also return pairwise turnover and 
### local and regional differentiation.
# =========================================


divHier <- function(X, hierLev="none", q.seq=seq(0,2,0.1), 
                    hier=NULL, relAb=NULL, pairwise = pairwise) {
  
  if (is.list(X) == TRUE) {
    comm <- X$comm
    
    if (hierLev != "none") {
      hier <- X$hier
    }
  } else {
    comm <- X
    if (!is.null(hier)) {
      hier <- hier
    }
  }
  
  if (!is.null(relAb)) {
    comm <- sweep(x = comm, MARGIN = 1, STATS = rowSums(comm), 
                  FUN = "/")
  }
  
  if(hierLev == "none") {
    Div <- noHier(comm=comm, q.seq=q.seq, pairwise=pairwise)
    return(Div)
  }
  
  if(hierLev == "region"){
    regionDiv <- regionDiversity(comm=comm, hier=hier, 
                                 q.seq=q.seq, pairwise=pairwise)
    return(regionDiv)
  }
  
  if(hierLev == "patch"){ 
    patchDiv <- patchDiversity(comm=comm, hier=hier, 
                               q.seq=q.seq, pairwise=pairwise)
    return(patchDiv)
  } 
  
  if(hierLev == "both"){
    regionDiv <- regionDiversity(comm=comm, hier=hier, 
                                 q.seq=q.seq, pairwise=pairwise)
    patchDiv <- patchDiversity(comm=comm, hier=hier, 
                               q.seq=q.seq, pairwise=pairwise)
    allDiv <- cbind(regionDiv, patchDiv)
    
    return(allDiv)
  }  	
}  


# =========================================
### diversify is a wrapper to calculate the
### diversities for communities with or without
### hierarchical structure for each iteration of the
### Bayesian posterior. It takes, as input, a list object
### outputted from thetaClean where each element of the
### list is a post-warmup community abundance matrix.
### Additionally, it requires the number of cores to use
### in parallel. Defaults to 4. The output of this
### function is a list where each element is a dataframe
### with order q in ascending order as rows and diversity
### metrics as columns. The list is the length of the
### number of post-warmup iterations.
# =========================================

diversify <- function(X, hierLev="none", q.seq=seq(0,2,0.1), 
                      cores = 4, relAb=NULL, pairwise=TRUE) {  
  
  divBayes <- mclapply(X, FUN=divHier, hierLev=hierLev, 
                       q.seq=q.seq, relAb=relAb, mc.cores=cores, 
                       pairwise=pairwise)
  
  return(divBayes)
}

# =========================================
### This function calculates summary statistics as needed
### across the giant list created in the source code.
# =========================================

bayesQuant <- function(postList, trueDiv) {
  
  quant <- function(divList) {
    postMat <- as.data.frame(do.call(rbind, divList$bayes))
    distr <- sapply(postMat, ecdf)
    divQuant <- mapply(distr, FUN=function(x, q) x(q), q=divList$true)
    return(divQuant)
  }
  
  bayesArray <- simplify2array(postList)
  temp <- vector("list", nrow(bayesArray))
  
  for(i in 1:nrow(bayesArray)) {
    temp[[i]] <- list(bayes =bayesArray[i,], true=trueDiv[,i])
  }
  
  outs <- sapply(temp, quant)
  colnames(outs) <- rownames(bayesArray) 
  return(outs) 
} 


sumStats <- function(X, q, stats = "mean", probs = NULL) {
  bayesArray <- simplify2array(X)
  temp <- vector("list", nrow(bayesArray))
  
  for(i in 1:nrow(bayesArray)) {
    temp[[i]] <- bayesArray[i,]
  }
  
  if(stats == "mean") {
    outs <- sapply(temp, FUN=function(x) colMeans(do.call(rbind,x)))
    colnames(outs) <- rownames(bayesArray) 
    return(outs) 
  }	
  
  if(stats == "min") {
    outs <- sapply(temp, FUN = function(x) sapply(as.data.frame(do.call(rbind,x)),min))
    colnames(outs) <- rownames(bayesArray) 
    return(outs)  	
  }
  
  if(stats == "max") {
    outs <- sapply(temp, FUN = function(x) sapply(as.data.frame(do.call(rbind,x)),max))
    colnames(outs) <- rownames(bayesArray) 
    return(outs) 
  }
  
  if(stats == "quantile") {
    outs <- sapply(temp, FUN = function(x) sapply(as.data.frame(do.call(rbind,x)),quantile,probs))
    colnames(outs) <- rownames(bayesArray) 
    return(outs) 
  }
}


###############
# For use with Entropart
##############


pairwiseFunc1 <- function(pairs, q.seq=q.seq, Correction=Correction, Biased=Biased){
  MC <- MetaCommunity(pairs)
  divProf <- DivProfile(MC=MC, q.seq=q.seq, Biased=Biased, Correction=Correction)[c("GammaDiversity","TotalAlphaDiversity")]
  beta <- divProf$GammaDiversity/divProf$TotalAlphaDiversity
  #beta <- ifelse(beta < 1, 1, beta)
  #beta <- ifelse(beta > 2, 2, beta)
  loc <- mapply(local,q=q.seq,beta=beta,N=2)
  reg <- mapply(regional,q=q.seq,beta=beta,N=2)
  # return(beta)
  data.frame(beta=beta,local=loc,region=reg)
}

noHier1 <- function(comm, q.seq=q.seq, Correction=Correction, 
                    Biased=Biased, cores=cores, pairwise=pairwise){
  MC <- MetaCommunity(comm)
  divProf <- DivProfile(MC=MC, q.seq=q.seq, Biased=Biased, 
                        Correction=Correction)[c("GammaDiversity", 
                                                 "TotalAlphaDiversity")]
  
  gamma <- divProf$GammaDiversity
  alpha <- divProf$TotalAlphaDiversity
  mBeta <- gamma/alpha
  
  if (pairwise == FALSE) {
    commDiv <- data.frame(alpha = alpha, mBeta = mBeta, 
                          gamma = gamma)
    return(commDiv)
  } else { 
    # -----if pairwise = TRUE------
    target <- combn(1:ncol(comm),2)
    pairs <- vector("list",length=ncol(target))
    
    for(i in 1:ncol(target)){
      pairs[[i]] <- cbind(comm[,target[1,i]],
                          comm[,target[2,i]])
    }
    
    betaPairs <- mclapply(pairs, FUN=pairwiseFunc1, 
                          q.seq=q.seq, mc.cores=cores, Correction=Correction, 
                          Biased=Biased)
    betaPairs <- simplify2array(betaPairs)
    
    pbeta <- rowMeans(do.call("cbind", betaPairs["beta",]), 
                      na.rm=TRUE)-1
    loc <- rowMeans(do.call("cbind", betaPairs["local",]), 
                    na.rm=TRUE)
    reg <- rowMeans(do.call("cbind", betaPairs["region",]), 
                    na.rm=TRUE)
    
    commDiv <- data.frame(alpha=alpha, mBeta = mBeta, 
                          gamma=gamma, pBeta = pbeta, lOver = loc, rOver = reg)
    return(commDiv)
  }
}

patchDiversity1 <- function(comm, hier=hier, q.seq=q.seq, 
                            Correction=Correction, Biased=Biased, cores=cores, 
                            pairwise=pairwise) {
  
  patchNum <- unique(hier)
  commDat <- vector("list", max(patchNum)) 
  for(p in 1:max(patchNum)){
    patchComm <- comm[,hier == patchNum[p]]
    commDat[[p]] <- patchComm  	  
  } 	
  
  patchList <- lapply(commDat, FUN=noHier1, q.seq=q.seq, 
                      Correction=Correction, Biased=Biased, cores=cores, 
                      pairwise=pairwise)
  patchDiv <- do.call(cbind, patchList)
  
  if(pairwise == FALSE) {
    names(patchDiv) <- paste(c("alphaP", "mBetaP", "gammaP"),
                             rep(patchNum, each = 3), sep = "")
  } else {
    names(patchDiv) <- paste(c("alphaP", "mBetaP", "gammaP", 
                               "pBetaP", "lDiffP", "rDiffP"), rep(patchNum, each = 6), 
                             sep = "") 
  }
  
  return(patchDiv)
}

regionDiversity1 <- function(comm, hier=hier, q.seq=q.seq, 
                             Correction=Correction, Biased=Biased, cores=cores, 
                             pairwise=pairwise) {
  
  MC <- MetaCommunity(comm)
  divProf <- DivProfile(MC=MC, q.seq=q.seq, Biased=Biased, 
                        Correction=Correction)[c("GammaDiversity",
                                                 "TotalAlphaDiversity")]
  
  gamma <- divProf$GammaDiversity
  alpha <- divProf$TotalAlphaDiversity
  mBeta <- gamma/alpha
  
  if (pairwise == FALSE) {
    regionDiv <- data.frame(alphaR = alpha, mBetaR = mBeta, 
                            gammaR = gamma)
  } else {
    #------Regional pairwise Beta------------
    nObs <- ncol(comm)
    a <- rep(hier, nObs)
    b <- rep(unique(hier), times = as.vector(table(hier)) * nObs)
    combs <- expand.grid(1:nObs,1:nObs)
    target <- subset(combs, Var1 < Var2 & a != b)
    pairs <- vector("list",length=nrow(target))
    
    for(i in 1:nrow(target)){
      pairs[[i]] <- cbind(comm[,target$Var1[i]],
                          comm[,target$Var2[i]])
    }
    
    betaPairs <- mclapply(pairs, FUN=pairwiseFunc1, q.seq=q.seq, 
                          mc.cores=cores, Correction=Correction, Biased=Biased)
    betaPairs <- simplify2array(betaPairs)
    
    pbeta <- rowMeans(do.call("cbind",betaPairs["beta",]), 
                      na.rm=TRUE) - 1
    loc <- rowMeans(do.call("cbind",betaPairs["local",]), 
                    na.rm=TRUE)
    reg <- rowMeans(do.call("cbind",betaPairs["region",]), 
                    na.rm=TRUE)
    
    regionDiv <- data.frame(alphaR = alpha, mBetaR = mBeta, 
                            gammaR = gamma,  pBetaR = pbeta, 
                            lDiffR = loc, rDiffR = reg) 
  }
  
  return(regionDiv)
}

divHier1 <- function(X, hierLev = "none", q.seq=seq(0,2,0.1),
                     Correction="None", Biased=FALSE, hier=NULL, cores=cores, 
                     pairwise=pairwise) {
  
  if (is.list(X) == TRUE) {
    comm <- t(X$comm)
    
    if (hierLev != "none") {
      hier <- X$hier
    }
  } else {
    comm <- t(X)
    if (!is.null(hier)) {
      hier <- hier
    }
  }
  
  if(hierLev == "none") {
    Div <- noHier1(comm=comm, q.seq=q.seq, Correction=Correction,
                   Biased=Biased, cores=cores, pairwise=pairwise)
    return(Div)
  }
  
  if(hierLev == "region"){
    regionDiv <- regionDiversity1(comm=comm, hier=hier, 
                                  q.seq=q.seq, Correction=Correction, Biased=Biased, 
                                  cores=cores, pairwise=pairwise)
    return(regionDiv)
  }
  
  if(hierLev == "patch"){ 
    patchDiv <- patchDiversity1(comm=comm, hier=hier, 
                                q.seq=q.seq, Correction=Correction, Biased=Biased, 
                                cores=cores, pairwise=pairwise)
    return(patchDiv)
  } 
  
  if(hierLev == "both"){
    regionDiv <- regionDiversity1(comm=comm, hier=hier, 
                                  q.seq=q.seq, Correction=Correction, Biased=Biased, 
                                  cores=cores, pairwise=pairwise)
    
    patchDiv <- patchDiversity1(comm=comm, hier=hier, 
                                q.seq=q.seq, Correction=Correction, Biased=Biased, 
                                cores=cores, pairwise=pairwise)
    
    allDiv <- cbind(regionDiv, patchDiv)
    
    return(allDiv)
  }  	
}  	 



