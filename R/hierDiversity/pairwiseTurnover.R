#' get_pairs is a function to generate 
#' @param comm A community dataset matrix
#' @param q.seq Level(s) of q that you want to calculate
#' @return The pairwise beta estimates
pairwiseTurnover<-function(comm, q.seq) {
  Diversity.seq <- sapply(q.seq, divPart, abundances = comm)
  gamma <- unlist(Diversity.seq["gamma", ])
  alpha <- unlist(Diversity.seq["alpha", ])
  mBeta <- gamma/alpha
  
 
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
    
  return(pbeta)
}