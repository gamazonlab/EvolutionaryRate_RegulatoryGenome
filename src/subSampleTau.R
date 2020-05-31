#subsamples a matrix of mean RPKM values across tissue and subsamples a certain number of 
#tissues and calculates tau for a specified number of replicates
#This uses taudata.txt as imput

subsampleTau <- function(matin, nsam, nreps, outputtable) {
  myresult <- matrix(nrow = nrow(matin), ncol = (nreps + 2))
  myresult[,1] <- matin[,1]
  myresult[,2] <- matin[,4]
  if(nreps == 1) {
    colnames(myresult) <- c("EnsemblID", "dN/dS", 1)
  } else {
    colnames(myresult) <- c("EnsemblID", "dN/dS", 1:nreps)
  }
  for(j in 3:(nreps+2)) {
    tsample <- sample(5:ncol(matin), nsam)
    tausamp <- matin[,tsample]
    for(i in 1:nrow(tausamp)) {
      maxrpkm <- max(tausamp[i,])
      tau <- sum(1 - (tausamp[i,]/maxrpkm))/(nsam - 1)
      myresult[i,j] <- tau
    }
  }
  write.table(myresult, file = outputtable, quote = FALSE, col.names = TRUE, row.names = FALSE, sep="\t")
}