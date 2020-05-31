getPeerCorrelation <- function(evofile, peerfile, outfile) {
    evodata <- read.table(evofile, header = TRUE, sep = "\t")
    peerdata <- read.table(peerfile, header = TRUE, sep = "\t", row.names = 1)
    proteins <- evodata[which(evodata$Type == "protein_coding"),]
    proteins <- proteins[order(proteins$GeneID),]
    proteins <- proteins[,c(1,14:19)]
    peerdata <- peerdata[order(rownames(peerdata)),]
    peerdata <- peerdata[which(rownames(peerdata) %in% proteins[,1]),]
    proteins <- proteins[which(proteins[,1] %in% rownames(peerdata)),]
    sdpeer <- apply(peerdata,1,sd)
    sdcor <- cor.test(proteins$Mouse_dN.dS,sdpeer, method = c("spearman"))
    outmat <- matrix(nrow=1, ncol=2)
    outmat[1,1] <- peerfile
    outmat[1,2] <- sdcor$estimate
    write.table(outmat,file=outfile, quote = FALSE, row.names = FALSE, col.names = FALSE, sep = "\t")
  }