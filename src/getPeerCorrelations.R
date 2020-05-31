getPeerCorrelation <- function(evofile, peerfle, outfile) {
    evodata <- read.table(evofile, header = TRUE, sep = "\t")
    peerdata <- read.table(peerfile, header = TRUE, sep = "\t", row.names = 1)
    proteins <- evodata[which(evodata$Type == "protein_coding"),]
    sortprotein <- proteins[order(proteins$GeneID),]
    proteins <- proteins[,c(1,14:19)]
    peerdata <- peerdata[order(rownames(peerdata)),]
    peerdata <- peerdata[which(rownames(peerdata) %in% proteins[,1]),]
    proteins <- evodata[which(proteins[,1] %in% rownames(peerdata)),]
    
  }