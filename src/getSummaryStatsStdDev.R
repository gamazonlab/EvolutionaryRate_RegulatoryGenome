getSummaryStats <- function(inputfile, outputfile) {
  fdata <- read.table(inputfile, sep = "\t", header = TRUE, stringsAsFactors = FALSE)
  resultmatrix <- matrix(nrow = ncol(fdata) - 8, ncol = 13)
  colnames(resultmatrix) <- c("Tissue", "Beta_dN", "Beta_dN_pvalue", "r_dN", "r_dN_pvalue", "Beta_dS", "Beta_dS_pvalue", "r_dS", "r_dS_pvalue", "Beta_dN/dS", "Beta_dN/dS_pvalue", "r_dN/dS", "r_dN/dS_pvalue")
  ntissue = 1
  for(i in 9:ncol(fdata)) {
    tissuestr <- strsplit(colnames(fdata)[i], split = "_Mean")
    print(tissuestr[[1]][1])
    resultmatrix[ntissue,1] <- tissuestr[[1]][1]
    if(length(which(fdata[,i] == 0))) {
      zerodev <- which(fdata[,i] == 0)
      fit_dN <- lm(fdata$Chimp_dN[-zerodev]~log2(fdata[-zerodev,i]))
      fit_dS <- lm(fdata$Chimp_dS[-zerodev]~log2(fdata[-zerodev,i]))
      fit_dNdS <- lm(fdata$Chimp_dN.dS[-zerodev]~log2(fdata[-zerodev,i]))
      beta_dN <- coef(summary(fit_dN))[2,1]
      p_dN <- coef(summary(fit_dN))[2,4]
      beta_dS <- coef(summary(fit_dS))[2,1]
      p_dS <- coef(summary(fit_dS))[2,4]
      beta_dNdS <- coef(summary(fit_dNdS))[2,1]
      p_dNdS <- coef(summary(fit_dNdS))[2,4]
      corfit_dN <- cor.test(fdata$Chimp_dN, fdata[,i], method=c("spearman"))
      corfit_dS <- cor.test(fdata$Chimp_dS, fdata[,i], method=c("spearman"))
      corfit_dNdS <- cor.test(fdata$Chimp_dN.dS, fdata[,i], method=c("spearman"))
      r_dN <- corfit_dN$estimate
      r_dS <- corfit_dS$estimate
      r_dNdS <- corfit_dNdS$estimate
      r_dN_pvalue <- corfit_dN$p.value
      r_dS_pvalue <- corfit_dS$p.value
      r_dNdS_pvalue <- corfit_dNdS$p.value
    } else {
      fit_dN <- lm(fdata$Chimp_dN~log2(fdata[,i]))
      fit_dS <- lm(fdata$Chimp_dS~log2(fdata[,i]))
      fit_dNdS <- lm(fdata$Chimp_dN.dS~log2(fdata[,i]))
      beta_dN <- coef(summary(fit_dN))[2,1]
      p_dN <- coef(summary(fit_dN))[2,4]
      beta_dS <- coef(summary(fit_dS))[2,1]
      p_dS <- coef(summary(fit_dS))[2,4]
      beta_dNdS <- coef(summary(fit_dNdS))[2,1]
      p_dNdS <- coef(summary(fit_dNdS))[2,4]
      corfit_dN <- cor.test(fdata$Chimp_dN, fdata[,i], method=c("spearman"))
      corfit_dS <- cor.test(fdata$Chimp_dS, fdata[,i], method=c("spearman"))
      corfit_dNdS <- cor.test(fdata$Chimp_dN.dS, fdata[,i], method=c("spearman"))
      r_dN <- corfit_dN$estimate
      r_dS <- corfit_dS$estimate
      r_dNdS <- corfit_dNdS$estimate
      r_dN_pvalue <- corfit_dN$p.value
      r_dS_pvalue <- corfit_dS$p.value
      r_dNdS_pvalue <- corfit_dNdS$p.value
    }
    resultmatrix[ntissue,2:ncol(resultmatrix)] <- c(beta_dN, p_dN, r_dN, r_dN_pvalue, beta_dS, p_dS, r_dS, r_dS_pvalue, beta_dNdS, p_dNdS, r_dNdS, r_dNdS_pvalue)
    ntissue <- ntissue + 1
  }
  write.table(resultmatrix, file=outputfile, quote = FALSE, col.names = TRUE, row.names = FALSE, sep="\t")
}