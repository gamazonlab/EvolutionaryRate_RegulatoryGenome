library(ggplot2)
library(reshape)
library(patchwork)
require(plotrix)

load("d://3_EXPRESSION_SELECTION/r2_GBM.RData")

# Figure 4a
f4.1 = ggplot(data.frame(r2), aes(r2)) + 
    geom_histogram(data = data.frame(r2), fill = "darkgreen", position="dodge", binwidth=0.01, alpha=0.5) + 
    xlab(expression("Adjusted R"^2)) + xlim(c(0.5, 0.7)) + ylab(expression("Count")) +
    theme(axis.text=element_text(face = "bold"), axis.title=element_text(face = "bold", size=22), text=element_text(size=22), panel.background = element_rect(fill="white"), panel.grid.major = element_line(color="white"), panel.grid.minor = element_line(color="white") ) + theme( axis.line = element_line(colour = "grey90", size = 1, linetype = "solid")) + labs(tag = "A")


# Figure 4b
mm.4 = read.table('j://3_EXPRESSION_SELECTION/InfluenceScore.100Perms.xls', header=F)
tmp.0 = mm.4[,c(1,2)]
colnames(tmp.0) = c("Variable", "Influence Score")
tmp.00 = mm.4[, c(3,4)]
colnames(tmp.00) = c("Variable", "Influence Score")
mm.4 = rbind(tmp.0, tmp.00)
colnames(mm.4) = c("Variable", "Influence.Score")

f4.2 = ggplot(data.frame(mm.4), aes(x=Influence.Score)) + 
   geom_histogram(data = data.frame(subset(mm.4, mm.4$Variable=="MaxVarianceRPKM")), fill = "darkgreen",  alpha=0.5)  + 
   geom_histogram(data = data.frame(subset(mm.4, mm.4$Variable=="MaxMeanRPKM")), fill = "darkgreen",  alpha=0.5)  + 
   theme(legend.position="top", plot.title = element_text(hjust = 0.5),  axis.title=element_text( size=22), text=element_text(size=22), panel.background = element_rect(fill="white"), panel.grid.major = element_line(color="white"), panel.grid.minor = element_line(color="white"))  + xlab("Influence Score") + ylab("Count") + theme( axis.line = element_line(colour = "grey90", size = 1, linetype = "solid")) + xlim(c(30, 60))  + labs(tag = "B") # + scale_fill_discrete(breaks=c("Mean","Variance")) 


# Figure 4c
evo_f = read.table("c://Users/egamazon/Dropbox/Papers/SelectionExpression/Paper_v1/dropbox_evolution_paper/EvoStats_Final2.txt", header=T)
evo=evo_f
protein = subset(evo_f, evo_f$Type=="protein_coding")
mend = subset(protein, protein$Mendelian==1)
ess = read.table("c://Users/egamazon/Dropbox/Papers/SelectionExpression/Paper_v1/dropbox_evolution_paper/essentialGenes.txt", header=T,sep="\t")
ess.1 = subset(evo, evo$GeneID %in% ess$Ensembl.gene.ID)

res = c()
for (i in 1:1000)
{
	tmp = protein[sample(1:nrow(protein), dim(mend)[1], replace=F),]
	r = cor.test(tmp$Mouse_dN.dS, tmp$MaxVarianceRPKM, method="spearman")
	res = rbind(res, r$estimate)
}

res1 = c()
for (i in 1:1000)
{
	tmp = protein[sample(1:nrow(protein), dim(ess.1)[1], replace=F),]
	r = cor.test(tmp$Mouse_dN.dS, tmp$MaxVarianceRPKM, method="spearman")
	res1 = rbind(res1, r$estimate)
}

   # observed values included
p1 = ggplot(data.frame(res), aes(rho, fill="darkgreen")) + 
    geom_histogram(data = data.frame(res), fill = "darkgreen", position="dodge", binwidth=0.01, alpha=0.5) + 
    xlab("") + ylab("Count") +
    theme(plot.title = element_text(hjust = 0.5, size=18), axis.title=element_text( size=22), text=element_text(size=22), panel.background = element_rect(fill="white"), panel.grid.major = element_line(color="white"), panel.grid.minor = element_line(color="white")) +
   geom_vline(xintercept = -0.070, col="orange", size=3) + labs(title="Mendelian Genes") + theme( axis.line = element_line(colour = "grey90", size = 1, linetype = "solid")) + labs(tag = "C")
p2 = ggplot(data.frame(res1), aes(rho, fill="darkgreen")) + 
    geom_histogram(data = data.frame(res1), fill = "darkgreen", position="dodge", binwidth=0.01, alpha=0.5) + 
    xlab(expression(paste("Spearman's ", rho, sep=""))) + ylab("Count") +
    theme(plot.title = element_text(hjust = 0.5, size=18), axis.title=element_text( size=22), text=element_text(size=22), panel.background = element_rect(fill="white"), panel.grid.major = element_line(color="white"), panel.grid.minor = element_line(color="white")) +
   geom_vline(xintercept = -0.076, col="orange", size=3) + labs(title="Essential Genes") + theme( axis.line = element_line(colour = "grey90", size = 1, linetype = "solid")) 


# PLOT
f4.1 + f4.2 + p1/p2

