library(ggplot2)
library(reshape)
library(patchwork)
require(plotrix)
library("pwr")


# Figure 5a
mm = read.table('j://3_EXPRESSION_SELECTION/Correlation_OutofSampleR2_EvoRate.txt', header=T)
mm$Tissue=gsub("-", " ", mm$Tissue)

f5.1 = ggplot(mm, aes(x=-log10(pvalue), y=Spearman )) + geom_point(size=4, col="darkgreen") + geom_label(aes(label=Tissue),hjust=1, vjust=1, size=4) + xlab(expression(-log[10]~(P))) + ylab(expression(paste("Spearman's ", rho, sep=""))) +
    theme(legend.position="none", plot.title = element_text(hjust = 0.5),  axis.title=element_text( size=22), text=element_text(size=22), panel.background = element_rect(fill="white"), panel.grid.major = element_line(color="white"), panel.grid.minor = element_line(color="white")) +
xlim(c(-2,18)) + geom_vline(xintercept = -log10(0.05/dim(mm)[1]), col="orange", size=2) + theme( axis.line = element_line(colour = "grey90", size = 1, linetype = "solid")) + labs(tag = "A")


# Figure 5b
r2 = read.table('j://3_BIOVU/getR2.sql.out', header=T)
r2.muscle = subset(r2, r2$tissue=="Muscle-Skeletal")

a = read.table('c://Users/egamazon/Dropbox/Papers/SelectionExpression/Paper_v1/dropbox_evolution_paper/EvoStats_final2.txt', header=T)
xx = merge(r2.muscle, a, by.x=1, by.y=3)

conserved = subset(xx, xx$Chimp_dN.dS<0.01)
dim(conserved)
fast = subset(xx, xx$Chimp_dN.dS>0.99)
dim(fast)
conserved.l = pwr.f2.test(u=1,v=999, conserved$R2/(1-conserved$R2), sig.level=0.05/dim(conserved)[1])
fast.l = pwr.f2.test(u=1,v=999, fast$R2/(1-fast$R2), sig.level=0.05/dim(fast)[1])

tmp.c = data.frame(conserved.l$power)
colnames(tmp.c) = c("Power")
tmp.f = data.frame(fast.l$power)
colnames(tmp.f) = c("Power")
erg = tmp.c
erg$geneclass = 'Conserved'
erg1 = tmp.f
erg1$geneclass = 'Fast-evolving'
erg2 = rbind(erg, erg1)


# Figure 5b
f5.2 = ggplot(erg2, aes(x = geneclass, y = Power)) + geom_violin(trim = TRUE, aes(fill=geneclass)) + geom_boxplot(width = 0.05) +
  scale_fill_manual(values = c("#00AFBB", "#E7B800")) + 
  theme(legend.position = "none")+
   theme(plot.title = element_text(hjust = 0.5),  axis.title=element_text( size=22), text=element_text(size=22), panel.background = element_rect(fill="white"), panel.grid.major = element_line(color="white"), panel.grid.minor = element_line(color="white")) + xlab("Gene Class") + theme( axis.line = element_line(colour = "grey90", size = 1, linetype = "solid")) + labs(tag = "B")


## PLOT
f5.1 + f5.2

