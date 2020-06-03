library(ggplot2)
library(reshape)
library(patchwork)
require(plotrix)

setwd("d://3_EXPRESSION_SELECTION/")
evo_f = read.table("c://Users/egamazon/Dropbox/Papers/SelectionExpression/Paper_v1/dropbox_evolution_paper/EvoStats_Final2.txt", header=T)
 
# Fig 2a
a = read.table('FUNCTIONAL_CATEGORIES.txt', header=T, sep="\t")

f2.1 = ggplot(data.frame(a), aes(x=Term, y=Fold_Enrichment) ) + geom_col(aes(fill=Tissue), col="black") +  
xlab(expression("Functional Category")) + ylab(expression("Fold Enrichment")) +
    theme(legend.position="top", plot.title = element_text(hjust = 0.5), axis.text.y=element_text( size=12),axis.text.x=element_text( size=22), axis.title=element_text( size=22), text=element_text(size=22), panel.background = element_rect(fill="white"), panel.grid.major = element_line(color="white"), panel.grid.minor = element_line(color="white")) +
     coord_flip()  + theme( axis.line = element_line(colour = "grey90", size = 1, linetype = "solid")) +
labs(tag = "A")


# Fig 2b
tauvals <- read.table("c://Users/egamazon/Dropbox/Papers/SelectionExpression/Paper_v1/dropbox_evolution_paper/FinalSubmission/Resubmission_CellReports/Cell_Reports_files/REVISION/PEERJ_SUBMISSION/Final_Submission/REVISED_FIGURES/figure2bdata.txt", sep = "\t", header = TRUE, row.names = 1)
colnames(tauvals) <- c(5,10,15,20,25,30,35,40)
tauvalsmelt <- melt(tauvals)

f2.2 = ggplot(tauvalsmelt, aes(x=variable, y=value)) + 
geom_boxplot(fill="green", color="black") +
    theme(legend.position="none", plot.title = element_text(hjust = 0.5), axis.text=element_text( size=22), axis.title=element_text( size=22), text=element_text(size=22), panel.background = element_rect(fill="white"), panel.grid.major = element_line(color="white"), panel.grid.minor = element_line(color="white"))  + labs(tag = "C") +  theme( axis.line = element_line(colour = "grey90", size = 1, linetype = "solid")) + ylab(expression("h"^2)) +
xlab("Number of Tissues Sampled") + ylab(expression(paste("S.D. of ", tau, " "))) + 
labs(tag = "B")


# Fig 2c
setwd("j://OLD_COMPUTER_I/2_EXPRESSION_SELECTION/")
h2 = read.table('../1_GTEX_h2/OUT_h2.xls', header=F)
h2.wb = subset(h2, h2$V3=="WholeBlood_TW")
h2.ms = subset(h2, h2$V3=="Muscle-Skeletal_TW")
h2.DGN = subset(h2, h2$V3=="DGN-WB")

protein = subset(evo_f, evo_f$Type=="protein_coding")

protein.copy = protein
protein.copy.1 = merge(protein.copy, h2.ms[,c(1,2,4,5)], by.x = c(3), by.y= c(1))
protein.copy.2 = merge(protein.copy.1, h2.DGN[,c(1,2,4,5)], by.x = c(1), by.y= c(1))
xx = protein.copy.2[,c(9,11:16,19,20,33)]  # changed to 33 from 30 before, changed 18 to 19, 25 to 20  July 5, 2017
colnames(xx) = c("gene length", "expression breadth", "expression level", "dN", "dS", "dN/dS", "branch", "node degree", "error-free misfolding", "heritability")

yy = subset(xx, !is.na(xx[,6]))
yy$CONSERVED = "NA"
for (i in 1:dim(yy)[1])
{
	if (yy[i,6] < 0.0615)
	{
		yy[i,11] = 'Conserved'

	}
	if (yy[i,6] > 0.80)

	{
		yy[i,11] = 'Fast-evolving'
	}
}

f2.3 = ggplot(yy[yy$CONSERVED!="NA",], aes(x=CONSERVED, y=heritability)) + geom_violin(draw_quantiles= c(0.25, 0.5, 0.75), aes(fill=CONSERVED)) + # theme(axis.title.x = element_text(size=22), axis.title.y = element_text(size=22), axis.text.x = element_text(angle=90, hjust=1, vjust=0.3, size=22), axis.text.y = element_text(angle=90, hjust=1, vjust=0.3, size=22),panel.background = element_rect(fill="white"), panel.grid.major = element_line(color="grey90"), panel.grid.minor = element_line(color="grey90")) + ylab(expression("h"^2)) + 
xlab(expression("Gene Class")) + ggtitle("") + labs(tag = "C") + theme(legend.position="none")  +  scale_fill_manual(values = c("green", "orange")) +
    theme(legend.position="none", plot.title = element_text(hjust = 0.5), axis.text=element_text( size=22), axis.title=element_text(size=22), text=element_text(size=22), panel.background = element_rect(fill="white"), panel.grid.major = element_line(color="white"), panel.grid.minor = element_line(color="white"))  + labs(tag = "C") +  theme( axis.line = element_line(colour = "grey90", size = 1, linetype = "solid")) + ylab(expression("h"^2)) +
ylim(c(0,0.30))


f2.1 | (f2.2 / f2.3)

