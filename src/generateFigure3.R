library(ggplot2)
library(reshape)
library(patchwork)
require(plotrix)
library(gbm)

evo_f = read.table("c://Users/egamazon/Dropbox/Papers/SelectionExpression/Paper_v1/dropbox_evolution_paper/EvoStats_Final2.txt", header=T)
protein = subset(evo_f, evo_f$Type=="protein_coding")
protein.copy = protein
xx = protein.copy[,c(9,11:13,17:19,22, 23)]  
colnames(xx) = c("Gene Length", "Expression Breadth", "Maximum Expression Level", "MaxVariance", "Mouse dN", "Mouse  dS", "Mouse dN/dS", "DGN Heritability", "No of Interactions")
corr1 = cor(xx, use="complete.obs", method="spearman")
corr1.m <- melt(corr1, id=rownames(corr1))
colnames(corr1.m) = c("Feature1", "Feature2" , "Value")
 

# Figure 3a
fig3.1 = ggplot(corr1.m, aes(Feature2, Feature1, fill= Value)) + 
  geom_tile() +
scale_fill_gradient(low="white", high="blue") +
  theme(text = element_text(size=18),axis.title=element_text( size=15), legend.title=element_text(size=17), panel.background = element_rect(fill="white", color="black"), panel.grid.major = element_line(color="white"), panel.grid.minor = element_line(color="white"), legend.key = element_rect(fill="white", size=.5, linetype="dotted"), plot.title = element_text( size=15)) +  
theme(axis.text.x=element_text(angle=15,hjust=1)) +
  scale_fill_gradient2() + xlab("Feature") + ylab("Feature")  + labs(tag = "A") 


# Figure 3b
load("c://Users/egamazon/Dropbox/Papers/SelectionExpression/Paper_v1/dropbox_evolution_paper/boostedgradientmodels.Rdata")

f1d.500 = ggplot(summary(boost.evo_1d_500trees), aes(x= reorder(var, rel.inf),y= rel.inf)) + geom_bar(stat = "identity", fill="darkgreen", col="black") + coord_flip() + theme(plot.title = element_text(hjust = 0.5, size=22), panel.background = element_rect(fill="white", color="black"), panel.grid.major = element_line(color="white"), panel.grid.minor = element_line(color="white"), axis.title = element_text(size = 17), axis.text = element_text(size = 15)) + ylab("") + xlab("Depth = 1") + labs(title="500 Trees")  +  labs(tag = "B")
f1d.1000 = ggplot(summary(boost.evo_1d_1000trees), aes(x= reorder(var, rel.inf),y= rel.inf)) + geom_bar(stat = "identity", fill="darkgreen", col="black") + coord_flip() + theme(plot.title = element_text(hjust = 0.5, size=22), panel.background = element_rect(fill="white", color="black"), panel.grid.major = element_line(color="white"), panel.grid.minor = element_line(color="white"), axis.title = element_text(size = 17), axis.text = element_text(size = 15)) + ylab("")  + xlab("") + labs(title="1,000 Trees") +  labs(tag = "C")
f1d.10k = ggplot(summary(boost.evo_1d), aes(x= reorder(var, rel.inf),y= rel.inf)) + geom_bar(stat = "identity", fill="darkgreen", col="black") + coord_flip() + theme(plot.title = element_text(hjust = 0.5, size=22), panel.background = element_rect(fill="white", color="black"), panel.grid.major = element_line(color="white"), panel.grid.minor = element_line(color="white"), axis.title = element_text(size = 17), axis.text = element_text(size = 15)) + ylab("") + xlab("") + labs(title="10,000 Trees") +  labs(tag = "D")

f4d.500 = ggplot(summary(boost.evo_4d_500trees), aes(x= reorder(var, rel.inf),y= rel.inf)) + geom_bar(stat = "identity", fill="darkgreen", col="black") + coord_flip() + theme(panel.background = element_rect(fill="white", color="black"), panel.grid.major = element_line(color="white"), panel.grid.minor = element_line(color="white"), axis.title = element_text(size = 17), axis.text = element_text(size = 15)) + ylab("")  + xlab("Depth = 4") +  labs(tag = "E") 
f4d.1000 = ggplot(summary(boost.evo_4d_1000trees), aes(x= reorder(var, rel.inf),y= rel.inf)) + geom_bar(stat = "identity", fill="darkgreen", col="black") + coord_flip() + theme(panel.background = element_rect(fill="white", color="black"), panel.grid.major = element_line(color="white"), panel.grid.minor = element_line(color="white"), axis.title = element_text(size = 17), axis.text = element_text(size = 15)) + ylab("")  + xlab("") +  labs(tag = "F")
f4d.10k = ggplot(summary(boost.evo_4d), aes(x= reorder(var, rel.inf),y= rel.inf)) + geom_bar(stat = "identity", fill="darkgreen", col="black") + coord_flip() + theme(panel.background = element_rect(fill="white", color="black"), panel.grid.major = element_line(color="white"), panel.grid.minor = element_line(color="white"), axis.title = element_text(size = 17), axis.text = element_text(size = 15)) + ylab("") + xlab("") +  labs(tag = "G")

f5d.500 = ggplot(summary(boost.evo_5d_500trees), aes(x= reorder(var, rel.inf),y= rel.inf)) + geom_bar(stat = "identity", fill="darkgreen", col="black") + coord_flip() + theme(panel.background = element_rect(fill="white", color="black"), panel.grid.major = element_line(color="white"), panel.grid.minor = element_line(color="white"), axis.title = element_text(size = 17), axis.text = element_text(size = 15)) + ylab("") + xlab("Depth = 5") +  labs(tag = "H")
f5d.1000 = ggplot(summary(boost.evo_5d_1000trees), aes(x= reorder(var, rel.inf),y= rel.inf)) + geom_bar(stat = "identity", fill="darkgreen", col="black") + coord_flip() + theme(panel.background = element_rect(fill="white", color="black"), panel.grid.major = element_line(color="white"), panel.grid.minor = element_line(color="white"), axis.title = element_text(size = 22), axis.text = element_text(size = 15)) + ylab("Relative Influence") + xlab("") +  labs(tag = "I")
f5d.10k = ggplot(summary(boost.evo_5d), aes(x= reorder(var, rel.inf),y= rel.inf)) + geom_bar(stat = "identity", fill="darkgreen", col="black") + coord_flip() + theme(panel.background = element_rect(fill="white", color="black"), panel.grid.major = element_line(color="white"), panel.grid.minor = element_line(color="white"), axis.title = element_text(size = 17), axis.text = element_text(size = 15)) + ylab("") + xlab("") +  labs(tag = "J")

#f1d.500 + f1d.1000 + f1d.10k + f4d.500 + f4d.1000 + f4d.10k + f5d.500 + f5d.1000 + f5d.10k  + plot_layout(ncol=3)
fig3.2 = f1d.500 + f1d.1000 + f1d.10k + f4d.500 + f4d.1000 + f4d.10k + f5d.500 + f5d.1000 + f5d.10k  + plot_layout(ncol=3) 

#( fig3.1 + plot_spacer() + plot_spacer() ) / (f1d.500 + f1d.1000 + f1d.10k + f4d.500 + f4d.1000 + f4d.10k + f5d.500 + f5d.1000 + f5d.10k  + plot_layout(ncol=3) ) + plot_layout(widths = unit(c(30, 30), c('cm', 'null')), heights = unit(c(1, 5), c('cm', 'null')))

( fig3.1 + plot_spacer() + plot_spacer() ) / (f1d.500 + f1d.1000 + f1d.10k + f4d.500 + f4d.1000 + f4d.10k + f5d.500 + f5d.1000 + f5d.10k  + plot_layout(ncol=3) ) + plot_layout( heights = unit(c(3.5, 3), c('cm', 'null')))
 
