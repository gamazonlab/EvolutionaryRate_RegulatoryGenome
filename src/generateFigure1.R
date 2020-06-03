
setwd("c://Users/egamazon/Dropbox/Papers/SelectionExpression/Paper_v1/")
setwd("dropbox_evolution_paper/")

library(ggplot2)
library(pheatmap)
library(reshape)

# devtools::install_github("thomasp85/patchwork")
library(patchwork)
require(plotrix)

#####################  FIGURE 1

m = read.table('FinalSubmission/June2018/AllTissue_COV_results.txt', header=T)

colnames(m) = sub("Spearman_rho_", "", colnames(m))
colnames(m) = sub("_", ", ", colnames(m))
colnames(m) = sub("_", ", ", colnames(m))
colnames(m) = sub("dnds", "dn/ds", colnames(m))
colnames(m) = sub("var.mean", "var/mean", colnames(m))
rownames(m) = m[,1]

masterM = data.frame(m[,c(1,2,4,6,8,10,12,14,16,18,20,22,24,26,28,30,32,34,36)])

mdata <- melt(masterM, id=c("Tissue"))
mdata$Tissue <- gsub("_", " ", mdata$Tissue)
mdata$variable <- gsub("\\.\\.", " ", mdata$variable)
mdata$variable <- gsub("\\.", "/", mdata$variable)
mdata$variable <- gsub("dn", "dN", mdata$variable)
mdata$variable <- gsub("ds", "dS", mdata$variable)
mdata$variable <- gsub("dN dS", "dN/dS", mdata$variable)


# Figure 1a
fig1.a = ggplot(mdata, aes(variable, Tissue, fill= value)) + 
  geom_tile() +
scale_fill_gradient(low="white", high="blue") +
  theme(text = element_text(size=15),axis.title=element_text( size=18), legend.title=element_text(size=22), panel.background = element_rect(fill="white", color="black"), panel.grid.major = element_line(color="white"), panel.grid.minor = element_line(color="white"), legend.key = element_rect(fill="white", size=.5, linetype="dotted"), plot.title = element_text( size=22)) +  
theme(axis.text.x=element_text(angle=30,hjust=1)) +
  scale_fill_gradient2() + labs(tag = "A") + xlab("Variable")

# Figure 1b
corvals <- read.table("c://users/egamazon/Dropbox/Papers/SelectionExpression/Paper_v1/dropbox_evolution_paper/FinalSubmission/June2018/AllTissue_COV_results.txt", sep = "\t", header = TRUE)
fig1.b = ggplot(corvals, aes(Spearman_rho_mean_dn_mouse, Spearman_rho_mean_dnds_chimp, col="dN mean")) + geom_point() + geom_point(aes(Spearman_rho_mean_dnds_mouse, Spearman_rho_mean_dnds_chimp, col="dN/dS mean")) + geom_point(aes(Spearman_rho_var_ds_mouse, Spearman_rho_mean_ds_chimp, col="dS mean")) + geom_point(aes(Spearman_rho_var.mean_dn_mouse, Spearman_rho_var.mean_dn_chimp, col="dN variance/mean")) + geom_point(aes(Spearman_rho_var.mean_dnds_mouse, Spearman_rho_var.mean_dnds_chimp, col="dN/dS variance/mean")) + geom_point(aes(Spearman_rho_var.mean_ds_mouse, Spearman_rho_var.mean_ds_chimp, col="dS variance/mean")) + theme(text = element_text(size=22),panel.background = element_rect(fill="white", color="black"), panel.grid.major = element_line(color="white"), panel.grid.minor = element_line(color="white"), legend.key = element_rect(fill="white", size=.5, linetype="dotted"),legend.title = element_text(size = 22), plot.title = element_text( size = 30)) + xlab(expression(paste("Spearman's ", rho, " (Mouse)", sep=""))) + ylab(expression(paste("Spearman's ", rho," (Chimp)", sep=""))) + scale_colour_manual(name= "Comparison", values = c("black", "red", "blue", "orange", "green", "violet")) +
 labs(tag = "B") +
 theme(text = element_text(size=22),axis.title=element_text( size=18), panel.background = element_rect(fill="white", color="black"), panel.grid.major = element_line(color="white"), panel.grid.minor = element_line(color="white"), legend.key = element_rect(fill="white", size=.5, linetype="dotted"), plot.title = element_text( size = 30)) 

# Figure 1c
peerdata <- read.table("c://Users/egamazon/Dropbox/Papers/SelectionExpression/Paper_v1/dropbox_evolution_paper/FinalSubmission/Resubmission_CellReports/Cell_Reports_files/REVISION/PEERJ_SUBMISSION/Final_Submission/REVISED_FIGURES/figure1cdata.txt", sep = "\t", header = TRUE, row.names = 1) 
fig1.c = ggplot(peerdata, aes(x=Spearman_rho_var_dnds_mouse, y=Spearman_rho_peer_var_dnds_mouse)) + geom_point(aes(size=2), show.legend="NA", col="darkgreen") + theme( axis.text.x = element_text( hjust=1, vjust=0.3), text=element_text(size=22), panel.background = element_rect(fill="white", color="black"), panel.grid.major = element_line(color="white"), panel.grid.minor = element_line(color="white"), axis.title = element_text( size = 18)) + xlab(expression(paste("Spearman's ",rho," (Unadjusted)"))) + ylab(expression(paste("Spearman's ",rho," (Adjusted)"))) +
 labs(tag = "C")


# Figure 1d
ggenes <- read.table("c://Users/egamazon/Dropbox/Papers/SelectionExpression/Paper_v1/dropbox_evolution_paper/FinalSubmission/Resubmission_CellReports/Cell_Reports_files/REVISION/PEERJ_SUBMISSION/Final_Submission/REVISED_FIGURES/figure1d_data_lowexp.txt", sep = "\t", header = TRUE, row.names = 1)
fig1.d = ggplot(ggenes, aes(x=log(Mean_RPKM,2), y=Mouse_dN.dS)) + geom_point(color="black") + theme(panel.background = element_rect(fill="white", color="black"), panel.grid.major = element_line(color="white"), panel.grid.minor = element_line(color="white"), text=element_text(size=18), axis.title = element_text(size = 18), axis.text = element_text(size = 18)) + xlab("log(Expression)") + ylab("dN/dS") +
labs(tag = "D") +
 scale_colour_manual(name= "Smooth Curve", values = c("LM" = "orange", "LOESS" = "red", "GAM" = "blue"), guide="legend") +
  guides(colour = guide_legend(override.aes = list(linetype=c(1,0), shape=c(NA, 16)))) +                                                      
stat_smooth(method = "lm", formula = y ~ x, size = 1, se = FALSE, colour = "orange") + 
    stat_smooth(method = "loess", formula = y ~ x, size = 1, se = FALSE, colour = "red") + 
    stat_smooth(method = "gam", formula = y ~ s(x), size = 1, se = FALSE, colour = "blue") # + stat_smooth(method = "gam",


# PLOT
fig1.a + (fig1.b + fig1.c + fig1.d + plot_layout(ncol=1)) + plot_layout(guides = 'keep')

