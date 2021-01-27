#########################################
### ADMIXTURE GRAPHS OUTLIERS HEATMAP ###
#########################################

# Loading required packages
library(wesanderson)
library(data.table)
library(ggplot2)
library(reshape)
library(scales)

sink(snakemake@log[["hm_log"]])
# Defining input using smk objects
Single_bb_df<-snakemake@input[["bb_outl_tab"]]
hm_out<-snakemake@output[["hm_s_out"]]
bb_file<-read.table(snakemake@input[["bb_list"]], header = F)
title=paste("Heat map, backbone: ",paste(bb_file$V1, sep="", collapse=" "), sep="", collapse=" ")
# Formatting data
tab<-read.table(Single_bb_df, header = T, sep="\t", dec = ".")
pal <- wes_palette("Zissou1",5, type = "continuous")
long_format_dat <- melt(tab)
colnames(long_format_dat) <- c("Target", "Model", "Outliers")

# Plotting and saving Heatmap
pdf(hm_out)
print(ggplot(long_format_dat, aes(x = Model, y = Target, fill = Outliers)) +
        ggtitle(title)+
        geom_tile() +
        scale_fill_gradientn(colours = pal ,values = rescale(c(0,0.01,0.1,0.25,0.5,1))) +
        theme(axis.text = element_text(size = 12), axis.title = element_blank())+
        scale_x_discrete(expand = c(0, 0)) +
        scale_y_discrete(expand = c(0, 0)))
dev.off()

sink()
