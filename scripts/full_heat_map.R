#########################################
### ADMIXTURE GRAPHS OUTLIERS HEATMAP ###
#########################################

library(wesanderson)
library(data.table)
library(ggplot2)
library(reshape)
library(scales)

# Creating the log file
sink(snakemake@log[["full_hm_log"]])

# Defining input/output variables using smk objects
All_bb_df<-snakemake@input[["outl_full_tab"]]
hm_out<-snakemake@output[["hm_f_out"]]

# Aggregating bacbone info and calculating for each target
# the average number of outliers across backbones under each model
tab<-read.table(All_bb_df, header = T, sep="\t",dec = ".")
dat<-aggregate(tab[,3:ncol(tab)], list(tab$sample), mean)

# Formatting data
long_format_dat <- melt(dat)
colnames(long_format_dat) <- c("Target", "Model", "Outliers")

# Setting colour palette for heat map
pal <- wes_palette("Zissou1",5, type = "continuous")

# Creating a Heatmap with values
pdf(hm_out)
print(ggplot(long_format_dat, aes(x = Model, y = Target, fill = Outliers)) +
  geom_tile() +
  scale_fill_gradientn(colours = pal ,values = rescale(c(0,0.01,0.1,0.25,0.5,1))) +
  theme(axis.text = element_text(size = 12), axis.title = element_blank())+
  scale_x_discrete(expand = c(0, 0)) +
  scale_y_discrete(expand = c(0, 0)) +
  geom_text(aes(label = Outliers), color = "black", size = 3))
dev.off()

source(snakemake@params[["models_schematics"]])

sink()
