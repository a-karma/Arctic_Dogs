#######################################################
### Testing NE/EU Ancestry Component in Arctic dogs ###
#######################################################

# This script is just a wrapper that prepares the input for the modelling script
# in a consistent manner. This should be sufficiently generalised to be independent
# from the actual graph that need to be fitted to the data.
# Unfortunately, changing completely the naming in the model fitting script
# will require to also modify this file

# loading required packages
library(admixturegraph)
library(ggplot2)

# Defining variables using snakemake objects
canids<- read.table(file = snakemake@input[["d_table"]], header =T)
target<- snakemake@params[["sample"]]
s=snakemake@params[["sigma"]]
bb=read.table(snakemake@input[["pop_names"]])
bb_prefix=snakemake@params[["bb_prefix"]]

# The canids dataset consists of 6 column:
# The first for are the sample names of the quadruplet considered in the d-stats
# The last two are the D-stat value and its Z-score (respectively)
# Becuse of the way in which C.I around D are calculated, a value of D exactly
# equal to zero will couse an error. To avoid this, we approximate those values
# by adding an extra digit when necessary
for (i in 1:length(canids$D)){if(canids$D[[i]]==0.0000){canids$D[[i]]=0.00001}}

# Creating a vector of named variables based on the second column
# of the input file containing the backbone
# this will be used to define the leaves of the admixture graph
sample_names<-as.vector(bb$V1)
pop_names<-as.vector(bb$V2)
for (i in 1:length(sample_names)){
  assign(pop_names[i],sample_names[i])}

# Initializing vector of the admixture graph outliers
outliers<-c(bb_prefix,target)

# Calling the script containing the definition of the admixture graphs to fit
# which will also produce the three requested output files
# while the logfile is available through the sink() command
sink(snakemake@log[["model_fitting_log"]])
source(snakemake@params[["models_to_fit"]])
sink()
