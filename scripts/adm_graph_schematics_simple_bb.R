
# Loading required library
library(admixturegraph)
library(ggplot2)

# Defining leaves with generic names
outgroup <- "Outgroup"
a_eu_sample<-"Ancient Europe"
a_ne_sample<-"Ancient Near East"
target<-"Target"
leaves <- c(a_ne_sample, a_eu_sample,"Zhokhov dog", target,outgroup)

###################
# Defining graphs #
###################

# Divergence
inner_nodes <- c("A_Arctic","Arctic","A_EU","NE","All_dogs","root")
edges <- parent_edges(c(edge(outgroup,"root"),
                        edge("All_dogs","root"),
                        edge("A_Arctic","All_dogs"),
                        edge("NE","All_dogs"),
                        edge("Arctic","A_Arctic"),
                        edge(a_ne_sample,"NE"),
                        edge(a_eu_sample,"A_EU"),
                        edge("Zhokhov dog","A_Arctic"),
                        edge(target,"Arctic"),
                        admixture_edge("A_EU","NE","Arctic","alpha")))
div_k9 <- agraph(leaves, inner_nodes, edges)

# Divergence variant
inner_nodes <- c("A_Arctic","Arctic","A_EU","NE","All_dogs","root")
edges <- parent_edges(c(edge(outgroup,"root"),
                        edge("All_dogs","root"),
                        edge("A_Arctic","All_dogs"),
                        edge("NE","All_dogs"),
                        edge("Arctic","A_Arctic"),
                        edge(a_ne_sample,"NE"),
                        edge(a_eu_sample,"A_EU"),
                        edge("Zhokhov dog","Arctic"),
                        edge(target,"Arctic"),
                        admixture_edge("A_EU","NE","A_Arctic","alpha")))
div_v_k9 <- agraph(leaves, inner_nodes, edges)

# Admixture Arctic/Europe
inner_nodes <- c("A_Arctic","ghost","Arctic","A_EU","NE","All_dogs","root")
edges <- parent_edges(c(edge(outgroup,"root"),
                        edge("All_dogs","root"),
                        edge("A_Arctic","All_dogs"),
                        edge("NE","All_dogs"),
                        edge(a_ne_sample,"NE"),
                        edge(a_eu_sample,"A_EU"),
                        edge("Arctic","A_Arctic"),
                        edge("Zhokhov dog","Arctic"),
                        edge(target,"ghost"),
                        admixture_edge("A_EU","NE","A_Arctic","alpha"),
                        admixture_edge("ghost","Arctic","A_EU","beta")))
adm_1_k9 <- agraph(leaves, inner_nodes, edges)

# Old admixture Arctic/Near East
inner_nodes <- c("A_Arctic","ghost","Arctic","A_EU","NE","All_dogs","root")
edges <- parent_edges(c(edge(outgroup,"root"),
                        edge("All_dogs","root"),
                        edge("A_Arctic","All_dogs"),
                        edge("NE","All_dogs"),
                        edge(a_ne_sample,"NE"),
                        edge(a_eu_sample,"A_EU"),
                        edge("Arctic","A_Arctic"),
                        edge("Zhokhov dog","Arctic"),
                        edge(target,"ghost"),
                        admixture_edge("ghost","NE","A_Arctic","alpha"),
                        admixture_edge("A_EU","Arctic","ghost","beta")))
adm_2_old_k9 <- agraph(leaves, inner_nodes, edges)

# Admixture Arctic/Near East
inner_nodes <- c("A_Arctic","ghost","Arctic","A_EU","NE","All_dogs","root")
edges <- parent_edges(c(edge(outgroup,"root"),
                        edge("All_dogs","root"),
                        edge("A_Arctic","All_dogs"),
                        edge("NE","All_dogs"),
                        edge(a_ne_sample,"NE"),
                        edge(a_eu_sample,"A_EU"),
                        edge("Arctic","A_Arctic"),
                        edge("Zhokhov dog","Arctic"),
                        edge(target,"ghost"),
                        admixture_edge("A_EU","NE","A_Arctic","alpha"),
                        admixture_edge("ghost","Arctic","NE","beta")))
adm_2_k9 <- agraph(leaves, inner_nodes, edges)

# Admixture Euorpe/Near East
inner_nodes <- c("ghost","Arctic","A_EU","A_NE","NE","All_dogs","root")
edges <- parent_edges(c(edge(outgroup,"root"),
                        edge("All_dogs","root"),
                        edge("A_NE","All_dogs"),
                        edge("Arctic","All_dogs"),
                        edge("NE","A_NE"),
                        edge(a_ne_sample,"NE"),
                        edge(a_eu_sample,"A_EU"),
                        edge("Zhokhov dog","Arctic"),
                        edge(target,"ghost"),
                        admixture_edge("A_EU","A_NE","Arctic","alpha"),
                        admixture_edge("ghost","A_EU","NE","beta")))
adm_3_k9 <- agraph(leaves, inner_nodes, edges)


#------------------------------------------------------------------------------------------------------------------------------------------------
# Printing graphs schematics
all_graphs<-snakemake@output[["graphs"]]

pdf(all_graphs)

name <- 'Divergence'
plot(div_k9, color='darkturquoise',draw_inner_nodes=T,
     inner_node_color='aquamarine', show_inner_node_labels = TRUE)
title(name, col.main="blue", cex.main=1.5)

name <- 'Div alternative'
plot(div_v_k9, color='darkturquoise',draw_inner_nodes=T,
     inner_node_color='aquamarine', show_inner_node_labels = TRUE)
title(name, col.main="blue", cex.main=1.5)

name <- 'Admixture EU/AR'
plot(adm_1_k9, color='darkturquoise',draw_inner_nodes=T,
     inner_node_color='aquamarine', show_inner_node_labels = TRUE)
title(name, col.main="blue", cex.main=1.5)

name <- 'Old admixture AR/NE'
plot(adm_2_old_k9, color='darkturquoise',draw_inner_nodes=T,
     inner_node_color='aquamarine', show_inner_node_labels = TRUE)
title(name, col.main="blue", cex.main=1.5)

name <- 'Admixture AR/NE'
plot(adm_2_k9, color='darkturquoise',draw_inner_nodes=T,
     inner_node_color='aquamarine', show_inner_node_labels = TRUE)
title(name, col.main="blue", cex.main=1.5)

name <- 'Admixture EU/NE'
plot(adm_3_k9, color='darkturquoise',draw_inner_nodes=T,
     inner_node_color='aquamarine', show_inner_node_labels = TRUE)
title(name, col.main="blue", cex.main=1.5)
dev.off()
