# Defining leaves
leaves<-c(outgroup, Zhokhov, a_ne_sample, a_eu_sample, target)


# Divergence models ===================================================================================

inner_nodes <- c("A_Arctic","Arctic","A_EU","NE","All_dogs","root")
edges <- parent_edges(c(edge(outgroup,"root"),
                        edge("All_dogs","root"),
                        edge("A_Arctic","All_dogs"),
                        edge("NE","All_dogs"),
                        edge("Arctic","A_Arctic"),
                        edge(a_ne_sample,"NE"),
                        edge(a_eu_sample,"A_EU"),
                        edge(Zhokhov,"A_Arctic"),
                        edge(target,"Arctic"),
                        admixture_edge("A_EU","NE","Arctic","alpha")))
div_k9 <- agraph(leaves, inner_nodes, edges)
#fitting the graph to our data
div_k9_fit <- fit_graph(canids, div_k9)
outliers<-append(outliers,nrow(poor_fits(div_k9_fit,sigma=s)))

inner_nodes <- c("A_Arctic","Arctic","A_EU","NE","All_dogs","root")
edges <- parent_edges(c(edge(outgroup,"root"),
                        edge("All_dogs","root"),
                        edge("A_Arctic","All_dogs"),
                        edge("NE","All_dogs"),
                        edge("Arctic","A_Arctic"),
                        edge(a_ne_sample,"NE"),
                        edge(a_eu_sample,"A_EU"),
                        edge(Zhokhov,"Arctic"),
                        edge(target,"Arctic"),
                        admixture_edge("A_EU","NE","A_Arctic","alpha")))
div_v_k9 <- agraph(leaves, inner_nodes, edges)
#fitting the graph to our data
div_v_k9_fit <- fit_graph(canids, div_v_k9)
outliers<-append(outliers,nrow(poor_fits(div_v_k9_fit,sigma=s)))

# Model 1: Target admix EU/Arctic ======================================================================
inner_nodes <- c("A_Arctic","ghost","Arctic","A_EU","NE","All_dogs","root")
edges <- parent_edges(c(edge(outgroup,"root"),
                        edge("All_dogs","root"),
                        edge("A_Arctic","All_dogs"),
                        edge("NE","All_dogs"),
                        edge(a_ne_sample,"NE"),
                        edge(a_eu_sample,"A_EU"),
                        edge("Arctic","A_Arctic"),
                        edge(Zhokhov,"Arctic"),
                        edge(target,"ghost"),
                        admixture_edge("A_EU","NE","A_Arctic","alpha"),
                        admixture_edge("ghost","Arctic","A_EU","beta")))
adm_1_k9 <- agraph(leaves, inner_nodes, edges)
#fitting the graph to our data
adm_1_k9_fit <- fit_graph(canids, adm_1_k9)
outliers<-append(outliers,nrow(poor_fits(adm_1_k9_fit,sigma=s)))

# Model 2: Target admix Arctic/NE =========================================================================
inner_nodes <- c("A_Arctic","ghost","Arctic","A_EU","NE","All_dogs","root")
edges <- parent_edges(c(edge(outgroup,"root"),
                        edge("All_dogs","root"),
                        edge("A_Arctic","All_dogs"),
                        edge("NE","All_dogs"),
                        edge(a_ne_sample,"NE"),
                        edge(a_eu_sample,"A_EU"),
                        edge("Arctic","A_Arctic"),
                        edge(Zhokhov,"Arctic"),
                        edge(target,"ghost"),
                        admixture_edge("ghost","NE","A_Arctic","alpha"),
                        admixture_edge("A_EU","Arctic","ghost","beta")))
adm_2_old_k9 <- agraph(leaves, inner_nodes, edges)
#fitting the graph to our data
adm_2_old_k9_fit <- fit_graph(canids, adm_2_old_k9)
outliers<-append(outliers,nrow(poor_fits(adm_2_old_k9_fit,sigma=s)))

inner_nodes <- c("A_Arctic","ghost","Arctic","A_EU","NE","All_dogs","root")
edges <- parent_edges(c(edge(outgroup,"root"),
                        edge("All_dogs","root"),
                        edge("A_Arctic","All_dogs"),
                        edge("NE","All_dogs"),
                        edge(a_ne_sample,"NE"),
                        edge(a_eu_sample,"A_EU"),
                        edge("Arctic","A_Arctic"),
                        edge(Zhokhov,"Arctic"),
                        edge(target,"ghost"),
                        admixture_edge("A_EU","NE","A_Arctic","alpha"),
                        admixture_edge("ghost","Arctic","NE","beta")))
adm_2_k9 <- agraph(leaves, inner_nodes, edges)
#fitting the graph to our data
adm_2_k9_fit <- fit_graph(canids, adm_2_k9)
outliers<-append(outliers,nrow(poor_fits(adm_2_k9_fit,sigma=s)))

# Model 3: Target admix EU/NE ============================================================================
inner_nodes <- c("ghost","Arctic","A_EU","A_NE","NE","All_dogs","root")
edges <- parent_edges(c(edge(outgroup,"root"),
                        edge("All_dogs","root"),
                        edge("A_NE","All_dogs"),
                        edge("Arctic","All_dogs"),
                        edge("NE","A_NE"),
                        edge(a_ne_sample,"NE"),
                        edge(a_eu_sample,"A_EU"),
                        edge(Zhokhov,"Arctic"),
                        edge(target,"ghost"),
                        admixture_edge("A_EU","A_NE","Arctic","alpha"),
                        admixture_edge("ghost","A_EU","NE","beta")))
adm_3_k9 <- agraph(leaves, inner_nodes, edges)
#fitting the  graph to our data
adm_3_k9_fit <- fit_graph(canids, adm_3_k9)
outliers<-append(outliers,nrow(poor_fits(adm_3_k9_fit,sigma=s)))

# Printing Summary (text file) =============================================================================
text_file<-snakemake@output[["mod_summary"]]
sink(text_file,append=TRUE)
cat('\n############\n MODEL: div \n############\n',file=text_file,append=TRUE)
summary(div_k9_fit)
cat('\n############\n MODEL: div variant \n############\n',file=text_file,append=TRUE)
summary(div_v_k9_fit)
cat('\n############\n MODEL: adm EU/AR \n############\n',file=text_file,append=TRUE)
summary(adm_1_k9_fit)
cat('\n############\n MODEL: adm AR/NE \n############\n',file=text_file,append=TRUE)
summary(adm_2_k9_fit)
cat('\n############\n MODEL: adm AR/NE old \n############\n',file=text_file,append=TRUE)
summary(adm_2_old_k9_fit)
cat('\n############\n MODEL: adm EU/NE \n############\n',file=text_file,append=TRUE)
summary(adm_3_k9_fit)
sink()


# Generating and savings plots ==============================================================================

pdf(snakemake@output[["graphs"]])

name <- 'Div'
res <- paste('residual:', round(sum_of_squared_errors(div_k9_fit), digits=3), sep = " ")
plot(div_k9, color='darkturquoise',draw_inner_nodes=T,
     inner_node_color='aquamarine', show_inner_node_labels = TRUE,
     title=paste('Model ',name,'; Target: ',target, sep = ""), sub= res)
print(plot(div_k9_fit,sigma=s) + ggtitle(paste('D-statistics: ',name,'_fit',' (',s,' sigma)', sep = "")))

name <- 'Div_variant'
res <- paste('residual:', round(sum_of_squared_errors(div_v_k9_fit), digits=3), sep = " ")
plot(div_v_k9, color='darkturquoise',draw_inner_nodes=T,
     inner_node_color='aquamarine', show_inner_node_labels = TRUE,
     title=paste('Model ',name,'; Target: ',target, sep = ""), sub= res)
print(plot(div_v_k9_fit,sigma=s) + ggtitle(paste('D-statistics: ',name,'_fit',' (',s,' sigma)', sep = "")))

name <- 'EU/AR'
res <- paste('residual:', round(sum_of_squared_errors(adm_1_k9_fit), digits=3), sep = " ")
plot(adm_1_k9, color='darkturquoise',draw_inner_nodes=T,
     inner_node_color='aquamarine', show_inner_node_labels = TRUE,
     title=paste('Model ',name,'; Target: ',target, sep = ""), sub= res)
sec_tit<-paste('D-statistics: ',name,'_fit',' (',s,' sigma)', sep = "")
print(plot(adm_1_k9_fit,sigma=s) + ggtitle(sec_tit))

name <- 'AR/NE old'
res <- paste('residual:', round(sum_of_squared_errors(adm_2_old_k9_fit), digits=3), sep = " ")
plot(adm_2_old_k9, color='darkturquoise',draw_inner_nodes=T,
     inner_node_color='aquamarine', show_inner_node_labels = TRUE,
     title=paste('Model ',name,'; Target: ',target, sep = ""), sub= res)
sec_tit<-paste('D-statistics: ',name,'_fit',' (',s,' sigma)', sep = "")
print(plot(adm_2_old_k9_fit,sigma=s) + ggtitle(sec_tit))

name <- 'AR/NE'
res <- paste('residual:', round(sum_of_squared_errors(adm_2_k9_fit), digits=3), sep = " ")
plot(adm_2_k9, color='darkturquoise',draw_inner_nodes=T,
     inner_node_color='aquamarine', show_inner_node_labels = TRUE,
     title=paste('Model ',name,'; Target: ',target, sep = ""), sub= res)
sec_tit<-paste('D-statistics: ',name,'_fit',' (',s,' sigma)', sep = "")
print(plot(adm_2_k9_fit,sigma=s) + ggtitle(sec_tit))

name <- 'EU/NE'
res <- paste('residual:', round(sum_of_squared_errors(adm_3_k9_fit), digits=3), sep = " ")
plot(adm_3_k9, color='darkturquoise',draw_inner_nodes=T,
     inner_node_color='aquamarine', show_inner_node_labels = TRUE,
     title=paste('Model ',name,'; Target: ',target, sep = ""), sub= res)
sec_tit<-paste('D-statistics: ',name,'_fit',' (',s,' sigma)', sep = "")
print(plot(adm_3_k9_fit,sigma=s) + ggtitle(sec_tit))
dev.off()


# Printing outliers file ==============================================================================
sink(snakemake@output[["outliers"]])
line<-paste(outliers,collapse="\t")
cat(paste(line,"\n",sep=""))
sink()

#####################################################################################
