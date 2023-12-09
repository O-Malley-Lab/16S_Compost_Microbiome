# Plotting NMDS from qiime2 data
# C. Ahern


if (!requireNamespace("devtools", quietly = TRUE)){install.packages("devtools")}
devtools::install_github("jbisanz/qiime2R")

library(qiime2R)

metadata <- read.delim("/Users/colleenahern/Documents/BASF/Analysis2/metadata.txt", header=T, stringsAsFactors=F, sep="\t")
metadata2 <- read.delim("/Users/colleenahern/Documents/BASF/Analysis2/metadata2.txt", header=T, stringsAsFactors=F, sep="\t")

weighted_dis <- read_qza("/Users/colleenahern/Documents/BASF/Analysis2/diversity/weighted_unifrac_distance_matrix.qza")
unweighted_dis <- read_qza("/Users/colleenahern/Documents/BASF/Analysis2/diversity/unweighted_unifrac_distance_matrix.qza")

setwd("/Users/colleenahern/Documents/BASF")

dsm <- qza_to_phyloseq(features="Analysis2/diversity/weighted_unifrac_distance_matrix.qza",tree="Analysis2/rooted-16S-tree-filteredSVs.qza", taxonomy = "Analysis2/16S-rep-seqs-taxonomy.qza", metadata= "Analysis2/metadata2.txt")

pswd <- qza_to_phyloseq(features="Analysis2/16S-table-noplant-rarefied-10000_filtered.qza",
                        tree="Analysis2/rooted-16S-tree-filteredSVs.qza", 
                        taxonomy = "Analysis2/16S-rep-seqs-taxonomy.qza", 
                        metadata= "Analysis2/metadata2.txt")

pswd.prop <- transform_sample_counts(pswd, function(otu) otu/sum(otu))

ord.nmds.bray_wd <- ordinate(pswd.prop, method="NMDS", distance="bray")

library(ggforce)
deg_plot <- plot_ordination(pswd.prop, ord.nmds.bray_wd, color="Substrate", shape = "Day", title="Bray NMDS") + scale_shape_binned(name = "Day", limits=c(0,90), breaks=c(0,45,90)) + geom_point(size = 4) 

deg_plot2 <- plot_ordination(pswd.prop, ord.nmds.bray_wd, color="Substrate", shape = "Day", title="Bray NMDS") + scale_shape_binned(name = "Day", limits=c(0,90), breaks=c(0,45,90)) + scale_color_manual(labels = c("1","2","3","4","5","6"), values = c("coral1","gold","green","turquoise","blue","deeppink2")) + geom_point(size = 4) 

devtools::install_github("cmartin/ggConvexHull")
library(ggConvexHull)

deg_plot + geom_convexhull(aes(fill = pswd.prop@sam_data$Group, 
                               colour = pswd.prop@sam_data$Substrate),
                           alpha = 0.5, show.legend = F)  

deg_plot2 + geom_convexhull(aes(fill = pswd.prop@sam_data$Group, 
                                colour = pswd.prop@sam_data$Substrate),
                            alpha = 0.5, show.legend = F)

deg_plot + ggforce::geom_mark_ellipse(aes(fill = pswd.prop@sam_data$Group,
                                          color = pswd.prop@sam_data$Substrate), show.legend = F)

deg_plot2 + ggforce::geom_mark_ellipse(aes(fill = pswd.prop@sam_data$Group,
                                           color = pswd.prop@sam_data$Substrate), show.legend = F)


