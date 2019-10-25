library("DESeq2")
library("dendextend")
library("phyloseq")
library("ggplot2")

count_tab <- read.table("example_data/sample_raw_counts.tsv", sep="\t", header=T, row.names=1)

sample_info_tab <- read.table("example_data/sample_info_tab.tsv", sep="\t", header=T, row.names=1)

# to make first we need to make a DESeq2 object
deseq_counts <- DESeqDataSetFromMatrix(count_tab, colData = sample_info_tab, design = ~treatment) 
# we have to include the "colData" and "design" arguments because they are 
# required, as they are needed for further downstream processing by DESeq2, 
# but for our purposes of simply transforming the data right now, they don't 
# matter

deseq_counts_vst <- varianceStabilizingTransformation(deseq_counts)
# NOTE: If you get this error here with your dataset: "Error in
# estimateSizeFactorsForMatrix(counts(object), locfunc =locfunc, : every
# gene contains at least one zero, cannot compute log geometric means", that
# can be because the count table is sparse with many zeroes, which is common
# with marker-gene surveys. In that case you'd need to use a specific
# function first that is equipped to deal with that. You could run:
# deseq_counts <- estimateSizeFactors(deseq_counts, type = "poscounts")
# now followed by the transformation function:
# deseq_counts_vst <- varianceStabilizingTransformation(deseq_counts)

# and here is pulling out our transformed table
vst_trans_count_tab <- assay(deseq_counts_vst)

# and calculating our Euclidean distance matrix
euc_dist <- dist(t(vst_trans_count_tab))

euc_clust <- hclust(euc_dist, method="ward.D2")

# hclust objects like this can be plotted with the generic plot() function
# plot(euc_clust) 
# but i like to change them to dendrograms for two reasons:
# 1) it's easier to color the dendrogram plot by groups
# 2) if wanted you can rotate clusters with the rotate() 
#    function of the dendextend package

euc_dend <- as.dendrogram(euc_clust, hang=0.1)
dend_cols <- as.character(sample_info_tab$color[order.dendrogram(euc_dend)])
labels_colors(euc_dend) <- dend_cols

# plot(euc_dend, ylab="VST Euc. dist.")

# writig out 
pdf("Hclust.pdf")
plot(euc_dend, ylab="VST Euc. dist.")
dev.off()


# making our phyloseq object with transformed table
vst_count_phy <- otu_table(vst_trans_count_tab, taxa_are_rows=T)
sample_info_tab_phy <- sample_data(sample_info_tab)
vst_physeq <- phyloseq(vst_count_phy, sample_info_tab_phy)

# generating and visualizing the PCoA with phyloseq
vst_pcoa <- ordinate(vst_physeq, method="MDS", distance="euclidean")
eigen_vals <- vst_pcoa$values$Eigenvalues # allows us to scale the axes according to their magnitude of separating apart the samples

# plot_ordination(vst_physeq, vst_pcoa, color="treatment") + 
#   geom_point(size=1) + labs(col="treatment") + 
#   geom_text(aes(label=rownames(sample_info_tab), hjust=0.3, vjust=-0.4)) + 
#   coord_fixed(sqrt(eigen_vals[2]/eigen_vals[1])) + ggtitle("PCoA") + 
#   scale_color_manual(values=unique(as.character(sample_info_tab$color)[order(sample_info_tab$treatment)])) + 
#   theme(legend.position="none")

# writing out
pdf("PCoA.pdf")
plot_ordination(vst_physeq, vst_pcoa, color="treatment") + 
  geom_point(size=1) + labs(col="treatment") + 
  geom_text(aes(label=rownames(sample_info_tab), hjust=0.5, vjust=-0.4)) + 
  coord_fixed(sqrt(eigen_vals[2]/eigen_vals[1])) + ggtitle("PCoA") + 
  scale_color_manual(values=unique(as.character(sample_info_tab$color)[order(sample_info_tab$treatment)])) + 
  theme(legend.position="none")
dev.off()
