library("DESeq2")
library("ggplot2")
library("apeglm")
library("pheatmap")

## Reading in data
count_tab <- read.table("example_data/sample_raw_counts.tsv", sep="\t", header=T, row.names=1)
sample_info_tab <- read.table("example_data/sample_info_tab.tsv", sep="\t", header=T, row.names=1)

## DESeq
  # making the deseq object
deseq <- DESeqDataSetFromMatrix(countData = count_tab, colData = sample_info_tab, design = ~treatment)
  # setting the baseline treatment as the "Low" treatment
deseq$treatment <- relevel(deseq$treatment, ref = "Low")
  # and running deseq standard analysis:
deseq <- DESeq(deseq)
  # pulling out our results table, we specify the object, the p-value we are going to use to filter our results, and what contrast we want to consider by first naming the column, then the two groups we care about (we don't necessarily need this here because in this case we set the base level and there are only 2, but this is how you would state which things you want to contrast)
high_vs_low_contrast <- results(deseq, alpha=0.01, contrast=c("treatment", "High", "Low"))
  # we can get a glimpse at what this table currently holds with the summary command
summary(high_vs_low_contrast)
    # this tells us out of ~20,000 CDSs, with adj-p < 0.01, there are 667 increased when comparing the high CO2 treatment to the low CO2 treatment, and 140 decreased
    # "decreased" in this case means at a lower expression level in the High CO2 treatment than in the Low CO2 treatment, and "increased" means greater expression in the High as compared to the Low

  # next let's stitch that together with KEGG annotations
anot_tab <- read.table("example_data/sample_annotation_classifications.tsv", header=T, row.names=1, sep="\t")[,1, drop=F]
  # reordering both for easier merging
anot_tab <- anot_tab[order(row.names(anot_tab)), , drop=F]
deseq_tab <- data.frame(high_vs_low_contrast)

# all.equal(row.names(anot_tab), row.names(deseq_tab)) # checking to make sure they are the same
deseq_res_with_KOs <- cbind(deseq_tab, anot_tab)

  # let's subset this table to only include these that pass our specified significance level
sigtab_high_vs_low <- deseq_res_with_KOs[which(deseq_res_with_KOs$padj < 0.01), ]

  # and now let's sort that table by the baseMean column
sigtab_high_vs_low <- sigtab_high_vs_low[order(sigtab_high_vs_low$baseMean, decreasing=T), ]

out_tab <- data.frame("CDS_ID"=row.names(sigtab_high_vs_low), sigtab_high_vs_low, row.names = NULL)
  # writing out table
write.table(out_tab, "DESeq_high_vs_low_contrast.tsv", sep="\t", quote=F, row.names=F)

## Some visualizations
  # visualizing top gene (based on adj. p-value)
topGene <- rownames(high_vs_low_contrast)[which.min(high_vs_low_contrast$padj)]
data <- plotCounts(deseq, gene=topGene, intgroup=c("treatment"), returnData = T)
top_gene_KO <- anot_tab[row.names(anot_tab) %in% topGene, ]

pdf("Most-sig-gene.pdf")
ggplot(data, aes(x=treatment, y=count, fill=treatment)) +
  scale_y_log10() + theme_bw() +
  geom_dotplot(binaxis="y", stackdir="center") +
  ggtitle(paste0(topGene, " (", top_gene_KO, ")"))
dev.off()

  # applying shrinkage to lessen the impact of those with very low expression and those highly variable for plotting
high_vs_low_shr <- lfcShrink(deseq, coef="treatment_High_vs_Low", type="apeglm")

  # plotting MA plots of both
pdf("DESeq_MA_plots.pdf")
par(mfrow=c(1,2))
# original logFC
plotMA(deseq, ylim=c(min(high_vs_low_contrast$log2FoldChange),max(high_vs_low_contrast$log2FoldChange)), alpha=0.01, main="No log-fold-change shrinkage")
# apeglm shrinkage logFC
plotMA(high_vs_low_shr, ylim=c(min(high_vs_low_shr$log2FoldChange),max(high_vs_low_shr$log2FoldChange)), alpha=0.01, main="With log-fold-change shrinkage")
dev.off()

  # getting variance stabilized transformed table
deseq_vst <- varianceStabilizingTransformation(deseq)
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
vst_tab <- assay(deseq_vst)

# making heatmap of top 20 most highly expressed genes
select <- order(rowMeans(counts(deseq, normalized=TRUE)),
                decreasing=TRUE)[1:20]

pdf("DESeq-highly-expressed-heatmap.pdf")
pheatmap(vst_tab[select, ], cluster_rows=FALSE, cluster_cols=FALSE, annotation_col=sample_info_tab[, 1, drop=F])
dev.off()

# making heatmap of sample distances
euc_dists <- dist(t(vst_tab))
euc_dist_mat <- as.matrix(euc_dists)

pdf("DESeq-sample-euc-dist-heatmap.pdf")
pheatmap(euc_dist_mat, clustering_distance_rows=euc_dists, clustering_distance_cols=euc_dists, clustering_method="ward.D2")
dev.off()





