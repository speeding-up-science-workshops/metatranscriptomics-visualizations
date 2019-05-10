library("DESeq2")
library("phyloseq")
library("ggplot2")

count_tab <- read.table("example_data/sample_raw_counts.tsv", sep="\t", header=T, row.names=1)

sample_info_tab <- read.table("example_data/sample_info_tab.tsv", sep="\t", header=T, row.names=1)

count_phy <- otu_table(count_tab, taxa_are_rows=T)
sample_info_phy <- sample_data(sample_info_tab)
physeq <- phyloseq(count_phy, sample_info_phy)

# now converting our phyloseq object to a deseq object
deseq <- phyloseq_to_deseq2(physeq, ~treatment)

# and running deseq standard analysis:
deseq <- DESeq(deseq)

# pulling out our results table, we specify the object, the p-value we are going to use to filter our results, and what contrast we want to consider by first naming the column, then the two groups we care about
deseq_res_treatment <- results(deseq, alpha=0.01, contrast=c("treatment", "Low", "High"))

# we can get a glimpse at what this table currently holds with the summary command
summary(deseq_res_treatment)
# this tells us out of ~20,000 CDSs, with adj-p < 0.01, there are 140 increased when comparing the low CO2 treatment to the high CO2 treatment, and about 667 decreased
# "decreased" in this case means at a lower expression level in the low CO2 treatment than in the high, and "increased" means greater expression

# next let's stitch that together with KEGG annotations
anot_tab <- read.table("example_data/sample_annotation_classifications.tsv", header=T, row.names=1, sep="\t")[,1, drop=F]
head(anot_tab)
# reordering both for easy merging
anot_tab <- anot_tab[order(row.names(anot_tab)), , drop=F]
deseq_tab <- data.frame(deseq_res_treatment)

# all.equal(row.names(anot_tab), row.names(deseq_tab)) # checking to make sure they are the same

deseq_res_with_KOs <- cbind(deseq_tab, anot_tab)


# let's subset this table to only include these that pass our specified significance level
sigtab_res_deseq <- deseq_res_with_KOs[which(deseq_res_with_KOs$padj < 0.01), ]

# now we can see this table only contains those we consider significantly differentially abundant
head(sigtab_res_deseq) 

# and now let's sort that table by the baseMean column
sigtab_deseq_altered_vs_glassy_with_tax[order(sigtab_deseq_altered_vs_glassy_with_tax$baseMean, decreasing=T), ]
sigtab_res_deseq <- sigtab_res_deseq[order(sigtab_res_deseq$baseMean, decreasing=T), ]

out_tab <- data.frame("CDS_ID"=row.names(sigtab_res_deseq), sigtab_res_deseq, row.names = NULL)

write.table(out_tab, "DESeq_contrast.tsv", sep="\t", quote=F)

topGene <- rownames(deseq_res_treatment)[which.min(deseq_res_treatment$padj)]
data <- plotCounts(deseq, gene=topGene, intgroup=c("treatment"), returnData = T)
top_gene_KO <- anot_tab[row.names(anot_tab) %in% topGene, ]


pdf("Most-sig-gene.pdf")
ggplot(data, aes(x=treatment, y=count, fill=treatment)) +
  scale_y_log10() + theme_bw() +
  geom_dotplot(binaxis="y", stackdir="center") +
  ggtitle(paste0(topGene, " (", top_gene_KO, ")"))
dev.off()

