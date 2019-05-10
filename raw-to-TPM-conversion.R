raw_tab <- read.table("example_data/sample_raw_counts.tsv", header=T, row.names=1, sep="\t")
lengths_bps <- read.table("example_data/Num_bps.tsv", header=T, row.names=1, sep="\t")

lengths_kbs <- lengths_bps / 1000

rpk_tab <- data.frame(apply(raw_tab, 2, function(x) ( x / lengths_kbs ) )) # normalizing by gene length
tpm_tab <- apply(rpk_tab, 2, function (x) ( x / sum(x)) * 1000000) # normalizing for sample depth

colnames(tpm_tab) <- colnames(raw_tab)
tpm_tab2 <- data.frame("CDS_ID"=row.names(tpm_tab), tpm_tab, row.names=NULL)

write.table(tpm_tab2, "example_data/sample_TPM.tsv")

