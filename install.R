pkgs = c("knitr", "rmarkdown", "ggplot2", "ggpubr", "reshape2","cowplot","superheat","plyr","dplyr", "vegan", "reshape", "devtools", "pheatmap")
ncores = parallel::detectCores()

install.packages(pkgs, Ncpus = ncores)


devtools::install_github("benjjneb/dada2")
devtools::install_github("benjjneb/decontam")
devtools::install_version("mvtnorm", version = "1.0-8", repos = "http://cran.us.r-project.org")
devtools::install_version("fpc", version = "2.1-11.1", repos = "http://cran.us.r-project.org")

BiocManager::install("phyloseq")
BiocManager::install("DESeq2")
BiocManager::install("ShortRead")
BiocManager::install("dendextend")
BiocManager::install("apeglm")
