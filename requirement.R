# requirements.R
# If an error occurs when using the install method or the BiocManager method directly,
#  use the github installation method provided by debtools
install.packages(c(
  "Seurat",
  "SingleR",
  "curl",
  "devtools",
  "usethis",
  "tidyverse",
  "SummarizedExperiment",
  "scuttle",
  "patchwork",
  "CellChat",
  "NMF",
  "spacexr",
  "Matrix",
  "ggplot2",
  "SeuratData",
  "copykat",
  "scCancer",
  "DropletUtils",
  "AnnotationDbi",
  "dplyr",
  "IRanges",
  "S4Vectors",
  "stats",
  "clusterProfiler",
  "topGO",
  "Rgraphviz",
  "pathview",
  "org.Hs.eg.db", # For mouse, replace with org.Mm.eg.db if needed
  "stringr",
  "graph",
  "ggnewscale",
  "msigdbr",
  "GSVA",
  "pheatmap",
  "gelnet",
  "Rmisc",
  "magrittr"
))

# Install Bioconductor packages
if (!requireNamespace("BiocManager", quietly = TRUE)) {
    install.packages("BiocManager")
}
BiocManager::install(c(
  "SingleR",
  "SummarizedExperiment",
  "scuttle",
  "DropletUtils",
  "AnnotationDbi",
  "IRanges",
  "S4Vectors",
  "clusterProfiler",
  "topGO",
  "Rgraphviz",
  "pathview",
  "org.Hs.eg.db" # For mouse, replace with org.Mm.eg.db if needed
))
