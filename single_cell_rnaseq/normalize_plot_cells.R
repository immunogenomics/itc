## Updated: 01/16/2019
## Normalize raw RNA-seq counts, scale genes, cosine normalize, and plot a PCA and UMAP.
## (As shown in figure 8)


# Libraries to import
library(data.table)
library(Matrix)
library(dplyr)
library(reticulate)
library(Rcpp)
##library(scales)
library(irlba)
library(Seurat)
library(ggplot2)

# Will use python for UMAP
use_python("/path/to/python")

## Install UMAP python package from Github:  https://github.com/lmcinnes/umap
umap <- import("umap")

# Additional functions to optomize data processing
sourceCpp("utils_ilya.cpp")
source("utils_ilya.R")

# Read gene counts and meta data
exprs_raw <- read.table("single_cell_rnaseq_gene_counts.txt", header = T, row.names = 1)
exprs_raw <- Matrix(as.matrix(exprs_raw), sparse = T)

meta_data <- read.table("single_cell_rnaseq_meta_data.txt", header = T)
# head(meta_data)
# all(meta_data$cell_id == colnames(exprs_raw))

# Log-normalize data
exprs_norm <- exprs_raw[,meta_data$cell_id] %>% normalizeData(method = "log")

# Scale data by gene for top 1000 variable genes within each donor
var_genes <- unique(data.table(findVariableGenes(exprs_norm, meta_data$donor))[, head(.SD, 1000), by = group][, symbol])
exprs_scaled <- exprs_norm[var_genes, ] %>% ScaleDataSeurat()

# Calculate top 20 principal components of cosine-normalized data
exprs_cosine <- exprs_scaled %>% cosine_normalize(2)
pca_res <- irlba::prcomp_irlba(t(exprs_cosine), 20)

# Calculate UMAP low-dimensional embedding of top 20 PCs
umap_res <- umap$UMAP(n_neighbors = 30L, metric = "correlation", min_dist = .1)$fit_transform(pca_res$x)


# Plot cells along top two PCs
ggplot(data = as.data.frame(pca_res$x), aes(x = pca_res$x[,1], y = pca_res$x[,2], col = meta_data$cell.type)) + 
    geom_point(size = 0.5, shape = 16) + 
    theme(axis.text.x = element_text(size = 8), 
          axis.text.y=element_text(size=8), 
          axis.title.x = element_text(size=8),
          axis.title.y = element_text(size=8)) +
    xlab("PC1") +
    ylab("PC2") +
    labs(col = "Cell type") +
    scale_color_manual(values=c("#0072B2", "#56B4E9", "#CC79A7", "#009E73", "#606060", "#E69F00", "#D55E00"))

# Plot cells in low-dimensional UMAP embedding
ggplot(data = as.data.frame(umap_res), aes(x = V1, y = V2, col = meta_data$cell.type)) + 
  geom_point(size = 0.5) + 
  theme(axis.text.x = element_text(size = 8), 
        axis.text.y=element_text(size=8), 
        axis.title.y = element_text(size=8),
        axis.title.x = element_text(size=8)) +
  xlab("UMAP 1") +
  ylab("UMAP 2") + 
  labs(col = "Cell type")+
  scale_color_manual(values=c("#0072B2", "#56B4E9", "#CC79A7", "#009E73", "#606060", "#E69F00", "#D55E00"))
