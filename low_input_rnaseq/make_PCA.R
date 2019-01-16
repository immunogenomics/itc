## Updated: 01/15/2019
## Plot a PCA from log2(tpm + 1) values
## (As shown in figure 2a-b)

#### SUMMARY
# Plot PC1 vs PC2 scores of PCA:
#   chosing expressed genes with sd > 1.4 (N = 1,007)
#   centering on samples

# Dependencies.
library(ggplot2)
library(RColorBrewer)

wd <- "/Users/mgutierr/Documents/work/ITC/manuscript/submission_NatCommunications_Mar2017/Revision/GEOsubmission/"
setwd(wd)

log2tpm <- read.table("low_input_rnaseq_gene_normalized.txt", header = T, stringsAsFactors = F, row.names = 1)
m <- read.table("low_input_rnaseq_meta_data.txt", header = T, stringsAsFactors = F)

reorder <- c(which(m$cell_type %in% "CD4"), which(m$cell_type %in% "CD8"), which(m$cell_type %in% "MAIT"), which(m$cell_type %in% "iNKT"), which(m$cell_type %in% "Vd1"), which(m$cell_type %in% "Vd2"), which(m$cell_type %in% "NK"))
m <- m[reorder,]
log2tpm <- log2tpm[,reorder]

# selec genes with a certain level of expression
filter <- apply(log2tpm, 1, function(x) length(x[x>2])>=10) 
filtered <- log2tpm[filter,]
genes <- rownames(filtered)[grep("ENS", rownames(filtered))]

# select 1,007 top variable genes based on standard deviation
row_sd <- apply(filtered, 1, sd)
idx_sd <- row_sd > 1.4

# scale genes
scaledGenes <- scale(t(filtered[idx_sd, ]))

# do PCA
pca1 <- prcomp(scaledGenes, scale = FALSE, center = FALSE)

# prepare data frame for plotting
rownames(m) <- m$sampleID
meta <- m[rownames(pca1$x),]
pca1_r <- cbind(pca1$x, meta)

pca1_r <- as.data.frame(pca1_r)
for(i in 1:ncol(pca1$x)){
  pca1_r[,i] <- as.numeric(as.character(pca1_r[,i]))
}

pca1_r$cell_type <- as.factor(pca1_r$cell_type)
pca1_r$cell_type <- factor(pca1_r$cell_type, levels = c("CD4" , "CD8" , "MAIT" ,  "iNKT"  , "Vd1" , "Vd2" , "NK"))

# colors used for each cell type
cbPalette <- c("#0072B2", "#56B4E9", "#009E73", "#CC79A7","#E69F00", "#D55E00",  "#606060")

# plot

## PC1 vs PC2 (figure 2a)
postscript(file = "PCs1vs2.ps", width=4, height=3, family = "Helvetica", paper = "special", horizontal = F)
ggplot(data = pca1_r) +
  geom_point(size = 2.5, shape = 1, stroke = 1, aes(PC1, PC2, color = cell_type)) +
  theme_bw(base_size = 8) +
  theme(axis.text.x=element_text(size=8), axis.text.y=element_text(size=8), legend.text=element_text(size=8)) +
  scale_colour_manual(values=cbPalette,
                      name="cell type",
                      breaks=c("CD4" , "CD8" , "MAIT" ,  "NKT"  , "Vd1" , "Vd2" , "NK"),
                      labels=c("CD4" , "CD8" , "MAIT" ,  "iNKT"  , expression(paste("V", delta, "1", sep = "")) , expression(paste("V", delta, "2", sep = "")) , "NK"))
dev.off()

## Boxplots for PC1, by cell type (figure 2b)

pdf("PC1_boxplots_perCellType.pdf", width=4, height=2, family = "Helvetica")
par(mai = c(.5,.6,0.2,.2), ps = 8)
boxplot(list(CD4=pca1_r$PC1[which(pca1_r$cell_type %in% "CD4")], CD8=pca1_r$PC1[which(pca1_r$cell_type %in% "CD8")], MAIT=pca1_r$PC1[which(pca1_r$cell_type %in% "MAIT")], NKT=pca1_r$PC1[which(pca1_r$cell_type %in% "NKT")], Vd1=pca1_r$PC1[which(pca1_r$cell_type %in% "Vd1")], Vd2=pca1_r$PC1[which(pca1_r$cell_type %in% "Vd2")], NK=pca1_r$PC1[which(pca1_r$cell_type %in% "NK")]), 
        col = cbPalette, 
        ylab = "", 
        main = "", 
        names = F, 
        yaxt = "n", 
        cex = 0.5)
text(seq(1,7,by=1), par("usr")[3]-5, 
     srt = 45, adj= 1, xpd = TRUE,
     labels = c("CD4+ T", "CD8+ T", "MAIT", "iNKT", expression(paste("V", delta, "1", sep = "")), expression(paste("V", delta, "2", sep = "")), "NK"), ps = 8)
axis(side = 2, tck = -.05, labels = NA)
axis(side = 2, lwd = 0, line = -.2, las = 1)
title(ylab="PC1, 16% of variance", line=1.7, ps = 8)
dev.off()
