## Updated: 01/16/2019
## Written by: Ilya Korsunsky, ilya.korsunsky@gmail.com
## Functions to more efficiently process raw single-cell RNA-seq data
## Used in normalize_plot_cells.R
## (Functions adapted from Seurat package [Butler et al., 2018])


# Normalize columns
normalizeData <- function(A, scaling_factor = 1e4, method) {
    if(!'dgCMatrix' %in% class(A)) A <- as(A, "dgCMatrix")

    if (method == "log") {
        A@x <- A@x / rep.int(Matrix::colSums(A), diff(A@p))
        A@x <- scaling_factor * A@x
        A@x <- log(1 + A@x)
    } else if (method == "fft") {
        A@x <- A@x / rep.int(Matrix::colSums(A), diff(A@p))
        A@x <- scaling_factor * A@x
        A@x <- sqrt(A@x) + sqrt(1 + A@x)
    } else if (method == "geneCLR") {
        A@x <- as.numeric(normalizeCLR_dgc(A@x, A@p, A@i, ncol(A), nrow(A), 1))
    } else if (method == "cellCLR") {
        A@x <- as.numeric(normalizeCLR_dgc(A@x, A@p, A@i, ncol(A), nrow(A), 2))
    } else {
        stop(sprintf("ERROR: method %s not implemented", method))
    }

        return(A)
}

# Find variable genes by group
findVariableGenes <- function (X, groups, min_expr = 0.1, max_expr = Inf, min_dispersion = 0, 
    max_dispersion = Inf, return_top_n = 0) 
{
    groups <- factor(groups)
    groups_int <- as.integer(factor(groups)) - 1
    groups_table <- table(groups_int)
    means_nonlog <- exp_mean(X@x, X@p, X@i, ncol(X), nrow(X), 
        groups_int, groups_table)
    colnames(means_nonlog) <- levels(groups)
    row.names(means_nonlog) <- row.names(X)
    vmr <- log_vmr(X@x, X@p, X@i, ncol(X), nrow(X), means_nonlog, 
        groups_int, groups_table)
    colnames(vmr) <- levels(groups)
    row.names(vmr) <- row.names(X)
    vargenes_df <- dplyr::inner_join(means_nonlog %>% log1p %>% 
        as.data.frame() %>% tibble::rownames_to_column("symbol") %>% 
        tidyr::gather(group, gene_mean, -symbol), vmr %>% as.data.frame() %>% 
        tibble::rownames_to_column("symbol") %>% tidyr::gather(group, 
        gene_dispersion, -symbol), by = c("symbol", "group")) %>% 
        dplyr::arrange(-gene_dispersion) %>% subset(gene_mean >= 
        min_expr & gene_mean <= max_expr) %>% subset(gene_dispersion >= 
        min_dispersion & gene_dispersion <= max_dispersion)
    ## For some reason, this isn't working
    if (return_top_n > 0) {
        vargenes_union <- unique(data.table(vargenes_df)[, head(.SD, 
            return_top_n), by = group][, symbol])
        return(vargenes_union)
    }
    else {
        return(vargenes_df)
    }
}

# Scale data by row
ScaleDataSeurat <- function (data.use, margin = 1, scale.max = 10,
                                block.size = 1000) {

    if (margin == 2) data.use %<>% t
    max.block <- ceiling(nrow(data.use)/block.size)

    ## Define data and functions to use in sparse and dense cases
    if (class(data.use) == "dgCMatrix" | class(data.use) == "dgTMatrix") {
        scale_fxn <- function(x) {
            FastSparseRowScale(mat = x, scale = TRUE, center = TRUE,
                               scale_max = scale.max, display_progress = FALSE)
        }
    } else {
        scale_fxn <- function(x) {
            FastRowScale(mat = x, scale = TRUE, center = TRUE,
                               scale_max = scale.max, display_progress = FALSE)
        }
        data.use <- as.matrix(data.use)
    }

    ## Do scaling, at once or in chunks
    if (max.block == 1) {
        scaled.data <- scale_fxn(data.use)
    } else {
        scaled.data <- matrix(NA, nrow(data.use), ncol(data.use))
        for (i in 1:max.block) {
            idx.min <- (block.size * (i - 1))
            idx.max <- min(nrow(data.use), (block.size * i - 1) + 1)
            my.inds <- idx.min:idx.max
            scaled.data[my.inds, ] <- scale_fxn(data.use[my.inds, , drop = F])
        }
    }

    colnames(scaled.data) <- colnames(data.use)
    row.names(scaled.data) <- row.names(data.use)
    scaled.data[is.na(scaled.data)] <- 0
    if (margin == 2) scaled.data %<>% t
    return(scaled.data)
}

environment(ScaleDataSeurat) <- asNamespace("Seurat")

# Cosine normalize values
cosine_normalize <- function(X, MARGIN = 1, do_safe = TRUE) {
    if (do_safe) {
        X <- sweep(X, MARGIN, apply(X, MARGIN, max), "/")
    }
    sweep(X, MARGIN, apply(X, MARGIN, function(x) sqrt(sum(x^2))), "/")
}
