library(plotly)
library(ComplexHeatmap)

# Function to check which variables are constant
check_constant_vars <- function(df, vars) {
  constant_vars <- sapply(vars, function(v) {
    length(unique(df[[v]])) == 1  # Check if variable has only one unique value
  })
  return(names(constant_vars)[constant_vars])  # Return names of constant variables
}

ensembl2symbol <- function(ensembl_ids, OrgDb) {
  # Load package if not already loaded
  if (!requireNamespace("AnnotationDbi", quietly = TRUE)) {
    stop("Package 'AnnotationDbi' is required but not installed.")
  }
  # Use AnnotationDbi::select to obtain the mapping
  mapping <- AnnotationDbi::select(
    OrgDb,
    keys = ensembl_ids,
    columns = c("SYMBOL"),
    keytype = "ENSEMBL"
  )
  # Remove rows with missing SYMBOL values
  mapping <- mapping[!is.na(mapping$SYMBOL), ]
  # Remove duplicates (optional, but useful for named vector)
  mapping <- mapping[!duplicated(mapping$ENSEMBL), ]
  # Create a named vector with SYMBOL as values and ENSEMBL IDs as names
  named_vector <- setNames(mapping$SYMBOL, mapping$ENSEMBL)
  return(named_vector)
}

mydeseq2 = function(counts, min_reads, sample_info, min_sample = 5, design_formula = ~ batch + condition, test = "Wald", ...) {
  # Filter lowly expressed genes
  counts = counts[rowSums(counts) > min_reads & rowSums(counts > 1) > min_sample, ]
  # Align sample_info rows with count matrix columns
  sample_info <- sample_info[colnames(counts), ]
  # Make condition names syntactically valid
  sample_info$condition <- make.names(sample_info$condition)
  # Create DESeq2 dataset object
  dds_current <- DESeqDataSetFromMatrix(
    countData = counts,
    colData = sample_info,
    design = design_formula
  )
  # Run DESeq2 differential expression analysis
  diff <- DESeq(dds_current)
  # Initialize lists to store results
  res_all <- list()
  df_all_genes <- list()
  # Get unique conditions
  conds <- unique(sample_info$condition)
  # Loop through all pairs of conditions for comparisons
  for (i in 1:(length(conds) - 1)) {
    for (j in (i + 1):length(conds)) {
      contrast <- c("condition", conds[i], conds[j])
      # Extract results for the contrast
      res <- results(diff, contrast = contrast, ...)
      # Convert results to data.frame
      res <- data.frame(res)
      # Add a column describing the comparison
      res$comparison_n_vs_d <- paste0(conds[i], "_vs_", conds[j])
      # Add gene identifiers as a column
      res$gene <- rownames(res)
      # Remove rows with NA values
      res <- na.omit(res)
      # Clean gene names by removing suffixes like ".1", ".2", etc.
      rownames(res) <- sub("\\.\\d+$", "", rownames(res))
      # Map Ensembl IDs to gene symbols using org.Hs.eg.db
      res$symbol <- mapIds(org.Hs.eg.db,
                           keys = rownames(res),
                           column = "SYMBOL",
                           keytype = "ENSEMBL",
                           multiVals = "first")
      # Store results in the lists using comparison name as key
      res_all[[paste0(conds[i], "_vs_", conds[j])]] <- res
      df_all_genes[[paste0(conds[i], "_vs_", conds[j])]] <- res
    }
  }
  # Perform regularized log transformation (rlog) on count data
  rld <- rlog(dds_current, blind = FALSE)
  # Retrieve available results names (not assigned here, just called)
  resultsNames(diff)
  # Return a list with combined results, rlog-transformed data, and sample info
  return(list(
    res = do.call(rbind, res_all),
    df_all_genes = do.call(rbind, df_all_genes),
    rld = rld,
    coldata = sample_info
  ))
}





specific_deseq2 = function(counts, min_reads, sample_info, exp, contr, min_sample = 5, design_formula = ~ batch + condition, ...) {
  # Filter lowly expressed genes
  counts = counts[rowSums(counts) > min_reads & rowSums(counts > 1) > min_sample, ]
  # Print dimensions of count matrix after filtering
  message("Counts dimensions after filtering: ", paste(dim(counts), collapse = " x "))
  # Stop if no genes pass filtering criteria
  if (nrow(counts) == 0) stop("No gene passed the filters. Reduce min_reads or min_sample.")
  # Align sample_info rows with count matrix columns
  sample_info <- sample_info[colnames(counts), ]
  # Make condition names syntactically valid
  sample_info$condition <- make.names(sample_info$condition)
  # Create DESeq2 dataset object
  dds_current <- DESeqDataSetFromMatrix(
    countData = counts,
    colData = sample_info,
    design = design_formula
  )
  # Run DESeq2 differential expression analysis with additional arguments
  diff <- DESeq(dds_current, ...)
  # Define contrast of interest for comparison
  contrast <- c("condition", exp, contr)
  # Extract results for specified contrast
  res <- results(diff, contrast = contrast)
  # Convert results to data.frame
  res <- data.frame(res)
  # Add a column describing the comparison
  res$comparison_exp_vs_contr <- paste0(exp, "_vs_", contr)
  # Add gene identifiers as a column
  res$gene <- rownames(res)
  # Remove rows with NA values
  res <- na.omit(res)
  # Clean gene names by removing suffixes like ".1", ".2", etc.
  rownames(res) <- sub("\\.\\d+$", "", rownames(res))
  # Map Ensembl IDs to gene symbols using org.Hs.eg.db
  res$symbol <- mapIds(org.Hs.eg.db,
                       keys = rownames(res),
                       column = "SYMBOL",
                       keytype = "ENSEMBL",
                       multiVals = "first")
  # Perform regularized log transformation (rlog) on count data
  rld <- rlog(dds_current, blind = FALSE)
  # Extract colData as data.frame
  coldata <- as.data.frame(colData(dds_current))
  # Subset rlog data and coldata to include only samples in the contrast
  rld_subset <- rld[, coldata$condition %in% c(exp, contr)]
  coldata_subset <- coldata[coldata$condition %in% c(exp, contr), ]
  # Retrieve available results names (called but not used)
  resultsNames(diff)
  # Return list with results, rlog data subset, and metadata subset
  return(list(
    res = res,
    df_all_genes = res,
    rld = rld_subset,
    coldata = coldata_subset
  ))
}





preprocessing_PCA = function(counts, min_reads, sample_info, normalization, min_sample = 5, design_formula = ~ batch + condition, ...) {
  # Filter lowly expressed genes
  counts = counts[rowSums(counts) > min_reads & rowSums(counts > 1) > min_sample, ]
  # Align sample_info rows with count matrix columns
  sample_info <- sample_info[colnames(counts), ]
  # Make condition names syntactically valid
  sample_info$condition <- make.names(sample_info$condition)
  # Create DESeq2 dataset object
  dds_current <- DESeqDataSetFromMatrix(
    countData = counts,
    colData = sample_info,
    design = design_formula
  )
  # Get unique condition names
  conds <- unique(sample_info$condition)
  # Perform normalization according to specified method
  if (normalization == "vst") {
    rld <- vst(dds_current, blind = FALSE)
  } else if (normalization == "rlog") {
    rld <- rlog(dds_current, blind = FALSE)
  }
  # Return normalized data and sample metadata
  return(list(
    rld = rld,
    coldata = sample_info
  ))
}




mypcaAnalysis <- function(title_1vs2, rld, intgroup1 = "condition_description", intgroup2 = NULL, ...) {
  # Automatically check if intgroup2 (e.g. batch) has more than one level
  if (!is.null(intgroup2)) {
    levels_intgroup2 <- unique(colData(rld)[[intgroup2]])
    if (length(levels_intgroup2) <= 1) {
      # Ignore intgroup2 if it has only one level or is NA
      intgroup2 <- NULL
    }
  }
  # Build intgroups vector to pass to plotPCA
  intgroups <- if (is.null(intgroup2)) intgroup1 else c(intgroup1, intgroup2)
  # Compute PCA data using DESeq2::plotPCA with returnData = TRUE
  pcaData <- DESeq2::plotPCA(rld, intgroup = intgroups, returnData = TRUE)
  # Add sample names and grouping columns for plotly visualization
  pcaData$name <- colnames(rld)
  pcaData$cond <- colData(rld)[[intgroup1]]
  if (!is.null(intgroup2)) {
    pcaData$group2 <- as.factor(colData(rld)[[intgroup2]])
  }
  # Calculate percentage variance explained by PC1 and PC2 using prcomp on assay data
  assayData <- SummarizedExperiment::assay(rld)
  pca <- prcomp(t(assayData))
  percentVar <- pca$sdev^2 / sum(pca$sdev^2)
  pc1_var <- round(100 * percentVar[1], 1)
  pc2_var <- round(100 * percentVar[2], 1)
  # Create plotly PCA scatter plot depending on presence of intgroup2
  if (is.null(intgroup2)) {
    p <- plotly::plot_ly(
      data = pcaData,
      x = ~PC1,
      y = ~PC2,
      text = ~name,
      color = ~cond,
      type = "scatter",
      mode = "markers",
      marker = list(size = 10, line = list(width = 1, color = 'black'))
    ) %>%
      layout(
        title = paste("PCA:", title_1vs2),
        xaxis = list(title = paste0("PC1 (", pc1_var, "%)")),
        yaxis = list(title = paste0("PC2 (", pc2_var, "%)")),
        legend = list(title = list(text = paste0("<b>", intgroup1, "</b>")))
      )
  } else {
    p <- plotly::plot_ly(
      data = pcaData,
      x = ~PC1,
      y = ~PC2,
      text = ~name,
      color = ~cond,
      symbol = ~group2,
      type = "scatter",
      mode = "markers",
      marker = list(size = 10, line = list(width = 1, color = 'black'))
    ) %>%
      layout(
        title = paste("PCA:", title_1vs2),
        xaxis = list(title = paste0("PC1 (", pc1_var, "%)")),
        yaxis = list(title = paste0("PC2 (", pc2_var, "%)")),
        legend = list(title = list(text = paste0("<b>", intgroup1, " / ", intgroup2, "</b>")))
      )
  }
  # Return the plotly PCA plot object
  return(p)
}







mypcaAnalysis_test <- function(title_1vs2, rld, intgroup1 = "condition_description", intgroup2 = NULL, ...) {
  # Check if intgroup2 (e.g., batch) has more than one unique level
  if (!is.null(intgroup2)) {
    levels_intgroup2 <- unique(colData(rld)[[intgroup2]])
    if (length(levels_intgroup2) <= 1) {
      intgroup2 <- NULL
    }
  }
  # Build vector of grouping variables to pass to plotPCA
  intgroups <- if (is.null(intgroup2)) intgroup1 else c(intgroup1, intgroup2)
  # Compute PCA data with DESeq2::plotPCA and get data frame
  pcaData <- DESeq2::plotPCA(rld, intgroup = intgroups, returnData = TRUE)
  # Add sample names and grouping columns for plotting
  pcaData$name <- colnames(rld)
  pcaData$cond <- colData(rld)[[intgroup1]]
  if (!is.null(intgroup2)) {
    pcaData$group2 <- as.factor(colData(rld)[[intgroup2]])
  }
  # Compute % variance explained by PC1 and PC2 using prcomp on assay data
  assayData <- SummarizedExperiment::assay(rld)
  pca <- prcomp(t(assayData))
  percentVar <- pca$sdev^2 / sum(pca$sdev^2)
  pc1_var <- round(100 * percentVar[1], 1)
  pc2_var <- round(100 * percentVar[2], 1)
  # Plot PCA with plotly depending on presence of intgroup2
  if (is.null(intgroup2)) {
    p <- plot_ly(
      data = pcaData,
      x = ~PC1,
      y = ~PC2,
      text = ~name,
      color = ~cond,
      type = "scatter",
      mode = "markers",
      marker = list(size = 10, line = list(width = 1, color = 'black'))
    ) %>%
      layout(
        title = paste("PCA:", title_1vs2),
        xaxis = list(title = paste0("PC1 (", pc1_var, "%)")),
        yaxis = list(title = paste0("PC2 (", pc2_var, "%)")),
        legend = list(title = list(text = paste0("<b>", intgroup1, "</b>")))
      )
  } else {
    # Define color palette for conditions using RColorBrewer or fallback to rainbow
    cond_levels <- unique(pcaData$cond)
    n_cond <- length(cond_levels)
    if (n_cond <= 8) {
      cond_colors <- brewer.pal(n_cond, "Set1")
    } else {
      cond_colors <- grDevices::rainbow(n_cond)
    }
    # Define symbols list for batch groups (group2)
    symbol_list <- c("circle", "cross", "diamond", "triangle-up", "triangle-down", "square", "x")
    group2_levels <- unique(pcaData$group2)
    n_batch <- length(group2_levels)
    # Use rainbow colors for batches
    batch_colors <- grDevices::rainbow(n_batch)
    # Initialize empty plotly object
    p <- plot_ly()
    # Add traces for conditions with fixed colors and circle markers
    for (i in seq_along(cond_levels)) {
      cond_level <- cond_levels[i]
      df_cond <- pcaData[pcaData$cond == cond_level, ]
      p <- p %>% add_markers(
        data = df_cond,
        x = ~PC1,
        y = ~PC2,
        text = ~name,
        marker = list(
          color = cond_colors[i],
          size = 12,
          symbol = "circle",
          line = list(width = 1, color = 'black')
        ),
        name = as.character(cond_level),
        showlegend = TRUE,
        legendgroup = "Condition"
      )
    }
    # Add traces for batch groups with different symbols and gray transparent colors
    for (j in seq_along(group2_levels)) {
      group2_level <- group2_levels[j]
      df_group2 <- pcaData[pcaData$group2 == group2_level, ]
      symb <- ifelse(j <= length(symbol_list), symbol_list[j], "circle") # fallback symbol
      p <- p %>% add_markers(
        data = df_group2,
        x = ~PC1,
        y = ~PC2,
        text = ~name,
        marker = list(
          color = batch_colors[j],
          size = 6,
          symbol = "circle",
          line = list(width = 1, color = 'black')
        ),
        name = as.character(group2_level),
        showlegend = TRUE,
        legendgroup = "Batch"
      )
    }
    # Layout settings for the plot
    p <- p %>% layout(
      title = paste("PCA:", title_1vs2),
      xaxis = list(title = paste0("PC1 (", pc1_var, "%)")),
      yaxis = list(title = paste0("PC2 (", pc2_var, "%)")),
      legend = list(title = list(text = paste0("<b>", intgroup1, " / ", intgroup2, "</b>")))
    )
  }
  # Return the plotly PCA plot object
  return(p)
}













check_batch_effect <- function(rld, batch_var = "batch", cond_var = "condition_description", var_threshold = 0.1) {
  # Extract assay data matrix from the SummarizedExperiment object
  assayData <- SummarizedExperiment::assay(rld)
  # Perform PCA on the transposed assay data (samples as rows)
  pca <- prcomp(t(assayData))
  # Create a dataframe with scores for the first two principal components (PC1 and PC2)
  scores <- as.data.frame(pca$x[, 1:2])
  # Add batch information from the colData metadata to the scores dataframe
  scores[[batch_var]] <- colData(rld)[[batch_var]]
  # Add condition information from the colData metadata to the scores dataframe
  scores[[cond_var]] <- colData(rld)[[cond_var]]
  # Perform ANOVA for PC1 using batch and condition as predictors
  model_PC1 <- aov(PC1 ~ get(batch_var) + get(cond_var), data = scores)
  # Perform ANOVA for PC2 using batch and condition as predictors
  model_PC2 <- aov(PC2 ~ get(batch_var) + get(cond_var), data = scores)
  # Extract the sum of squares for each factor from the ANOVA summary for PC1
  ss_PC1 <- summary(model_PC1)[[1]][, "Sum Sq"]
  # Extract the sum of squares for each factor from the ANOVA summary for PC2
  ss_PC2 <- summary(model_PC2)[[1]][, "Sum Sq"]
  # Calculate the proportion of variance in PC1 explained by the batch variable
  perc_batch_PC1 <- ss_PC1[1] / sum(ss_PC1)
  # Calculate the proportion of variance in PC2 explained by the batch variable
  perc_batch_PC2 <- ss_PC2[1] / sum(ss_PC2)
  # Calculate the average variance explained by batch across PC1 and PC2
  mean_batch_var <- mean(c(perc_batch_PC1, perc_batch_PC2))
  # Determine if the batch effect is significant based on the variance threshold
  batch_significant <- mean_batch_var > var_threshold
  # Return a list with the significance flag and the variance explained by batch
  return(list(
    batch_significant = batch_significant,
    batch_var_exp = mean_batch_var
  ))
}


check_batch_condition_confounding <- function(sample_info, batch_var = "batch", cond_var = "condition_description", pval_threshold = 0.05, alpha = 0.5) {
  # Create contingency table between batch and condition variables
  tab <- table(sample_info[[batch_var]], sample_info[[cond_var]])
  # Calculate parameters for low count threshold
  n_batches <- nrow(tab)
  n_conditions <- ncol(tab)
  total_samples <- sum(tab)
  expected_mean <- total_samples / (n_batches * n_conditions)
  threshold <- alpha * expected_mean
  # Identify cells in contingency table with low counts (<= threshold)
  low_count_cells <- which(tab <= threshold, arr.ind = TRUE)
  # Initialize warning message variable
  warning_msg <- NULL
  # Check for zero counts indicating strong confounding
  if(any(tab == 0)) {
    warning_msg <- "WARNING: Some batch-condition combinations have zero samples, indicating potential strong confounding. <br>This means some batches and conditions never occur together, making it impossible to disentangle batch effects from biological conditions."
  } 
  # Otherwise, warn if there are low count cells under the threshold
  else if(nrow(low_count_cells) > 0) {
    warning_msg <- paste0(
      "WARNING: Some batch-condition combinations have low sample counts (<= ", round(threshold, 2), "), which may increase confounding risk.<br><br>",
      "The threshold to define low sample counts is calculated as <code><span style='white-space: nowrap;'>alpha * [z / (n * m)]</span></code>, where:<br>",
      " - <b>z</b>: total number of samples, in this case ", total_samples, "<br>",
      " - <b>n</b>: number of batches, in this case ", n_batches, "<br>",
      " - <b>m</b>: number of conditions, in this case ", n_conditions, "<br>",
      " - <b>alpha</b>: stringency parameter to define how low the count must be to consider it problematic, in this case ", alpha, "<br><br>",
      "Low counts reduce statistical power and may bias batch effect correction.<br>",
      "Therefore, applying batch correction in the presence of strong confounding or very uneven sample distribution <b>may remove true biological signals or introduce artifacts</b>.<br>",
      "Batch correction should be applied cautiously in such scenarios."
    )
  }
  # Perform chi-squared test of independence on contingency table
  test <- suppressWarnings(chisq.test(tab))
  # Determine if confounding is significant based on p-value threshold
  confounded <- test$p.value < pval_threshold
  # Return results as a list
  return(list(
    confounded = confounded,            # TRUE if significant confounding detected
    p_value = test$p.value,             # p-value from chi-squared test
    contingency_table = tab,            # The contingency table used for testing
    warning = warning_msg,              # Warning message if any
    low_count_cells = low_count_cells,  # Positions of low count cells in table
    threshold = threshold               # Threshold used to define low counts
  ))
}










my_MA_andVolcano_plot = function(title_1vs2, res, qvalue, logfc, ...) {
  # Identify significant genes based on adjusted p-value and log2 fold change thresholds
  select <- which(res$padj < qvalue & abs(res$log2FoldChange) > logfc)
  cat("Total number of significant genes:", length(select), "\n")
  # Add a 'significant' column to indicate significant vs non-significant genes
  res$significant <- "not_significant"
  res$significant[select] <- "significant"
  # Create labels combining gene symbol and gene ID
  res$label <- paste0(res$symbol, " (", res$gene, ")")
  # Create MA plot with plotly
  ma_plot <- plotly::plot_ly(
    data = res,
    x = ~log10(baseMean + 1),
    y = ~log2FoldChange,
    color = ~significant,
    text = ~label,
    type = "scatter",
    mode = "markers"
  ) %>%
    layout(title = paste("MA Plot: ", title_1vs2))
  # Create Volcano plot with plotly
  volcano_plot <- plotly::plot_ly(
    data = res,
    x = ~log2FoldChange,
    y = ~-log10(padj),
    color = ~significant,
    text = ~label,
    type = "scatter",
    mode = "markers"
  ) %>%
    layout(title = paste("Volcano Plot: ", title_1vs2))
  # Return a list containing both plots
  return(list(ma = ma_plot, volcano = volcano_plot))
}



my_genetable = function(res, title_1vs2, qvalue, logfc, ...) {
  # Remove duplicate rows based on rownames (genes)
  res <- res[!duplicated(rownames(res)), ]
  # Add FoldChange column calculated from log2FoldChange
  res$FoldChange <- 2^res$log2FoldChange
  # Filter significant genes by padj and log2FoldChange thresholds
  deg <- res[which(res$padj < qvalue & abs(res$log2FoldChange) > logfc), ]
  # Initialize differential expression column as "no"
  res$differentially_expressed <- rep("no", nrow(res))
  # Find common genes between full results and filtered DEGs
  common_genes <- intersect(rownames(res), rownames(deg))
  # Mark common genes as "yes" for differential expression
  res[rownames(res) %in% common_genes, "differentially_expressed"] <- "yes"
  # Create interactive datatable for all genes with export buttons
  res_dt <- datatable(
    res,
    extensions = 'Buttons',
    options = list(
      dom = 'Blfrtip',
      buttons = c('copy', 'csv', 'excel', 'pdf', 'print'),
      lengthMenu = list(c(10, 25, 50, -1),
                        c(10, 25, 50, "All"))
    ),
    rownames = FALSE,
    caption = paste("All genes: ", title_1vs2)
  )
  # Create interactive datatable for differentially expressed genes (DEGs)
  deg_dt <- datatable(
    deg,
    extensions = 'Buttons',
    options = list(
      dom = 'Blfrtip',
      buttons = c('copy', 'csv', 'excel', 'pdf', 'print'),
      lengthMenu = list(c(10, 25, 50, -1),
                        c(10, 25, 50, "All"))
    ),
    rownames = FALSE,
    caption = paste("DEGs: ", title_1vs2)
  )
  # Prepare underscored title string for export filenames (not used here)
  title_1vs2_underscore <- gsub(" ", "_", title_1vs2)
  # (commented) Export DEGs to CSV file
  # write.table(deg, file = paste("./output/gse_results/DEG_", title_1vs2_underscore, ".csv"), sep = ",", row.names = FALSE, quote = FALSE)
  # Return list with both datatables and original and filtered results
  return(list(res_dt = res_dt, deg_dt = deg_dt, res = res, deg = deg))
}




my_heatmaps <- function(deg, rld, title_1vs2) {
  # Create named vector mapping genes to their symbols
  gene_symbols <- setNames(deg$symbol, deg$gene)
  # Select top 20 genes sorted by adjusted p-value
  top20_genes <- deg %>% as.data.frame() %>% arrange(padj) %>% head(20) %>% dplyr::pull(gene)
  # Ensure top 20 genes are present in the assay matrix
  top20_genes <- intersect(top20_genes, rownames(assay(rld)))
  # Warn if no top genes found, set heatmap variable to NULL
  if (length(top20_genes) == 0) {
    warning("No top 20 genes found in the assay matrix.")
    h_20 <- NULL
  } else {
    # Extract expression matrix for top 20 genes
    expr_top20 <- assay(rld)[top20_genes, , drop = FALSE]
    # Replace rownames with gene symbols where available
    rownames(expr_top20) <- coalesce(gene_symbols[rownames(expr_top20)], rownames(expr_top20))
    # Scale expression data by rows (genes)
    expr_top20_scaled <- t(scale(t(expr_top20)))
    # Remove rows containing non-finite values after scaling
    expr_top20_scaled <- expr_top20_scaled[apply(expr_top20_scaled, 1, function(x) all(is.finite(x))), , drop = FALSE]
    # Warn if scaled matrix is empty, set heatmap variable to NULL
    if (nrow(expr_top20_scaled) == 0) {
      warning("Top 20 genes expression matrix is empty after scaling.")
      h_20 <- NULL
    } else {
      # Generate heatmap for top 20 genes without clustering
      h_20 <- pheatmap(
        expr_top20_scaled,
        cluster_rows = FALSE,
        cluster_cols = FALSE,
        show_rownames = TRUE,
        show_colnames = TRUE,
        color = colorRampPalette(c("blue", "white", "red"))(50),
        main = paste("Top 20 DE Genes:", title_1vs2),
        cellwidth = 12,
        cellheight = 12,
        legend = TRUE,
        legend_breaks = c(-2, 0, 2),
        legend_labels = c("-2 (Z-score)", "0", "2 (Z-score)")
      )
    }
  }
  # Extract all DE genes
  all_de_genes <- deg %>% dplyr::pull(gene)
  # Keep only those present in assay matrix
  all_de_genes <- intersect(all_de_genes, rownames(assay(rld)))
  # Warn and set heatmap variable to NULL if none found
  if (length(all_de_genes) == 0) {
    warning("No DE genes found in the assay matrix.")
    h_all <- NULL
  } else {
    # Extract expression matrix for all DE genes
    expr_de <- assay(rld)[all_de_genes, , drop = FALSE]
    # Replace rownames with gene symbols when available
    rownames(expr_de) <- coalesce(gene_symbols[rownames(expr_de)], rownames(expr_de))
    # Scale expression data by rows (genes)
    expr_de_scaled <- t(scale(t(expr_de)))
    # Remove rows with non-finite values
    expr_de_scaled <- expr_de_scaled[apply(expr_de_scaled, 1, function(x) all(is.finite(x))), , drop = FALSE]
    # Warn and set heatmap variable to NULL if empty matrix
    if (nrow(expr_de_scaled) == 0) {
      warning("All DE genes expression matrix is empty after scaling.")
      h_all <- NULL
    } else {
      # Generate clustered heatmap for all DE genes
      h_all <- pheatmap(
        expr_de_scaled,
        cluster_rows = TRUE,
        cluster_cols = TRUE,
        show_rownames = FALSE,
        show_colnames = TRUE,
        color = colorRampPalette(c("blue", "white", "red"))(50),
        main = paste("All DE Genes:", title_1vs2),
        cellwidth = 15,
        legend = TRUE,
        legend_breaks = c(-2, 0, 2),
        legend_labels = c("-2 (Z-score)", "0", "2 (Z-score)")
      )
    }
  }
  # Return list of heatmaps for top 20 and all DE genes
  return(list(h_20 = h_20, h_all = h_all))
}






run_gsea_ensembl = function(genelist, contrast_label, OrgDb, ont, toType = "ENSEMBL", pvalueCutoff = 0.2, ...) {
  gse <- gseGO(geneList=sort(genelist, decreasing = T), 
                  ont = ont, 
                  keyType = toType, 
                  minGSSize = 3, 
                  maxGSSize = 800, 
                  pvalueCutoff = pvalueCutoff, 
                  verbose = TRUE, 
                  OrgDb = OrgDb, 
                  pAdjustMethod = "BH")
  return(gse)
}


run_gsea_symbol = function(genelist, contrast_label, OrgDb, ont, toType = "SYMBOL", pvalueCutoff = 0.2, ...) {
  gse <- gseGO(geneList=sort(genelist, decreasing = T), 
                  ont = ont, 
                  keyType = toType, 
                  minGSSize = 3, 
                  maxGSSize = 800, 
                  pvalueCutoff = pvalueCutoff, 
                  verbose = TRUE, 
                  OrgDb = OrgDb, 
                  pAdjustMethod = "BH")
  return(gse)
}

goenrichment = function(markers, down = FALSE, n_genes = 100, logfc_col = "avg_log2FC", pvalueCutoff = 0.05, qvalueCutoff = 0.05){
  # If down is TRUE, invert the sign of log fold changes
  if(down){
    markers[[logfc_col]] = -markers[[logfc_col]]
  }
  # Select markers with positive log fold change
  up = markers[markers[[logfc_col]] > 0, ]
  # Order markers by log fold change decreasingly
  up = up[order(up[[logfc_col]], decreasing = T), ]
  # Number of available markers after filtering
  n_available <- nrow(up)
  # Select top n_genes or all available if fewer
  if(n_genes >= n_available){
    marker_genes <- up$gene
  } else {
    marker_genes <- head(up$gene, n_genes)
  }
  # Map gene symbols (ENSEMBL IDs) to Entrez IDs using org.Hs.eg.db
  entrez_ids <- mapIds(
    org.Hs.eg.db,
    keys = marker_genes,
    column = "ENTREZID",
    keytype = "ENSEMBL",
    multiVals = "first"
  )
  # Remove NA values from Entrez IDs
  entrez_ids <- na.omit(entrez_ids)
  # Perform GO enrichment analysis with enrichGO function
  go_enrichment <- enrichGO(
    gene = entrez_ids,
    OrgDb = org.Hs.eg.db,
    keyType = "ENTREZID",
    ont = "ALL",
    pAdjustMethod = "BH",
    pvalueCutoff = pvalueCutoff,
    qvalueCutoff = qvalueCutoff,
    readable = TRUE
  )
  # Filter enrichment results based on adjusted p-value and q-value cutoffs
  go_enrichment@result = go_enrichment@result[go_enrichment@result$p.adjust < go_enrichment@pvalueCutoff &
                                                go_enrichment@result$qvalue < go_enrichment@qvalueCutoff, ]
  # Return the filtered GO enrichment object
  return(go_enrichment)
}



goenrichment_reactome = function(markers, logfc_col = "avg_log2FC", down = FALSE, n_genes = 100, pvalueCutoff = 0.05, qvalueCutoff = 0.05){
  # If down is TRUE, invert the sign of log fold changes
  if(down){
    markers[[logfc_col]] = -markers[[logfc_col]]
  }
  # Select markers with positive log fold change
  up = markers[markers[[logfc_col]] > 0, ]
  # Order markers by log fold change decreasingly
  up = up[order(up[[logfc_col]], decreasing = TRUE), ]
  # Number of available markers after filtering
  n_available <- nrow(up)
  # Check if n_genes is "all" (case-insensitive) and select all genes if so
  if(is.character(n_genes) && tolower(n_genes) == "all"){
    marker_genes <- up$gene
  # Otherwise, check if n_genes is numeric and positive
  } else if (is.numeric(n_genes)) {
    if(n_genes <= 0){
      stop("n_genes must be a positive number or 'all'")
    }
    # Select top n_genes or all available if fewer
    if(n_genes >= n_available){
      marker_genes <- up$gene
    } else {
      marker_genes <- head(up$gene, n_genes)
    }
  # Stop with error if n_genes is not valid
  } else {
    stop("n_genes must be a positive number or 'all'")
  }
  # Map gene symbols (ENSEMBL IDs) to Entrez IDs using org.Hs.eg.db
  entrez_ids <- mapIds(
    org.Hs.eg.db,
    keys = marker_genes,
    column = "ENTREZID",
    keytype = "ENSEMBL",
    multiVals = "first"
  )
  # Remove NA values from Entrez IDs
  entrez_ids <- na.omit(entrez_ids)
  # Perform Reactome pathway enrichment analysis with enrichPathway function
  go_enrichment <- enrichPathway(
    gene = entrez_ids,
    organism = "human",
    pAdjustMethod = "BH",
    pvalueCutoff = pvalueCutoff,
    qvalueCutoff = qvalueCutoff,
    readable = TRUE
  )
  # Filter enrichment results based on adjusted p-value and q-value cutoffs
  go_enrichment@result = go_enrichment@result[
    go_enrichment@result$p.adjust < go_enrichment@pvalueCutoff &
    go_enrichment@result$qvalue < go_enrichment@qvalueCutoff, ]
  # Return the filtered Reactome enrichment object
  return(go_enrichment)
}






my_WGCNA_analysis <- function(sample_info, dds, dds_subset = NULL) {
  # Normalize and correct batch effect if subset is provided
  if (!is.null(dds_subset)) {
    vst_counts <- assay(vst(dds_subset))
    vst_bc <- removeBatchEffect(vst_counts, batch = dds_subset$batch)
    final_counts <- vst_bc
  } else {
    final_counts <- assay(dds)
  }
  # Identify good samples and genes to remove outliers
  gsg <- goodSamplesGenes(t(final_counts), verbose = 3)
  # If not all samples and genes are good, filter out bad ones
  if (!gsg$allOK) {
    final_counts <- final_counts[gsg$goodSamples, ]
    keep <- rowSums(counts(dds)[gsg$goodSamples, ] >= 10) >= 6
  } else {
    keep <- rowSums(counts(dds) >= 10) >= 6
  }
  # Filter counts based on genes passing the threshold
  filtered_counts <- t(final_counts[keep, ])
  # Subset final counts to good genes only
  wgcna_data <- final_counts[gsg$goodGenes == TRUE,]
  # Perform PCA on transposed WGCNA data
  pca <- prcomp(t(wgcna_data))
  pca.dat <- pca$x
  # Calculate variance from PCA standard deviations
  pca.var <- pca$sdev^2
  # Calculate percentage variance explained by each PC
  pca.var.percent <- round(pca.var/sum(pca.var)*100, digits = 2)
  pca.dat <- as.data.frame(pca.dat)
  # Plot the first two principal components with labels
  ggplot_PCA <- ggplot(pca.dat, aes(PC1, PC2)) +
    geom_point() +
    geom_text(label = rownames(pca.dat)) +
    labs(x = paste0('PC1: ', pca.var.percent[1], ' %'),
         y = paste0('PC2: ', pca.var.percent[2], ' %'))
  # Define candidate powers for soft-thresholding
  powers <- c(1:10, seq(12, 50, 2))
  # Pick soft threshold power based on scale-free topology criteria
  sft <- pickSoftThreshold(filtered_counts, powerVector = powers, networkType = "signed", verbose = 5)
  soft_power <- 14
  sft.data <- sft$fitIndices
  # Plot scale-free topology fit index versus power
  a1 <- ggplot(sft.data, aes(Power, SFT.R.sq, label = Power)) +
    geom_point() +
    geom_text(nudge_y = 0.1) +
    geom_hline(yintercept = 0.8, color = 'red') +
    labs(x = 'Power', y = 'Scale free topology model fit, signed R^2') +
    theme_classic()
  # Plot mean connectivity versus power
  a2 <- ggplot(sft.data, aes(Power, mean.k., label = Power)) +
    geom_point() +
    geom_text(nudge_y = 0.1) +
    labs(x = 'Power', y = 'Mean Connectivity') +
    theme_classic()
  power_plots <- grid.arrange(a1, a2, nrow = 2)
  # Ensure filtered_counts are numeric for network construction
  filtered_counts[] <- sapply(filtered_counts, as.numeric)
  # Construct the signed weighted gene co-expression network
  net <- blockwiseModules(filtered_counts,
                          power = soft_power,
                          TOMType = "signed",
                          maxBlockSize = 5000,
                          numericLabels = FALSE,
                          randomSeed = 1234,
                          verbose = 3)
  # Extract module eigengenes from network results
  MEs <- net$MEs
  # Encode conditions as factors
  sample_info$condition <- as.factor(sample_info$condition)
  # Create design matrix with no intercept for multi-group comparisons
  design <- model.matrix(~ 0 + condition, data = sample_info)
  colnames(design) <- levels(sample_info$condition)
  # Fit linear model on module eigengenes versus condition design
  fit <- lmFit(t(MEs), design)
  fit <- eBayes(fit)
  # Extract statistics for all contrasts adjusted by FDR
  stats_df <- topTable(fit, number = Inf, adjust = "fdr") %>%
    rownames_to_column("module")
  # Compute correlations between module eigengenes and condition traits
  traits <- model.matrix(~ 0 + condition, data = sample_info)
  mod_trait_cor <- cor(MEs, traits, use = "p")
  # Calculate p-values for correlation significance
  pvals <- corPvalueStudent(mod_trait_cor, nrow(filtered_counts))
  # Prepare data for heatmap combining eigengenes and traits
  heatmap_data <- cbind(MEs, traits)
  # Create heatmap of module-trait relationships
  CLP <- CorLevelPlot(heatmap_data,
               x = colnames(traits),
               y = names(MEs),
               col = c("blue1", "skyblue", "white", "pink", "red"),
               rotTitleX = 90
               )
  # Return all key outputs in a list
  return(list(
    ggplot_PCA = ggplot_PCA,
    power_plots = power_plots,
    module_colors = net$colors,
    module_eigengenes = MEs,
    stats = stats_df,
    module_trait_correlation = mod_trait_cor,
    module_trait_pvals = pvals,
    CLP = CLP
    )
  )
}




