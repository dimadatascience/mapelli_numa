library(plotly)
library(ComplexHeatmap)
# Funzione per verificare quali variabili sono costanti
check_constant_vars <- function(df, vars) {
  constant_vars <- sapply(vars, function(v) {
    length(unique(df[[v]])) == 1
  })
  return(names(constant_vars)[constant_vars])
}

ensembl2symbol <- function(ensembl_ids, OrgDb) {
  # Carica pacchetto se non già caricato
  if (!requireNamespace("AnnotationDbi", quietly = TRUE)) {
    stop("Package 'AnnotationDbi' is required but not installed.")
  }
  
  # Usa AnnotationDbi::select per ottenere il mapping
  mapping <- AnnotationDbi::select(
    OrgDb,
    keys = ensembl_ids,
    columns = c("SYMBOL"),
    keytype = "ENSEMBL"
  )
  
  # Rimuovi eventuali righe con SYMBOL mancante
  mapping <- mapping[!is.na(mapping$SYMBOL), ]
  
  # Rimuovi duplicati (opzionale, ma utile per named vector)
  mapping <- mapping[!duplicated(mapping$ENSEMBL), ]
  
  # Crea named vector
  named_vector <- setNames(mapping$SYMBOL, mapping$ENSEMBL)
  
  return(named_vector)
}

mydeseq2 = function(counts, min_reads, sample_info, min_sample = 5, design_formula = ~ batch + condition, ...) {
  # Filtro geni poco espressi
  counts = counts[rowSums(counts) > min_reads & rowSums(counts > 1) > min_sample, ]
  
  # Allineamento sample_info
  sample_info <- sample_info[colnames(counts), ]
  sample_info$condition <- make.names(sample_info$condition)

  # Crea oggetto DESeq2
  dds_current <- DESeqDataSetFromMatrix(
    countData = counts,
    colData = sample_info,
    design = design_formula
  )
  
  # DESeq2
  diff <- DESeq(dds_current, ...)
  
  res_all <- list()
  df_all_genes <- list()
  
  conds <- unique(sample_info$condition)
  
  for (i in 1:(length(conds) - 1)) {
    for (j in (i + 1):length(conds)) {
      contrast <- c("condition", conds[i], conds[j])
      res <- results(diff, contrast = contrast)
      res <- data.frame(res)
      res$comparison_n_vs_d <- paste0(conds[i], "_vs_", conds[j])
      res$gene <- rownames(res)
      res <- na.omit(res)
      rownames(res) <- sub("\\.\\d+$", "", rownames(res))  # rimuove ".1", ".2", etc
      res$symbol <- mapIds(org.Hs.eg.db,
                           keys = rownames(res),
                           column = "SYMBOL",
                           keytype = "ENSEMBL",
                           multiVals = "first")
      res_all[[paste0(conds[i], "_vs_", conds[j])]] <- res
      df_all_genes[[paste0(conds[i], "_vs_", conds[j])]] <- res
    }
  }

  rld <- rlog(dds_current, blind = FALSE)
  resultsNames(diff)

  return(list(
    res = do.call(rbind, res_all),
    df_all_genes = do.call(rbind, df_all_genes),
    rld = rld,
    coldata = sample_info
  ))
}

preprocessing_PCA = function(counts, min_reads, sample_info, normalization, min_sample = 5, design_formula = ~ batch + condition, ...) {
  # Filtro geni poco espressi
  counts = counts[rowSums(counts) > min_reads & rowSums(counts > 1) > min_sample, ]
  
  # Allineamento sample_info
  sample_info <- sample_info[colnames(counts), ]
  sample_info$condition <- make.names(sample_info$condition)

  # Crea oggetto DESeq2
  dds_current <- DESeqDataSetFromMatrix(
    countData = counts,
    colData = sample_info,
    design = design_formula
  )
  
  conds <- unique(sample_info$condition)

  if (normalization == "vst") {
    rld <- vst(dds_current, blind = FALSE)
  } else if (normalization == "rlog") {
    rld <- rlog(dds_current, blind = FALSE)
  }

  return(list(
    rld = rld,
    coldata = sample_info
  ))
}

specific_deseq2 = function(counts, min_reads, sample_info, exp, contr, min_sample = 5, design_formula = ~ batch + condition, ...) {
  # Filtro geni poco espressi
  counts = counts[rowSums(counts) > min_reads & rowSums(counts > 1) > min_sample, ]

  message("Dimensioni counts dopo filtro: ", paste(dim(counts), collapse = " x "))
  if (nrow(counts) == 0) stop("Nessun gene passato il filtro. Riduci min_reads o min_sample.")
  
  # Allineamento sample_info
  sample_info <- sample_info[colnames(counts), ]
  sample_info$condition <- make.names(sample_info$condition)

  # Crea oggetto DESeq2
  dds_current <- DESeqDataSetFromMatrix(
    countData = counts,
    colData = sample_info,
    design = design_formula
  )
  
  # DESeq2
  diff <- DESeq(dds_current, ...)
  
  res_all <- list()
  df_all_genes <- list()
  
  contrast <- c("condition", exp, contr)
  res <- results(diff, contrast = contrast)
  res <- data.frame(res)
  res$comparison_exp_vs_contr <- paste0(exp, "_vs_", contr)
  res$gene <- rownames(res)
  res <- na.omit(res)
  rownames(res) <- sub("\\.\\d+$", "", rownames(res))  # rimuove ".1", ".2", etc
  res$symbol <- mapIds(org.Hs.eg.db,
                        keys = rownames(res),
                        column = "SYMBOL",
                        keytype = "ENSEMBL",
                        multiVals = "first")
  res_all[[paste0(exp, "_vs_", contr)]] <- res
  df_all_genes[[paste0(exp, "_vs_", contr)]] <- res

  rld <- rlog(dds_current, blind = FALSE)
  resultsNames(diff)

  return(list(
    res = do.call(rbind, res_all),
    df_all_genes = do.call(rbind, df_all_genes),
    rld = rld,
    coldata = sample_info
  ))
}


mypcaAnalysis <- function(title_1vs2, rld, intgroup, ...) {
  # Calcola PCA
  pcaData <- DESeq2::plotPCA(rld, intgroup = intgroup, returnData = TRUE)
  
  # Estrai la colonna da colData
  plot_info <- colData(rld)[[intgroup]]
  
  # PCA per percentuali
  assayData <- SummarizedExperiment::assay(rld)
  pca <- prcomp(t(assayData))
  percentVar <- pca$sdev^2 / sum(pca$sdev^2)
  pc1_var <- round(100 * percentVar[1], 1)
  pc2_var <- round(100 * percentVar[2], 1)
  
  # Plot
  plotly::plot_ly(
    data = pcaData,
    x = ~PC1,
    y = ~PC2,
    text = ~name,
    color = ~plot_info,
    type = "scatter",
    mode = "markers"
  ) %>%
    layout(
      title = paste("PCA:", title_1vs2),
      xaxis = list(title = paste0("PC1 (", pc1_var, "%)")),
      yaxis = list(title = paste0("PC2 (", pc2_var, "%)"))
    )
}




my_MA_andVolcano_plot = function(title_1vs2, res, qvalue, logfc, ...){
  select <- which(res$padj < qvalue & abs(res$log2FoldChange) > logfc)
  cat("Total number of significant genes:", length(select), "\n")

  res$significant <- "not_significant"
  res$significant[select] <- "significant"
  res$label <- paste0(res$symbol, " (", res$gene, ")")

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

  return(list(ma = ma_plot, volcano = volcano_plot))
}

my_genetable = function(res, title_1vs2, qvalue, logfc, ...){
  res <- res[!duplicated(rownames(res)), ]
  # Aggiungi FoldChange prima di tutto
  res$FoldChange <- 2^res$log2FoldChange
  # Filtro per geni significativi
  deg <- res[which(res$padj < qvalue & abs(res$log2FoldChange) > logfc), ]
  # Inizializzo la colonna con "no"
  res$differentially_expressed <- rep("no", nrow(res))
  # Trovo geni comuni (intersezione)
  common_genes <- intersect(rownames(res), rownames(deg))

  # Segno come "yes" solo i geni comuni
  res[rownames(res) %in% common_genes, "differentially_expressed"] <- "yes"
  # Tabella interattiva di tutti i geni (se vuoi mostrarla)
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
  # Tabella interattiva solo dei geni differenzialmente espressi
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
  title_1vs2_underscore <- gsub(" ", "_", title_1vs2)
  # Esportazione CSV
  write.table(deg, file = paste("./output/gse_results/DEG_", title_1vs2_underscore, ".csv"), sep = ",", row.names = FALSE, quote = FALSE)
  return(list(res_dt = res_dt, deg_dt = deg_dt, res = res, deg = deg))
}



my_heatmaps = function(deg, rld, title_1vs2) {
  if (!requireNamespace("ComplexHeatmap", quietly = TRUE)) stop("Install ComplexHeatmap package.")
  gene_symbols <- setNames(deg$symbol, deg$gene)

  # Seleziona i top 20 geni significativi ordinati per padj
  top20_genes <- deg %>%
    as.data.frame() %>%
    arrange(padj) %>%
    head(20) %>%
    pull(gene)

  expr_top20 <- assay(rld)[top20_genes, ]

  new_rownames <- coalesce(gene_symbols[rownames(expr_top20)], rownames(expr_top20))
  rownames(expr_top20) <- new_rownames

  expr_top20_scaled <- t(scale(t(expr_top20)))

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

  # Tutti i DE genes
  all_de_genes <- deg %>% pull(gene)
  expr_de <- assay(rld)[all_de_genes, ]

  new_rownames <- coalesce(gene_symbols[rownames(expr_de)], rownames(expr_de))
  rownames(expr_de) <- new_rownames

  expr_de_scaled <- t(scale(t(expr_de)))

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
  return(list(h_20 = h_20, h_all = h_all))
}




# Funzione riutilizzabile per GSEA e salvataggio
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

# Funzione riutilizzabile per GSEA e salvataggio
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


goenrichment = function(markers, logfc_col = "avg_log2FC", down = FALSE, n_genes = 100){
  if(down){
    markers[[logfc_col]] = -markers[[logfc_col]]
  }
  
  up = markers[markers[[logfc_col]] > 0, ]
  up = up[order(up[[logfc_col]], decreasing = TRUE), ]
  
  n_available <- nrow(up)
  
  if(is.character(n_genes) && tolower(n_genes) == "all"){
    marker_genes <- up$gene
  } else if (is.numeric(n_genes)) {
    if(n_genes <= 0){
      stop("n_genes deve essere un numero positivo o 'all'")
    }
    if(n_genes >= n_available){
      marker_genes <- up$gene
    } else {
      marker_genes <- head(up$gene, n_genes)
    }
  } else {
    stop("n_genes deve essere un numero positivo o 'all'")
  }
  
  entrez_ids <- mapIds(
    org.Hs.eg.db,
    keys = marker_genes,
    column = "ENTREZID",
    keytype = "ENSEMBL",
    multiVals = "first"
  )
  entrez_ids <- na.omit(entrez_ids)
  
  go_enrichment <- enrichPathway(
    gene = entrez_ids,
    organism = "human",
    pAdjustMethod = "BH",
    pvalueCutoff = 0.05,
    qvalueCutoff = 0.05,
    readable = TRUE
  )
  
  go_enrichment@result = go_enrichment@result[
    go_enrichment@result$p.adjust < go_enrichment@pvalueCutoff &
    go_enrichment@result$qvalue < go_enrichment@qvalueCutoff, ]
  
  return(go_enrichment)
}





my_WGCNA_analysis <- function(sample_info, dds, dds_subset = NULL) {

  # Normalizzazione e batch effect
  if (!is.null(dds_subset)) {
    vst_counts <- assay(vst(dds_subset))
    vst_bc <- removeBatchEffect(vst_counts, batch = dds_subset$batch)
    final_counts <- vst_bc
  } else {
    final_counts <- assay(dds)
  }

  # Rimozione outlier genetici
  gsg <- goodSamplesGenes(t(final_counts), verbose = 3)
  if (!gsg$allOK) {
    final_counts <- final_counts[gsg$goodSamples, ]  # Filtra righe (geni)
    keep <- rowSums(counts(dds)[gsg$goodSamples, ] >= 10) >= 6
  } else {
    keep <- rowSums(counts(dds) >= 10) >= 6
  }
  filtered_counts <- t(final_counts[keep, ])

  wgcna_data <- final_counts[gsg$goodGenes == TRUE,]

  pca <- prcomp(t(wgcna_data))
  pca.dat <- pca$x
  #Obtain the variance from standard deviation 
  pca.var <- pca$sdev^2 #standard deviations of the principal components
  #Obtain the percentage of the variance
  pca.var.percent <- round(pca.var/sum(pca.var)*100, digits = 2)
  pca.dat <- as.data.frame(pca.dat)
  # Plot the first two principal components
  ggplot_PCA <- ggplot(pca.dat, aes(PC1, PC2)) +
    geom_point() +
    geom_text(label = rownames(pca.dat)) +
    labs(x = paste0('PC1: ', pca.var.percent[1], ' %'), #the percentage of variance explained by PC1
        y = paste0('PC2: ', pca.var.percent[2], ' %')) #the percentage of variance explained by PC2


  # Soft thresholding
  powers <- c(1:10, seq(12, 50, 2))
  sft <- pickSoftThreshold(filtered_counts, powerVector = powers, networkType = "signed", verbose = 5)
  soft_power <- 14

  sft.data <- sft$fitIndices

  # visualization to pick power

  a1 <- ggplot(sft.data, aes(Power, SFT.R.sq, label = Power)) +
    geom_point() +
    geom_text(nudge_y = 0.1) +
    geom_hline(yintercept = 0.8, color = 'red') +
    labs(x = 'Power', y = 'Scale free topology model fit, signed R^2') +
    theme_classic()


  a2 <- ggplot(sft.data, aes(Power, mean.k., label = Power)) +
    geom_point() +
    geom_text(nudge_y = 0.1) +
    labs(x = 'Power', y = 'Mean Connectivity') +
    theme_classic()

  power_plots <- grid.arrange(a1, a2, nrow = 2)

  # Costruzione rete
  filtered_counts[] <- sapply(filtered_counts, as.numeric)
  net <- blockwiseModules(filtered_counts,
                          power = soft_power,
                          TOMType = "signed",
                          maxBlockSize = 5000,
                          numericLabels = FALSE,
                          randomSeed = 1234,
                          verbose = 3)

  # Calcolo eigengenes
  MEs <- net$MEs

  # Regressione con tutti i livelli della condizione (multi-gruppo)
  sample_info$condition <- as.factor(sample_info$condition)
  design <- model.matrix(~ 0 + condition, data = sample_info)  # evita riferimento implicito
  colnames(design) <- levels(sample_info$condition)

  fit <- lmFit(t(MEs), design)
  fit <- eBayes(fit)

  # Statistiche per tutti i contrasti (comparazioni)
  stats_df <- topTable(fit, number = Inf, adjust = "fdr") %>%
    rownames_to_column("module")

  # Calcolo correlazione moduli-condizione (codificata come dummy matrix)
  traits <- model.matrix(~ 0 + condition, data = sample_info)
  mod_trait_cor <- cor(MEs, traits, use = "p")
  pvals <- corPvalueStudent(mod_trait_cor, nrow(filtered_counts))

  # Heatmap
  heatmap_data <- cbind(MEs, traits)

  CLP <- CorLevelPlot(heatmap_data,
               x = colnames(traits),
               y = names(MEs),
               col = c("blue1", "skyblue", "white", "pink", "red"),
               rotTitleX = 90
               )

  # Output
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



