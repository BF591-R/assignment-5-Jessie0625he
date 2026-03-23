library('tidyverse')
library('SummarizedExperiment')
library('DESeq2')
library('biomaRt')
library('testthat')
library('fgsea')

#' Function to generate a SummarizedExperiment object with counts and coldata
#' to use in DESeq2
#'
#' @param csv_path (str): path to the file verse_counts.tsv
#' @param metafile (str): path to the metadata sample_metadata.csv
#' @param selected_times (list): list of sample timepoints to use
#' @return SummarizedExperiment object with subsetted counts matrix
#'   and sample data. Ensure that the timepoints column used as input 
#'   to the model design has 'vP0' set as the reference factor level. Your 
#'   colData dataframe should have columns named samplename and timepoint.
#' @export
#'
#' @examples se <- make_se('verse_counts.tsv', 'sample_metadata.csv', c('vP0', 'vAd'))
make_se <- function(counts_csv, metafile_csv, selected_times) {
  # read counts table
  counts_df <- readr::read_delim(counts_csv, delim = "\t", show_col_types = FALSE)
  
  # read metadata file
  meta_data <- read.csv(metafile_csv, header = TRUE, sep = ",")
  
  # subset metadata using selected_times
  colData <- meta_data[meta_data$timepoint %in% selected_times, c("samplename", "timepoint")]
  
  # set vP0 as reference level
  colData$timepoint <- factor(
    colData$timepoint,
    levels = c("vP0", setdiff(selected_times, "vP0"))
  )
  
  # subset counts using sample names from metadata
  filtered_counts <- counts_df[, c("gene", colData$samplename)]
  
  # create rowData
  rowData <- filtered_counts["gene"]
  
  # remove gene column and assign row names
  counts_mat <- as.data.frame(filtered_counts[, -1])
  rownames(counts_mat) <- rowData$gene
  
  # check sample order matches
  stopifnot(identical(colnames(counts_mat), colData$samplename))
  
  # create SummarizedExperiment
  se <- SummarizedExperiment(
    assays = list(counts = as.matrix(counts_mat)),
    colData = colData,
    rowData = rowData
  )
  return(se)
}

#' Function that runs DESeq2 and returns a named list containing the DESeq2
#' results as a dataframe and the dds object returned by DESeq2
#'
#' @param se (obj): SummarizedExperiment object containing counts matrix and
#' coldata
#' @param design: the design formula to be used in DESeq2
#'
#' @return list with DESeqDataSet object after running DESeq2 and results from
#'   DESeq2 as a dataframe
#' @export
#'
#' @examples results <- return_deseq_res(se, ~ timepoint)
return_deseq_res <- function(se, design) {
  dds <- DESeqDataSet(se, design = design)
  dds <- DESeq(dds)
  deseq2_res <- as.data.frame(DESeq2::results(dds))
  
  return(list(
    dds = dds,
    results = deseq2_res
  ))
}

#' Function that takes the DESeq2 results dataframe, converts it to a tibble and
#' adds a column to denote plotting status in volcano plot. Column should denote
#' whether gene is either 1. Significant at padj < .10 and has a positive log
#' fold change, 2. Significant at padj < .10 and has a negative log fold change,
#' 3. Not significant at padj < .10. Have the values for these labels be UP,
#' DOWN, NS, respectively. The column should be named `volc_plot_status`. Ensure
#' that the column name for your rownames is called "genes". 
#'
#' @param deseq2_res (df): results from DESeq2 
#' @param padj_threshold (float): threshold for considering significance (padj)
#'
#' @return Tibble with all columns from DESeq2 results and one additional column
#'   labeling genes by significant and up-regulated, significant and
#'   downregulated, and not significant at padj < .10.
#'   
#' @export
#'
#' @examples labeled_results <- label_res(res, .10)
label_res <- function(deseq2_res, padj_threshold) {
  res_tbl <- as.data.frame(deseq2_res) %>%
    tibble::rownames_to_column(var = "genes") %>%
    tibble::as_tibble() %>%
    dplyr::mutate(
      volc_plot_status = dplyr::case_when(
        !is.na(padj) & padj < padj_threshold & log2FoldChange > 0 ~ "UP",
        !is.na(padj) & padj < padj_threshold & log2FoldChange < 0 ~ "DOWN",
        TRUE ~ "NS"
      )
    )
  
  return(res_tbl)
}

#' Function to plot the unadjusted p-values as a histogram
#'
#' @param labeled_results (tibble): Tibble with DESeq2 results and one additional
#' column denoting status in volcano plot
#'
#' @return ggplot: a histogram of the raw p-values from the DESeq2 results
#' @export
#'
#' @examples pval_plot <- plot_pvals(labeled_results)
plot_pvals <- function(labeled_results) {
    pval_plot <- ggplot(labeled_results, ggplot2::aes(x = pvalue)) +
      geom_histogram(bins = 30, color = "black", fill = "skyblue") +
      labs(
        x = "p-value",
        y = "Count",
        title = "Histogram of DESeq2 raw p-values (vP0 vs vAd)"
      ) +
      ggplot2::theme_minimal()
    return(pval_plot)
}

#' Function to plot the log2foldchange from DESeq2 results in a histogram
#'
#' @param labeled_results (tibble): Tibble with DESeq2 results and one additional
#' column denoting status in volcano plot
#' @param padj_threshold (float): threshold for considering significance (padj)
#'
#' @return ggplot: a histogram of log2FC values from genes significant at padj 
#' threshold of 0.1
#' @export
#'
#' @examples log2fc_plot <- plot_log2fc(labeled_results, .10)
plot_log2fc <- function(labeled_results, padj_threshold) {
  filtered_res <- labeled_results %>% 
    filter(!is.na(padj) & padj < padj_threshold)
  
  log2fc_plot <- ggplot(filtered_res, aes(x = log2FoldChange)) +
    geom_histogram(bins = 30, color = "black", fill = "salmon") +
    labs(
      x = "log2FoldChange",
      y = "Count",
      title = "Histogram of log2FC for DE genes (vP0 vs vAd)"
    ) +
    theme_minimal()
  
  return(log2fc_plot)
}

#' Function to make scatter plot of normalized counts for top ten genes ranked
#' by ascending padj
#'
#' @param labeled_results (tibble): Tibble with DESeq2 results and one
#'   additional column denoting status in volcano plot
#' @param dds_obj (obj): The object returned by running DESeq (dds) containing
#' the updated DESeqDataSet object with test results
#' @param num_genes (int): Number of genes to plot
#'
#' @return ggplot: a scatter plot with the normalized counts for each sample for
#' each of the top ten genes ranked by ascending padj
#' @export
#'
#' @examples norm_counts_plot <- scatter_norm_counts(labeled_results, dds, 10)
scatter_norm_counts <- function(labeled_results, dds_obj, num_genes){
  if (is.null(dds_obj)) {
    stop("dds_obj is NULL. Check return_deseq_res().")
  }
  
  top_genes <- labeled_results %>%
    filter(!is.na(padj)) %>% 
    arrange(padj) %>%
    slice_head(n = num_genes) %>% 
    pull(genes)
  
  # extract normalized counts from DESeq2 object
  norm_counts <- counts(dds_obj, normalized = TRUE)
  
  # subset to top genes
  norm_counts_df <- as.data.frame(norm_counts[top_genes, , drop = FALSE])
  norm_counts_df$genes <- rownames(norm_counts_df)
  
  
  # convert to long format
  plot_df <- norm_counts_df %>%
      pivot_longer(
        cols = -genes,
        names_to = "samplenames",
        values_to = "norm_counts"
    ) %>% 
      mutate(log10_norm_counts = log10(norm_counts + 1))
  
  norm_counts_plot <- ggplot(plot_df, aes(x = genes, y = log10_norm_counts, color = samplenames)) +
    geom_point(size = 3, position = ggplot2::position_jitter(width = 0.1)) +
    labs(
      title = paste("Plot of Log10(normalized counts) for top", num_genes, "DE genes"),
      x = NULL,
      y = "log10(norm_counts)",
      color = "samplenames"
    ) +
    ggplot2::theme_bw() +
    ggplot2::theme(
      axis.text.x = ggplot2::element_text(angle = 90, hjust = 1)
    )
  
  return(norm_counts_plot)
}

#' Function to generate volcano plot from DESeq2 results
#'
#' @param labeled_results (tibble): Tibble with DESeq2 results and one
#'   additional column denoting status in volcano plot
#'
#' @return ggplot: a scatterplot (volcano plot) that displays log2foldchange vs
#'   -log10(padj) and labeled by status
#' @export
#'
#' @examples volcano_plot <- plot_volcano(labeled_results)
#' 
plot_volcano <- function(labeled_results) {
  volcano_plot <- labeled_results %>%
    filter(!is.na(padj), !is.na(log2FoldChange)) %>%
    ggplot(aes(x = log2FoldChange,y = -log10(padj),color = volc_plot_status)) +
    geom_point() +
    labs(
      x = "log2 Fold Change",
      y = "-log10 adjusted p-value",
      color = "volc_plot_status"
    ) +
    theme_minimal()
  return(volcano_plot)
}

#' Function to generate a named vector ranked by log2FC descending
#'
#' @param labeled_results (tibble): Tibble with DESeq2 results and one
#'   additional column denoting status in volcano plot
#' @param id2gene_path (str): Path to the file containing the mapping of
#' ensembl IDs to MGI symbols
#'
#' @return Named vector with gene symbols as names, and log2FoldChange as values
#' ranked in descending order
#' @export
#'
#' @examples rnk_list <- make_ranked_log2fc(labeled_results, 'data/id2gene.txt')

make_ranked_log2fc <- function(labeled_results, id2gene_path) {
  id2gene <- read.csv(id2gene_path, header = FALSE, sep = "\t", stringsAsFactors = FALSE)
  colnames(id2gene) <- c("genes", "mgi_symbol")
  
  ranked_df <- as.data.frame(labeled_results)[, c("genes", "log2FoldChange")] %>%
    dplyr::left_join(id2gene, by = "genes") %>%
    dplyr::filter(!is.na(mgi_symbol), mgi_symbol != "", !is.na(log2FoldChange)) %>%
    dplyr::arrange(dplyr::desc(log2FoldChange)) %>%
    dplyr::distinct(mgi_symbol, .keep_all = TRUE)
  
  rnk_list <- ranked_df$log2FoldChange
  names(rnk_list) <- ranked_df$mgi_symbol
  return(rnk_list)
}


#' Function to run fgsea with arguments for min and max gene set size
#'
#' @param gmt_file_path (str): Path to the gene sets of interest in GMT format
#' @param rnk_list (named vector): Named vector generated previously with gene 
#' symbols and log2Fold Change values in descending order
#' @param min_size (int): Minimum number of genes in gene sets to be allowed
#' @param max_size (int): Maximum number of genes in gene sets to be allowed
#'
#' @return Tibble of results from running fgsea
#' @export
#'
#' @examples fgsea_results <- run_fgsea('data/m2.cp.v2023.1.Mm.symbols.gmt', rnk_list, 15, 500)

run_fgsea <- function(gmt_file_path, rnk_list, min_size, max_size) {
  pathways <- fgsea::gmtPathways(gmt_file_path)
  
  fgsea_results <- fgsea::fgsea(
    pathways = pathways,
    stats = rnk_list,
    minSize = min_size,
    maxSize = max_size
  )
  
  fgsea_results <- tibble::as_tibble(fgsea_results)
  return(fgsea_results)
}
  

#' Function to plot top ten positive NES and top ten negative NES pathways
#' in a barchart
#'
#' @param fgsea_results (tibble): the fgsea results in tibble format returned by
#'   the previous function
#' @param num_paths (int): the number of pathways for each direction (top or
#'   down) to include in the plot. Set this at 10.
#'
#' @return ggplot with a barchart showing the top twenty pathways ranked by positive
#' and negative NES
#' @export
#'
#' @examples fgsea_plot <- top_pathways(fgsea_results, 10)
top_pathways <- function(fgsea_results, num_paths) {
  top_pos <- fgsea_results %>%
    dplyr::filter(!is.na(NES), NES > 0) %>%
    dplyr::arrange(dplyr::desc(NES)) %>%
    dplyr::slice_head(n = num_paths)
  
  top_neg <- fgsea_results %>%
    dplyr::filter(!is.na(NES), NES < 0) %>%
    dplyr::arrange(NES) %>%
    dplyr::slice_head(n = num_paths)
  
  plot_df <- dplyr::bind_rows(top_pos, top_neg) %>%
    dplyr::mutate(
      NES_sign = ifelse(NES > 0, "positive", "negative"),
      pathway = factor(pathway, levels = pathway[order(NES)])
    )
  
  fgsea_plot <- ggplot2::ggplot(
    plot_df,
    ggplot2::aes(x = pathway, y = NES, fill = NES_sign)
  ) +
    ggplot2::geom_col(width = 0.8) +
    ggplot2::coord_flip() +
    ggplot2::scale_fill_manual(values = c("positive" = "firebrick", "negative" = "steelblue")) +
    ggplot2::labs(
      x = "Pathway",
      y = "NES",
      fill = "NES sign",
      title = paste("Top", num_paths, "positive and negative pathways")
    ) +
    ggplot2::theme_bw(base_size = 10)
  
  return(fgsea_plot)
}

