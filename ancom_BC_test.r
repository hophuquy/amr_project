### load library
library(EnhancedVolcano)
library(dplyr)
library(ANCOMBC)
library(tidyverse)
library(ggplot2)

### Function to perform the analysis
run_ancombc_analysis <- function(data_file, metadata, output_dir) {

  # Load the data
  df <- read.csv(data_file, sep="\t")
  d <- nrow(df)
  df2 <- df[, 2:ncol(df)]
  rownames(df2) <- paste0("T", seq_len(d))
  colnames(df2) <- paste0("S", seq_len(ncol(df2)))

  matrix <- as.matrix(df2)

  # Prepare metadata
  rownames(metadata) <- paste0("S", seq_len(ncol(df2)))

  # Create TreeSummarizedExperiment object
  assays <- S4Vectors::SimpleList(counts = matrix)
  smd <- S4Vectors::DataFrame(metadata)
  tse <- TreeSummarizedExperiment::TreeSummarizedExperiment(assays = assays, colData = smd)

  # Run ANCOM-BC
  set.seed(123)
  output <- ancombc2(data = tse, assay_name = "counts", tax_level = NULL,
                     fix_formula = "group", rand_formula = NULL,
                     p_adj_method = "holm", pseudo_sens = TRUE,
                     prv_cut = 0.10, lib_cut = 1000, s0_perc = 0.05,
                     group = "group", struc_zero = FALSE, neg_lb = FALSE,
                     alpha = 0.05, n_cl = 2, verbose = TRUE,
                     global = TRUE, pairwise = TRUE,
                     dunnet = TRUE, trend = FALSE,
                     iter_control = list(tol = 1e-2, max_iter = 20,
                                         verbose = FALSE),
                     em_control = list(tol = 1e-5, max_iter = 100),
                     lme_control = NULL,
                     mdfdr_control = list(fwer_ctrl_method = "holm", B = 100),
                     trend_control = NULL)

  res_prim <- output$res
  res_pair <- output$res_pair
  res_global <- output$res_global

  # Save results to CSV files
  write.csv(res_pair, file = file.path(output_dir, "res_pair.csv"), row.names = FALSE)

  # Return important results as a list
  return(list(
    res_prim = res_prim,
    res_pair = res_pair,
    res_global = res_global,
    t_test_results = t_test_results_df,
    plot = p
  ))
}

# draw Volcano_plot
draw_EnhancedVolcano_plot <- function(res_pair,output_path = ".") {
  ### Compare groups
  compare_groups <- c(
      "lfc_groupHCVS.phase.2_groupHC.phase.3",
      "lfc_groupHCVS.phase.3_groupHC.phase.3",
      "lfc_groupNon.fertilize.phase.1_groupHC.phase.3",
      "lfc_groupVC.phase.2_groupHC.phase.3",
      "lfc_groupVC.phase.3_groupHC.phase.3",
      "lfc_groupHCVS.phase.3_groupHCVS.phase.2",
      "lfc_groupNon.fertilize.phase.1_groupHCVS.phase.2",
      "lfc_groupVC.phase.2_groupHCVS.phase.2",
      "lfc_groupVC.phase.3_groupHCVS.phase.2",
      "lfc_groupNon.fertilize.phase.1_groupHCVS.phase.3",
      "lfc_groupVC.phase.2_groupHCVS.phase.3",
      "lfc_groupVC.phase.3_groupHCVS.phase.3",
      "lfc_groupVC.phase.2_groupNon.fertilize.phase.1",
      "lfc_groupVC.phase.3_groupNon.fertilize.phase.1",
      "lfc_groupVC.phase.3_groupVC.phase.2"
    )

  # Replace 'lfc_group' with 'p_group' in each element
  compare_groups_replaced <- gsub('lfc_group', 'p_group', compare_groups)
  print('pass')
  for (i in 1:length(compare_groups)) {
    filted_res_pair <- res_pair[, c(compare_groups[i], compare_groups_replaced[i])]
    rownames(filted_res_pair) <- res_pair$taxon
    colnames(filted_res_pair) <- c('log2FoldChange', 'padj')

    # Create the EnhancedVolcano plot
    volcano_plot <- EnhancedVolcano(filted_res_pair,
      lab = rownames(filted_res_pair),
      title = paste("Volcano Plot:", i),
      subtitle = 'log2FoldChange cutoff: ±1, adjusted P-value cutoff: 0.05',
      x = 'log2FoldChange',
      y = 'padj',
      xlab = 'Log2 Fold Change',
      ylab = '-Log10 Adjusted P-value',
      xlim = c(min(filted_res_pair$log2FoldChange), 2),
      ylim = c(0, -log10(1e-9)),
      pCutoff = 0.05,
      FCcutoff = 1.0,
      pointSize = 2,
      labSize = 4,
      legendLabels = c('Non-significant', 'Log2FC', 'Adjusted P-value', 'Adjusted P-value & Log2FC'),
      caption = paste("Comparison:", compare_groups[i]),
      legendPosition = 'top',
      boxedLabels = TRUE,
      drawConnectors = TRUE
    )

    # Construct the output file path
    safe_filename <- gsub("[^[:alnum:]_]", "_", compare_groups[i])  # Ensure valid filenames
    full_output_path <- file.path(output_path, paste0("volcano_plot_", safe_filename, ".png"))

    # Save the plot to the specified output path
    ggsave(full_output_path, plot = volcano_plot, width = 8, height = 6, dpi = 300)
  }
}

#### Function to perform the analysis
abuncance_analysis <- function(feature_table_path, bracken_list_path, metadata, output_dir) {
  print('----###########----')
  print('----RUN ANCOMBC----')
  ### run ancombc
  # ancom_result = run_ancombc_analysis(feature_table_path, metadata, output_dir)

  ### draw Volcano_plot
  # load data
  # res_pair <- ancom_result['res_pair']
  print('pass1')
  res_pair <- read.csv('/content/res_pair.csv')
  res_pair <- na.omit(res_pair)
  print('pass2')


  df <- read.csv(feature_table_path, sep="\t")
  df$taxon <- paste0("T", seq_len(nrow(df)))
  print('pass3')
  res_pair_merged <- merge(res_pair, df[, c('tax_id', 'taxon')], by = 'taxon', all = FALSE)
   print('pass4')
  # Load braken taxonomy level
  tax_F_level <- read.csv(bracken_list_path[1], sep="\t")
  tax_G_level <- read.csv(bracken_list_path[2], sep="\t")
  tax_S_level <- read.csv(bracken_list_path[3], sep="\t")
   print('pass5')
  # Filter in each level
  res_pair_F <- res_pair_merged %>% filter(tax_id %in% tax_F_level$taxonomy_id)
  res_pair_G <- res_pair_merged %>% filter(tax_id %in% tax_G_level$taxonomy_id)
  res_pair_S <- res_pair_merged %>% filter(tax_id %in% tax_S_level$taxonomy_id)
  list_path <- c('family', 'genus', 'species')
  res_pair_tax_list <- list(res_pair_F, res_pair_G, res_pair_S)

  print('----###################----')
  print('----RUN EnhancedVolcano----')

  for (i in 1:length(list_path)) {
    print('pass6')
    output_path <- paste0(output_dir, "/", list_path[i])
    res_pair_element <- res_pair_tax_list[[i]]

    draw_EnhancedVolcano_plot(res_pair_element, output_path)
  }
}

# TEST RUNING
# feature_table_path <- 'C:/Users/Admin/OneDrive/Documents/R_tutorial/kraken_table_no_header.tsv'
# bracken_list_path <- c('C:/Users/Admin/OneDrive/Documents/R_tutorial/braken_output/merged_bracken_F.txt', 'C:/Users/Admin/OneDrive/Documents/R_tutorial/braken_output/merged_bracken_G.txt', 'C:/Users/Admin/OneDrive/Documents/R_tutorial/braken_output/merged_bracken_S.txt')
# metadata <- data.frame('group' = c('Non-fertilize phase 1','Non-fertilize phase 1','Non-fertilize phase 1',
#                                   'VC phase 2','VC phase 2','VC phase 2',
#                                   'HC phase 2','HC phase 2','HC phase 2',
#                                   'HCVS phase 2','HCVS phase 2','HCVS phase 2',
#                                   'VC phase 3','VC phase 3','VC phase 3',
#                                   'HC phase 3','HC phase 3','HC phase 3',
#                                   'HCVS phase 3','HCVS phase 3','HCVS phase 3'))
# output_dir <- 'C:/Users/Admin/OneDrive/Documents/R_tutorial/test'

# abuncance_analysis(feature_table_path, bracken_list_path, metadata, output_dir)
