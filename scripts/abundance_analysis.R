#@title 3.2 Abundace different analysis
#create folder for analysis
# File: abundance_analysis.R

# Create folders for analysis
dir.create("/content/drive/MyDrive/output_funcscan/work_dir/abundance_analysis/vocano_plot/family", recursive = TRUE)
dir.create("/content/drive/MyDrive/output_funcscan/work_dir/abundance_analysis/vocano_plot/genus", recursive = TRUE)
dir.create("/content/drive/MyDrive/output_funcscan/work_dir/abundance_analysis/vocano_plot/species", recursive = TRUE)

# Load required R script for analysis
source('/content/ancom_BC_test.r')

# Set paths
feature_table_path <- '/content/drive/MyDrive/output_funcscan/work_dir/kraken_visualize/kraken_feature_no_header_table.tsv'

bracken_list_path <- c('/content/drive/MyDrive/output_funcscan/work_dir/bracken_output/merged/merged_bracken_F.txt',
                       '/content/drive/MyDrive/output_funcscan/work_dir/bracken_output/merged/merged_bracken_G.txt',
                       '/content/drive/MyDrive/output_funcscan/work_dir/bracken_output/merged/merged_bracken_S.txt')

# Metadata with groups
metadata <- data.frame('group' = c('Non-fertilize phase 1','Non-fertilize phase 1','Non-fertilize phase 1',
                                  'VC phase 2','VC phase 2','VC phase 2',
                                  'HC phase 2','HC phase 2','HC phase 2',
                                  'HCVS phase 2','HCVS phase 2','HCVS phase 2',
                                  'VC phase 3','VC phase 3','VC phase 3',
                                  'HC phase 3','HC phase 3','HC phase 3',
                                  'HCVS phase 3','HCVS phase 3','HCVS phase 3'))

# Output directory for results
output_dir <- '/content/drive/MyDrive/output_funcscan/work_dir/abundance_analysis/vocano_plot'

# Run the abundance analysis function
abuncance_analysis(feature_table_path, bracken_list_path, metadata, output_dir)

