import pandas as pd
import numpy as np
import statsmodels.api as sm
from scipy.stats import skew
from statsmodels.formula.api import ols
from statsmodels.stats.multicomp import pairwise_tukeyhsd
from scipy.stats import yeojohnson, boxcox
from skbio.diversity import alpha_diversity, beta_diversity
from skbio.stats.ordination import pcoa
import matplotlib.pyplot as plt
import seaborn as sns
from scipy.spatial.distance import pdist, squareform
from skbio.stats.distance import permanova
import skbio
import subprocess
import os
from rpy2.robjects import r, globalenv
from rpy2.robjects.packages import importr
from rpy2.robjects import pandas2ri
from bs4 import BeautifulSoup
import json
from ete3 import NCBITaxa

###### Statistic test ######

# Function to apply Yeojohnson or Boxcox transformation and normalization to a column
def normalize(column, trans_type='yeojohnso'):
    # Handle zero and negative values before applying Boxcox
    if trans_type == 'boxcox':
        # Shift values to be positive if any are zero or negative
        min_value = column.min()
        if min_value <= 0:
            column = column - min_value + 1e-6  # Add a small positive constant

    if trans_type == 'yeojohnso':
        transformed_data, _ = yeojohnson(column)
    elif trans_type == 'boxcox':
        transformed_data, _ = boxcox(column)
    else:
        raise ValueError("Invalid transformation type. Use 'yeojohnso' or 'boxcox'.")
    normalized_data = (transformed_data - np.min(transformed_data)) / (np.max(transformed_data) - np.min(transformed_data))
    return normalized_data

def anova_test(alpha_div, metadata):

  # convert to long data frame
  df_alpha_simpson = pd.DataFrame({'sample':metadata,
                                   'shannon':[float(i) for i in alpha_div[0]],
                                   'pielou_e':[float(i) for i in alpha_div[1]],
                                   'chao1':[float(i) for i in alpha_div[2]],
                                   'inv_simpson':[float(i) for i in alpha_div[3]]})
  skewness_results_before = df_alpha_simpson.set_index('sample').apply(skew)

  # normalize data if not standard distribution
  normalized_alpha_table = df_alpha_simpson.set_index('sample').apply(lambda col: normalize(col, trans_type='yeojohnso'), axis=0)
  skewness_results_after = normalized_alpha_table.apply(skew)


  df = normalized_alpha_table
  df['sample'] = metadata
  # Step 2: Perform ANOVA
  index = [i for i in df.columns.to_list() if i != 'sample']
  tukey_results = {}
  anova_tables = {}
  for i in index:
    # print(i)
    model = ols(f'{i} ~ sample', data=df).fit()
    anova_table = sm.stats.anova_lm(model, typ=2)

    # Step 3: Perform Tukey's HSD test
    tukey_result = pairwise_tukeyhsd(df[i], df['sample'], alpha=0.05)
    tukey_results[i] = tukey_result
    anova_tables[i] = anova_table

  return tukey_results, anova_tables, skewness_results_before, skewness_results_after
def permanova_test(df):

  # Calculate the distance matrix
  coords = df[['PC1', 'PC2', 'PC3']]
  distance_matrix = squareform(pdist(coords, metric='euclidean'))

  # Convert to skbio DistanceMatrix
  distance_matrix = skbio.DistanceMatrix(distance_matrix, ids=df.index)

  # Perform PERMANOVA
  group = df['group']
  result_permanova = permanova(distance_matrix, group, permutations=100000)
  return result_permanova


# pair_wise_list = [filter1, filter2, filter3]
# result_permanova_pair = []
# for i in range(3):
#   print(i)
#   if i == 0:
#     result_permanova_pair.append(permanova_test(pca_df[filter1]))
#   elif i == 1:
#     result_permanova_pair.append(permanova_test(pca_df[filter2]))
#   else:
#     result_permanova_pair.append(permanova_test(pca_df[filter3]))



###### Diversity ######

def diversity_calculation(file_path):
  # Read count table into a Pandas DataFrame
  count_table = pd.read_csv(file_path, sep='\t')
  count_table = count_table.iloc[:,1:]
  # Transpose the DataFrame to have samples as rows and features as columns
  # Scikit-bio only accepts this format
  transposed_table = count_table.T

  # Convert count values to numeric
  transposed_table = transposed_table.apply(pd.to_numeric, errors='coerce')
  # Remove rows with NaN values (if any)
  transposed_table = transposed_table.dropna()

  # Calculate alpha diversity (Shannon entropy in this example)

  alpha_shannon = alpha_diversity('shannon', transposed_table, ids=transposed_table.index)
  alpha_pielou_e = alpha_diversity('pielou_e', transposed_table, ids=transposed_table.index)
  alpha_chao1 = alpha_diversity('chao1', transposed_table, ids=transposed_table.index)
  alpha_inv_simpson = alpha_diversity('inv_simpson', transposed_table, ids=transposed_table.index)
  alpha_div =  (alpha_shannon, alpha_pielou_e, alpha_chao1, alpha_inv_simpson)
  # calculate betadiversity

  beta_jaccard = beta_diversity('jaccard', transposed_table)
  beta_braycurtis = beta_diversity('braycurtis', transposed_table)
  beta_div = {'jaccard':beta_jaccard, 'braycurtis':beta_braycurtis}
  pcoa_results = {}
  # Perform PCoA
  for i in  beta_div:
    pcoa_results[i] = pcoa( beta_div[i], number_of_dimensions=3)
  return  alpha_div, beta_div, pcoa_results

def alpha_plot(alpha_div, metadata, output_path):

  data = {
    # 'Group': ['SM1032307A01',	'SM1032307A02',	'SM1032307A03',	'SM1032308A01',	'SM1032308A02',	'SM1032308A03',	'SM1032308A04',	'SM1032308A05',	'SM1032308A06',	'SM1032308A07',	'SM1032308A08',	'SM1032308A09',	'SM1032309A01',	'SM1032309A02',	'SM1032309A03',	'SM1032309A04',	'SM1032309A05',	'SM1032309A06',	'SM1032309A07', 'SM1032309A08','SM1032309A09'],
    'Group': metadata,
    'shannon': alpha_div[0],
    'pielou_e': alpha_div[1],
    'chao1': alpha_div[2],
    'inv_simpson': alpha_div[3]
}
  df = pd.DataFrame(data)

  tukey_results, anova_tables, skewness_results_before, skewness_results_after = anova_test(alpha_div, metadata)
  # Set the aesthetic style of the plots
  sns.set(style="whitegrid")

  # Create a figure with one subplot (you had 1, 1 which is equivalent to one subplot)_
  index = [i for i in df.columns.to_list() if i != 'Group']

  for i in index:
    # Apply the classic style
    plt.style.use('ggplot')

    # Chao1 Richness Plot
    fig, ax = plt.subplots(1, 1, figsize=(14, 7))
    sns.scatterplot(x='Group', y= i, hue='Group', data= df, s=100, ax=ax)
    if i == 'shannon'  or i == 'inv_simpson':
      ax.set_title(f"{i} Index of diversity")
      ax.set_ylabel("Diversity")
      ax.set_xlabel("")
    elif i == 'chao1':
      ax.set_title(f"{i} Index of richness")
      ax.set_ylabel("Richness")
      ax.set_xlabel("")
    else:
      ax.set_title(f"{i} Index of eveness")
      ax.set_ylabel("Eveness")
      ax.set_xlabel("")

    # Customize the legend to be outside the chart
    legend = ax.legend(title='Group', bbox_to_anchor=(1.05, 1), loc='upper left', fontsize='medium')
    legend.get_title().set_fontsize('medium')  # Set legend title font size


    # Add ANOVA results as text in the plot
    anova_text = f"ANOVA:\nF={round(anova_tables[i].iloc[0, 2],2)}, p={round(anova_tables[i].iloc[0, 3],2)}"
    ax.text(0.95, 0.05, anova_text, horizontalalignment='right', verticalalignment='bottom', transform=ax.transAxes, bbox=dict(facecolor='white', alpha=0.5))

    # save plot
    plt.savefig(f'{output_path}/{i}_plot_image.png')

def beta_plot(pcoa_results, metadata, output_path):
  for method in pcoa_results:
    pcoa_df = pcoa_results[method].samples
    pcoa_df['group'] = metadata

    # Apply the classic style
    plt.style.use('ggplot')

    plt.figure(figsize=(16, 8))
    ax = sns.scatterplot(data=pcoa_df, x='PC1', y='PC2', hue='group', s=100)
    ax.set_title("Beta diversity analysis")
    ax.set_ylabel(f"PCoA2({round(pcoa_results[method].proportion_explained[1]*100,2)}%)")
    ax.set_xlabel(f"PCoA1({round(pcoa_results[method].proportion_explained[0]*100,2)}%)")
    # plt.scatter(centroids[:, 0], centroids[:, 1], s=300, c='red')
    # Customize the legend to be outside the chart
    legend = ax.legend(title='Group', bbox_to_anchor=(1.05, 1), loc='upper left', fontsize='medium')
    legend.get_title().set_fontsize('medium')  # Set legend title font size

    result_permanova = permanova_test(pcoa_df)

    # Add ANOVA results as text in the plot
    anova_text = f"PERMANOVA:\ntest statistic={round(result_permanova['test statistic'],2)}, p={round(result_permanova['p-value'],2)}"
    ax.text(0.95, 0.05, anova_text, horizontalalignment='right', verticalalignment='bottom', transform=ax.transAxes, bbox=dict(facecolor='white', alpha=0.5))

    # save plot
    plt.savefig(f'{output_path}/beta_{method}_plot_image.png')

###### Detecting Pathogens - using Kraken2 report #######

def read_hamronization(path):
  with open(path, 'r') as file:
    lines = file.readlines()
  amr_contigs_collection = {}
  for line in lines:
    infor = line.strip().split('\t')
    amr_contigs_collection[infor[20]] = infor[1]
  return amr_contigs_collection

def read_table(path, index):
  with open(path, 'r') as file:
    lines = file.readlines()
  result = {}
  for line in lines:
    infor = line.strip().split('\t')
    if infor[0] in result:
      result[infor[index]].append(infor)
    else:
      result[infor[index]] = [infor]
  return result


def read_kraken_output(file_path):
    """
    Reads a custom Kraken-like report file and parses its content.

    Parameters:
        file_path (str): The path to the report file.

    Returns:
        list: A list of dictionaries where each dictionary represents a line in the report.
    """
    report_data = {}

    with open(file_path, 'r') as file:
        for line in file:
            line = line.strip()
            if not line:
                continue  # Skip empty lines

            parts = line.split('\t')
            if len(parts) < 5:
                continue  # Skip lines that don't have enough columns

            status = parts[0]
            read_id = parts[1]
            tax_id = int(parts[2])
            length = parts[3]
            extra_info = parts[4]
            if status == "C":
              report_data[read_id] = {
                  'status': status,
                  'tax_id': tax_id,
                  'length': length,
                  'extra_info': extra_info
              }

    return report_data

def read_kraken2_report(file_path, mode='species'):
    """
    Reads a Kraken2 report file and parses its content.

    Parameters:
        file_path (str): The path to the Kraken2 report file.

    Returns:
        list: A list of dictionaries where each dictionary represents a line in the report.
    """
    print('mode: ',mode)
    report_data = {}
    with open(file_path, 'r') as file:
        for line in file:
            line = line.strip()
            if not line:
                continue  # Skip empty lines

            parts = line.split('\t')
            if len(parts) < 6:
                continue  # Skip lines that don't have enough columns

            percentage = float(parts[0])
            num_reads = int(parts[1])
            num_reads_direct = int(parts[2])
            rank_code = parts[5]
            taxon_id = int(parts[6])
            taxon_name = parts[7].strip()
            # print('pass')
            if percentage >= 0.05:
              # print('pass1')
              if rank_code == 'S':
                # print('pass2')
                if mode == 'species':
                  report_data[taxon_name] = {
                      'percentage': percentage,
                      'num_reads': num_reads,
                      'num_reads_direct': num_reads_direct,
                      'rank_code': rank_code,
                      'taxon_id': taxon_id,
                      'taxon_name': taxon_name
                      }
                elif mode == 'strain':
                  report_data[taxon_id] = {
                      'percentage': percentage,
                      'num_reads': num_reads,
                      'num_reads_direct': num_reads_direct,
                      'rank_code': rank_code,
                      'taxon_id': taxon_id,
                      'taxon_name': taxon_name
                }
    return report_data

def retrive_tanomomy_from_amr_contigs(kraken2_path, kraken2_report_path, hamronization_path):
  hamronization_data = read_hamronization(hamronization_path)
  # print(len(hamronization_data))
  kraken2_data = read_kraken_output(kraken2_path)
  # print('kraken2:', len(kraken2_data))
  kraken2_report_data = read_kraken2_report(kraken2_report_path,  mode='strain')
  # print('kraken2 report: ', len(kraken2_report_data.keys()))
  result = {}
  not_match = []
  for protein_id in hamronization_data:
    contig_id = '_'.join(protein_id.split('_')[:2])
    if contig_id in kraken2_data:
      tax_id = kraken2_data[contig_id]['tax_id']
      # print('pass1')
      if tax_id in kraken2_report_data:
        # print('pass2')
        result[kraken2_report_data[tax_id]['taxon_id']] = [kraken2_report_data[tax_id]['percentage'], tax_id, kraken2_report_data[tax_id]['num_reads'], kraken2_report_data[tax_id]['taxon_name'], [protein_id, hamronization_data[protein_id]]]
    else:
      not_match.append(protein_id)
  return result, not_match

def read_table_pathogen_db(path, index):
  with open(path, 'r') as file:
    lines = file.readlines()
  result = {}
  header = True
  for line in lines:
    if not header:
      infor = line.strip().split(',')
      if len(infor[index]) != 0:
          result[infor[index]] = infor
    else:
      header = False
  return result
# display pathogen using plot
def pathogen_table(pathogen_detected, pathogen_taxo_id, output_dir):
  # create pathogen table
  get_pathogen_isolates = []
  for tax_id in pathogen_detected:
    get_pathogen_isolates.append(pathogen_taxo_id[tax_id])

  get_pathogen_infor_usage = {}
  list_columns = ["Genome Name", "NCBI Taxon ID", "Taxonomy Rank",'Assembly Accession',
                  'GenBank Accessions', "Genome Status", "Strain", "MLST", "Host Name"]
  for column in list_columns:
    get_pathogen_infor_usage[column] = []
  for i in get_pathogen_isolates:
    get_pathogen_infor_usage["Genome Name"].append(i[1]) # Access elements of i directly
    get_pathogen_infor_usage["NCBI Taxon ID"].append(i[3])
    get_pathogen_infor_usage["Taxonomy Rank"].append("|".join([i[6]] + i[8:14]))
    get_pathogen_infor_usage["Genome Status"].append(i[14])
    get_pathogen_infor_usage["Strain"].append(i[15])
    get_pathogen_infor_usage["MLST"].append(i[19])
    get_pathogen_infor_usage['Assembly Accession'].append(i[41])
    get_pathogen_infor_usage["GenBank Accessions"].append(i[43])
    get_pathogen_infor_usage["Host Name"].append(i[75])
  df = pd.DataFrame(get_pathogen_infor_usage)
  df.to_csv(f'{output_dir}/pathogen_detected_infor.csv', sep=',', index=False)

def display_pathogen_boxplot(all_sample_taxonomy, pathogen_detected, sample_names, output_dir):

  pathogen_names = []
  percentage_times = {}
  for sample in sample_names:
    result_taxonomy = all_sample_taxonomy[sample]
    for i in result_taxonomy:
      if str(result_taxonomy[i][1]) in pathogen_detected:
        pathogen_names.append(result_taxonomy[i][3])
        if result_taxonomy[i][3] not in percentage_times:
          percentage_times[result_taxonomy[i][3]] = [[sample, result_taxonomy[i][0]]]
        else:
          percentage_times[result_taxonomy[i][3]].append([sample, result_taxonomy[i][0]])
  set_pathogen_names = set(pathogen_names)

  # create boxplot
  pathogen_percentage_times = {'Name':[], 'Sample': [], 'Percentage in sample':[], 'Frequency':[]}
  for i in percentage_times:
    for j in range(len(percentage_times[i])):
      pathogen_percentage_times['Name'].append(i)
      pathogen_percentage_times['Sample'].append(percentage_times[i][j][0])
      pathogen_percentage_times['Percentage in sample'].append(float(percentage_times[i][j][1]))
      pathogen_percentage_times['Frequency'].append(len(percentage_times[i]))
  df_pathogen_percentage_times = pd.DataFrame(pathogen_percentage_times)
  df_pathogen_percentage_times.head(50)

  df = df_pathogen_percentage_times

  # Plot
  fig, ax = plt.subplots(figsize=(16, 8))
  bars = sns.boxplot(x='Name', y='Percentage in sample', hue='Name',data=df, ax=ax)

  plt.title('percentage and frequency in all sample', fontsize=16)
  plt.xlabel('Pathogen Name', fontsize=14)
  plt.ylabel('Percentage in Samples', fontsize=14)
  plt.xticks(rotation=-25, fontsize= 10)  # Rotate x-axis labels if needed

  # Setting custom x-tick intervals
  max_value = df['Percentage in sample'].max()
  interval = 0.02  # Interval for x-ticks
  ticks = np.arange(0, max_value + interval, interval)
  plt.yticks(ticks)  # Adjusting y-ticks for better readability
  plt.savefig(f'{output_dir}/pathogen_boxplot_image.png')

  # create heatmap
  # Create an empty DataFrame with species as index and samples as columns
  df = pd.DataFrame(0, index=percentage_times.keys(), columns=sample_names)

  # Fill the DataFrame with the provided values
  for species, values in percentage_times.items():
      for sample, abundance in values:
          df.at[species, sample] = abundance

  # Create the heatmap
  plt.figure(figsize=(16, 12))
  sns.heatmap(df, cmap="YlGnBu", annot=True, fmt=".2f", linewidths=.5)
  plt.title("Heatmap of Abundance Values")
  plt.ylabel("Species")
  plt.xlabel("Samples")
  plt.xticks(rotation=45, ha="right")
  plt.tight_layout()

  # Show the plot
  plt.savefig(f'{output_dir}/pathogen_heatmap_image.png')



def pathogen_detecting(nonpathogen, kraken_dir, hamronization_path , sample_names, pathogen_db_path, output_dir):
  all_sample_taxonomy = {}
  for sample_id in sample_names:
    result_taxonomy, not_match = retrive_tanomomy_from_amr_contigs(f'{kraken_dir}/{sample_id}.kraken', f'{kraken_dir}/{sample_id}.report', hamronization_path)
    all_sample_taxonomy[sample_id] = result_taxonomy

  pathogen_taxo_id = read_table_pathogen_db(pathogen_db_path,3)


  tax_ids = []
  for sample in sample_names:
    result_taxonomy = all_sample_taxonomy[sample]
    for i in result_taxonomy:
      tax_ids.append(str(result_taxonomy[i][1]))

  pathogen_detected = set(pathogen_taxo_id.keys()).intersection(set(tax_ids))
  print('before', pathogen_detected)
  pathogen_detected = set([i for i in pathogen_detected if i not in nonpathogen])
  print('after', pathogen_detected)
  # create pathogen table
  pathogen_table(pathogen_detected, pathogen_taxo_id, output_dir)
  # display using heatmap and boxplot
  display_pathogen_boxplot(all_sample_taxonomy, pathogen_detected, sample_names, output_dir)

def convert_kraken_report_to_feature_table(kraken_dir, output_dir):
    # Bash script to convert Kraken report files to feature table
    command = f'''
    for i in {kraken_dir}/*.report.txt; do
        awk -F "\\t" '{{print $1 "\\t"  $2 "\\t"  $3 "\\t"  $6 "\\t"  $7 "\\t"  $8}}' "$i" > "${{i%.report.txt}}_6line_report.txt";
    done;
    ls {kraken_dir}/*_6line_report.txt | xargs kraken-biom --fmt tsv -o {output_dir}/kraken_feature_table.tsv
    '''

    # Run the command with bash
    result = subprocess.run(command, shell=True, executable='/bin/bash', capture_output=True, text=True)

    # Check for errors
    if result.returncode != 0:
        print("Error:", result.stderr)
    else:
        print("Success:", result.stdout)



# ######@@ Abundace diferent ananlysis @@######

from rpy2.robjects import r, pandas2ri
from rpy2 import robjects
from rpy2.robjects.packages import importr
from collections import OrderedDict
import os  # Use os for creating directories
import pandas as pd


def run_ancom_BC(feature_table_path, braken_list,metadata, output_dir):

    # Define the folder names
    folders = ["family", "genus", "species"]

    # Loop through the folder names and create each folder
    for folder in folders:
        directory_path = os.path.join(output_dir, folder)
        try:
            os.makedirs(directory_path, exist_ok=True)  # Use os.makedirs for creating directories
            print(f"Directory creation successful for {folder}: {directory_path}")
        except Exception as e:
            print(f"Error creating {folder} directory:", str(e))

    # Load the R script containing the ANCOM-BC function
    r.source('/content/ancom_BC_test.r')

    # Access the R function defined in the script
    abuncance_analysis = r['abuncance_analysis']

    # Run the R function with arguments
    result = abuncance_analysis(feature_table_path, braken_list, metadata, output_dir)

    return result

# # Example usage of the function
# ancombc, enhancedvolcano, dplyr = install_and_load_r_packages()


# # Example usage of the function
# run_ancom_BC("/path/to/feature_table.csv", "/path/to/metadata.csv", "/output/directory")


###### Amr gene analysis #######



### AMR intersection deduplication

def create_table(genes):
  sample_amr = {}
  set_sample_amr = {}
  amr_collection = {}

  for amr in genes:
    sample_name = amr[0].split('.')[0]
    amrfinder_tools =  amr[0].split('.')[-1]
    amr[1] = amr[1].upper()
    if sample_name in sample_amr:
      if amrfinder_tools in sample_amr[sample_name]:
        sample_amr[sample_name][amrfinder_tools].append((amr[1],amr[-1]))
        amr_collection[sample_name][amrfinder_tools].append(amr[1])
      else:
        sample_amr[sample_name][amrfinder_tools] = [(amr[1],amr[-1])]
        amr_collection[sample_name][amrfinder_tools] = [amr[1]]
    else:
      sample_amr[sample_name] = {amrfinder_tools:[(amr[1],amr[-1])]}
      amr_collection[sample_name] = {amrfinder_tools:[amr[1]]}
  return sample_amr,  amr_collection

def retrive_intersection(result2):
  print(len(result2))
  compare_list = {}
  for sample in result2:
    collection = []
    for amr_finder in result2[sample]:
      collection += result2[sample][amr_finder]
    compare_list[sample] = set(collection)
  lis_amrs = list(compare_list.values())
  for i in range(len(lis_amrs)):
      lis_amrs[i] = set( i.upper() for i in lis_amrs[i])
      print(lis_amrs[i])
  intersection_result = lis_amrs[0] & lis_amrs[1] & lis_amrs[2] & lis_amrs[3] & lis_amrs[4] & lis_amrs[5] & lis_amrs[6] & lis_amrs[7] & lis_amrs[8] & lis_amrs[9] & lis_amrs[10] & lis_amrs[11] & lis_amrs[12] & lis_amrs[13] & lis_amrs[14] & lis_amrs[15] & lis_amrs[16] & lis_amrs[17] & lis_amrs[18] & lis_amrs[19] & lis_amrs[20]
  return intersection_result

def display_table(intersection_result_genes, genes):
  loi_format = {}
  for gene_name in intersection_result_genes:
    loi_format[gene_name] = {}

  for sample in genes:
    for amr_finder in genes[sample]:
      if amr_finder == 'deeparg':
        if amr_finder in genes[sample]:
            for element in genes[sample][amr_finder]:

              if element[0] in loi_format:
                if sample in loi_format[element[0]]:
                  loi_format[element[0]][sample] += 1
                else:
                  loi_format[element[0]][sample] = 1

  return loi_format

def draw_heat_map(pivot_table, output_dir):
    import pandas as pd
    import matplotlib.pyplot as plt
    import seaborn as sns

    # Set 'Genes' column as the index
    df.set_index('Genes', inplace=True)

    # Plot heatmap
    plt.figure(figsize=(24, 10))
    # sns.heatmap(df, annot=True, cmap="YlGnBu")
    sns.heatmap(df, annot=True, cmap="coolwarm")
    plt.title('Gene Duplication in Each Sample')
    plt.xlabel('Samples')
    plt.ylabel('Genes')
    plt.savefig(f'{output_dir}/heat_map.png')

def AMR_intersection_deduplication(hamronization_path, output_dir):
  genes_75_result1, genes_75_result2 = create_table(genes_75)
  intersection_result_genes_75 = retrive_intersection(genes_75_result2)
  filter_75 = df_75['gene_symbol'].isin(intersection_result_genes_75)
  df_75[filter_75].to_csv('/content/amr_final_annalysis/hamronization_combined_only_intersection_report_Nhu_75.tsv', sep='\t', index=False)
  result_75 = display_table(intersection_result_genes_75, genes_75_result1)

  # Convert the dictionary to a pandas DataFrame
  df_loi_format_75 = pd.DataFrame(result_75)
  # Transpose the DataFrame to have ARM genes as rows and samples as columns
  df_loi_format_75 = df_loi_format_75.transpose()
  # Fill NaN values with 0
  df_loi_format_75 = df_loi_format_75.fillna(0)
  df_loi_format_75['Genes'] = df_loi_format_75.index.to_list()
  # Display the DataFrame
  df_loi_format_75.to_csv('/content/amr_final_annalysis/loi_genes_75.tsv', sep='\t', index=False)
  draw_heat_map(df_loi_format_75, output_dir)

### AMR ontology

def get_accession(ontology_genes, keyword):
  for key in ontology_genes.keys():
    if keyword in key:
      return ontology_genes[key][0]

def read_table(path):
  if path.split('.')[-1] == 'tsv':
    sep='\t'
  else:
     sep=','
  with open(path,'r') as file:
    raw = file.readlines()
  result = {}
  for element in raw:
    info =  element.rstrip('\n').split(sep)
    result[info[4].lower()] = info
  return result

def read_json (file_path):
    with open(file_path, 'r') as file:
        data = json.load(file)
    return data

def parse_html_to_data (file_path):
    # Reading the HTML file content
    with open(file_path, 'r') as file:
        html_content = file.read()

    # Using BeautifulSoup to parse the HTML content
    soup = BeautifulSoup(html_content, 'html.parser')

    # Extracting the relevant information from the HTML content
    data = {
        "Accession": [],
        "CARD Short Name": [],
        "Definition": [],
        "AMR Gene Family": [],
        "Drug Class": [],
        "Resistance Mechanism": [],
        # "Classification": [],
        "Parent Term(s)": [],
        "Publications": []
    }

    # Extracting the data from the HTML content with error handling
    keys = data.keys()
    for key in keys:
        element = soup.find(string=key)
        if element:
            value = element.find_next('td').text.strip()
        else:
            value = 'N/A'  # Handle the case where the key is not found
        data[key].append(value)
    return data

def find_cvterm_id(ontology_json, accession):
  for i in ontology_json:
    if i['accession'] == accession:
      return i['cvterm_id']

def find_accession_using_gene_name(amr_genes, ontology_genes):
  Accession_soil_amr = []
  notmatch = []
  for i in amr_genes:
    # print(i)
    try:
      # print(ontology_genes[i.lower()])
      Accession_soil_amr.append(ontology_genes[i.lower()][0])
    except:
      accession = get_accession(ontology_genes, i.lower())
      # print(accession)
      if accession == None:
        notmatch.append(i)
        continue
      else:
        Accession_soil_amr.append(accession)
  return Accession_soil_amr, notmatch

def download_file(path, cvterm_id, output_dir):
  import requests
  # Define the URL of the file to download
  url = f'https://card.mcmaster.ca/ontology/{cvterm_id}'

  # Send a GET request to the URL
  response = requests.get(url)

  # Check if the request was successful
  if response.status_code == 200:
      # Open a file in write-binary mode
      with open(f'{ output_dir}/card_ontology_html/{cvterm_id}', 'wb') as file:
          # Write the content of the response to the file
          file.write(response.content)
      print("File downloaded successfully.")
      return f'{ output_dir}/card_ontology_html/{cvterm_id}'
  else:
      print(f"Failed to download file. Status code: {response.status_code}")
      return None

def retrieve_amr_ontology(amr_genes_path, ontology_genes_path_list, ontology_dir):

  import pandas as pd
  df = pd.read_csv(amr_genes_path, sep='\t')
  amr_genes = df['Genes'].to_list()
  amr_genes[3] = 'OTR(A)'
  amr_genes = amr_genes[:len(amr_genes) -1]
  amr_genes.append('LMRA')

  ontology_genes = read_table(ontology_genes_path_list[0])
  ontology_json = read_json(ontology_genes_path_list[1])

  accession_list, notmatch = find_accession_using_gene_name(amr_genes, ontology_genes)
  cvterm_ids = ['36713','36210'] # for some genes no match by name with database
  data_ontology = {
          "Accession": [],
          "CARD Short Name": [],
          "Definition": [],
          "AMR Gene Family": [],
          "Drug Class": [],
          "Resistance Mechanism": [],
          # "Classification": [],
          "Parent Term(s)": [],
          "Publications": []
      }

  for i, accession in enumerate(accession_list):
    if i in [1,2]:
      cvterm_id = cvterm_ids[i-1]
    else:
      cvterm_id = find_cvterm_id(ontology_json, accession)
    path = download_file(f'{ontology_dir}/card_ontology_html/{cvterm_id}.html', cvterm_id, ontology_dir)
    if path != None:
      data = parse_html_to_data(path)
      for field in data_ontology:
        if field == 'Accession':
          data_ontology[field].append(accession)
        else:
          data_ontology[field].append(', '.join(data[field]))

  df_ontology = pd.DataFrame(data_ontology)
  df_ontology.to_csv(f'{ontology_dir}/ontology_amr_genes.csv', index=False)

### Drug tpm abundace table

def retrieve_amr_ontology_for_drug(amr_genes_path, ontology_genes_path_list, ontology_dir):

  import pandas as pd

  df = pd.read_csv(amr_genes_path, sep='\t')
  amr_genes= []
  for gene in list(set(df['gene_symbol'].to_list())):
     if gene == 'MULTIDRUG_ABC_TRANSPORTER':
      amr_genes.append('LMRA')
     elif  gene == 'OTRA':
      amr_genes.append('OTR(A)')
     else:
      amr_genes.append(gene)

  ontology_genes = read_table(ontology_genes_path_list[0])
  ontology_json = read_json(ontology_genes_path_list[1])

  # retrieve acccestion list and filter it

  accession_list, notmatch = find_accession_using_gene_name(amr_genes, ontology_genes)
  convert_notmatch = notmatch
  uniprot_genes = []
  uniprot_notmatch = []
  for i in notmatch:
    try:
      uniprot_genes.append(uniport_search(i))

    except:
      uniprot_notmatch.append(i)
  accession_list2, notmatch = find_accession_using_gene_name(uniprot_genes, ontology_genes)
  # print(notmatch)
  accession_list += accession_list2
  notmatch += uniprot_notmatch

  card_gene_ids_from_notmatch = [('36306','tetC'),('36712', 'tet43'),('36329', 'tetO'), ('36325', 'tetM'),('39364', 'vanR-O'), ('40381', 'vanX gene in vanI cluster'),('36320', 'TETV'), ('40381',  'VANXI')]
  card_gene_from_notmatch = [i[1] for i in card_gene_ids_from_notmatch]
  card_id_from_notmatch = [i[0] for i in card_gene_ids_from_notmatch]
  notmatch = set(notmatch).difference(set(card_gene_from_notmatch))
  accession_list += card_id_from_notmatch
  convert_notmatch = list(set(convert_notmatch).difference(set(notmatch)))
  data_ontology = {
    "Accession": [],
    "CARD Short Name": [],
    "Definition": [],
    "AMR Gene Family": [],
    "Drug Class": [],
    "Resistance Mechanism": [],
    # "Classification": [],
    "Parent Term(s)": [],
    "Publications": []
}

  for i in accession_list:
    if 'ARO' in i:
      cvterm_id = find_cvterm_id(ontology_json, i)
    else:
      cvterm_id = i
    path = download_file(f'{ontology_dir}/card_ontology_html/{cvterm_id}.html', cvterm_id, ontology_dir)
    # print(path)
    if path != None:
      data = parse_html_to_data(path)
      for field in data_ontology:
        if field == 'Accession':
          data_ontology[field].append(i)
        else:
          data_ontology[field].append(', '.join(data[field]))
    else:
        print('not find cvterm_id:', cvterm_id)

  df_ontology = pd.DataFrame(data_ontology)
  df_ontology.to_csv(f'{ontology_dir}/ontology_gff75_amr_genes.tsv', index=False)
  return  df_ontology

def uniport_search(query):
  import requests
  print(query)
  base_url = "https://rest.uniprot.org/uniprotkb/search"
  params = {
      "query": query,  # Example UniProt ID
      "format": "json"
  }
  gene= []
  response = requests.get(base_url, params=params)

  if response.status_code == 200:
      data = response.json()
      if "results" in data:

              entry = data["results"][0]
              if entry :
                # print(entry['genes'][0]['geneName']['value'])
                # print(entry['proteinDescription']['submissionNames'][0]['fullName']['value'])
                return (entry['genes'][0]['geneName']['value'])
      else:
          print("No results found")
          return None
  else:
      print(f"Failed to retrieve data: {response.status_code}")
      return None

def create_drug_class_abundace_tpm_table(amr_genes_path, ontology_genes_path_list, ontology_dir, HAMRONIZATION_DIR, BUILD_FEATURE_TABLE_DIR):

  df_ontology = retrieve_amr_ontology_for_drug(amr_genes_path, ontology_genes_path_list, ontology_dir)
  Drug_dict = {}
  Resistance_Mechanism_dict = {}
  Drug_Class = df_ontology['Drug Class'].to_list()
  CARD_Short_Name =  df_ontology['CARD Short Name'].to_list()
  Resistance_Mechanism = df_ontology['Resistance Mechanism'].to_list()

  for i, list_dug in enumerate(Drug_Class):
    drugs = list_dug.split(', ')
    if Resistance_Mechanism[i] in Resistance_Mechanism_dict:
      Resistance_Mechanism_dict[Resistance_Mechanism[i]].append(CARD_Short_Name[i])
    else:
      Resistance_Mechanism_dict[Resistance_Mechanism[i]] = [CARD_Short_Name[i]]
    for drug in drugs:
      if drug in Drug_dict:
        Drug_dict[drug].append(CARD_Short_Name[i].upper())
      else:
        Drug_dict[drug] = [CARD_Short_Name[i].upper()]

  convert_notmatch_dict = {'vanR-O': 'vanR_in_vanO_cl', 'TET43':'tet(43)', 'RIFAMPIN_MONOOXYGENASE': 'iri',
                         'ADP-RIBOSYLATING_TRANSFERASE_ARR':'arr-1', 'MULTIDRUG_ABC_TRANSPORTER':'LMRA',
                         'vanX gene in vanI cluster': 'vanX_in_vanI_cl', 'CHLORAMPHENICOL_EXPORTER':'MdtK',
                         'TETA(48)':'tet(C)', 'TETO':'tet(O)', 'VANXI':'vanX_in_vanI_cl', 'OTRA':'otr(A)S.rim',
                         'TRANSCRIPTIONAL_REGULATORY_PROTEIN_CPXR_CPXR':'Paer_CpxR_MULT', 'Mycobacterium tuberculosis rpsL mutations conferring resistance to Streptomycin': 'Mtub_rpsL_STR',
                          'TETM':'tet(M)', 'TETV':'tet(V)'}

  df_counts = pd.read_csv(f"{BUILD_FEATURE_TABLE_DIR}/amr_genes_tpm_result.tsv", sep='\t')
  df_counts_dict =  df_counts.to_dict(orient='list')
  chr_counts_dict = {}
  for i in range(len(df_counts_dict['Sequence_id'])):
    chr_counts_dict[df_counts.iloc[i,:][0]]=[round(i,2)for i in df_counts.iloc[i,1:].to_list()]

  import numpy as np

  df_hamronization_combined = pd.read_csv(f"{HAMRONIZATION_DIR}/hamronization_combined_report_Nhu_75.tsv", sep='\t')
  amr_names = df_hamronization_combined ['gene_symbol'].to_list()
  input_sequence_id	 = df_hamronization_combined ['input_sequence_id'].to_list()
  Drug_sequence_id_dict = {}
  for i, amr_name in enumerate(amr_names):
    if amr_name in convert_notmatch_dict:
      amr_name = convert_notmatch_dict[amr_name]
    for item in Drug_dict.items():
      if amr_name in item[1]:
        if item[0] in Drug_sequence_id_dict:
          array1 = np.array(Drug_sequence_id_dict[item[0]])
          array2 = np.array(chr_counts_dict[input_sequence_id[i]])
          result = array1 + array2
          Drug_sequence_id_dict[item[0]] = result.tolist()
        else:
          Drug_sequence_id_dict[item[0]] = chr_counts_dict[input_sequence_id[i]]
  Drug_tpm_dict = {'Drug_names':[], 'SM1032307A01': [],'SM1032307A02':[],
 'SM1032307A03':[], 'SM1032308A01':[], 'SM1032308A02':[], 'SM1032308A03':[],
 'SM1032308A04':[], 'SM1032308A05':[], 'SM1032308A06':[],
 'SM1032308A07':[], 'SM1032308A08':[], 'SM1032308A09':[],
 'SM1032309A01':[], 'SM1032309A02':[], 'SM1032309A03':[],
 'SM1032309A04':[], 'SM1032309A05':[], 'SM1032309A06':[],
 'SM1032309A07':[], 'SM1032309A08':[], 'SM1032309A09':[]}

  for i in Drug_sequence_id_dict:
    Drug_tpm_dict['Drug_names'].append(i)
    for i, tpm_score in enumerate(Drug_sequence_id_dict[i]):
      Drug_tpm_dict[sample_names[i]].append(tpm_score)

  df_drug_tpm = pd.DataFrame(Drug_tpm_dict)
  df_drug_tpm.to_csv(f'{ontology_dir}/drug_tpm_result.tsv', sep='\t', index=False)


####### Build Phylogeny tree #######
# link itol: https://itol.embl.de/

def build_phylogeny_tree(tax_ids, output_dir):

  ncbi = NCBITaxa()
  # Get the scientific names for the specified taxa
  taxid2name = ncbi.get_taxid_translator(tax_ids)
  tree = ncbi.get_topology(tax_ids)
  # Annotate the tree nodes with scientific names
  for node in tree.traverse():
      if node.is_leaf():
          node.name = taxid2name[int(node.name)]
  newick_format = tree.write(format=1)
  with open(f'{output_dir}/phylogeny_tree.nwk', 'w') as file:
      file.write(newick_format)


###### Kraken2 display bars plot ######
# link github: https://github.com/acvill/bracken_plot?tab=readme-ov-file

def convert_kraken_report_to_bracken(KRAKEN_DIR, KRAKEN_DB, output_dir):
    script = f"""
    #!/bin/bash

    # Create output directory
    mkdir -p {output_dir}

    # Define levels
    levels="P,C,O,F,G,S,S1"

    # Loop over each Kraken2 report file in the specified directory
    for sample in {KRAKEN_DIR}/*kraken.report.txt; do
      # Extract sample ID (this assumes filenames have a pattern like SMxxx_...)
      sample_id=$(echo $sample | grep -oP 'SM.*(?=\\.)')
      echo "Processing sample: $sample_id"

      # Loop over each taxonomic level
      for level in $(echo $levels | sed "s/,/ /g"); do
        # Run Bracken
        bracken \\
          -d {KRAKEN_DB} \\
          -i $sample \\
          -o "{output_dir}/${{sample_id}}.bracken_${{level}}.txt" \\
          -r 75 \\
          -l $level
      done
    done

    # Create merged output directory
    mkdir -p {output_dir}/merged

    # Combine Bracken outputs for each level
    for level in $(echo $levels | sed "s/,/ /g"); do
      combine_bracken_outputs.py \\
        --files {output_dir}/*.kraken.report.bracken_${{level}}.txt \\
        --output {output_dir}/merged/merged_bracken_${{level}}.txt
      sed -i "s/.kraken.report.bracken_${{level}}.txt//g" {output_dir}/merged/merged_bracken_${{level}}.txt
    done
    """

    result = subprocess.run(script, shell=True, capture_output=True, text=True, executable='/bin/bash')

    if result.returncode != 0:
        print("Error:", result.stderr)
    else:
        print("Success:", result.stdout)
