# @title 3.1 Down stream analysis
# @markdown Please execute this cell by pressing the *Play* button on
# @markdown the left.

### Steps:
1. Data Preprocessing
2. Quality Control
3. Statistical Analysis
4. Visualization


try:
  with tqdm.notebook.tqdm(total=100, bar_format=TQDM_BAR_FORMAT) as pbar:
    with io.capture_output() as captured:

      #### variable site #####

      KRAKEN_DIR = f'{WORK_DIR}/kraken_visualize'
      ARM_ONTOLOGY_DIR = f'{WORK_DIR}/amr_ontology'
      DIVERSITY_DIR = f'{WORK_DIR}/diverisity'
      PATHOGEN_DIR = f'{WORK_DIR}/pathogen_detected'
      metadata_sample = ['Time of sampling phase 1']*3 + ['Time of sampling phase 2']*9 +  ['Time of sampling phase 3']*9
      metadata_fertilized_sample = ['Non-fertilize phase 1']*3 + ['VC phase 2']*3 + ['HC phase 2']*3+ ['HCVS phase 2']*3 + ['VC phase 3']*3 + ['HC phase 3']*3+ ['HCVS phase 3']*3
      pbar.update(10)

      #### import site #####
      from function_for_part_tes_3 import *

      #### main run #####

      ## alpha beta diversity analysis

      !mkdir -p $DIVERSITY_DIR
      alpha_div, beta_div, pcoa_results = diversity_calculation('/content/drive/MyDrive/amr_final_annalysis/amr_genes_tpm_result.tsv')
      alpha_plot(alpha_div, metadata_fertilized_sample, DIVERSITY_DIR)
      beta_plot( pcoa_results, metadata_sample, DIVERSITY_DIR)
      pbar.update(10)
      ### abundance analysis using kraken2 report (ancom bc)

      ### visualize using krona plot and braken plot bar
      !mkdir $KRAKEN_DIR

      # Example usage:
      convert_kraken_report_to_bracken(
          "/content/drive/MyDrive/output_funcscan/kraken_output",
          "/content/drive/MyDrive/kraken2_db",
          "/content/drive/MyDrive/output_funcscan/bracken_output"
      )

      conver kraken2 report to feature table using kraken-biom (https://github.com/smdabdoub/kraken-biom)
      convert_kraken_report_to_feature_table('/content/drive/MyDrive/output_funcscan/kraken_output', KRAKEN_DIR)
      ! awk -F "\t" 'NR!=1{OFS="\t"; if (NR==2){gsub(/#OTU ID/, "tax_id", $1); gsub(/.kraken_6line_report/, "")} ; print}' /content/drive/MyDrive/output_funcscan/work_dir/kraken_visualize/kraken_feature_table.tsv > /content/drive/MyDrive/output_funcscan/work_dir/kraken_visualize/kraken_feature_no_header_table.tsv
      pbar.update(10)

      ### detect duplicate  amr gennes in all sample identity 75%, 90%

      ## detect pathogen and their amr genes
      Detecting Pathogens - using Kraken2 report
      !mkdir -p $PATHOGEN_DIR
      nonpathogen = ['1404', '2775496', '459858']
      pathogen_detecting(nonpathogen, '/content/drive/MyDrive/output_funcscan/kraken_output', '/content/drive/MyDrive/amr_final_annalysis/hamronization_combined_report_Nhu_75.tsv', sample_names, '/content/drive/MyDrive/BVBRC_genome_human.csv', PATHOGEN_DIR)
      pbar.update(10)
      ## get amr ontology using card database
      # download card database
      result = subprocess.run(f'''
            bash -c '
            mkdir -p {ARM_ONTOLOGY_DIR}
            mkdir -p {ARM_ONTOLOGY_DIR}/card_ontology
            wget -O {ARM_ONTOLOGY_DIR}/card_ontology/ontology-v3.2.9.tar.bz2 https://card.mcmaster.ca/download/5/ontology-v3.2.9.tar.bz2
            tar -xvjf {ARM_ONTOLOGY_DIR}/card_ontology/ontology-v3.2.9.tar.bz2 -C {ARM_ONTOLOGY_DIR}/card_ontology
            mkdir -p {ARM_ONTOLOGY_DIR}/card_ontology_html
            '
            ''', shell=True, capture_output=True, text=True)

      if result.returncode != 0:
          print("Error:", result.stderr)
      else:
          print("Success:", result.stdout)
      amr_intersection_path = '/content/drive/MyDrive/amr_final_annalysis/loi_genes_75.tsv'
      ontology_genes_path_list = [f'{ARM_ONTOLOGY_DIR}/card_ontology/aro.tsv', f'{ARM_ONTOLOGY_DIR}/card_ontology/aro.json']
      retrieve_amr_ontology(amr_intersection_path, ontology_genes_path_list, ARM_ONTOLOGY_DIR)
      pbar.update(10)

      ## get amr drug tpm abundace table
      amr_genes_path = '/content/drive/MyDrive/amr_final_annalysis/hamronization_combined_report_Nhu_75.tsv'
      ontology_genes_path_list = [f'{ARM_ONTOLOGY_DIR}/card_ontology/aro.tsv', f'{ARM_ONTOLOGY_DIR}/card_ontology/aro.json']
      ontology_dir = ARM_ONTOLOGY_DIR
      HAMRONIZATION_DIR = '/content/drive/MyDrive/amr_final_annalysis'
      BUILD_FEATURE_TABLE_DIR = '/content/drive/MyDrive/amr_final_annalysis'
      create_drug_class_abundace_tpm_table(amr_genes_path, ontology_genes_path_list, ontology_dir, HAMRONIZATION_DIR, BUILD_FEATURE_TABLE_DIR)
      pbar.update(10)

      ## build phylogeny tree of pathogen detected
      pathogen_tax_id  = pd.read_csv(f'{PATHOGEN_DIR}/pathogen_detected_infor.csv')['NCBI Taxon ID'].to_list()
      # after build file newick we use itol tool to visualize
      build_phylogeny_tree(pathogen_tax_id, PATHOGEN_DIR)
      pbar.update(40)
except subprocess.CalledProcessError:
  print(captured)
  raise
