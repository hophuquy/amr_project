#@title 2 Create feature table

# @markdown Please execute this cell by pressing the *Play* button on
# @markdown the left.

from IPython.utils import io
import pandas as pd
import os
import subprocess
import tqdm.notebook
from function_for_part_2 import *

TQDM_BAR_FORMAT = '{l_bar}{bar}| {n_fmt}/{total_fmt} [elapsed: {elapsed} remaining: {remaining}]'

#### set PATH for analysis
WORK_DIR = "/content/drive/MyDrive/output_funcscan/work_dir"
# DATASET = "/content/work_dir"
# FASTQ_PATH = '/content/drive/MyDrive/amr_final_fastq'
FUNCSCAN_DIR = '/content/drive/MyDrive/output_funcscan'
HAMRONIZATION_DIR = f'{WORK_DIR}/hamronization_summary'
AMR_CONTIGS_DIR = f'{WORK_DIR}/amr_cds_filter'
BUILD_FEATURE_TABLE_DIR = f'{WORK_DIR}/build_feature_table'
import os
os.environ['HAMRONIZATION_DIR'] = f'{WORK_DIR}/hamronization_summary'

add_compare = True

gff_files  = [f'{FUNCSCAN_DIR}/output_funcscan/annotation/prodigal/SM1032307A01/SM1032307A01.gff.gz',
           f'{FUNCSCAN_DIR}/output_funcscan/annotation/prodigal/SM1032307A02/SM1032307A02.gff.gz',
           f'{FUNCSCAN_DIR}/output_funcscan/annotation/prodigal/SM1032307A03/SM1032307A03.gff.gz',
           f'{FUNCSCAN_DIR}/output_funcscan/annotation/prodigal/SM1032308A01/SM1032308A01.gff.gz',
           f'{FUNCSCAN_DIR}/output_funcscan/annotation/prodigal/SM1032308A02/SM1032308A02.gff.gz',
           f'{FUNCSCAN_DIR}/output_funcscan/annotation/prodigal/SM1032308A03/SM1032308A03.gff.gz',
           f'{FUNCSCAN_DIR}/output_funcscan/annotation/prodigal/SM1032308A04/SM1032308A04.gff.gz',
           f'{FUNCSCAN_DIR}/output_funcscan/annotation/prodigal/SM1032308A05/SM1032308A05.gff.gz',
           f'{FUNCSCAN_DIR}/output_funcscan/annotation/prodigal/SM1032308A06/SM1032308A06.gff.gz',
           f'{FUNCSCAN_DIR}/output_funcscan/annotation/prodigal/SM1032308A07/SM1032308A07.gff.gz',
           f'{FUNCSCAN_DIR}/output_funcscan/annotation/prodigal/SM1032308A08/SM1032308A08.gff.gz',
           f'{FUNCSCAN_DIR}/output_funcscan/annotation/prodigal/SM1032308A09/SM1032308A09.gff.gz',
           f'{FUNCSCAN_DIR}/output_funcscan/annotation/prodigal/SM1032309A01/SM1032309A01.gff.gz',
           f'{FUNCSCAN_DIR}/output_funcscan/annotation/prodigal/SM1032309A02/SM1032309A02.gff.gz',
           f'{FUNCSCAN_DIR}/output_funcscan/annotation/prodigal/SM1032309A03/SM1032309A03.gff.gz',
           f'{FUNCSCAN_DIR}/output_funcscan/annotation/prodigal/SM1032309A04/SM1032309A04.gff.gz',
           f'{FUNCSCAN_DIR}/output_funcscan/annotation/prodigal/SM1032309A05/SM1032309A05.gff.gz',
           f'{FUNCSCAN_DIR}/output_funcscan/annotation/prodigal/SM1032309A06/SM1032309A06.gff.gz',
           f'{FUNCSCAN_DIR}/output_funcscan/annotation/prodigal/SM1032309A07/SM1032309A07.gff.gz',
           f'{FUNCSCAN_DIR}/output_funcscan/annotation/prodigal/SM1032309A08/SM1032309A08.gff.gz',
           f'{FUNCSCAN_DIR}/output_funcscan/annotation/prodigal/SM1032309A09/SM1032309A09.gff.gz',
           f'{FUNCSCAN_DIR}/compare_soil_sample/ERR436057.gff.gz',
           f'{FUNCSCAN_DIR}/compare_soil_sample/SRR22278226.gff.gz',
           f'{FUNCSCAN_DIR}/compare_soil_sample/SRR5207566.gff.gz',
           f'{FUNCSCAN_DIR}/compare_soil_sample/SRR5829540.gff.gz',
           ]

protein_files = [file.replace('.gff.gz', '.faa.gz') for file in gff_files]
cds_files = [file.replace('.gff.gz', '.fna.gz') for file in gff_files]

# set sample names
sample_names = []
for i, path  in enumerate(gff_files):
 if i <= len(gff_files)-5:
   sample_name = path.replace(f'{FUNCSCAN_DIR}/output_funcscan/annotation/prodigal/', '').replace('.gff.gz', '').split('/')[0]
 else:
   sample_name = path.replace(f'{FUNCSCAN_DIR}/compare_soil_sample/', '').replace('.gff.gz', '')
 sample_names.append(sample_name)

try:
 with tqdm.notebook.tqdm(total=100, bar_format=TQDM_BAR_FORMAT) as pbar:
   with io.capture_output() as captured:
     pass
     #### Filter and convert result of funscan to reference
     ### create folder for analysis
     %shell mkdir -p $WORK_DIR

     ### hamronization filter
     %shell mkdir $HAMRONIZATION_DIR

     if add_compare:
       %shell cat /content/drive/MyDrive/output_funcscan/output_funcscan/reports/hamronization_summarize/hamronization_combined_report.tsv > "$HAMRONIZATION_DIR/final_hamronization_combined_report.tsv"
       %shell awk -F "\t" 'NR!=1{print $0}' /content/drive/MyDrive/output_funcscan/compare_soil_sample/hamronization_combined_report.tsv >> "$HAMRONIZATION_DIR/final_hamronization_combined_report.tsv"

     ! awk -F "\t" '{if($36 >= 75){print $0}}' "$HAMRONIZATION_DIR/final_hamronization_combined_report.tsv"   > "$HAMRONIZATION_DIR/hamronization_combined_report_Nhu_75.tsv"
     ! awk -F "\t" '{if($36 >= 90){print $0}}' "$HAMRONIZATION_DIR/final_hamronization_combined_report.tsv"   > "$HAMRONIZATION_DIR/hamronization_combined_report_Nhu_90.tsv"
     ### retrieve contigs from hamronization filter file
     !mkdir -p $AMR_CONTIGS_DIR

     identity_75 = retrieve_taxid_from_hamronization_file(f"{HAMRONIZATION_DIR}/hamronization_combined_report_Nhu_75.tsv")
     identity_90 = retrieve_taxid_from_hamronization_file(f"{HAMRONIZATION_DIR}/hamronization_combined_report_Nhu_90.tsv")
     identity_75.pop('input_file_name')
     identity_90.pop('input_file_name')

     identities = [identity_75, identity_90]
     identity_precentage = ['identity_75', 'identity_90']

     for i, precentage in enumerate(identities):
       for j, sample_identity in enumerate(precentage.values()):
         retrive_amr_contigs(sample_identity,cds_files[j].replace('.gz', ''), f'{AMR_CONTIGS_DIR}/amr_{list(precentage.keys())[j]}_{identity_precentage[i]}.fna')

     list_amr_contigs_75 = [f'{AMR_CONTIGS_DIR}/amr_{i}_identity_75.fna' for i in sample_names]
     list_amr_contigs_90 = [f'{AMR_CONTIGS_DIR}/amr_{i}_identity_90.fna' for i in sample_names]

     combine_amr_contigs(list_amr_contigs_75, f'{AMR_CONTIGS_DIR}/combine_amr_contigs_75.fna')
     combine_amr_contigs(list_amr_contigs_90, f'{AMR_CONTIGS_DIR}/combine_amr_contigs_90.fna')
     pbar.update(25)

     ### Mapping and create feature table
     %shell mkdir -p  $BUILD_FEATURE_TABLE_DIR

     # create gff file

     %shell prodigal -f gff -i "$AMR_CONTIGS_DIR/combine_amr_contigs_75.fna" -o "$BUILD_FEATURE_TABLE_DIR/combine_amr_contigs_75.gff" -p meta
     # build index

     %shell bowtie2-build -q "$AMR_CONTIGS_DIR/combine_amr_contigs_75.fna" "$BUILD_FEATURE_TABLE_DIR/combine_amr_contigs_75"
     # mapping using bowties2

     result1 = subprocess.run(f'''
     bash -c '
     for i in /content/drive/MyDrive/amr_final_fastq/mapped/*1.fq; do
         id=$(echo $i | grep -oP "(SM|SRR).*?(?=_)");
         s2=${{i/1.fq/2.fq}};
         echo "Processing $i and $s2 with ID $id";
         bowtie2 -x {BUILD_FEATURE_TABLE_DIR}/combine_amr_contigs_75 -1 $i -2 $s2 -S {BUILD_FEATURE_TABLE_DIR}/$id.sam -p 8;
     done
     '
     ''', shell=True, capture_output=True, text=True)

     if result1.returncode != 0:
         print("Error:", result1.stderr)
     else:
         print("Success:", result1.stdout)

     retrieve only mapped read

     result2 = subprocess.run(f'''
     for id in $(ls "{BUILD_FEATURE_TABLE_DIR}"/*.sam | cut -d'.' -f1); do
         samtools view -bS -F 4 ${{id}}.sam > ${{id}}.bam;
         echo 'pass1';
         samtools sort ${{id}}.bam -@ 8 -o ${{id}}_sorted.bam;
         echo 'pass2';
     done
     ''', shell=True, capture_output=True, text=True)

     if result2.returncode != 0:
         print("Error:", result2.stderr)
     else:
         print("Success:", result2.stdout)

     # pbar.update(25)

     ## Create and Normalize feature table using using tpm normalization.

     result3 = subprocess.run(f'''
     # create feature count table
     featureCounts -F GFF -T 4 -p --countReadPairs -t 'CDS' -g 'ID' -a {BUILD_FEATURE_TABLE_DIR}/combine_amr_contigs_75.gff -o {BUILD_FEATURE_TABLE_DIR}/counts.txt {BUILD_FEATURE_TABLE_DIR}/*_sorted.bam

     # remove header
     awk -F "\t" 'NR!=1{{print $0}}' {BUILD_FEATURE_TABLE_DIR}/counts.txt > {BUILD_FEATURE_TABLE_DIR}/counts_no_header.tsv
     ''', shell=True, capture_output=True, text=True)

     if result3.returncode != 0:
         print("Error:", result3.stderr)
     else:
         print("Success:", result3.stdout)

     # rewrite columns
     columns = ['Geneid', 'Chr',   'Start', 'End', 'Strand',   'Length'] + sample_names # Remove the extra brackets around sample_names
     pd.read_csv(f"{BUILD_FEATURE_TABLE_DIR}/counts_no_header.tsv", sep='\t', header=None, names = columns).iloc[1:, :].to_csv(f"{BUILD_FEATURE_TABLE_DIR}/counts_no_header.tsv", sep='\t', index=False)
     tpm normalization
     calculate_tpm(f"{BUILD_FEATURE_TABLE_DIR}/counts_no_header.tsv", f"{BUILD_FEATURE_TABLE_DIR}/amr_genes_tpm_result.tsv")
     pbar.update(50)
except subprocess.CalledProcessError:
 print(captured)
 raise
