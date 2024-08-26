import pandas as pd

def retrieve_taxid_from_hamronization_file(file_path):
    try:
        with open(file_path, 'r') as file:
            lines = file.readlines()
    except FileNotFoundError:
        print(f"Error: The file at {file_path} was not found.")
        return {}
    except Exception as e:
        print(f"Error: An error occurred while reading the file: {e}")
        return {}

    tax_id = {}

    for line in lines:
        try:
            fields = line.strip().split('\t')

            if len(fields) < 21:
                print(f"Warning: Line does not have enough fields: {line.strip()}")
                continue

            sample_id = fields[0].split('.')[0]
            taxonomic_id = fields[20]

            if sample_id not in tax_id:
                tax_id[sample_id] = [taxonomic_id]
            else:
                tax_id[sample_id].append(taxonomic_id)

        except Exception as e:
            print(f"Error: An error occurred while processing the line: {line.strip()}")
            print(f"Exception: {e}")

    return tax_id

def retrive_amr_contigs (contigs_id, path_contigs, path_amr):
  from Bio import SeqIO
  # Read the FASTA file
  seqs = {}
  for record in SeqIO.parse(path_contigs, "fasta"):
      if record.id in contigs_id:
        seqs[record.id] = [record.description, record.seq]

  with open(path_amr, 'w') as file:
    for contig in seqs.values():
      for i, part in enumerate(contig):
        if i == 0:
          file.write(">" + part.split(' ')[0] + "\n")
        else:
          file.write(str(part) + "\n")

def combine_amr_contigs( path_contigs, path_conbine_amr):
  from Bio import SeqIO
  # Read the FASTA file
  seqs = {}
  for path_contig in path_contigs:
    for record in SeqIO.parse(path_contig, "fasta"):
        if record.id not in seqs:
          seqs[record.id] = [record.description, record.seq]

  with open(path_conbine_amr, 'w') as file:
    for contig in seqs.values():
      for i, part in enumerate(contig):
        if i == 0:
          file.write(">" + part.split(' ')[0] + "\n")
        else:
          file.write(str(part) + "\n")

def calculate_tpm(featurecounts_file, output_file):

    """
    Calculate TPM (Transcripts Per Million) from a featureCounts result file.

    Parameters:
    - featurecounts_file: Path to the input featureCounts file (TSV format).
    - output_file: Path to save the TPM result file (TSV format).

    The function assumes the first 6 columns are metadata and subsequent columns are read counts for each sample.
    """

    # Load featureCounts result
    df = pd.read_csv(featurecounts_file, sep='\t', on_bad_lines='skip')

    # Extract gene lengths and counts
    gene_lengths = df['Length']
    counts = df.iloc[:, 6:]  # Read counts start from the 7th column

    # Calculate RPK (Reads Per Kilobase)
    rpk = counts.div(gene_lengths, axis=0) * 1e3

    # Calculate per-sample scaling factor
    scaling_factors = rpk.sum(axis=0) / 1e6

    # Calculate TPM (Transcripts Per Million)
    tpm = rpk.div(scaling_factors, axis=1)

    # Add gene IDs back to the TPM DataFrame
    tpm.insert(0, 'Sequence_id', df['Chr'])

    # Save the TPM result to a new file
    tpm.to_csv(output_file, sep='\t', index=False)

    # print(f"TPM calculation completed. Results saved to '{output_file}'")
