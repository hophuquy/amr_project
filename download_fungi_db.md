# Ensembl Fungi Release 60, accessed on November 9, 2024 : 1505 fasta
wget -r -np -nH --cut-dirs=5 -A "*.dna.toplevel.fa.gz" ftp://ftp.ensemblgenomes.org/pub/fungi/release-60/fasta/

# FungiDB Release 68, Genome + Popset, Genetic Variation + Sequence: 424 fasta
https://fungidb.org/fungidb/app/downloads

# NCBI, accessed on November 9, 2024 
micromamba install -c bioconda ncbi-datasets-cli
datasets download genome taxon fungi --include genome --filename fungi_genomes_dna.zip
