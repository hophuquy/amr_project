#!/bin/bash

# Bacteria
## Download the Kraken2 "pluspf" (plus, protein, fungi) database archive for bacteria
wget https://genome-idx.s3.amazonaws.com/kraken/k2_pluspf_20240605.tar.gz

## Extract the contents of the bacteria database archive to the directory named "bacteria_db"
tar -xvf k2_pluspf_20240605.tar.gz -C bacteria_db

# Virus
## Download the Kraken2 viral database archive from the web
wget https://genome-idx.s3.amazonaws.com/kraken/k2_viral_20240605.tar.gz

## Extract the contents of the viral database archive to the directory named "virus_db"
tar -xvf k2_viral_20240605.tar.gz -C virus_db

# Fungi
# Download the NCBI taxonomy data required for the Kraken2 database construction 
# and store it in the "fungi_db" directory
kraken2-build --download-taxonomy --db fungi_db

# Download the fungal library sequences from NCBI and add them to the "fungi_db" directory
kraken2-build --download-library fungi --db fungi_db

# Build the Kraken2 database for fungi using the downloaded taxonomy and library data
kraken2-build --build --db fungi_db
