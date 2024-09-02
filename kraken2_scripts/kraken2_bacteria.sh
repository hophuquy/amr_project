#!/bin/bash

# Define paths and parameters
KRAKEN2_DB="/path/to/bacteria_db"  # Path to the Kraken2 database for bacterial sequences
INPUT_DIR="/path/to/contigs"       # Directory containing input contigs (.fa.gz files)
OUTPUT_DIR="/path/to/kraken2"      # Directory where Kraken2 output will be stored
SAMPLE_LIST="$INPUT_DIR/sample_list.txt"  # File to store the list of sample names (now located in INPUT_DIR)

# Remove the sample list file if it already exists to avoid appending to an old list
rm -f "$SAMPLE_LIST"

# Extract unique base names from .fa.gz files in INPUT_DIR and write them to the sample list
for fa_gz_file in "$INPUT_DIR"/*.fa.gz; do
    if [[ -f "$fa_gz_file" ]]; then
        # Extract the base name by removing the .contigs.fa.gz suffix
        base_name=$(basename "$fa_gz_file" | sed -e 's/\.contigs\.fa\.gz$//')
        echo "$base_name" >> "$SAMPLE_LIST"  # Append the base name to the sample list
    fi
done

# Remove duplicate entries and sort the sample list alphabetically
sort -u "$SAMPLE_LIST" -o "$SAMPLE_LIST"
echo "Sample list created: $SAMPLE_LIST"

# Create output directory for bacteria-related files if it does not exist
mkdir -p "$OUTPUT_DIR/bacteria"

# Process each sample listed in SAMPLE_LIST
while IFS= read -r SAMPLE; do
    # Construct the full path for the corresponding .fa.gz file for the sample
    fa_file="$INPUT_DIR/${SAMPLE}.contigs.fa.gz"

    if [[ -f "$fa_file" ]]; then
        # Define the paths for classified output, Kraken2 output, and report files
        CLASSIFIED_FILE="$OUTPUT_DIR/bacteria/${SAMPLE}.classified-bacteria.fa"
        OUTPUT_FILE="$OUTPUT_DIR/bacteria/${SAMPLE}.bacteria.kraken2"
        REPORT_FILE="$OUTPUT_DIR/bacteria/${SAMPLE}.bacteria.kraken2.report"

        echo "Running Kraken2 on $fa_file ..."
        
        # Run Kraken2 with specified options to classify the sequences in each sample
        kraken2 --db "$KRAKEN2_DB" \  # Specify the path to the Kraken2 database to be used for classification
                --classified-out "$CLASSIFIED_FILE" \  # Output the sequences that are classified by Kraken2 to this file
                --output "$OUTPUT_FILE" \  # Save the detailed classification results (sequence ID, taxonomic ID, etc.) to this file
                --confidence 0.001 \  # Set the confidence threshold for classification; Kraken2 requires at least 0.1% confidence to make a classification
                --report "$REPORT_FILE" \  # Generate a report summarizing the classification results at various taxonomic levels (e.g., species, genus)
                --memory-mapping \  # Enable memory mapping of the Kraken2 database to reduce memory usage by accessing parts of the database from disk as needed
                --gzip-compressed "$fa_file"  # Indicate that the input file is gzip-compressed; Kraken2 will decompress it on the fly during processing

    else
        # If the .fa.gz file does not exist, skip processing for that sample
        echo "File $fa_file does not exist. Skipping..."
    fi
done < "$SAMPLE_LIST"

# Step 3: Process the Kraken2 report using kreport2mpa.py for each sample
cd "$OUTPUT_DIR/bacteria" || { echo "Cannot change directory to $OUTPUT_DIR/bacteria"; exit 1; }
while IFS= read -r SAMPLE; do
    echo "Processing kreport2mpa for $SAMPLE ..."
    # Run kreport2mpa.py to generate an MPA-style report from the Kraken2 report
    kreport2mpa.py -r "${SAMPLE}.bacteria.kraken2.report" --display-header -o "${SAMPLE}.bacteria.kraken2.mpa"
    
    # Extract counts from the MPA report, format them, and save to a count file
    tail -n+2 "${SAMPLE}.bacteria.kraken2.mpa" | cut -f 2 | sed "1 s/^/${SAMPLE} /" > "${SAMPLE}.bacteria.kraken2.count"
done < "$SAMPLE_LIST"
