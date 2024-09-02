#!/bin/bash

# Define paths and parameters
KRAKEN2_DB="/path/to/fungi_db"
INPUT_DIR="/path/to/contigs"
OUTPUT_DIR="/path/to/fungi"
SAMPLE_LIST="$INPUT_DIR/sample_list.txt"  # Đã thay đổi để nằm trong INPUT_DIR

# Remove the sample list file if it already exists
rm -f "$SAMPLE_LIST"

# Extract unique base names from .fa.gz files in INPUT_DIR and write to the sample list
for fa_gz_file in "$INPUT_DIR"/*.fa.gz; do
    if [[ -f "$fa_gz_file" ]]; then
        # Extract base name (remove the .contigs.fa.gz part)
        base_name=$(basename "$fa_gz_file" | sed -e 's/\.contigs\.fa\.gz$//')
        echo "$base_name" >> "$SAMPLE_LIST"
    fi
done

# Remove duplicate entries and sort the sample list
sort -u "$SAMPLE_LIST" -o "$SAMPLE_LIST"
echo "Sample list created: $SAMPLE_LIST"

# Create output directory if it does not exist
mkdir -p "$OUTPUT_DIR/kraken2_fungi"

# Process each sample listed in SAMPLE_LIST
while IFS= read -r SAMPLE; do
    fa_file="$INPUT_DIR/${SAMPLE}.contigs.fa.gz"  # Construct the full path for the .fa.gz file

    if [[ -f "$fa_file" ]]; then
        # Define output and report file paths
        CLASSIFIED_FILE="$OUTPUT_DIR/kraken2_fungi/${SAMPLE}.classified-fungi.fa"
        OUTPUT_FILE="$OUTPUT_DIR/kraken2_fungi/${SAMPLE}.fungi.kraken2"
        REPORT_FILE="$OUTPUT_DIR/kraken2_fungi/${SAMPLE}.fungi.kraken2.report"

        echo "Running Kraken2 on $fa_file ..."
        
        # Run Kraken2 with specified options
        kraken2 --db "$KRAKEN2_DB" \
                --classified-out "$CLASSIFIED_FILE" \
                --output "$OUTPUT_FILE" \
                --confidence 0.001 \
                --report "$REPORT_FILE" \
                --memory-mapping \
                --gzip-compressed "$fa_file"
    else
        echo "File $fa_file does not exist. Skipping..."
    fi
done < "$SAMPLE_LIST"

# Step 3: Process kreport2mpa.py for each sample
cd "$OUTPUT_DIR/kraken2_fungi" || { echo "Cannot change directory to $OUTPUT_DIR/kraken2_fungi"; exit 1; }
while IFS= read -r SAMPLE; do
    echo "Processing kreport2mpa for $SAMPLE ..."
    /home/lucianhu/AMR/fungi/kreport2mpa.py -r "${SAMPLE}.fungi.kraken2.report" --display-header -o "${SAMPLE}.fungi.kraken2.mpa"
    tail -n+2 "${SAMPLE}.fungi.kraken2.mpa" | cut -f 2 | sed "1 s/^/${SAMPLE} /" > "${SAMPLE}.fungi.kraken2.count"
done < "$SAMPLE_LIST"
