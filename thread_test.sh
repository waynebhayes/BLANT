#!/bin/bash

# Output file
OUTPUT_FILE="output_threads8.txt"

# Clear the output file before writing
> "$OUTPUT_FILE"

# List of k values
K_VALUES=(3 4 5 6 7 8)

# Loop through each k value and execute the command
for k in "${K_VALUES[@]}"; do
    CMD="./blant -n10000000 -k$k -s NBE networks/syeast.el"
    echo "Running: $CMD" | tee -a "$OUTPUT_FILE"
    eval "$CMD" >> "$OUTPUT_FILE" 2>&1
    echo -e "\n---\n" >> "$OUTPUT_FILE"
done


# manually change _MAX_THREADS in blant.c to 1 or 8, make it, then manually change the output file, then  run this script