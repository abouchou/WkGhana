#!/bin/bash
input_file="../data/rvf_africa.tsv"
output_file="../data/newcomplete.tsv"
awk -F'\t' '$2 == "Complete" {print $5}' "$input_file" > "$output_file"
echo "Filtered file generated: $output_file"