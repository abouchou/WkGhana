---
title: "Genetic Analysis of Rift Valley Fever Virus (RVFV) in Africa"
author: "Your Name"
date: "`r Sys.Date()`"
output:
  html_document:
    toc: true
    toc_depth: 2
    theme: united
    highlight: tango
---

```{r setup, include=FALSE}
knitr::opts_chunk$set(echo = TRUE, warning = FALSE, message = FALSE)
```

## Introduction

This report presents a genomic analysis of the Rift Valley Fever Virus (RVFV) using data from African countries. It includes steps for metadata preparation (Bash), GenBank data retrieval (Python/Biopython), data cleaning and analysis (R/tidyverse), and result visualization. The report is reproducible and includes the versions of the tools used.

## Tool Versions

Below are the simulated versions of the tools used in this project (adjust according to your environment):

- **Bash**: 5.1.16
- **Python**: 3.9.5
- **Biopython**: 1.79
- **R**: 4.3.0
- **tidyverse**: 2.0.0
- **lubridate**: 1.9.2
- **ggplot2**: 3.4.2


## Task 1: Environment Setup

### Commands and Setup

We set up a Conda environment to manage the required tools for this project, ensuring reproducibility across different systems. The following steps were performed in the terminal:

#### Step 1: Create the Conda Environment

A Conda environment named `quarto-bio` was created with the necessary tools (Bash, Python, Biopython, R, tidyverse, etc.):

```bash
# Create a new Conda environment with Python and R
conda create -n quarto-bio python=3.9 r-base=4.3 -y

# Activate the environment
conda activate quarto-bio

# Install required Python packages
pip install biopython==1.79 pandas

# Install required R packages
R -e "install.packages(c('tidyverse', 'lubridate', 'ggplot2'), repos='https://cloud.r-project.org')"
```

#### Step 2: Export the Environment to a `.yml` File

The Conda environment was exported to a `.yml` file to ensure reproducibility:

```bash
# Export the environment to a YAML file
conda env export > rvfv_env.yml
```


## Task 2: Metadata Preparation with Bash

### Commands and Results

We used Bash commands to inspect and analyze the `rvf_africa.tsv` file.

#### Directory Organization

Files are organized in a `data/` directory for raw data and an `output/` directory for results.

#### File Inspection

Here are the commands used to understand the content of `rvf_africa.tsv`:

- **Number of unique GenBank accession numbers**:
  ```bash
  cut -f 5 rvf_africa.tsv | sort | uniq | wc -l
  ```
  **Result**: 1435 unique accession numbers.

- **Number of complete and partial sequences**:
  ```bash
  grep -i "complete" data/rvf_africa.tsv | wc -l
  ```
  **Result**: 833 complete sequences.
  ```bash
  grep -i "partial" rvf_africa.tsv | wc -l
  ```
  **Result**: 603 partial sequences.

- **Number of S, L, and M segments**:
  ```bash
  cut -f 4 rvf_africa.tsv | sort | uniq -c
  ```
  **Result**: 368 L segments, 476 M segments, 454 S segments.

- **Number of complete S, L, and M segments**:
  ```bash
  awk -F'\t' '$4 == "S" && $2 == "Complet" {print}' rvf_africa.tsv | wc -l
  ```
  **Result**: 319 complete S segments.
  ```bash
  awk -F'\t' '$4 == "L" && $2 == "Complet" {print}' rvf_africa.tsv | wc -l
  ```
  **Result**: 257 complete L segments.
  ```bash
  awk -F'\t' '$4 == "M" && $2 == "Complet" {print}' rvf_africa.tsv | wc -l
  ```
  **Result**: 254 complete M segments.

#### Interpretation

The file contains 1435 unique sequences, with 833 complete and 603 partial sequences. The segments are distributed as follows: 454 S, 368 L, and 476 M. Among them, 319 S, 257 L, and 254 M segments are complete, indicating a good representation of complete sequences for further analysis.

## Task 3: GenBank Data Retrieval with Python/Biopython

### Scripts and Results

#### Filtering Complete L Segment Sequences

A Bash script was used to filter complete genomes sequences:

```bash
#!/bin/bash
input_file="data/rvf_africa.tsv"
output_file="data/newcomplete.tsv"
awk -F'\t' '$2 == "Complete" {print $5}' "$input_file" > "$output_file"
echo "Filtered file generated: $output_file"
```

A Bash script was used to filter complete L segment sequences:

```bash
#!/bin/bash
input_file="data/rvf_africa.tsv"
output_file="data/newcompleteL.tsv"
awk -F'\t' '$2 == "Complete" && $4 == "L" {print $5}' "$input_file" > "$output_file"
echo "Filtered file generated: $output_file"
```


#### Downloading genomes

```python
#!/usr/bin/env python3

import argparse
from pathlib import Path
from Bio import Entrez, SeqIO


def download_genbank_records(accession_numbers, output_folder):
    Entrez.email = "christianpaulako@gmail.com"  # Update if needed

    output_path = Path(output_folder)
    output_path.mkdir(parents=True, exist_ok=True)  # Create output folder if not exists

    for accession in accession_numbers:
        try:
            handle = Entrez.efetch(db="nucleotide", id=accession, rettype="gb", retmode="text")
            record = SeqIO.read(handle, "genbank")
            filename = output_path / f"{accession}.gb"
            with open(filename, "w") as output_file:
                SeqIO.write(record, output_file, "genbank")
            print(f"[INFO] Downloaded {accession} → {filename}")
        except Exception as e:
            print(f"[ERROR] Could not download {accession}: {e}")


def main():
    parser = argparse.ArgumentParser(
        description="""
Download GenBank records using accession numbers from a file.
Each accession should be on a separate line.
The script will fetch .gb files and save them in the specified output directory.
""",
        epilog="""
Example:
  python download_genbank_records.py -acc accessions.txt -out ./genbank_files

Required:
  - A text file with accession numbers (one per line)
  - Internet connection (uses NCBI Entrez API)
"""
    )

    parser.add_argument("--accession", "-acc", required=True,
                        help="File containing GenBank accession numbers (one per line)")
    parser.add_argument("--output", "-out", required=True,
                        help="Output directory to save downloaded .gb files")

    args = parser.parse_args()

    try:
        with open(args.accession) as f:
            accession_numbers = [line.strip() for line in f if line.strip()]
    except FileNotFoundError:
        print(f"[ERROR] Accession file not found: {args.accession}")
        return

    try:
        download_genbank_records(accession_numbers, args.output)
        print(f"[DONE] All downloads saved in: {args.output}")
    except Exception as e:
        print(f"[FATAL] Script failed: {e}")


if __name__ == "__main__":
    main()
```



#### Downloading and Analyzing FASTA Files

A Python script using Biopython was used to downlod  sequences:

```python
#!/usr/bin/env python3
import argparse
from Bio import Entrez, SeqIO

def download_sequences(accession_numbers, output_folder):
    # Set your email address for Entrez
    Entrez.email = "christianpaulako@gmail.com"
    
    for accession in accession_numbers:
        try:
            handle = Entrez.efetch(db="nucleotide", id=accession, rettype="fasta", retmode="text")
            record = SeqIO.read(handle, "fasta")
            filename = f"{output_folder}/{accession}.fasta"
            with open(filename, "w") as output_file:
                SeqIO.write(record, output_file, "fasta")
            print(f"Downloaded {accession} and saved to {filename}")
        except Exception as e:
            print(f"Error downloading {accession}: {e}")

def main():
    description = (
        "Download FASTA sequences from GenBank using accession numbers."
        " Requires a file with accession numbers and an output folder."
    )
    parser = argparse.ArgumentParser(description=description)
    
    parser.add_argument("--accession", "-acc", required=True, help="Path to the file containing accession numbers.")
    parser.add_argument("--output", "-out", required=True, help="Path to the output folder where FASTA files will be saved.")
    
    args = parser.parse_args()

    try:
        with open(args.accession) as accession_file:
            accession_numbers = [line.strip() for line in accession_file.readlines()]
    except FileNotFoundError:
        print(f"Error: Accession file '{args.accession}' not found.")
        return

    try:
        download_sequences(accession_numbers, args.output)
        print("All sequences downloaded successfully.")
    except Exception as e:
        print(f"An error occurred: {e}")

if __name__ == "__main__":
    main()
```

A Python script using Biopython was used to analyze the L segment sequences:

```python
from Bio import SeqIO
import os
import glob

folder = "./L_fast"  # adjust to your directory path
fasta_files = glob.glob(os.path.join(folder, "*.fasta"))

all_lengths = []
for file in fasta_files:
    for record in SeqIO.parse(file, "fasta"):
        all_lengths.append(len(record.seq))

if all_lengths:
    mean_length = sum(all_lengths) / len(all_lengths)
    print(f"Total sequences: {len(all_lengths)}")
    print(f"Mean sequence length: {mean_length:.2f}")
else:
    print("No sequences found.")
```

**Result**: Total sequences: 257, Mean sequence length: 6399.75 nucleotides.

#### Interpretation

We retrieved 257 complete L segment sequences, with an average length of 6399.75 nucleotides, consistent with the expected size of RVFV L segments.

## Task 4: Data Cleaning and Analysis with R/tidyverse

### R Script and Results

Here is the R script used for data cleaning and analysis:

```{r task4, results='asis'}
# Load necessary packages
library(tidyverse)
library(dplyr)
library(ggplot2)
library(lubridate)

# Load the data
data <- read_tsv("rvf_africa.tsv")

# 1. Remove entries missing key fields
filtered_data <- data %>%
  filter(Segment != "" & Segment != "NA" & !is.na(Segment))
write_tsv(filtered_data, "filtered_rvf_data.tsv")

# 2. Group and summarize by country, year, and segment
summary_by_country_year_segment <- filtered_data %>%
  group_by(IsolationCountry, CollectionYear, Segment) %>%
  summarise(Count = n(), .groups = 'drop')

# 3. Create derived columns
filtered_data <- filtered_data %>%
  mutate(Region = case_when(
    IsolationCountry %in% c("Kenya", "Uganda", "Tanzania") ~ "East Africa",
    IsolationCountry %in% c("Nigeria", "Ghana") ~ "West Africa",
    TRUE ~ "Other"
  ))

# 4. Count sequences by host
host_summary <- filtered_data %>%
  group_by(HostCommonName) %>%
  summarise(Count = n(), .groups = 'drop') %>%
  arrange(desc(Count))

# 5. Filter specific records
filtered_specific <- filtered_data %>%
  filter(IsolationCountry %in% c("Kenya", "Tanzania"),
         CollectionYear >= 2007,
         Segment %in% c("S", "L"))

# 6. Rename collection_year to year
filtered_data <- filtered_data %>%
  rename(Year = CollectionYear)

# 7. Group and summarize by country, year, segment
summary_by_country_year_segment_renamed <- filtered_data %>%
  group_by(IsolationCountry, Year, Segment) %>%
  summarise(Count = n(), .groups = 'drop')

# 8. Complete sequences by country
complete_sequences_by_country <- filtered_data %>%
  filter(GenomeStatus == "Complete") %>%
  group_by(IsolationCountry) %>%
  summarise(Count = n(), .groups = 'drop') %>%
  arrange(desc(Count))

# 9. Top occurring hosts
top_hosts <- filtered_data %>%
  group_by(HostCommonName) %>%
  summarise(Count = n(), .groups = 'drop') %>%
  arrange(desc(Count)) %>%
  head(10)

# Descriptive analysis
isolates_by_country <- filtered_data %>%
  group_by(IsolationCountry) %>%
  summarise(Count = n(), .groups = 'drop') %>%
  arrange(desc(Count))

isolates_by_year <- filtered_data %>%
  group_by(Year) %>%
  summarise(Count = n(), .groups = 'drop') %>%
  arrange(desc(Count))

segments_distribution_by_country <- filtered_data %>%
  group_by(IsolationCountry, Segment) %>%
  summarise(Count = n(), .groups = 'drop') %>%
  pivot_wider(names_from = Segment, values_from = Count, values_fill = 0)

# Display results
cat("### Results\n")
cat("#### Summary by Country, Year, and Segment\n")
print(summary_by_country_year_segment)
cat("\n#### Host Summary\n")
print(host_summary)
cat("\n#### Filtered Specific Data (Kenya/Tanzania, >= 2007, S/L)\n")
print(filtered_specific)
cat("\n#### Summary After Renaming\n")
print(summary_by_country_year_segment_renamed)
cat("\n#### Complete Sequences by Country\n")
print(complete_sequences_by_country)
cat("\n#### Top 10 Hosts\n")
print(top_hosts)
cat("\n#### Isolates by Country\n")
print(isolates_by_country)
cat("\n#### Isolates by Year\n")
print(isolates_by_year)
cat("\n#### Segment Distribution by Country\n")
print(segments_distribution_by_country)
```

#### Interpretation

- **Data Cleaning**: Entries missing the `Segment` field were removed, ensuring data quality for subsequent analyses.
- **Summary by Country, Year, and Segment**: The data shows varied distribution across countries and years, potentially indicating specific outbreak foci.
- **Main Hosts**: The top 10 hosts reveal the species most affected by RVFV, crucial for epidemiological studies.
- **Geographic Distribution**: Isolates by country and segment distribution (S, M, L) highlight the predominance of certain segments in specific countries, possibly reflecting viral strain differences.

## Task 5: Data Visualization

### R Script and Results

Here is the R script for the visualizations:

```{r task5, fig.width=10, fig.height=6}
# 📅 Step 1: Clean and format dates
filtered_data <- filtered_data %>%
  mutate(CollectionDate = parse_date_time(
    CollectionDate,
    orders = c("ymd", "Ymd", "dmy", "dmY", "Y-m", "Y"),
    tz = "UTC"
  ))

# 🎁 1. Boxplot: Distribution of Collection Dates by Country
boxplot_date_by_country <- ggplot(filtered_data, aes(x = IsolationCountry, y = CollectionDate)) +
  geom_boxplot(fill = "#69b3a2", color = "black", outlier.shape = NA) +
  theme_minimal() +
  labs(title = "Distribution of Collection Dates by Country", x = "Country", y = "Collection Date") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
print(boxplot_date_by_country)
ggsave("boxplot_date_by_country.png", plot = boxplot_date_by_country, width = 10, height = 6)

# 🎁 2. Boxplot: Distribution of Collection Years by Segment
boxplot_year_by_segment <- ggplot(filtered_data, aes(x = Segment, y = Year)) +
  geom_boxplot(fill = "#f8766d", color = "black") +
  theme_minimal() +
  labs(title = "Distribution of Collection Years by Segment", x = "Segment", y = "Year")
print(boxplot_year_by_segment)
ggsave("boxplot_year_by_segment.png", plot = boxplot_year_by_segment, width = 8, height = 6)

# 📦 3. Barplot: Number of Isolates by Country
barplot_by_country <- ggplot(isolates_by_country, aes(x = reorder(IsolationCountry, -Count), y = Count)) +
  geom_bar(stat = "identity", fill = "#00BFC4") +
  theme_minimal() +
  labs(title = "Number of Isolates by Country", x = "Country", y = "Number of Isolates") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
print(barplot_by_country)
ggsave("barplot_isolates_by_country.png", plot = barplot_by_country, width = 10, height = 6)

# 📆 4. Barplot: Number of Isolates by Year
barplot_by_year <- ggplot(isolates_by_year, aes(x = Year, y = Count)) +
  geom_bar(stat = "identity", fill = "#C77CFF") +
  theme_minimal() +
  labs(title = "Number of Isolates by Year", x = "Year", y = "Number of Isolates")
print(barplot_by_year)
ggsave("barplot_isolates_by_year.png", plot = boxplot_year_by_segment, width = 8, height = 5)
```

#### Interpretation

- **Collection Dates by Country**: The boxplot shows temporal variations in collections by country, potentially indicating specific epidemic periods.
- **Collection Years by Segment**: The S, M, and L segments have similar temporal distributions, suggesting uniform collection of different segments over time.
- **Isolates by Country and Year**: The barplots highlight the countries and years with the most isolates, useful for identifying high-risk areas and periods.

## Conclusion

This report analyzed RVFV genomic data from Africa, covering metadata preparation, sequence retrieval, data cleaning, analysis, and visualization. The findings highlight important geographic and temporal trends for epidemiological surveillance.