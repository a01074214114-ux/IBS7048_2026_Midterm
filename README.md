# IBS7048_2026_Midterm
# IBS7048 Midterm Exam Repository

## Description
This repository contains code and output files for Problems 1–6 of the IBS7048 midterm exam.

## Uploaded Files
- `중간.R`: R script containing pipelines for sequence sorting, DADA2 processing, NJ tree generation, and metagenomic calculations (Problems 2, 5, 6).
- `top5_ASVs.nwk`: Newick format tree output file for top 5 ASVs (Problem 5c).

## How to Run the Code

### 1. File Sorting (Problem 1 & 2)
- **R sorting**: Run the sorting section in `중간.R` to process files into `sorted_R/` and cross-validate.

### 2. DADA2 Analysis & Tree Building (Problem 5)
- Run the DADA2 pipeline section in `중간.R` to generate the sequence table, perform chimera removal, and export `top5_ASVs.nwk`.

### 3. Metagenomic Calculations (Problem 6)
- Execute the calculation section in `중간.R` to obtain microbial read counts and fold changes.
