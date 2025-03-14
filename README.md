# Genome Assembly as Shortest Superstring

## Overview
This project is a solution for Rosalind's **"Genome Assembly as Shortest Superstring"** problem. It was developed as part of my **Fall 2021 Bioinformatics Class final project**.

## Problem Statement
Given a collection of DNA strings, the goal is to find the shortest possible superstring that contains all given sequences as substrings. This problem models the reconstruction of a chromosome from sequencing reads.

### Input:
- A set of **at most 50 DNA strings** in **FASTA format**.
- Each string has a length of approximately **1 kbp**.
- The dataset guarantees a unique way to reconstruct the chromosome by overlapping pairs of reads by **more than half their length**.

### Output:
- A **shortest superstring** containing all input DNA strings.

## Implementation
The solution:
- Uses **Python** and the **Biopython** library for parsing FASTA sequences.
- Implements a **recursive algorithm** to iteratively merge overlapping sequences into a single superstring.

## Dependencies
To run this project, you need:
- Python (>=3.x)
- Biopython (`conda install -c conda-forge biopython` or `pip install biopython`)

## Usage
Run the script using:
```sh
python genome_assembly.py
```

## Example
### Input:
```
>Rosalind_56
ATTAGACCTG
>Rosalind_57
CCTGCCGGAA
>Rosalind_58
AGACCTGCCG
>Rosalind_59
GCCGGAATAC
```
### Output:
```
ATTAGACCTGCCGGAATAC
```

## Author
Developed for **Fall 2021 Bioinformatics Class** final project.

