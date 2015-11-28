reduce_redundancy
=================
## Purpose
The goal of this script is to reduce the redundancy in a set of sequences.

## Input
- A FASTA file of sequences.  FASTA entry descriptions are limited to less than 128 characters.
- Distance threshold.

## Process
Iterate through each transcript in the dictionary and:
 1. Generate a distance matrix from the sequences in the input file using Clustal Omega.
 2. Parse distance matrix to dictionary.
 3. For each FASTA entry, identify as a cluster those other sequences with a distance less than given threshold.
 4. For each sequence in each cluster, generate a total combined distance score to the other members in the cluster.
 5. Keep the entry with the lowest distance to the other members of the cluster as a representative sequence.  Ignore the other members of the cluster.
 6. Keep all other entries which have not included in any cluster.

## Output
- A FASTA file of all cluster representative sequences and non-clustered sequences.

## Dependencies
- Currently only tested on Linux
- Python 2.7
- Biopython (http://biopython.org/wiki/Download#Installation_Instructions)
- clustalo `sudo apt-get install clustalo`

## Usage
- Edit the input_seqs_filename variable of the 'Initiating Variables' section to point to the input sequences.
- 12,610 related protein sequences (retrieved by iterative BLAST search) reduced to 2,048 sequences, using a distance threshold of 0.10, in 20 minutes on an i7 processor.
