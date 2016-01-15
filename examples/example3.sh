#! /bin/sh

### Run the poaconsv3.py script to produce th econsensus, and evaluate it with dnadiff.

contigs="../tests/ecoli/contigs_fast.fasta"
reads="../tests/ecoli/reads/"
results="../tests/ecoli/results/poacons-ra.fasta"
reference="../tests/ecoli/escherichia_coli.fa"

/usr/bin/time ../src/poaconsv3.py $contigs $reads $results

mkdir -p results/dnadiff-ra
dnadiff -p results/dnadiff-ra/consensus-ra $reference $results
dnadiff -p results/dnadiff-ra/raw-ra $reference $contigs
