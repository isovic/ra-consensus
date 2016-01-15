#! /bin/sh

### Run the poaconsv3.py script to produce th econsensus, and evaluate it with dnadiff.

contigs="../tests/layout_20151114_221431/contigs_fast.fasta"
reads="../tests/sample-dataset/reads-lambda-R73-without_problematic.fasta"
results="results/poacons-ra.fasta"
reference="../tests/sample-dataset/NC_001416.fa"

/usr/bin/time ../src/poaconsv3.py $contigs $reads $results

mkdir -p results/dnadiff-ra
dnadiff -p results/dnadiff-ra/consensus-ra $reference $results
dnadiff -p results/dnadiff-ra/raw-ra $reference $contigs
