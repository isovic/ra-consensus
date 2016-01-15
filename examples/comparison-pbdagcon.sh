#! /bin/sh

### Run Pbdagcon to produce th econsensus, and evaluate it with dnadiff.

contigs="../tests/layout_20151114_221431/contigs_fast.fasta"
reads="../tests/sample-dataset/reads-lambda-R73-without_problematic.fasta"
results="results/consensus-pbdagcon.fasta"
reference="../tests/sample-dataset/NC_001416.fa"
mappings="results/temp/mapped-pbdagcon.m5"

/usr/bin/time ../comparison/blasr/alignment/bin/blasr $reads $contigs -bestn 1 -m 5 -out $mappings
/usr/bin/time ../comparison/pbdagcon/src/cpp/pbdagcon $mappings -j 1 > $results

mkdir -p results/dnadiff-ra
dnadiff -p results/dnadiff-ra/consensus-pbdagcon $reference $results
dnadiff -p results/dnadiff-ra/raw-ra $reference $contigs
