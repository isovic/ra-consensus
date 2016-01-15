#! /bin/sh

../src/poaconsv3.py ../tests/miniasm/layout.fasta ../tests/sample-dataset/reads-lambda-R73-without_problematic.fasta results/poacons-miniasm.fasta

mkdir -p results/dnadiff-miniasm
dnadiff -p results/dnadiff-miniasm/consensus-miniasm ../tests/sample-dataset/NC_001416.fa results/poacons-miniasm.fasta
dnadiff -p results/dnadiff-miniasm/raw-miniasm ../tests/sample-dataset/NC_001416.fa ../tests/miniasm/layout.fasta
