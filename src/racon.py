#! /usr/bin/python

import os;
import sys;
import subprocess;

SCRIPT_PATH = os.path.dirname(os.path.realpath(__file__));

class Tle:
	def __init__(self):
		self.clr_s = 0;
		self.clr_e = 0;
		self.off = 0;
		self.src = 0;
		self.rvc = 0;

def parse_afg(afg_path):
	try:
		fp = open(afg_path, 'r');
	except:
		sys.stderr.write('ERROR: Could not open file "%s" for reading.\n' % (afg_path));
		exit(1);

	lines = fp.readlines();
	fp.close();

	reads = [];

	for line in lines:
		line = line.strip();
		if ('src' in line):
			qid = int(line.split(':')[-1]);
			reads.append(qid);
			print qid;
	# print '\n'.join(str(reads);

# ../samscripts/src/fastqfilter.py readid temp/ids.txt tests/sample-dataset/reads-lambda-R73.fasta > temp/reads_for_cons.fasta
# tools/graphmap/bin/graphmap-not_release -a anchor -b 3 -r tests/layout_20151114_221431/contigs_fast.fasta -d temp/reads_for_cons.fasta -o temp/aligned.sam
# tools/samscripts/src/consensus.py tests/layout_20151114_221431/contigs_fast.fasta 0 temp/consensus temp/aligned.sam
# src/consfromvcf.py  tests/layout_20151114_221431/contigs_fast.fasta temp/consensus-cov_0.variant.vcf temp/polished.fa
# dnadiff -p temp/dnadiff-ref_vs_unpolished/out tests/sample-dataset/NC_001416.fa tests/layout_20151114_221431/contigs_fast.fasta
# dnadiff -p temp/dnadiff-ref_vs_polished/out tests/sample-dataset/NC_001416.fa temp/polished.fa



def main():
	parse_afg('/home/isovic/work/eclipse-workspace/git/ra-consensus/tests/layout_20151114_221431/contigs.afg');

if __name__ == "__main__":
	main();
