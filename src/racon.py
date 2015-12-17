#! /usr/bin/python

# Copyright Ivan Sovic, 2015. www.sovic.org

import os;
import sys;
import subprocess;

SCRIPT_PATH = os.path.dirname(os.path.realpath(__file__));
TOOLS_PATH = '%s/../tools/' % (SCRIPT_PATH);
SAMSCRIPTS_PATH = '%s/../tools/samscripts/src' % (SCRIPT_PATH);

sys.path.append(SAMSCRIPTS_PATH);

import fastqparser;

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
	qids = [];

	for line in lines:
		line = line.strip();
		if ('src' in line):
			qid = int(line.split(':')[-1]);
			reads.append(qid);
			# print qid;
			qids.append(qid);
	# print '\n'.join(str(reads);
	return qids;

# ../samscripts/src/fastqfilter.py readid temp/ids.txt tests/sample-dataset/reads-lambda-R73.fasta > temp/reads_for_cons.fasta
# tools/graphmap/bin/graphmap-not_release -a anchor -b 3 -r tests/layout_20151114_221431/contigs_fast.fasta -d temp/reads_for_cons.fasta -o temp/aligned.sam
# tools/samscripts/src/consensus.py tests/layout_20151114_221431/contigs_fast.fasta 0 temp/consensus temp/aligned.sam
# src/consfromvcf.py  tests/layout_20151114_221431/contigs_fast.fasta temp/consensus-cov_0.variant.vcf temp/polished.fa
# dnadiff -p temp/dnadiff-ref_vs_unpolished/out tests/sample-dataset/NC_001416.fa tests/layout_20151114_221431/contigs_fast.fasta
# dnadiff -p temp/dnadiff-ref_vs_polished/out tests/sample-dataset/NC_001416.fa temp/polished.fa


def execute_command_w_dryrun(dry_run, command):
	sys.stderr.write('[Executing] "%s"\n' % command);
	if (dry_run == False):
		subprocess.call(command, shell=True);
	sys.stderr.write('\n');

def execute_command_with_ret(dry_run, command):
	sys.stderr.write('Executing command: "%s"\n' % command);
	if (dry_run == False):
		p = subprocess.Popen(command, shell=True, stdin=subprocess.PIPE, stdout=subprocess.PIPE, stderr=subprocess.PIPE);
	[output, err] = p.communicate()
	rc = p.returncode
	sys.stderr.write('\n');
	return [rc, output, err];

def consensus_from_afg(contigs_file, afg_file, reads_file, out_file, ref_file=''):
	path_qids = '%s/tmp.pathreads.qid' % (os.path.dirname(out_file));
	path_reads = '%s/tmp.pathreads.fasta' % (os.path.dirname(out_file));
	path_aligns = '%s/tmp.pathreads.sam' % (os.path.dirname(out_file));
	path_variants = '%s/tmp.pathreads.vcf' % (os.path.dirname(out_file));
	path_dnadiff_raw = '%s/tmp.dnadiff/raw' % (os.path.dirname(out_file));
	path_dnadiff_polished = '%s/tmp.dnadiff/polished' % (os.path.dirname(out_file));

	if (not os.path.exists(os.path.dirname(path_dnadiff_raw))):
		os.makedirs(os.path.dirname(path_dnadiff_raw));

	reads_on_path = parse_afg(afg_file);

	try:
		fp = open(path_qids, 'w');
		fp.write('\n'.join([str(val) for val in reads_on_path]));
		fp.close();
	except:
		sys.stderr.write('ERROR: Could not open file "%s" for writing! Exiting.\n' % (path_qids));
		exit(1);

	dry_run = False;
	execute_command_w_dryrun(dry_run, '%s/fastqfilter.py readid %s %s > %s' % (SAMSCRIPTS_PATH, path_qids, reads_file, path_reads));
	execute_command_w_dryrun(dry_run, '%s/graphmap/bin/graphmap-not_release -a anchor -b 3 -r %s -d %s -o %s' % (TOOLS_PATH, contigs_file, path_reads, path_aligns));
	execute_command_w_dryrun(dry_run, '%s/denovoconsensus.py %s 0 %s %s' % (SCRIPT_PATH, contigs_file, path_variants, path_aligns));
	execute_command_w_dryrun(dry_run, '%s/consfromvcf.py %s %s %s' % (SCRIPT_PATH, contigs_file, path_variants, out_file));

	if (ref_file != ''):
		execute_command_w_dryrun(dry_run, 'dnadiff -p %s %s %s' % (path_dnadiff_raw, ref_file, contigs_file));
		execute_command_w_dryrun(dry_run, 'dnadiff -p %s %s %s' % (path_dnadiff_polished, ref_file, out_file));

def main():
	dry_run = False;

	# parse_afg('/home/isovic/work/eclipse-workspace/git/ra-consensus/tests/layout_20151114_221431/contigs.afg');

	# execute_command_with_ret(dry_run, '../samscripts/src/fastqfilter.py readid temp/ids.txt tests/sample-dataset/reads-lambda-R73.fasta > temp/reads_for_cons.fasta');
	# execute_command_with_ret(dry_run, 'tools/graphmap/bin/graphmap-not_release -a anchor -b 3 -r tests/layout_20151114_221431/contigs_fast.fasta -d temp/reads_for_cons.fasta -o temp/aligned.sam');
	# execute_command_with_ret(dry_run, 'tools/samscripts/src/consensus.py tests/layout_20151114_221431/contigs_fast.fasta 0 temp/consensus temp/aligned.sam');
	# execute_command_with_ret(dry_run, 'src/consfromvcf.py  tests/layout_20151114_221431/contigs_fast.fasta temp/consensus-cov_0.variant.vcf temp/polished.fa');
	# execute_command_with_ret(dry_run, 'dnadiff -p temp/dnadiff-ref_vs_unpolished/out tests/sample-dataset/NC_001416.fa tests/layout_20151114_221431/contigs_fast.fasta');
	# execute_command_with_ret(dry_run, 'dnadiff -p temp/dnadiff-ref_vs_polished/out tests/sample-dataset/NC_001416.fa temp/polished.fa');

	###

	# execute_command_with_ret(dry_run, 'tools/graphmap/bin/graphmap-not_release -a anchor -b 3 -r tests/layout_20151114_221431/contigs_fast.fasta -d tests/sample-dataset/reads-lambda-R73.fasta -o temp/aligned-all.sam');
	# execute_command_with_ret(dry_run, 'tools/graphmap/bin/graphmap-not_release -a anchorgotoh -b 3 -r tests/layout_20151114_221431/contigs_fast.fasta -d tests/sample-dataset/reads.fasta -o temp/aligned-all.sam -z 0 -c 40');
	# execute_command_with_ret(dry_run, 'tools/samscripts/src/consensus.py tests/layout_20151114_221431/contigs_fast.fasta 0 temp/consensus-all temp/aligned-all-anchorgotoh.sam');
	# execute_command_w_dryrun(dry_run, '%s/denovoconsensus.py tests/layout_20151114_221431/contigs_fast.fasta 0 temp/consensus-all temp/aligned-all-anchorgotoh.sam' % (SCRIPT_PATH));
	# execute_command_w_dryrun(dry_run, 'src/consfromvcf.py tests/layout_20151114_221431/contigs_fast.fasta temp/consensus-all-cov_0.variant.vcf temp/polished-all.fa');

	# execute_command_with_ret(dry_run, 'dnadiff -p temp/dnadiff-ref_vs_unpolished/out tests/sample-dataset/NC_001416.fa tests/layout_20151114_221431/contigs_fast.fasta');
	# execute_command_with_ret(dry_run, 'dnadiff -p temp/dnadiff-ref_vs_polished/out tests/sample-dataset/NC_001416.fa temp/polished.fa');
	# execute_command_w_dryrun(dry_run, 'mkdir -p temp/dnadiff-ref_vs_polished-all; dnadiff -p temp/dnadiff-ref_vs_polished-all/out tests/sample-dataset/NC_001416.fa temp/polished-all.fa');

	consensus_from_afg('tests/layout_20151114_221431/contigs_fast.fasta', 'tests/layout_20151114_221431/contigs.afg', 'tests/sample-dataset/reads-lambda-R73.fasta', 'temp/polished-on_path.fa', 'tests/sample-dataset/NC_001416.fa');

if __name__ == "__main__":
	main();
