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

def execute_command(command):
	sys.stderr.write('[Executing] "%s"\n' % (command));
	subprocess.call(command, shell=True);
	sys.stderr.write('\n');

### Mutates the given reference, and writes the mutations in a vcf file.
### Requires Mutatrix - it needs to be placed in the TOOLS_PATH.
def generate_mutated_reference(reference_path, snp_rate, indel_rate, out_file): # , out_path):
	out_file = os.path.abspath(out_file);
	reference_path = os.path.abspath(reference_path);
	out_path = os.path.abspath(os.path.dirname(out_file));
	if (not os.path.exists(out_path)):
		os.makedirs(out_path);

	# out_prefix = '%s/mutated_%s_snp%f_indel%f' % (out_path, os.path.splitext(os.path.basename(reference_path))[0], snp_rate, indel_rate);
	# out_vcf = os.path.abspath('%s.vcf' % (out_prefix));
	out_vcf = '%s.vcf' % (os.path.splitext(out_file)[0]);
	out_rev_vcf = '%s.rev.vcf' % (os.path.splitext(out_file)[0]);
	ref_ext = os.path.splitext(reference_path)[-1];
	out_ref_file = out_file; # '%s%s' % (out_prefix, ref_ext);

	sys.stderr.write('Mutating the reference using Mutatrix, output VCF file: "%s".\n' % (out_vcf));
	if (indel_rate != 0):
		execute_command('cd %s; %s/mutatrix/mutatrix --snp-rate %f --population-size 1 --microsat-min-len 0 --mnp-ratio 0 --indel-rate %f --indel-max 10 %s > %s' % (out_path, TOOLS_PATH, snp_rate, indel_rate, reference_path, out_vcf));
	else:
		execute_command('cd %s; %s/mutatrix/mutatrix --snp-rate %f --population-size 1 --microsat-min-len 0 --mnp-ratio 0 --indel-rate 0 --indel-max 0 %s > %s' % (out_path, TOOLS_PATH, snp_rate, reference_path, out_vcf));

	sys.stderr.write('Reversing the SNP bases in the VCF file, output VCF file: "%s".\n' % (out_rev_vcf));
	execute_command(r"cat %s | awk -F '\t' 'BEGIN {OFS = FS} {if ($0 == /^#.*/) print ; else {a=$4; $4=$5; $5=a; print } }' > %s" % (out_vcf, out_rev_vcf));

	sys.stderr.write('Compressing and indexing the VCF file.\n');
	execute_command('bgzip -c %s > %s.gz' % (out_rev_vcf, out_rev_vcf));
	execute_command('tabix -p vcf %s.gz' % (out_rev_vcf));

	### Mutatrix splits all reference sequences into separate files. This part of code joins them back into one FASTA file.
	[headers, lengths] = fastqparser.get_headers_and_lengths(reference_path);
	print headers;
	all_files = ['"%s/1:%s:0%s"' % (out_path, header.split(' ')[0], ref_ext) for header in headers];
	# all_files = ['"%s/1:%s:0%s"' % (out_path, header, ref_ext) for header in headers];
	if (os.path.exists(out_ref_file)):
		os.rename(out_ref_file, '%s.bak' % (out_ref_file));
	for ref_file in all_files:
		### Take care of the special characters.
		# escaped_ref_file = ref_file.replace('|', '\|');
		escaped_ref_file = ref_file;
		print escaped_ref_file;
		execute_command('cat %s >> %s' % (escaped_ref_file, out_ref_file));
		if (len(ref_file) > 0 and ('*' in ref_file) == False):
			print 'Removing file: "%s".' % (ref_file);
#			os.remove(escaped_ref_file);
			execute_command('rm %s' % (escaped_ref_file));

def TEST_SIMULATE():
	for i in xrange(0, 10):
		# generate_mutated_reference(('%s/../reference-genomes/escherichia_coli.fa' % SCRIPT_PATH), 0.0006, 0.0067, 'data/test-msa/ecoli-%d.fa' % (i));
		# seqs_for_mafft
		if (i == 0):
			execute_command('cat %s > %s' % ('data/test-msa/ecoli-%d.fa' % (i), 'data/test-msa/ecoli-all.fa'));
		else:
			execute_command('cat %s >> %s' % ('data/test-msa/ecoli-%d.fa' % (i), 'data/test-msa/ecoli-all.fa'));

if __name__ == "__main__":
	TEST_SIMULATE();

	# if (len(sys.argv) != 4):
	# 	sys.stderr.write('Uses GATK to construct a mutated reference from the original reference given a VCF file with mutations.\n');
	# 	sys.stderr.write('Usage:\n');
	# 	sys.stderr.write('\t%s <reference.fa> <vcf_file> <out_mutated_reference.fa>\n' % (sys.argv[0]));
	# 	exit(1);

	# make_consensus_reference_from_vcf(sys.argv[1], sys.argv[2], sys.argv[3]);
