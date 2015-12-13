#! /usr/bin/python

import os;
import sys;
import subprocess;

SCRIPT_PATH = os.path.dirname(os.path.realpath(__file__));
TOOLS_PATH = '%s/../tools/' % (SCRIPT_PATH);

def execute_command(command):
	sys.stderr.write('[Executing] "%s"\n' % (command));
	subprocess.call(command, shell=True);
	sys.stderr.write('\n');

def make_consensus_reference_from_vcf(reference_file, vcf_file, out_consensus_sequence_file):
	if (not os.path.exists(os.path.dirname(out_consensus_sequence_file))):
		sys.stderr.write('Creating a folder on path: "%s".\n' % (os.path.dirname(out_consensus_sequence_file)));
		os.makedirs(os.path.dirname(out_consensus_sequence_file));

	sys.stderr.write('%s\n' % (reference_file));
	sys.stderr.write('%s\n' % (vcf_file));
	sys.stderr.write('%s\n' % (out_consensus_sequence_file));
	sys.stderr.write('\n');

	fp = open(vcf_file, 'r');
	fp_temp = open('.temp.vcf', 'w');
	for line in fp:
		line = line.strip();
		if (len(line) == 0 or (len(line) > 0 and line[0] == '#')):
			fp_temp.write(line + '\n');
			continue;
		split_line = line.split('\t');
		if (split_line[3] == 'N' and split_line[4] == 'N'):
			continue;
		fp_temp.write(line + '\n');
	fp.close();
	fp_temp.close();
	os.rename(vcf_file, '%s.bak' % (vcf_file));
	os.rename('.temp.vcf', vcf_file);

	sys.stderr.write('Making a Picard dictionary of the reference.\n');
	execute_command('java -jar %s/picard-tools-1.138/picard.jar CreateSequenceDictionary R= %s O= %s.dict' % (TOOLS_PATH, reference_file, os.path.splitext(reference_file)[0]));
	execute_command('samtools faidx %s' % (reference_file));

	sys.stderr.write('Applying the VCF file to the reference to generate the alternate (consensus) sequence.\n');
	execute_command('java -jar %s/GenomeAnalysisTK-3.4-46/GenomeAnalysisTK.jar -T FastaAlternateReferenceMaker -R %s -o %s -V %s' % (TOOLS_PATH, reference_file, out_consensus_sequence_file, vcf_file));

if __name__ == "__main__":
	if (len(sys.argv) != 4):
		sys.stderr.write('Uses GATK to construct a mutated reference from the original reference given a VCF file with mutations.\n');
		sys.stderr.write('Usage:\n');
		sys.stderr.write('\t%s <reference.fa> <vcf_file> <out_mutated_reference.fa>\n' % (sys.argv[0]));
		exit(1);

	make_consensus_reference_from_vcf(sys.argv[1], sys.argv[2], sys.argv[3]);
