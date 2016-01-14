#! /usr/bin/python

# Copyright Ivan Sovic, 2015. www.sovic.org

import os;
import sys;
import operator;
import subprocess;
from time import gmtime, strftime

SCRIPT_PATH = os.path.dirname(os.path.realpath(__file__));
TOOLS_PATH = '%s/../tools/' % (SCRIPT_PATH);
SAMSCRIPTS_PATH = '%s/../tools/samscripts/src' % (SCRIPT_PATH);

sys.path.append(SAMSCRIPTS_PATH);

import fastqparser;
import fastqfilter;
import utility_sam;

import altctg;

DEBUG_VERBOSE = False;
DEBUG_VERBOSE = True;

def execute_command(command):
	sys.stderr.write('[Executing] "%s"\n' % (command));
	subprocess.call(command, shell=True);
	sys.stderr.write('\n');

def align(ref, reads, out_sam):
	out_sam = os.path.abspath(out_sam);
	out_sam_basename = '%s/%s' % (os.path.dirname(out_sam), os.path.splitext(os.path.basename(out_sam))[0]);
	path_aligns = ('%s.unsorted.sam' % (out_sam_basename)) if (out_sam_basename.endswith('.sorted') == False) else ('%s.unsorted.sam' % (out_sam_basename.split('.sorted')[0]));

	execute_command('%s/graphmap/bin/Linux-x64/graphmap -a anchor -b 3 -r %s -d %s -o %s' % (TOOLS_PATH, ref, reads, path_aligns));
	execute_command('samtools view -Sb %s | samtools sort - %s && samtools view -h %s.bam > %s' % (path_aligns, out_sam_basename, out_sam_basename, out_sam));

def fragment_contig(contig_raw, breakpoints):
	ret = [];
	if (len(breakpoints) == 0): return [];
	ret.append(contig_raw[0:breakpoints[0]]);
	for i in xrange(1, len(breakpoints)):
		ret.append(contig_raw[breakpoints[i-1]:breakpoints[i]]);
	ret.append(contig_raw[breakpoints[-1]:-1]);
	return ret;

def consensus_from_pir(pir_path):
	[headers, seqs, quals] = fastqparser.read_fastq(pir_path);

	cons_seq = '';
	for i in xrange(0, len(seqs[0])):
		base_counts = {'A': 0, 'C': 0, 'T': 0, 'G': 0, '.': 0};
		for j in xrange(0, len(seqs)):
			base_counts[seqs[j][i]] += 1;
		sorted_base_counts = sorted(base_counts.items(), key=operator.itemgetter(1));
		# print sorted_base_counts;
		if (sorted_base_counts[-1][0] != '.'):
			cons_seq += sorted_base_counts[-1][0]

	return cons_seq;

def run(raw_contig_path, reads_path, out_consensus_path):
	out_folder = os.path.dirname(out_file);
	temp_folder = '%s/temp' % (out_folder);
	out_basename = os.path.splitext(os.path.basename(out_file))[0];

	if (not os.path.exists(out_folder)):
		os.path.makedirs(os.path.dirname(out_folder));
	if (not os.path.exists(temp_folder)):
		os.path.makedirs(os.path.dirname(temp_folder));

	# ### Input paths
	# raw_contig = '/home/isovic/work/eclipse-workspace/git/ra-consensus/tests/miniasm/layout.fasta';
	# reads_path = '/home/isovic/work/eclipse-workspace/git/ra-consensus/tests/sample-dataset/reads-lambda-R73.fasta';
	# ### Output path.
	# out_consensus_path = '/home/isovic/work/eclipse-workspace/git/ra-consensus/tests/out/temp-new1/consensus.fasta';

	### Paths to temporary files, useful for debug purposes.
	sorted_alns = '%s/tmp.allreads.sorted.sam' % (temp_folder);
	alt_ctg_path = '%s/alt_ctgs.fasta' % (temp_folder);
	frags_path = '%s/frags_for_msa.fasta' % (temp_folder);
	pir_path = '%s/frags_for_msa.pir' % (temp_folder);

	align(raw_contig_path, reads_path, sorted_alns);

	[ctg_headers, ctg_seqs, ctg_quals] = fastqparser.read_fastq(raw_contig_path);
	for i in xrange(0, len(ctg_headers)):
		contig_name = ctg_headers[i];
		sys.stderr.write('Testing for contig_name = "%s", contig_length = %d\n' % (contig_name, len(ctg_seqs[0])));

		### Extract only alignments for the current contig.
		all_contig_sams = altctg.extract_alns_for_contig(contig_name, sorted_alns, 0.70, 0.40);
		sys.stderr.write('Total %d alignments for contig "%s".\n' % (len(all_contig_sams), contig_name));

		### Determine new alternate contigs from given alignments, in the form of an array of SAM lines.
		alt_contigs = altctg.extract_alternate_contigs(all_contig_sams, 0.80, len(ctg_seqs[0]), 0.01);

		### Process alternate contigs to create actual sequences and corresponding coordinates for subsampling.
		all_alt_contig_frags = [];
		prev_sample_pos_keys = [];
		fp_ctg = open(alt_ctg_path, 'w');
		for i in xrange(0, len(alt_contigs)):
		# for i in xrange(5, 6):
			if (DEBUG_VERBOSE == True): sys.stderr.write('Alternate contig %d subsampled starting positions:\n' % (i));
			alt_contig_seq = altctg.construct_contig_from_sams(ctg_seqs[0], alt_contigs[i]);
			ctg_sample_pos = altctg.subsample_ref_positions(ctg_seqs[0], alt_contigs[i], 5000);

			### Sanity check that the reference sample positions are the same for every contig.
			sample_pos_keys = sorted(ctg_sample_pos.keys());
			if (i > 0 and sample_pos_keys != prev_sample_pos_keys):
				sys.stderr.write('ERROR: Sample positions calculated wrongly - missing a sample key for contig #%d!\n' % (i));
			prev_sample_pos_keys = sample_pos_keys;

			### Split the contig into fragments between the specified sampling positions, and store them for later.
			sample_pos_on_ctg = [ctg_sample_pos[val] for val in sample_pos_keys];
			alt_contig_frags = fragment_contig(alt_contig_seq, sample_pos_on_ctg);
			all_alt_contig_frags.append(alt_contig_frags);

			### Write contig to disk for debug purposes.
			fp_ctg.write('>Alternate contig %d, len: %d\n%s\n' % (i, len(alt_contig_seq), alt_contig_seq));

			### Just verbose.
			if (DEBUG_VERBOSE == True):
				for sample_pos_key in sorted(ctg_sample_pos.keys()):
					sys.stderr.write('[%d, %d]\n' % (sample_pos_key, ctg_sample_pos[sample_pos_key]));
				sys.stderr.write('[%d, %d] (end)\n' % (len(ctg_seqs[0]), len(alt_contig_seq)));
				sys.stderr.write('\n');
		fp_ctg.close();

		### Process all fragments with a MSA tool. First, for every fragment, all alternate sequences need to be output.
		fp_consensus = open(out_consensus_path, 'w');
		timestamp = strftime("%Y/%m/%d %H:%M:%S", gmtime());
		fp_consensus.write('>Consensus_with_POA %s\n' % (timestamp));
		num_fragments = (0) if (len(all_alt_contig_frags) == 0) else (len(all_alt_contig_frags[0]));
		for i in xrange(0, num_fragments):
			if (DEBUG_VERBOSE == True): sys.stderr.write('Writing fragment #%d for alternate contigs.\n' % (i));
			fp_frags = open(frags_path, 'w');
			for j in xrange(0, len(all_alt_contig_frags)):
				fp_frags.write('>Fragment %d for alt contig %d\n' % (i, j));
				fp_frags.write('%s\n' % (all_alt_contig_frags[j][i]));
			fp_frags.close();
			if (DEBUG_VERBOSE == True): sys.stderr.write('Fragments written, proceeding to MSA.\n');

			execute_command('%s/poaV2/poa -do_global -do_progressive -read_fasta %s -pir %s %s/poaV2/blosum80.mat' % (TOOLS_PATH, frags_path, pir_path, TOOLS_PATH));
			# execute_command('%s/poaV2/poa -do_global -do_progressive -read_fasta %s -pir %s %s/poaV2/all1.mat' % (TOOLS_PATH, temp_subseq_file, temp_msa_file, TOOLS_PATH));
			frag_consensus = consensus_from_pir(pir_path);
			fp_consensus.write('%s' % (frag_consensus));
		fp_consensus.close();

		# execute_command('dnadiff -p dnadiff/consensus ../../sample-dataset/NC_001416.fa consensus.fasta');
		# execute_command('dnadiff -p dnadiff/raw-miniasm ../../sample-dataset/NC_001416.fa ../../miniasm/layout.fasta');

def main():
	# ### Input paths
	# raw_contig = '/home/isovic/work/eclipse-workspace/git/ra-consensus/tests/miniasm/layout.fasta';
	# reads_path = '/home/isovic/work/eclipse-workspace/git/ra-consensus/tests/sample-dataset/reads-lambda-R73.fasta';
	# ### Output path.
	# out_consensus_path = '/home/isovic/work/eclipse-workspace/git/ra-consensus/tests/out/temp-new1/consensus.fasta';

	if (len(sys.argv) != 4):
		sys.stderr.write('Proof-of-Concept Consensus for de novo genome assemblies.\n');
		sys.stderr.write('This consensus tool first creates several instances of the contig through alignment (e.g. 10x; these instances have similar error rate to the original data), and then applies POA on chunks (windows) of alternate contigs to produce consensus.\n');
		sys.stderr.write('Usage:\n');
		sys.stderr.write('\t%s <contigs.fasta> <reads.fasta> <out_polished.fasta>\n' % (sys.argv[0]));
		exit(1);

	contigs_file = os.path.abspath(sys.argv[1]);
	reads_file = os.path.abspath(sys.argv[2]);
	out_file = os.path.abspath(sys.argv[3]);

	run(contigs_file, reads_file, out_file);

if __name__ == "__main__":
	main();
