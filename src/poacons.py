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

def run_poa_sequentially(seq_path):
	temp_subseq_file = '%s/tmp.subseq.fasta' % (os.path.dirname(seq_path));
	temp_msa_file = '%s/tmp.subseq.fasta.pir' % (os.path.dirname(seq_path));
	out_consensus_file = '%s/consensus-poa.fasta' % (os.path.dirname(seq_path));
	out_consensus_file_chunks = '%s/consensus-poa-chunks.fasta' % (os.path.dirname(seq_path));

	fp_out_all = open(out_consensus_file, 'a');
	fp_out_chunks = open(out_consensus_file_chunks, 'a');

	timestamp = strftime("%Y/%m/%d %H:%M:%S", gmtime());
	fp_out_all.write('>Consensus_with_POA all %s\n' % (timestamp));

	[ret_string, num_seqs, total_seq_len, average_seq_len, max_seq_len] = fastqparser.count_seq_length(seq_path);

	window_len = 1000;

	start_coord = 0;
	while (start_coord < max_seq_len):
		end_coord = start_coord + window_len;
		if (end_coord > (max_seq_len - window_len)):
			end_coord = max_seq_len;

		sys.stderr.write('Window: start = %d, end = %d\n' % (start_coord, end_coord));
		execute_command('%s/fastqfilter.py subseqs %s %d %d %s' % (SAMSCRIPTS_PATH, seq_path, start_coord, end_coord, temp_subseq_file));

		execute_command('%s/poaV2/poa -read_fasta %s -hb -pir %s %s/poaV2/blosum80.mat' % (TOOLS_PATH, temp_subseq_file, temp_msa_file, TOOLS_PATH));

		timestamp = strftime("%Y/%m/%d %H:%M:%S", gmtime());
		fp_out_chunks.write('>Consensus_with_POA %d-%d %s\n' % (start_coord, end_coord, timestamp));
		[headers, seqs, quals] = fastqparser.read_fastq(temp_msa_file);
		# print temp_subseq_file;
		# print headers;
		i = 0;
		while (i < len(headers)):
			if ('consensus' in headers[i]):
				# print seqs[i];
				# print seqs[i].replace('.', '');
				chunk_seq = seqs[i].replace('.', '');
				fp_out_all.write('%s' % (chunk_seq));
				fp_out_chunks.write('%s\n' % (chunk_seq));
				break;
			i += 1;

		# break;
		start_coord = end_coord;

	fp_out_all.write('\n');
	fp_out_all.close();
	fp_out_chunks.close();

def run_poa_sequentially_v2(seq_path):
	temp_subseq_file = '%s/tmp.subseq.fasta' % (os.path.dirname(seq_path));
	temp_msa_file = '%s/tmp.subseq.fasta.pir' % (os.path.dirname(seq_path));
	out_consensus_file = '%s/consensus-poa.fasta' % (os.path.dirname(seq_path));
	out_consensus_file_chunks = '%s/consensus-poa-chunks.fasta' % (os.path.dirname(seq_path));

	fp_out_all = open(out_consensus_file, 'a');
	fp_out_chunks = open(out_consensus_file_chunks, 'a');

	timestamp = strftime("%Y/%m/%d %H:%M:%S", gmtime());
	fp_out_all.write('>Consensus_with_POA all %s\n' % (timestamp));

	[ret_string, num_seqs, total_seq_len, average_seq_len, max_seq_len] = fastqparser.count_seq_length(seq_path);

	window_len = 5000;
	# window_len = 1000;
	# window_len = max_seq_len;

	start_coord = 0;
	while (start_coord < max_seq_len):
		end_coord = start_coord + window_len;
		if (end_coord > (max_seq_len - window_len)):
			end_coord = max_seq_len;

		sys.stderr.write('Window: start = %d, end = %d\n' % (start_coord, end_coord));
		execute_command('%s/fastqfilter.py subseqs %s %d %d %s' % (SAMSCRIPTS_PATH, seq_path, start_coord, end_coord, temp_subseq_file));

		if (start_coord == 0 or end_coord == max_seq_len):
			# execute_command('%s/poaV2/poa -do_progressive -read_fasta %s -pir %s %s/poaV2/blosum80.mat' % (TOOLS_PATH, temp_subseq_file, temp_msa_file, TOOLS_PATH));
			execute_command('%s/poaV2/poa -do_progressive -read_fasta %s -pir %s %s/poaV2/all1.mat' % (TOOLS_PATH, temp_subseq_file, temp_msa_file, TOOLS_PATH));
		else:
			# execute_command('%s/poaV2/poa -do_global -do_progressive -read_fasta %s -pir %s %s/poaV2/blosum80.mat' % (TOOLS_PATH, temp_subseq_file, temp_msa_file, TOOLS_PATH));
			execute_command('%s/poaV2/poa -do_global -do_progressive -read_fasta %s -pir %s %s/poaV2/all1.mat' % (TOOLS_PATH, temp_subseq_file, temp_msa_file, TOOLS_PATH));

		timestamp = strftime("%Y/%m/%d %H:%M:%S", gmtime());
		fp_out_chunks.write('>Consensus_with_POA %d-%d %s\n' % (start_coord, end_coord, timestamp));
		[headers, seqs, quals] = fastqparser.read_fastq(temp_msa_file);

		cons_seq = '';
		for i in xrange(0, len(seqs[0])):
			base_counts = {'A': 0, 'C': 0, 'T': 0, 'G': 0, '.': 0};
			for j in xrange(0, len(seqs)):
				base_counts[seqs[j][i]] += 1;
			sorted_base_counts = sorted(base_counts.items(), key=operator.itemgetter(1));
			# print sorted_base_counts;
			if (sorted_base_counts[-1][0] != '.'):
				cons_seq += sorted_base_counts[-1][0]

		fp_out_all.write('%s' % (cons_seq));
		fp_out_chunks.write('%s\n' % (cons_seq));

		# # print temp_subseq_file;
		# # print headers;
		# i = 0;
		# while (i < len(headers)):
		# 	if ('consensus' in headers[i]):
		# 		# print seqs[i];
		# 		# print seqs[i].replace('.', '');
		# 		chunk_seq = seqs[i].replace('.', '');
		# 		fp_out_all.write('%s' % (chunk_seq));
		# 		fp_out_chunks.write('%s\n' % (chunk_seq));
		# 		break;
		# 	i += 1;

		# break;
		start_coord = end_coord;

	fp_out_all.write('\n');
	fp_out_all.close();
	fp_out_chunks.close();



def check_overlap(sam_line1, sam_line2):
	start1 = sam_line1.pos;
	end1 = start1 + len(sam_line1.seq) - sam_line1.clip_count_front - sam_line1.clip_count_back;

	start2 = sam_line2.pos;
	end2 = start2 + len(sam_line2.seq) - sam_line2.clip_count_front - sam_line2.clip_count_back;

	if (start2 > end1 or end2 < start1):
		# print 'start1 = %d, end1 = %d, start2 = %d, end2 = %d' % (start1, end1, start2, end2);
		return 0;
	if ((start2 >= start1 and end2 <= end1) or (start1 >= start2 and end1 <= end2)):
		# print 'start1 = %d, end1 = %d, start2 = %d, end2 = %d' % (start1, end1, start2, end2);
		return -1;

	return (end2 - start1);

def construct_contig_from_overlapping_sams(ctg_seqs, contig_sams):
	new_contig = '';

	overhang_before = contig_sams[0].pos;
	# new_contig += ctg_seqs[0][0:overhang_before];

	non_clipped_len = 0;

	i = 0;
	while (i < len(contig_sams)):
		sam_line = contig_sams[i];
		start_pos = sam_line.clip_count_front;
		end_pos = (len(sam_line.seq) - sam_line.clip_count_back) if ((i + 1) == len(contig_sams)) else (contig_sams[i+1].pos - sam_line.pos);
		new_chunk = sam_line.seq[start_pos:end_pos];
		non_clipped_len += len(new_chunk);
		new_contig += new_chunk;
		i += 1;

	overhang_after = contig_sams[-1].pos + len(contig_sams[-1].seq) - contig_sams[-1].clip_count_front - contig_sams[-1].clip_count_back;
	# new_contig += ctg_seqs[0][overhang_after:-1];

	return [new_contig, non_clipped_len];

def extract_same_contigs(contigs_file, reads_file, out_sam_file, out_alt_ctg_file, ref_file=''):
	path_aligns_basename = '%s/temp/tmp.allreads' % (os.path.dirname(out_alt_ctg_file));
	path_aligns = '%s/temp/tmp.allreads.sam' % (os.path.dirname(out_alt_ctg_file));
	path_aligns_sorted = '%s/temp/tmp.allreads.sorted.sam' % (os.path.dirname(out_alt_ctg_file));

	# if (not os.path.exists(os.path.dirname(path_aligns))):
	# 	os.makedirs(os.path.dirname(path_aligns));

	# execute_command('%s/graphmap/bin/graphmap-not_release -a anchor -z 0 -c 40 -b 3 -r %s -d %s -o %s' % (TOOLS_PATH, contigs_file, reads_file, path_aligns));

	# execute_command('%s/graphmap/bin/Linux-x64/graphmap -a anchor -z 0 -c 40 -b 3 -r %s -d %s -o %s' % (TOOLS_PATH, contigs_file, reads_file, path_aligns));
	# execute_command('%s/graphmap/bin/Linux-x64/graphmap -a anchor -b 3 -r %s -d %s -o %s' % (TOOLS_PATH, contigs_file, reads_file, path_aligns));
	execute_command('/home/isovic/work/eclipse-workspace/graphmap/bin/graphmap-not_release -a anchor -b 3 -r %s -d %s -o %s' % (contigs_file, reads_file, path_aligns));
	exit(1);
	execute_command('samtools view -Sb %s | samtools sort - %s.sorted && samtools view %s.sorted.bam > %s.sorted.sam' % (path_aligns, path_aligns_basename, path_aligns_basename, path_aligns_basename));

	[ctg_headers, ctg_seqs, ctg_quals] = fastqparser.read_fastq(contigs_file);

	[headers, all_sam_lines] = utility_sam.LoadSAM(path_aligns_sorted);

	# sam_lines = all_sam_lines;

	print 'len(all_sam_lines) = %d' % (len(all_sam_lines));

	sam_lines = [];
	for sam_line in all_sam_lines:
		cigop_counts = sam_line.CountAlignmentOps();
		print cigop_counts;
		try:
			matches = cigop_counts['='];
			errors = cigop_counts['X'] + cigop_counts['D'] + cigop_counts['I'];
		except:
			continue;
		seq_len = len(sam_line.seq) - sam_line.clip_count_front - sam_line.clip_count_back;
		# print 'matches = %d (%f)' % (matches, float(matches) / float(seq_len));
		# print 'errors = %d (%f)' % (errors, float(errors) / float(seq_len));
		# if (sam_line.qname == 'channel_416_read_32_twodirections'):
		# 	print sam_line.qname;
		# 	print cigop_counts;
		# 	print seq_len;
		# 	exit(1);

		if (float(matches) / float(seq_len) >= 0.70 and float(errors) / float(seq_len) < 0.40):
			sam_lines.append(sam_line);

	print 'len(sam_lines) = %d' % (len(sam_lines));
	print headers;

	fp_out_alt_ctg = open(out_alt_ctg_file, 'w');
	fp_out = open(out_sam_file, 'w');
	fp_out.write('\n'.join(headers) + '\n');

	sams_to_process = sam_lines;
	coverage = 0;
	while (coverage < 100 and len(sams_to_process) > 0):
		coverage += 1;
		contig_sams = [];
		unused_sams = [];
		i = 0;
		candidate_read = i;
		contig_sams.append(sams_to_process[candidate_read]);
		while ((candidate_read + 1) < len(sams_to_process)):
			max_overlap_len = 0;
			max_overlap_id = -1;
			j = candidate_read + 1;
			while (j < len(sams_to_process)):
				overlap_len = check_overlap(sams_to_process[candidate_read], sams_to_process[j]);
				if (overlap_len == 0):
					print 'break 1';
					print '  j = %d' % (j);
					break;
				elif (overlap_len == -1):
					j += 1;
					continue;

				if (max_overlap_id == -1 or overlap_len >= max_overlap_len):
					max_overlap_len = overlap_len;
					max_overlap_id = j;
				j += 1;
			if (max_overlap_id > 0):
				print '  candidate_read = %d' % (max_overlap_id);
				print '  unused reads: %d - %d' % ((candidate_read + 1), max_overlap_id);
				unused_sams += sams_to_process[(candidate_read + 1):max_overlap_id];
				candidate_read = max_overlap_id;
				contig_sams.append(sams_to_process[candidate_read]);
			else:
				print 'break 2';
				break;

		print '  unused reads: %d - %d' % ((candidate_read + 1), len(sams_to_process));
		unused_sams += sams_to_process[(candidate_read + 1):len(sams_to_process)];
		sams_to_process = unused_sams + [];

		# if ((candidate_read + 1) == len(sam_lines)):
		# 	break;

		# i += 1;
		# max_overlap_len = 0;
		# max_overlap_id = -1;
		# while (i < len(sam_lines)):
		# 	overlap_len = check_overlap(sam_lines[candidate_read], sam_lines[i + 1]);
		# 	if ((i + 1) >= len(sam_lines) or overlap_len <= 0):
		# 		break;
		# 	else:
		# 		unused_sams.append(sam_lines[i]);
		# 		overlap_len = check_overlap(sam_lines[candidate_read], sam_lines[i]);
		# 		if (overlap_len >= max_overlap_len):
		# 			max_overlap_len = overlap_len;
		# 			max_overlap_id = i;
		# 	i += 1;
		# contig_sams.append(sam_lines[candidate_read]);
		# # candidate_read = i;
		# if (max_overlap_id > 0):
		# 	candidate_read = max_overlap_id;
		# else:
		# 	break;

		# i += 1;

		print len(sams_to_process);
		print len(contig_sams);
		print len(unused_sams);


		[new_contig, non_clipped_len] = construct_contig_from_overlapping_sams(ctg_seqs, contig_sams);

		print 'len(new_contig) = %d, non_clipped_len = %d' % (len(new_contig), non_clipped_len);

		if (float(non_clipped_len) < 0.99*float(len(ctg_seqs[0]))):
			# print 'Tu sam!';
			# exit(1);
			continue;

		fp_out_alt_ctg.write('>%s %d\n' % (ctg_headers[0], coverage));
		fp_out_alt_ctg.write('%s\n' % (new_contig));

		for sam_line in contig_sams:
			fp_out.write(sam_line.original_line + '\n');
		fp_out.write('\n');

	fp_out.close();
	fp_out_alt_ctg.close();



if __name__ == "__main__":
	# TEST_SIMULATE();
	# run_poa_sequentially('data/test-msa/ecoli-all.fa');
	contigs_file = 'tests/layout_20151114_221431/contigs_fast.fasta';

	contigs_file = 'tests/miniasm/layout.fasta';

	reads_file = 'tests/sample-dataset/reads-lambda-R73.fasta';
	out_file = 'tests/out/contigs.same.fasta';
	out_sam_file = 'tests/out/contigs.same.sam';
	out_alt_ctg_file = 'tests/out/contigs.same.fasta';

	extract_same_contigs(contigs_file, reads_file, out_sam_file, out_alt_ctg_file);

	run_poa_sequentially_v2(out_alt_ctg_file);

	# if (len(sys.argv) != 4):
	# 	sys.stderr.write('Uses GATK to construct a mutated reference from the original reference given a VCF file with mutations.\n');
	# 	sys.stderr.write('Usage:\n');
	# 	sys.stderr.write('\t%s <reference.fa> <vcf_file> <out_mutated_reference.fa>\n' % (sys.argv[0]));
	# 	exit(1);

	# make_consensus_reference_from_vcf(sys.argv[1], sys.argv[2], sys.argv[3]);
