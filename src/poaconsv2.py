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

DEBUG_VERBOSE = False;
DEBUG_VERBOSE = True;

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

def run_poa_sequentially_v2(seq_path, out_consensus_file):
	temp_subseq_file = '%s/tmp.subseq.fasta' % (os.path.dirname(out_consensus_file));
	temp_msa_file = '%s/tmp.subseq.fasta.pir' % (os.path.dirname(out_consensus_file));
	# out_consensus_file = '%s/consensus-poa.fasta' % (os.path.dirname(seq_path));
	out_consensus_file_chunks = '%s/tmp.consensus.chunks.fasta' % (os.path.dirname(out_consensus_file));

	fp_out_all = open(out_consensus_file, 'w');
	fp_out_chunks = open(out_consensus_file_chunks, 'w');

	timestamp = strftime("%Y/%m/%d %H:%M:%S", gmtime());
	fp_out_all.write('>Consensus_with_POA all %s\n' % (timestamp));

	print 'seq_path = "%s"' % (seq_path);

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

		# if (start_coord == 0 or end_coord == max_seq_len):
		# 	execute_command('%s/poaV2/poa -do_progressive -read_fasta %s -pir %s %s/poaV2/blosum80.mat' % (TOOLS_PATH, temp_subseq_file, temp_msa_file, TOOLS_PATH));
			# execute_command('%s/poaV2/poa -do_progressive -read_fasta %s -pir %s %s/poaV2/all1.mat' % (TOOLS_PATH, temp_subseq_file, temp_msa_file, TOOLS_PATH));
		# else:
		execute_command('%s/poaV2/poa -do_global -do_progressive -read_fasta %s -pir %s %s/poaV2/blosum80.mat' % (TOOLS_PATH, temp_subseq_file, temp_msa_file, TOOLS_PATH));
			# execute_command('%s/poaV2/poa -do_global -do_progressive -read_fasta %s -pir %s %s/poaV2/all1.mat' % (TOOLS_PATH, temp_subseq_file, temp_msa_file, TOOLS_PATH));

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



### SAM lines are expected to be sorted. sam_line1 should come before sam_line2 (that is sam_line1.pos should be <= sam_line2.pos).
### min_overlap_percent specifies the minimum percentage of each read that should be covered by the overlap. This should help avoid introduction of long indels.
def check_overlap(sam_line1, sam_line2, min_overlap_percent):
	start1 = sam_line1.pos - 1;
	# end1 = start1 + len(sam_line1.seq) - sam_line1.clip_count_front - sam_line1.clip_count_back;
	end1 = start1 + sam_line1.CalcReferenceLengthFromCigar();

	start2 = sam_line2.pos - 1;
	# end2 = start2 + len(sam_line2.seq) - sam_line2.clip_count_front - sam_line2.clip_count_back;
	end2 = start2 + sam_line2.CalcReferenceLengthFromCigar();

	### In this case, reads do not overlap.
	if (start2 > end1 or end2 < start1):
		return 0;

	### One read is contained in the other.
	if ((start2 >= start1 and end2 <= end1) or (start1 >= start2 and end1 <= end2)):
		return -1;

	# overlap_len = (end2 - start1);
	overlap_len = (end1 - start2);

	### The reads *do* overlap, but the length of the overlap is too short, so dismiss the overlap.
	# if (overlap_len < min_overlap_percent*(end1 - start1) or overlap_len < min_overlap_percent*(end2 - start2)):
	if (overlap_len < min_overlap_percent):

		print 'overlap_len = ', overlap_len;
		print 'start1 = %d, end1 = %d' % (start1, end1);
		print 'start2 = %d, end2 = %d' % (start2, end2);
		print 'min_overlap_percent*(end1 - start1) = ', min_overlap_percent*(end1 - start1);
		print 'min_overlap_percent*(end2 - start2) = ', min_overlap_percent*(end2 - start2);
		print 'len(sam_line1.seq) = %d' % (len(sam_line1.seq));
		print 'len(sam_line2.seq) = %d' % (len(sam_line2.seq));
		print '';

		# print 'Tu sam 1!';
		# exit(1);
		# return -2;
		return 0;

	# return overlap_len;

	joined_length = (end2 - start1);

	return joined_length;

def construct_contig_from_overlapping_sams(ctg_seqs, contig_sams):
	new_contig = '';
	new_contig_cigar = [];

	overhang_before = contig_sams[0].pos - 1;
	if (DEBUG_VERBOSE == True):
		print 'overhang_before = %d' % (contig_sams[0].pos - 1);

	new_contig += ctg_seqs[0][0:overhang_before];
	if (overhang_before > 0):
		new_contig_cigar.append([overhang_before, '=']);

	non_clipped_len = 0;

	i = 0;
	while (i < len(contig_sams)):
		sam_line = contig_sams[i];
		start_pos = sam_line.clip_count_front;
		end_pos = (len(sam_line.seq) - sam_line.clip_count_back) if ((i + 1) == len(contig_sams)) else sam_line.FindBasePositionOnRead(contig_sams[i+1].pos - 1); # (contig_sams[i+1].pos - sam_line.pos);

# Stao sam ovdje, duljina reada procitana iz CIGAR stringa nije jednaka duljini seq polja.
# Potencijalno je to zbog krivih koordinata na CIGAR stringu.
# Trebao bih jos jednu stvar napraviti u utility_sam, a to je da funkcija CalcCigarStartingPositions moze razlomiti sve eventove na pojedinacne baze, a ne samo matcheve.

		if (DEBUG_VERBOSE == True):
			print 'start_pos = %d' % (start_pos);
			print 'end_pos = %d' % (end_pos);
			print 'len(new_contig) = %d' % (len(new_contig));
			print 'contig_sams[i].pos - 1 = %d' % (contig_sams[i].pos - 1);
			print 'contig_sams[i].CalcReferenceLengthFromCigar() = %d' % (contig_sams[i].CalcReferenceLengthFromCigar());
			print '(contig_sams[i+1].pos - 1) = %d' % (-1 if ((i+1) == len(contig_sams)) else (contig_sams[i+1].pos - 1));
			print '(contig_sams[i].pos - 1) = %d' % ((contig_sams[i].pos - 1));
			print 'sam_line.CalcReferenceLengthFromCigar() = %d' % sam_line.CalcReferenceLengthFromCigar();
			print 'contig_sams[i].pos - 1 + sam_line.CalcReferenceLengthFromCigar() = %d' % (contig_sams[i].pos - 1 + sam_line.CalcReferenceLengthFromCigar());
			print 'len(sam_line.seq) = %d' % (len(sam_line.seq));
			print 'sam_line.CalcReadLengthFromCigar() = %d' % (sam_line.CalcReadLengthFromCigar() - sam_line.clip_count_front - sam_line.clip_count_back);

		new_chunk = sam_line.seq[start_pos:end_pos];
		non_clipped_len += len(new_chunk);
		new_contig += new_chunk;
		new_contig_cigar += sam_line.GetCigarBetweenBases(start_pos, end_pos+1);
		# new_contig_cigar.append([1, '\n']);
		i += 1;

		if (DEBUG_VERBOSE == True):
			print 'len(new_contig) = %d' % (len(new_contig));
			print '';

	# overhang_after = contig_sams[-1].pos + len(contig_sams[-1].seq) - contig_sams[-1].clip_count_front - contig_sams[-1].clip_count_back;
	overhang_after = contig_sams[-1].pos - 1 + contig_sams[-1].CalcReferenceLengthFromCigar();
	if (DEBUG_VERBOSE == True):
		print 'len(new_contig) = %d' % (len(new_contig));
		print 'overhang_after = %d' % (overhang_after);
		print 'contig_sams[-1].pos - 1 = %d' % (contig_sams[-1].pos - 1);
		print 'contig_sams[-1].CalcReferenceLengthFromCigar() = %d' % (contig_sams[-1].CalcReferenceLengthFromCigar());
	new_contig += ctg_seqs[0][overhang_after:-1];
	if ((len(ctg_seqs[0]) - overhang_after) > 0):
		new_contig_cigar.append([(len(ctg_seqs[0]) - overhang_after), '=']);

	if (DEBUG_VERBOSE == True):
		print 'len(new_contig) after adding the overhang = %d' % (len(new_contig));
		print 'len(ctg_seqs[0]) = %d' % (len(ctg_seqs[0]));
	
	# for val in new_contig_cigar:
	# 	print val;
	new_contig_cigar_string = ''.join(('%d%s' % (val[0], val[1])) for val in new_contig_cigar);
	# fp1 = open('temp.txt', 'w');
	# fp1.write('%s\n' % new_contig_cigar_string);
	# fp1.close();
	# exit(1);

	return [new_contig, non_clipped_len, new_contig_cigar_string];

### Parameter 'single_contig_file' is the path to a file containing only single contig sequence. If the original contig file was a multifasta, then a single contig
### from that multifasta needs to be extracted to a separate file, which is the file provided through this parameter.
### Parameter 'out_alt_ctg_file' will contain all the alternate contigs for the given input contig.
def extract_alternate_contigs(single_contig_file, reads_file, out_alt_ctg_file, ref_file=''):
	### Generate file paths for some temporary files.
	path_aligns_basename = '%s/tmp.allreads' % (os.path.dirname(out_alt_ctg_file));
	path_aligns = '%s.sam' % (path_aligns_basename);
	path_aligns_sorted_basename = '%s.sorted' % (path_aligns_basename);
	path_aligns_sorted_sam = '%s.sam' % (path_aligns_sorted_basename);
	path_alt_contig_sams = '%s.altctgs.sam' % (path_aligns_basename);

	if (not os.path.exists(os.path.dirname(out_alt_ctg_file))):
		os.path.makedirs(os.path.dirname(out_alt_ctg_file));

	### Generate alignments.
	# execute_command('%s/graphmap/bin/Linux-x64/graphmap -a anchor -b 3 -r %s -d %s -o %s' % (TOOLS_PATH, single_contig_file, reads_file, path_aligns));
	# execute_command('samtools view -Sb %s | samtools sort - %s && samtools view -h %s.bam > %s' % (path_aligns, path_aligns_sorted_basename, path_aligns_sorted_basename, path_aligns_sorted_sam));

	[ctg_headers, ctg_seqs, ctg_quals] = fastqparser.read_fastq(single_contig_file);
	[headers, all_sam_lines] = utility_sam.LoadSAM(path_aligns_sorted_sam);

	sys.stderr.write('Number of lines in the original SAM file: %d\n' % (len(all_sam_lines)));

	sam_lines = [];
	for sam_line in all_sam_lines:
		if (sam_line.IsMapped() == False):
			continue;
		seq_len = len(sam_line.seq) - sam_line.clip_count_front - sam_line.clip_count_back;
		cigop_counts = sam_line.CountAlignmentOps();
		### Check if the CIGAR string is actually in the extended format.
		if ('M' in cigop_counts):
			sys.stderr.write('Warning: alignment does not contain the *extended* CIGAR format! Skipping alignment.\n');
			exit(1);
		else:
			matches = cigop_counts['='];
			errors = cigop_counts['X'] + cigop_counts['D'] + cigop_counts['I'];

		if (float(matches) / float(seq_len) >= 0.70 and float(errors) / float(seq_len) < 0.40):
			sam_lines.append(sam_line);

	sys.stderr.write('Number of filtered SAM lines (only mapped and with errors below threshold): %d\n' % (len(sam_lines)));

	fp_out_alt_ctg = open(out_alt_ctg_file, 'w');
	fp_out_alt_ctg_sams = open(path_alt_contig_sams, 'w');
	fp_out_alt_ctg_sams.write('\n'.join(headers) + '\n');

	### Find alternate contigs from alignments.
	sams_to_process = sam_lines;
	coverage = 0;
	while (coverage < 100 and len(sams_to_process) > 0):
		coverage += 1;
		print '---------------------------------------';
		print 'Coverage = %d' % (coverage);
		sys.stderr.write('Number of alignments in pool: %d\n' % (len(sams_to_process)));
		contig_sams = [];
		unused_sams = [];
		i = 0;
		candidate_read = i;
		contig_sams.append(sams_to_process[candidate_read]);
		# for candidate_read in xrange((i+1), len(sams_to_process)):
		start1 = sams_to_process[candidate_read].pos - 1;
		end1 = start1 + sams_to_process[candidate_read].CalcReferenceLengthFromCigar();
		print 'candidate: start = %d, end = %d' % (start1, end1);

		while ((candidate_read + 1) < len(sams_to_process)):
			max_overlap_len = 0;
			max_overlap_id = -1;
			# j = candidate_read + 1;
			# while (j < len(sams_to_process)):
			for j in xrange(candidate_read + 1, len(sams_to_process)):
				overlap_len = check_overlap(sams_to_process[candidate_read], sams_to_process[j], 0);
				if (overlap_len == 0):
					print 'break 1';
					print '  j = %d (in the range of %d to %d)' % (j, candidate_read + 1, len(sams_to_process));
					break;
				elif (overlap_len == -1 or overlap_len == -2):	### -1 is for contained sequences, and -2 is for overlaps which are below the threshold.
					# j += 1;
					continue;

				if (max_overlap_id == -1 or overlap_len >= max_overlap_len):
					max_overlap_len = overlap_len;
					max_overlap_id = j;
				# j += 1;

			if (max_overlap_id > 0):
				print '  starting read = %d' % (candidate_read);
				print '  candidate_read = %d' % (max_overlap_id);
				print '  max_overlap_len = %d' % (max_overlap_len);
				print '  unused overlapping reads: %d - %d' % ((candidate_read + 1), max_overlap_id);

				start1 = sams_to_process[max_overlap_id].pos - 1;
				end1 = start1 + sams_to_process[max_overlap_id].CalcReferenceLengthFromCigar();
				print '  candidate: start = %d, end = %d' % (start1, end1);
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

		print '  after coverage %d:' % (coverage);
		print '    len(sams_to_process) = %d' % (len(sams_to_process));
		print '    len(contig_sams) = %d' % len(contig_sams);
		print '    len(unused_sams) = %d' % len(unused_sams);


		[new_contig, non_clipped_len, new_contig_cigar] = construct_contig_from_overlapping_sams(ctg_seqs, contig_sams);

		test_sam_line = utility_sam.SAMLine();
		test_sam_line.seq = new_contig;
		test_sam_line.cigar = new_contig_cigar;

		print 'test_sam_line.CalcReadLengthFromCigar() = %d' % (test_sam_line.CalcReadLengthFromCigar());
		print 'test_sam_line.CalcReferenceLengthFromCigar() = %d' % (test_sam_line.CalcReferenceLengthFromCigar());
		print 'len(test_sam_line.seq) = %d' % (len(test_sam_line.seq));

		print '*********************    len(new_contig) = %d, non_clipped_len = %d' % (len(new_contig), non_clipped_len);
		exit(1);

		if (float(non_clipped_len) < 0.85*float(len(ctg_seqs[0]))):
			# print 'Tu sam!';
			# exit(1);
			continue;
		else:
			print '++++++++++++++++++++++++++++++++++++++++';


		fp_out_alt_ctg.write('>%s %d\n' % (ctg_headers[0], coverage));
		fp_out_alt_ctg.write('%s\n' % (new_contig));

		for sam_line in contig_sams:
			fp_out_alt_ctg_sams.write(sam_line.original_line + '\n');
		fp_out_alt_ctg_sams.write('\n');

	fp_out_alt_ctg_sams.close();
	fp_out_alt_ctg.close();




def TEST_SIMULATE():
	for i in xrange(0, 10):
		# generate_mutated_reference(('%s/../reference-genomes/escherichia_coli.fa' % SCRIPT_PATH), 0.0006, 0.0067, 'data/test-msa/ecoli-%d.fa' % (i));
		# seqs_for_mafft
		if (i == 0):
			execute_command('cat %s > %s' % ('data/test-msa/ecoli-%d.fa' % (i), 'data/test-msa/ecoli-all.fa'));
		else:
			execute_command('cat %s >> %s' % ('data/test-msa/ecoli-%d.fa' % (i), 'data/test-msa/ecoli-all.fa'));

def TEST_SAM_TO_CONTIG(single_contig_file, contig_sam, output_alt_contig_fasta):
	[ctg_headers, ctg_seqs, ctg_quals] = fastqparser.read_fastq(single_contig_file);
	[headers, contig_sams] = utility_sam.LoadSAM(contig_sam);
	[new_contig, non_clipped_len, new_contig_cigar] = construct_contig_from_overlapping_sams(ctg_seqs, contig_sams);

	fp = open(output_alt_contig_fasta, 'w');
	fp.write('>Alternate contig\n');
	fp.write('%s\n' % new_contig);
	fp.close();

### This function opens the contigs file and separates contigs one-by-one to a temporary file, to which reads will be aligned to.
def process_contigs(contigs_file, reads_file, out_file):
	try:
		fp_in = open(contigs_file, 'r');
	except IOError:
		sys.stderr.write('ERROR: Could not open file "%s" for reading!\n' % contigs_file);
		return;

	single_contig_file = '%s/temp/tmp.singlecontig.raw.fasta' % (os.path.dirname(out_file));
	out_alt_ctg_file = '%s/temp/tmp.singlecontig.alt.fasta' % (os.path.dirname(out_file));
	out_consensus_single_contig = '%s/temp/tmp.singlecontig.cons.fasta' % (os.path.dirname(out_file));

	try:
		fp_out_single_contig = open(single_contig_file, 'w');
		fp_out_single_contig.close();
	except IOError:
		sys.stderr.write('ERROR: Could not open file "%s" for writing!\n' % single_contig_file);
		return;

	try:
		fp_out = open(out_file, 'a');
		fp_out.close();
	except IOError:
		sys.stderr.write('ERROR: Could not open file "%s" for writing!\n' % out_file);
		return;

	while True:
		[header, read] = fastqparser.get_single_read(fp_in);
		if (len(header) == 0):
			break;
		seq = read[1];

		fp_out_single_contig = open(single_contig_file, 'w');
		fp_out_single_contig.write('>%s\n' % header);
		fp_out_single_contig.write('%s\n' % seq);
		fp_out_single_contig.close();
		
		extract_alternate_contigs(single_contig_file, reads_file, out_alt_ctg_file);
		run_poa_sequentially_v2(out_alt_ctg_file, out_consensus_single_contig);

		execute_command('cat %s >> %s' % (out_consensus_single_contig, out_file));

	fp_in.close();

def backup_old_polish_results(polished_fasta, backup_fasta):
	polished_fasta = os.path.abspath(polished_fasta);
	backup_fasta = os.path.abspath(backup_fasta);
	if (os.path.exists(polished_fasta)):
		execute_command('cat %s >> %s' % (polished_fasta, backup_fasta));
		os.remove(polished_fasta);

def RUN_TEST_SAM_TO_CONTIG():
	single_contig_file = '/home/isovic/work/eclipse-workspace/git/ra-consensus/tests/out/temp/tmp.singlecontig.raw.fasta';
	contig_sam = '/home/isovic/work/eclipse-workspace/git/ra-consensus/tests/out/temp/third.sam';
	output_alt_contig_fasta = '/home/isovic/work/eclipse-workspace/git/ra-consensus/tests/out/temp/third.fasta';
	TEST_SAM_TO_CONTIG(single_contig_file, contig_sam, output_alt_contig_fasta);

	exit(1);

def TEST_SAMPLE_DATA():
	# run_poa_sequentially('data/test-msa/ecoli-all.fa');
	# contigs_file = 'tests/layout_20151114_221431/contigs_fast.fasta';
	contigs_file = '/home/isovic/work/eclipse-workspace/git/ra-consensus/tests/miniasm/layout.fasta';
	reads_file = '/home/isovic/work/eclipse-workspace/git/ra-consensus/tests/sample-dataset/reads-lambda-R73.fasta';
	reads_file = '/home/isovic/work/eclipse-workspace/git/ra-consensus/tests/sample-dataset/reads-lambda-R73-without_problematic.fasta';
	out_file = '/home/isovic/work/eclipse-workspace/git/ra-consensus/tests/out/contigs.consensus.fasta';
	# out_sam_file = 'tests/out/contigs.same.sam';

	backup_old_polish_results(out_file, '%s.bak' % (out_file));
	process_contigs(contigs_file, reads_file, out_file);

	execute_command('dnadiff -p tests/out/temp/dnadiff/raw tests/sample-dataset/NC_001416.fa tests/miniasm/layout.fasta');
	execute_command('dnadiff -p tests/out/temp/dnadiff/polished tests/sample-dataset/NC_001416.fa tests/out/contigs.consensus.fasta');

	# make_consensus_reference_from_vcf(sys.argv[1], sys.argv[2], sys.argv[3]);
	exit(1);

def main():
	# TEST_SIMULATE();
	# RUN_TEST_SAM_TO_CONTIG();
	TEST_SAMPLE_DATA();

	# if (len(sys.argv) == 1):
	# 	TEST_SAMPLE_DATA();

	if (len(sys.argv) != 4):
		sys.stderr.write('Proof-of-Concept Consensus for de novo genome assemblies.\n');
		sys.stderr.write('This consensus tool first creates several instances of the contig through alignment (e.g. 10x; these instances have similar error rate to the original data), and then applies POA on chunks (windows) alternate contigs to produce consensus.\n');
		sys.stderr.write('Usage:\n');
		sys.stderr.write('\t%s <reads.fasta> <contigs.fasta> <out_polished.fasta>\n' % (sys.argv[0]));
		exit(1);

	reads_file = os.path.abspath(sys.argv[1]);
	contigs_file = os.path.abspath(sys.argv[2]);
	out_file = os.path.abspath(sys.argv[3]);

	out_folder = os.path.dirname(out_file);
	out_basename = os.path.splitext(os.path.basename(out_file))[0];
	out_sam_file = '%s/tmp-%s/%s.all.sam' % (out_folder, out_basename, out_basename);

if __name__ == "__main__":
	main();
