#! /usr/bin/python

### Module for constructing alternate contigs from a given aligned SAM file of reads to that specific contig.

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

DEBUG_VERBOSE = True;

def construct_contig_from_sams(contig_raw, contig_sams):
	new_contig = '';
	# new_contig_cigar = [];

	if (DEBUG_VERBOSE == True):	print 'Constructing the alternate contig sequence.';

	prev_end_ref = 0;
	for i in xrange(0, len(contig_sams)):
		sam_line = contig_sams[i];

		start_ref = (prev_end_ref) if ((sam_line.pos - 1) < prev_end_ref) else (sam_line.pos - 1);
		end_ref = (sam_line.pos - 1) + sam_line.CalcReferenceLengthFromCigar() - 1;

		start_seq = sam_line.FindBasePositionOnRead(start_ref);
		end_seq = sam_line.FindBasePositionOnRead(end_ref);

		### If there was a gap between the current and the previous read, fill it with the corresponding part of the given contig_raw.
		### Also, if the first sam_line did not start at position 0, the leading chunk will be filled in.
		if (start_ref > prev_end_ref):
			if (DEBUG_VERBOSE == True):	print '\t[%d] Filling gap: start_ref = %d, end_ref = %d, prev_end_ref = %d, start_seq = %d, end_seq = %d, new_contig_len (before) = %d, new_contig_len (after) = %d' % (i, start_ref, end_ref, prev_end_ref, start_seq, end_seq, len(new_contig), (len(new_contig) + (end_seq - start_seq)));
			new_contig += contig_raw[prev_end_ref:start_ref];

		### Copy the new sequence into the alternate contig.
		if (DEBUG_VERBOSE == True):	print '\t[%d] Adding SAM: start_ref = %d, end_ref = %d, prev_end_ref = %d, start_seq = %d, end_seq = %d, new_contig_len (before) = %d, new_contig_len (after) = %d, len(ref) = %d, len(seq) = %d' % (i, start_ref, end_ref, prev_end_ref, start_seq, end_seq, len(new_contig), (len(new_contig) + (end_seq - start_seq)), len(contig_raw), len(sam_line.seq));
		new_contig += sam_line.seq[start_seq:end_seq];

		### Check if this was the last sam_line. If it was, and it didn't end on the end of the reference, fill in the missing part on the back.
		if ((i + 1) == len(contig_sams)):
			if (end_ref < len(contig_raw)):
				if (DEBUG_VERBOSE == True):	print '\t[%d] Fixing the end: (I) with contig_raw: start_ref = %d, end_ref = %d, prev_end_ref = %d, new_contig_len (before) = %d, new_contig_len (after) = %d' % (i, start_ref, end_ref, prev_end_ref, len(new_contig), (len(new_contig) + (len(contig_raw) - end_ref)));
				new_contig += contig_raw[end_ref:-1];
			else:
				if (DEBUG_VERBOSE == True):	print '\t[%d] Fixing the end: (II) with contig_raw: start_ref = %d, end_ref = %d, prev_end_ref = %d, new_contig_len (before) = %d, new_contig_len (after) = %d' % (i, start_ref, end_ref, prev_end_ref, len(new_contig), len(new_contig) + (len(sam_line.seq) - end_seq));
				new_contig += sam_line.seq[end_seq:-1];

		if (DEBUG_VERBOSE == True):	print '';

		prev_end_ref = end_ref;

	if (DEBUG_VERBOSE == True): print 'len(new_contig) = %d' % (len(new_contig));
	if (DEBUG_VERBOSE == True): print 'Finished constructing the alternate contig sequence.';

	return new_contig;

### Determines the coodinate of given positions on the contig_raw sequence to the corresponding position on the newly-generated alternate sequence.
def subsample_ref_positions(contig_raw, contig_sams, sample_step):
	ref_sample_pos = sample_step;
	new_contig_len = 0;
	samples_on_alt_ctg = [];

	prev_end_ref = 0;
	for i in xrange(0, len(contig_sams)):
		sam_line = contig_sams[i];
		if (DEBUG_VERBOSE == True):	print 'sam_line.pos - 1 = %d' % (sam_line.pos - 1);
		if (DEBUG_VERBOSE == True):	print 'sam_line.CalcReferenceLengthFromCigar() = %d' % (sam_line.CalcReferenceLengthFromCigar());

		start_ref = (prev_end_ref) if ((sam_line.pos - 1) < prev_end_ref) else (sam_line.pos - 1);
		end_ref = (sam_line.pos - 1) + sam_line.CalcReferenceLengthFromCigar() - 1;

		start_seq = sam_line.FindBasePositionOnRead(start_ref);
		end_seq = sam_line.FindBasePositionOnRead(end_ref);

		### If there was a gap between the current and the previous read, fill it with the corresponding part of the given contig_raw.
		### Also, if the first sam_line did not start at position 0, the leading chunk will be filled in.
		if (start_ref > prev_end_ref):
			while (ref_sample_pos < start_ref):
				samples_on_alt_ctg.append([ref_sample_pos, new_contig_len + (ref_sample_pos - prev_end_ref + 1)]); ### ref_sample_pos-prev_end_ref makes sense bekause we dont move down with new_contig_len until after this loop is done.
				# samples_on_alt_ctg.append(new_contig_len + (ref_sample_pos - ???start_ref + 1));
				ref_sample_pos += sample_step;
			if (DEBUG_VERBOSE == True):	print '[%d] Filling gap: start_ref = %d, end_ref = %d, prev_end_ref = %d, start_seq = %d, end_seq = %d, new_contig_len (before) = %d, new_contig_len (after) = %d' % (i, start_ref, end_ref, prev_end_ref, start_seq, end_seq, new_contig_len, (new_contig_len + (end_seq - start_seq)));
			new_contig_len += start_ref - prev_end_ref; # len(contig_raw[prev_end_ref:start_ref]);

		### Copy the new sequence into the alternate contig.
		while (ref_sample_pos < end_ref):
			sample_pos_on_seq = sam_line.FindBasePositionOnRead(ref_sample_pos);
			samples_on_alt_ctg.append([ref_sample_pos, new_contig_len + (sample_pos_on_seq - start_seq + 1)]);
			ref_sample_pos += sample_step;
		if (DEBUG_VERBOSE == True):	print '[%d] Adding SAM: start_ref = %d, end_ref = %d, prev_end_ref = %d, start_seq = %d, end_seq = %d, new_contig_len (before) = %d, new_contig_len (after) = %d, len(ref) = %d, len(seq) = %d' % (i, start_ref, end_ref, prev_end_ref, start_seq, end_seq, new_contig_len, (new_contig_len + (end_seq - start_seq)), len(contig_raw), len(sam_line.seq));
		new_contig_len += end_seq - start_seq; # len(sam_line.seq[start_seq:end_seq]);

		### Check if this was the last sam_line. If it was, and it didn't end on the end of the reference, fill in the missing part on the back.
		if ((i + 1) == len(contig_sams)):
#			new_contig += contig_raw[end_ref:-1];
			while (ref_sample_pos < len(contig_raw)):
				samples_on_alt_ctg.append([ref_sample_pos, new_contig_len + (ref_sample_pos - end_ref + 1)]);
				ref_sample_pos += sample_step;
			if (DEBUG_VERBOSE == True):	print '[%d] Fixing the end: start_ref = %d, prev_end_ref = %d, new_contig_len (before) = %d, new_contig_len (after) = %d' % (i, start_ref, prev_end_ref, new_contig_len, (new_contig_len + (start_ref - prev_end_ref)));
			new_contig_len += len(contig_raw) - end_ref; # len(contig_raw[end_ref:-1]);

		if (DEBUG_VERBOSE == True):	print '';
		prev_end_ref = end_ref;

	if (DEBUG_VERBOSE == True):	print 'new_contig_len = %d' % (new_contig_len);

	samples_on_alt_ctg_dict = {};
	for sample in samples_on_alt_ctg:
		samples_on_alt_ctg_dict[sample[0]] = sample[1];

	# return samples_on_alt_ctg;
	return samples_on_alt_ctg_dict;

def extract_alternate_contigs(contig_sams, coverage_threshold, raw_ctg_len, percent_overlap):
	# sorted_lines_by_quality = sorted(sam_lines, reverse=True, key=lambda sam_line: sam_line.chosen_quality);

	### Find alternate contigs from alignments.
	alt_contigs = [];
	sams_to_process = contig_sams;
	coverage = 0;
	while (coverage < 100 and len(sams_to_process) > 0):
		coverage += 1;
		contig_sams = [];
		unused_sams = [];

		### Get alternate contigs from non-overlapping SAM alignments. Parts in-between will be filled with the original raw contig sequence.
		previous_candidate = 0;
		sam_ref_len = sams_to_process[0].CalcReferenceLengthFromCigar();
		new_ctg_bases = sam_ref_len;	### Number of bases on the reference (original raw contig) covered by new alignments.
		old_ctg_bases = sams_to_process[0].pos - 1;					### Number of bases in between alignments (or, those comming from the original raw contig).
		previous_start = sams_to_process[0].pos - 1;
		previous_end = previous_start + sam_ref_len;
		prev_ref_len = sam_ref_len;
		contig_sams.append(sams_to_process[0]);
		for i in xrange(1, len(sams_to_process)):
			candidate_sam = sams_to_process[i];
			sam_ref_len = sams_to_process[i].CalcReferenceLengthFromCigar();
			if ((candidate_sam.pos - 1) >= (previous_end - prev_ref_len * percent_overlap) and ((candidate_sam.pos - 1) + sam_ref_len) > previous_end):
				# new_ctg_bases += sam_ref_len;
				new_ctg_bases += (sam_ref_len) if ((candidate_sam.pos - 1) >= previous_end) else (sam_ref_len - (previous_end - candidate_sam.pos + 1));
				# old_ctg_bases += (sams_to_process[i].pos - 1 - previous_end);
				old_ctg_bases += (sams_to_process[i].pos - 1 - previous_end) if ((candidate_sam.pos - 1) >= previous_end) else 0; # (previous_end - (sams_to_process[i].pos - 1));

				unused_sams += sams_to_process[(previous_candidate + 1):i];
				previous_candidate = i;
				previous_start = sams_to_process[i].pos - 1;
				previous_end = previous_start + sam_ref_len;
				prev_ref_len = sam_ref_len;

				contig_sams.append(sams_to_process[i]);

		unused_sams += sams_to_process[(previous_candidate + 1):-1];
		old_ctg_bases += raw_ctg_len - previous_end;

		### Check the threshold for the covered bases. If above threshold, accept the contig.
		if ((float(new_ctg_bases) / float(new_ctg_bases + old_ctg_bases)) < coverage_threshold):
			### Just verbose.
			if (DEBUG_VERBOSE == True):
				sys.stderr.write('Failed alternate contig #%d, coverage = %d, coverage = %f, num_alignments = %d, new_ctg_bases = %d, old_ctg_bases = %d\n' % (len(alt_contigs), coverage, (float(new_ctg_bases) / float(new_ctg_bases + old_ctg_bases)), len(contig_sams), new_ctg_bases, old_ctg_bases));
				for i in xrange(0, len(contig_sams)):
					sys.stderr.write('  -- [%d] start = %d, end = %d, len(seq) = %d, len_on_ref = %d\n' % (i, contig_sams[i].pos - 1, contig_sams[i].pos - 1 + contig_sams[i].CalcReferenceLengthFromCigar(), len(contig_sams[i].seq), contig_sams[i].CalcReferenceLengthFromCigar()));
			unused_sams = sams_to_process[1:-1];
		else:
			alt_contigs.append(contig_sams);
			### Just verbose.
			if (DEBUG_VERBOSE == True):
				sys.stderr.write('Generated alternate contig #%d, coverage = %d, coverage = %f, num_alignments = %d, new_ctg_bases = %d, old_ctg_bases = %d\n' % (len(alt_contigs), coverage, (float(new_ctg_bases) / float(new_ctg_bases + old_ctg_bases)), len(contig_sams), new_ctg_bases, old_ctg_bases));
				for i in xrange(0, len(contig_sams)):
					sys.stderr.write('  [%d] start = %d, end = %d, len(seq) = %d, len_on_ref = %d, qname = "%s"\n' % (i, contig_sams[i].pos - 1, contig_sams[i].pos - 1 + contig_sams[i].CalcReferenceLengthFromCigar(), len(contig_sams[i].seq), contig_sams[i].CalcReferenceLengthFromCigar(), contig_sams[i].qname));

		sams_to_process = unused_sams;

	return alt_contigs;

### The alignments_path should point to a sorted SAM file.
def extract_alns_for_contig(contig_name, alignments_path, match_threshold, error_threshold):
	### Extract alignments which correspond only to the given contig name.
	contig_name_trimmed = contig_name.split()[0];
	[headers, all_sam_lines] = utility_sam.LoadSAM(alignments_path);
	all_contig_sams = [];
	for sam_line in all_sam_lines:
		if (sam_line.IsMapped() == False):
			continue;
		if (sam_line.rname.split()[0] == contig_name_trimmed):
			### Count CIGAR operations to determine the error rate of the alignments.
			cigop_counts = sam_line.CountAlignmentOps();
			seq_len = len(sam_line.seq) - sam_line.clip_count_front - sam_line.clip_count_back;

			### Check if the CIGAR string is actually in the extended format.
			if ('M' in cigop_counts):
				sys.stderr.write('Warning: alignment does not contain the *extended* CIGAR format! Skipping alignment.\n');
				exit(1);
			else:
				matches = cigop_counts['='];
				errors = cigop_counts['X'] + cigop_counts['D'] + cigop_counts['I'];

			if ((float(matches) / float(seq_len) >= match_threshold) and
				(float(errors) / float(seq_len)) < error_threshold):
				all_contig_sams.append(sam_line);

	return all_contig_sams;



if __name__ == "__main__":
	raw_contig = '/home/isovic/work/eclipse-workspace/git/ra-consensus/tests/miniasm/layout.fasta';
	sorted_alns = '/home/isovic/work/eclipse-workspace/git/ra-consensus/tests/out/temp/tmp.allreads.sorted.sam';

	[ctg_headers, ctg_seqs, ctg_quals] = fastqparser.read_fastq(raw_contig);	
	contig_name = ctg_headers[0];
	sys.stderr.write('Testing for contig_name = "%s", contig_length = %d\n' % (contig_name, len(ctg_seqs[0])));

	all_contig_sams = extract_alns_for_contig(contig_name, sorted_alns, 0.70, 0.40);
	sys.stderr.write('Total %d alignments for contig "%s".\n' % (len(all_contig_sams), contig_name));
	alt_contigs = extract_alternate_contigs(all_contig_sams, 0.80, len(ctg_seqs[0]), 0.01);

	fp_ctg = open('/home/isovic/work/eclipse-workspace/git/ra-consensus/tests/out/temp/alt_ctgs.fasta', 'w');
	for i in xrange(0, len(alt_contigs)):
	# for i in xrange(5, 6):
		print 'Alternate contig %d subsampled starting positions:' % (i);
		alt_contig_seq = construct_contig_from_sams(ctg_seqs[0], alt_contigs[i]);
		ctg_sample_pos = subsample_ref_positions(ctg_seqs[0], alt_contigs[i], 5000);
		fp_ctg.write('>Alternate contig %d, len: %d\n%s\n' % (i, len(alt_contig_seq), alt_contig_seq));
		for sample_pos_key in sorted(ctg_sample_pos.keys()):
			print '[%d, %d]' % (sample_pos_key, ctg_sample_pos[sample_pos_key]);
		print '[%d, %d] (end)' % (len(ctg_seqs[0]), len(alt_contig_seq));
		print '';
	fp_ctg.close();
