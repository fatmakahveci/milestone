#!/usr/bin/env python3

'''
author: @rafael and @fatmakhv
date: April 30s, 2021
aim: Compare sample's variants to assign allele ID using CDS alleles of reference genome
'''


import argparse, pysam, re

from Bio import SeqIO
from collections import Counter
from copy import deepcopy
from itertools import groupby

class VCF:

	def __init__(self, vcf_line):

		fields = vcf_line.split('\t')

		self.cds = get_cds_name_from_allele_name(fields[0])
		self.pos = fields[1]
		self.id = fields[2]
		self.ref = fields[3]

		alt_list = fields[4].split(',')

		if len(alt_list) > 1:
			self.alt = alt_list

		else:
			self.alt = fields[4]

		self.qual = fields[5]
		self.filter = fields[6]
		self.info = fields[7]
		self.format = fields[8]
		self.sample = fields[9]

	def __repr__(self):
		return f"cds: {self.cds}\tpos: {self.pos}\tid: {self.id}\tref: {self.ref}\talt: {self.alt}\tqual: {self.qual}\tfilter: {self.filter}\tinfo: {self.info}\tformat: {self.format}\tsample: {self.sample}"


class REFERENCE_INFO:

	def __init__(self, info_line):

		fields = info_line.strip('\n').split(' ')

		self.cds = get_cds_name_from_allele_name(fields[0])
		self.allele_id = get_allele_id_from_allele_name(fields[0])

		pos_list, alt_list, qual_list = [], [], []
		for variation in fields[1].split(','):
			pos_list.append(variation[:variation.index('*')])
			alt_list.append(variation[variation.index('>')+1:variation.index('-')])
			qual_list.append(variation[variation.index('-')+1:])

		self.pos_list = pos_list
		self.alt_list = alt_list
		self.qual_list = qual_list

	def __repr__(self):
		return f'cds: {self.cds} allele_id: {self.allele_id} positions: {self.pos_list} quals: {self.qual_list}'


class SAMPLE_INFO:

	def __init__(self, vcf_line):
	
		self.cds = vcf_line.cds
		self.pos_list = [vcf_line.pos]

		alt_list = []
		if len(vcf_line.alt) > 1:
			alt_list.append(vcf_line.alt)
		else:
			alt_list = vcf_line.alt.split(',')

		self.alt_list = alt_list
		self.qual_list = vcf_line.qual.split(',')

	def __repr__(self):
		return f'cds: {self.cds} positions: {self.pos_list} alts: {self.alt_list} quals: {self.qual_list}'


class SAM:

	def __init__(self, sam_line):

		fields = sam_line.strip('\n').split('\t')

		self.query_name = fields[0]
		self.read_info = fields[1]
		self.target_sequence_name = fields[2]
		self.position = int(fields[3])
		self.mapping_quality = fields[4]
		self.cigar = fields[5]
		self.rnext = fields[6]
		self.pnext = fields[7]
		self.tlen = fields[8]
		self.sequence = fields[9]
		self.read_quality = fields[10]

	def __repr__(self):
		return f'QNAME: {self.query_name} FLAG: {self.read_info} RNAME: {self.target_sequence_name} POS: {self.position} MAPQ: {self.mapping_quality} CIGAR: {self.cigar} RNEXT: {self.rnext} PNEXT: {self.pnext} TLEN: {self.tlen} SEQ: {self.sequence} QUAL: {self.read_quality}'


class PAF:

	def __init__(self, paf_line):

		self.query_name = paf_line[0]
		self.query_length = int(paf_line[1])
		self.query_start = int(paf_line[2])
		self.query_end = int(paf_line[3])
		self.strand = str(paf_line[4])
		self.target_sequence_name = str(paf_line[5])
		self.target_sequence_length = int(paf_line[6])
		self.target_start_on_original_strand = int(paf_line[7]) # 0-based
		self.target_end_on_original_strand = int(paf_line[8]) # 0-based
		self.number_of_sequence_matches = int(paf_line[9]) # required
		self.alignment_block_length = int(paf_line[10]) # required
		self.mapping_quality = int(paf_line[11]) # 0-255; 255 for missing
		self.identity = float(paf_line[12])
		self.cigar = str(paf_line[13])
		self.sequence = str(paf_line[14])
		self.cs = str(paf_line[15])
		self.matching = paf_line[16]
		self.single_position_coverage = paf_line[17]

	def __repr__(self):
		return f'QNAME: {self.query_name} QUERY LENGTH: {self.query_length} QUERY START: {self.query_start} QUERY END: {self.query_end} STRAND: {self.strand} TARGET SEQUENCE NAME: {self.target_sequence_name} TARGET SEQUENCE LENGTH: {self.target_sequence_length} TARGET START ON ORIGINAL STRAND: {self.target_start_on_original_strand} TARGET END ON ORIGINAL STRAND: {self.target_end_on_original_strand} NUMBER OF SEQUENCE MATCHES: {self.number_of_sequence_matches} ALIGNMENT BLOCK LENGTH: {self.alignment_block_length} MAPPING QUALITY: {self.mapping_quality} IDENTITY: {self.identity} CIGAR: {self.cigar} SEQUENCE: {self.sequence} CS: {self.cs} MATCHING: {self.matching}'


def get_allele_id_from_allele_name(allele_name: str) -> str:

	return allele_name.strip("\n").split("_")[-1]


def read_vcf_file(vcf_file: str) -> dict:

	sample_variant_dict = {}

	with open(vcf_file, 'r') as file:

		for line in file.readlines():
			
			if not line.startswith('#'):

				vcf_line = VCF(line)

				sample_info = SAMPLE_INFO(vcf_line)

				if not vcf_line.cds in sample_variant_dict.keys():

					sample_variant_dict[vcf_line.cds] = sample_info

				else:

					sample_variant_dict[vcf_line.cds].pos_list.append(vcf_line.pos)
					sample_variant_dict[vcf_line.cds].alt_list.append(vcf_line.alt)
					sample_variant_dict[vcf_line.cds].qual_list.append(vcf_line.qual)

		file.close()

	return sample_variant_dict


def get_cds_name_from_allele_name(allele_name: str) -> str:

	return allele_name.strip("\n").split("_")[0]


def read_info_txt_file(info_txt_file: str) -> dict:

	info_line_dict = {}

	with open(info_txt_file, 'r') as file:

		for line in file.readlines():

			reference_info = REFERENCE_INFO(line)
			
			info_line_dict[reference_info.cds] = reference_info

		file.close()

	return info_line_dict


def assign_allele_id(sample_variant_dict, reference_variant_dict):

	for sample_cds, sample_info in sample_variant_dict.items():

		if sample_cds in reference_variant_dict.keys():

			is_known_allele = True

			for reference_cds, reference_info in reference_variant_dict.items():

				if reference_info.pos_list == sample_info.pos_list:
					print('a')

			print(sample_variant_dict[sample_cds], reference_variant_dict[sample_cds])

		else:
			print(len(sample_variant_dict[sample_cds].pos_list), sample_variant_dict[sample_cds])


def merge_intervals(start_end_coverage_list: list) -> list:
	"""
	Merges overlapping intervals.
	
	Parameters
	----------
	start_end_coverage_list : [[<start>, <end>, {base_at_start: cov, ..., base_at_end: cov}], [...], ...]

	Returns
	-------
	merged : [ merge of overlapping intervals of which coverage data is updated and incremented for positions in common. ]
	"""

	merged = [ deepcopy(start_end_coverage_list[0]) ] # merged = [ first interval ]
	
	for current in start_end_coverage_list[1:]: # intervals except the first interval
	
		previous = merged[-1]
	
		# current and previous intervals intersect
		if current[0] <= previous[1]:
	
			# determine top position
			previous[1] = max( previous[1], current[1] )
	
			# merge coverage dictionaries
			previous_cov = previous[2]
			current_cov = current[2]
	
			for k, v in current_cov.items():
	
				if k not in previous_cov:
					previous_cov[k] = v
	
				else:
					previous_cov[k] += v
	
			previous[2] = previous_cov
	
		# current and previous intervals do not intersect
		else:
			merged.append( deepcopy(current) )
	
	return merged


def select_variants( variant_dict: dict, interval_dict: dict ) -> list[ dict, dict ]:
	"""
	Identifies variants with high frequency.
	
	Parameters
	----------
	variant_dict : {<cds> : [ identified variants, variant]} 
	interval_dict : {<cds> : [[<start>, <end>, {base_at_start: cov, ..., base_at_end: cov}], [...], ...]}

	Returns
	-------
	[
	selected_variant_dict: {<cds> : a dictionary with position:variant tuples as keys and with a position_coverage:variant_coverage tuple as value. Stores variants with frequency equal or greater than value passed to minimum_frequency.
	high_frequency : Same structure as previous returned variable but only stores variants with frequency value equal or greater than value passed to variant_frequency.
	]
	"""

	# select indels or SNPs based on frequency
	high_frequency = {}
	selected_variant_dict = {}
	
	# for each locus and detected indels/SNPs
	for loci_id, v in variant_dict.items():
	
		counts = v[1]
		selected = {}
		probable = {}
	
		# for each detected indel/SNP
		for p, c in counts.items():
	
			# get coverage values for locus positions
			coverage_values = interval_dict[loci_id]
	
			# determine which region the indel/SNP is in
			for region in coverage_values:
	
				if p[0] >= region[0] and p[0] <= region[1]:
	
					pos_cov = region[2][p[0]]
	
					# store indel/SNP info if coverage value for the position
					# is 0 or the indel/SNP has high frequency
					frequency = ( c / ( c + pos_cov ) )
	
					if pos_cov == 0 or frequency >= args.minimum_frequency:
	
						# key has position and variant (ref to reads)
						# value has position coverage and indel/SNP frequency
						selected[p] = ( pos_cov, c )
	
					if pos_cov == 0 or frequency >= args.variant_frequency:
	
						probable[p] = ( pos_cov, c )
	
		selected_variant_dict[loci_id] = selected
		high_frequency[loci_id] = probable
	
	return [ selected_variant_dict, high_frequency ]


def determine_breadth_coverage( cds: str, start_end_coverage_list: list ) -> list[ float, int ]:
	"""
	Determines the percentage and total number of covered bases according to a set of coverage intervals.

	Parameters
	----------
	start_end_coverage_list : [[<start>, <end>, {base_at_start: cov, ..., base_at_end: cov}], [...], ...]
	cds_length : Total number of bases in CDS.
	
	Returns
	-------
	[
	breadth_coverage : Percentage of covered bases.
	covered_bases : Total number of covered bases.
	]
	"""

	covered_bases = 0

	for start_end_coverage in start_end_coverage_list:

		coverage_dict = start_end_coverage[2]
		covered_bases += sum([ 1 for position, coverage in coverage_dict.items() if coverage > 0 ])
	
	breadth_coverage = covered_bases / cds_length_dict[cds]
	
	return [ breadth_coverage, covered_bases ]


def determine_missing_intervals( cds: str, start_end_coverage_list: list ) -> list[ list, int ]:
	"""
	Determines sequence intervals that are not covered by any reads.
	
	Parameters
	----------
	cds : cds name
	start_end_coverage_list : [[<start>, <end>, {base_at_start: cov, ..., base_at_end: cov}], [...], ...]
	
	Returns
	-------
	[
	missing_region_list : [ [start and stop positions for a sequence interval that is not covered by reads], [...],...]
	not_covered_length : total number of bases not covered by reads.
	]
	"""

	start, not_covered_length, missing_region_list = 0, 0, []

	for start_end_coverage in start_end_coverage_list:
	
		difference = start_end_coverage[0] - start
	
		if difference > 0:
			missing_region_list.append([ start, start + difference ])
			not_covered_length += difference
			start += difference
	
		# create groups of equal values
		values_groups = [ list(j) for i, j in groupby( start_end_coverage[2].values() ) ]
	
		for g in values_groups:

			g_length = len(g)

			if g[0] == 0:

				missing_region_list.append([ start, start + g_length ])
				not_covered_length += g_length
				start += g_length
	
			else:
				start += g_length

	cds_length = cds_length_dict[cds]

	# add terminal region
	if start != cds_length:
	
		missing_region_list.append([ start, cds_length ])
		not_covered_length += cds_length - start
	
	return [ missing_region_list, not_covered_length ]


def determine_depth_coverage( cds: str, start_end_coverage_list: list ) -> dict:
	"""
	Determine depth of coverage for a sequence.

	Parameters
	----------
	start_end_coverage_list : [ [<start>, <end>, {base_at_start: cov, ..., base_at_end: cov}], [...], ... ]
	cds_length : total length of the sequence.
	
	Returns
	-------
	count_dict : { <coverage value> : total number of positions with that coverage value> }
	"""
	
	# create dictionary to add coverage for all positions
	position_depth_dict = { position : 0 for position in list(range(cds_length_dict[cds])) }

	for start_end_coverage in start_end_coverage_list:

		# increment coverage values based on intervals
		for position, coverage in start_end_coverage[2].items(): # start_end_coverage[2] := coverage_dict
			position_depth_dict[position] += coverage

	# determine coverage distribution
	count_dict = sorted( Counter( position_depth_dict.values() ).most_common(), key=lambda x:x[0] )

	count_dict = { count[0]: count[1] for count in count_dict }

	return count_dict


def process_paf() -> list[dict, list, dict, dict]:
	"""
	Processes results in a PAF file in order to determine coverage
	statistics and identify variants.

	Returns
	-------
	[
	merged_interval_dict : {<cds> : [[<start>, <end>, {base_at_start: cov, ..., base_at_end: cov}], [...], ...]}
	 : {<cds> : a dictionary with position:variant tuples as keys and with a position_coverage:variant_coverage tuple as value}
	 : {<cds> : a list as value. Each list has the following elements: a list with the breadth of coverage value and the number of covered bases. A list with a dictionary with cds as key and a list with missing positions as value and the number of missing positions. Mean depth of coverage. Number of positions with 0 coverage. Number of positions below minimum coverage threshold.}
	]
	"""

	# get indels and SNPs for each locus
	loci_miss = {}
	for paf in paf_list:
		loci_miss.setdefault( paf.target_sequence_name, [] ).extend( paf.single_position_coverage[1] )

	# indels and SNPs counts
	loci_miss = { cds: [coverage, Counter(coverage)] for cds, coverage in loci_miss.items() }

	# identify subsequences that are well covered by reads
	covered_interval_dict = {}
	for paf in paf_list:
		covered_interval_dict.setdefault( paf.target_sequence_name, [] ).append( [ paf.target_start_on_original_strand, paf.target_end_on_original_strand, paf.single_position_coverage[0] ] )

	# sort covered intervals
	covered_interval_dict_sorted = { cds: sorted(start_end_coverage, key=lambda start: start[0]) for cds, start_end_coverage in covered_interval_dict.items() }

	# merge overlapping intervals
	# deepcopy to avoid altering original intervals
	merged_interval_dict = { cds: merge_intervals(start_end_coverage) for cds, start_end_coverage in covered_interval_dict_sorted.items() }

	# identify variants and highly probable variants
	selected_miss_dict, probable_variant_dict = select_variants( loci_miss, merged_interval_dict )

	coverage_stats = { cds: [] for cds in merged_interval_dict }

	for cds, start_end_coverage in merged_interval_dict.items():

		# breadth of coverage and number of covered bases
		coverage = determine_breadth_coverage( cds, start_end_coverage )

		# determine subsequences that are not covered
		missing = determine_missing_intervals( cds, start_end_coverage )
		
		# determine mean depth of coverage
		depth_count_dict = determine_depth_coverage( cds, start_end_coverage )

		# get number of positions with 0 coverage
		zero_coverage = depth_count_dict.get(0, 0)

		# get number of positions below low coverage threshold
		below_coverage = sum([ v for k, v in depth_count_dict.items() if k < args.low_coverage ])

		# get mean depth of coverage
		depth_sum = sum([ k*v for k, v in depth_count_dict.items() ])
		mean_depth = round( depth_sum / cds_length_dict[cds], 4 )

		coverage_stats[cds].extend( [ coverage, missing, mean_depth, zero_coverage, below_coverage ] )

	return [ merged_interval_dict, selected_miss_dict, coverage_stats, probable_variant_dict ]


def get_single_position_coverage(coverage_info: list, start: int) -> list[dict, list]:
	"""
	Determine if positions in a subsequence are covered based on information in the cs field of PAF.

	Parameters
	----------
	coverage_info : [subsequent operations extracted from PAF's cs field]
	start : Subsequence start position in the complete sequence.

	Returns
	-------
	coverage : {<sequence position> : coverage for each position}
	mismatches : [mismatch regions]
	"""

	coverage = {}
	mismatches = []

	for m in coverage_info:

		# subsequence part with exact matches
		if m[0] == ':':

			# create dictionary entries with coverage = 1
			new_cov = { i: 1 for i in range( start, start + int(m[1:]) ) }
			coverage = { **coverage, **new_cov }

			# increment start position
			start = start + int(m[1:])

		# position with substitution
		elif m[0] == '*':

			coverage[start] = 0
			mismatches.extend( [ (start, m) ] )
			start += 1

		# position with deletion
		elif m[0] == '-':

			# coverage 0 for missing bases
			new_cov = { i: 0 for i in range( start, start + len(m[1:]) ) }
			coverage = { **coverage, **new_cov }
			mismatches.extend([ (start, m) ])
			start = start + len(m[1:])

		# insertion
		elif m[0] == '+':

			# do not add coverage values for positions because
			# insertion does not exist in reference
			mismatches.extend([ (start, m) ])

	return [ coverage, mismatches ]


def convert_sam_into_paf(sam_file_name: str) -> list[PAF]:
	"""
	Convert SAM-formatted file into PAF list.

	Parameters
	----------
	sam_file_name : Name of the SAM file to be converted.

	Returns
	-------
	paf_list : List of PAF items for the further analysis. 
	"""

	paf_list = []

	with open(sam_file_name, "r") as sam_file:

		cds_length_dict = {}

		for sam_idx, sam_line in enumerate(sam_file.readlines()):

			# skip header lines except the line containing CDS length information
			if sam_line.startswith("@"): # header starts with @
				
				# ['@SQ', 'SN:<cds_name>', 'LN:<cds_length>']
				if sam_line.startswith("@SQ"):

					sq_fields = sam_line.strip('\n').split('\t')
					cds_name, cds_length = sq_fields[1][3:], int(sq_fields[2][3:])
					cds_length_dict[cds_name] = cds_length
					
				continue
			
			sam = SAM(sam_line)
			
			if (sam.sequence != "*" and sam.read_quality != "*" and len(sam.sequence) != len(sam.read_quality)):
				raise ValueError("ERROR at line " + str(sam_idx) + ":inconsistent SEQ and QUAL lengths - " + str(len(sam.sequence)) + " != " + str(len(sam.read_quality)))
			
			# read_info == 4 := unmapped read; read_info == 0x100 := more than one location for mapped read
			if sam.target_sequence_name == '*' or sam.read_info == '4' or sam.read_info == '0x100':
				continue

			# Get the reference length for this alignment
			if sam.target_sequence_name in cds_length_dict:
				target_sequence_length = cds_length_dict[sam.target_sequence_name]

			# nn := no_ambiguous_bases_in_aln
			match = re.findall(r"\tnn:i:(\d+)", sam_line)
			no_ambiguous_bases_in_aln = [ int(match[0]) if match else 0 ][0]

			# NM := Total number of mismatches and gaps in the alignment
			match = re.findall(r"\tNM:i:(\d+)", sam_line)
			mismatch_and_gaps = [ int(match[0])if match else 0 ][0]  + no_ambiguous_bases_in_aln
			have_mismatch_and_gaps = [ True if match else False ][0]

			# See sequana.cigar for more information
			clip = [0, 0]
			I = [0, 0]  # Insertion
			D = [0, 0]  # Deletion
			M, N = 0, 0 # Matches
			query_length, target_length, mm = 0, 0, 0
			ext_cigar = False
			n_cigar = 0

			cigar, cs = "", ""
			for no_letter, letter in re.findall(r"(\d+)([MIDSHNX=])", sam.cigar):

				no_letter = int(no_letter)

				# for alignment MATCH
				if (letter == 'M'):
					M += no_letter
					query_length += no_letter
					target_length += no_letter
					ext_cigar = False
					cigar += f"{no_letter}M"
					cs += f":{no_letter}"

				# for Insertion to the reference
				elif (letter == 'I'):
					I[0] += 1
					I[1] += no_letter
					query_length += no_letter
					cigar += f"{no_letter}I"
					cs += f"+{no_letter}"

				# for deletion from the reference
				elif (letter == 'D'):
					D[0] += 1
					D[1] += no_letter
					target_length += no_letter
					cigar += f"{no_letter}D"
					cs += f"-{no_letter}"

				# for skipped region from the reference
				elif (letter == 'N'):
					N += no_letter
					target_length += no_letter

				# clipped sequence present in seq
				elif (letter == 'S'):
					clip[ 0 if M == 0 else 1 ] = no_letter
					query_length += no_letter

				# clipped sequence NOT present in seq
				elif (letter == 'H'):
					clip[0 if M == 0 else 1] = no_letter

				# for equal
				elif (letter == '='):
					M += no_letter
					query_length += no_letter
					target_length += no_letter
					ext_cigar = True

				# for diff (sequence mismatched)
				elif (letter == 'X'):
					M += no_letter
					query_length += no_letter
					target_length += no_letter
					mm += no_letter
					ext_cigar = True
					cs += f"*{no_letter}"

				n_cigar += 1

			# CIGAR operations or alignment end position larger than ref length or SEQ length inconsistent with CIGAR(" + str(len(sequence)) + " != " + str(query_length) + ")
			if (n_cigar > 65535) or (target_length + sam.position - 1 > target_sequence_length) or (sam.sequence != '*' and len(sam.sequence) != query_length):
				continue

			if (have_mismatch_and_gaps is False or ext_cigar):
				mismatch_and_gaps = I[1] + D[1] + mm

			# mismatch_and_gaps is less than the total number of gaps
			if (mismatch_and_gaps < I[1] + D[1] + mm):
				mismatch_and_gaps = I[1] + D[1] + mm

			query_length = M + I[1] + clip[0] + clip[1]

			if sam.read_info == '16': # if read is from reverse strand
				query_start = clip[1]
				query_end = query_length - clip[0]
				strand = '-'

			else:  # if read is not from reverse strand
				query_start = clip[0]
				query_end = query_length - clip[1]
				strand = '+'

			target_start = sam.position - 1
			target_end = target_start + M + D[1] + N

			number_of_sequence_matches = ( M - ( mismatch_and_gaps - I[1] - D[1] ) ) + no_ambiguous_bases_in_aln # match + no_ambiguous_bases_in_aln
			alignment_block_length = ( M + I[1] + D[1] ) - no_ambiguous_bases_in_aln # base length - no_ambiguous_bases_in_aln

			# to filter out alignments below defined identity
			identity = float(number_of_sequence_matches / alignment_block_length)

			# match alignment string with regex
			matching = re.findall(r':[0-9]+|\*[a-z][a-z]|\+[a-z]+|-[a-z]+', cs)

			# get information about positions that match to determine coverage
			single_position_coverage = get_single_position_coverage(matching, target_start)

			paf = PAF([sam.query_name, query_length, query_start, query_end, strand, sam.target_sequence_name, target_sequence_length, target_start, target_end, number_of_sequence_matches, alignment_block_length, sam.mapping_quality, identity, sam.cigar, sam.sequence, cs, matching, single_position_coverage])

			# Quality check for mapping and identity
			if (args.identity <= paf.identity) or (args.mapping_quality <= paf.mapping_quality):
				paf_list.append(paf)

	return [ cds_length_dict, paf_list ]


if __name__ == "__main__":

	parser = argparse.ArgumentParser()

	parser.add_argument('--vcf', type=str, required=True,
						help='Sample\'s VCF file to be compared.')

	parser.add_argument('--info_txt', type=str, required=True,
						help='<reference>_info.txt file to be compared.')
	
	parser.add_argument('--identity', type=float, required=False, default=0.95,
						help='Minimum identity percentage for mapped reads.')

	parser.add_argument('--mapping_quality', type=int, required=False, default=30,
						help='Minimum quality for read mapping.')

	parser.add_argument('--minimum_frequency', type=float, required=False, default=0.10,
						help='Minimum relative frequency of indels or SNPs that are reported.')

	parser.add_argument('--low_coverage', type=int, required=False, default=10,
						help='Minimum coverage value. Positions with smaller values are reported as low coverage positions.')

	parser.add_argument('--variant_frequency', type=float, required=False, default=0.30,
						help='Indels or SNPs with a relative frequency that is equal or greater that this value are reported as highly probable.')

	# @todo add to workflow convert <sample.bam> into <sample.sam> before running this
	parser.add_argument('--sam', type=str, required=True,
						help='Create SAM file with mapping results.')

	parser.add_argument('--sample_fasta', type=str, required=True,
						help='Path to sample\'s FASTA file containing CDS sequences')

	args = parser.parse_args()
	
	assign_allele_id( sample_variant_dict = read_vcf_file(args.vcf), reference_variant_dict = read_info_txt_file(args.info_txt) )
	
	# cds_length_dict, paf_list = convert_sam_into_paf(args.sam)

	# merged_interval_dict, selected_miss_dict, coverage_stats, probable_variant_dict = process_paf()


