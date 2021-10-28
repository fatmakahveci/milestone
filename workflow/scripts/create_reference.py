#!/usr/bin/env python3

######################################################
## Author: @fatmakhv @rafaelmamede                  ##
## The latest update: 27/10/2021                    ##
## Aim: Reference FASTA, VCF and INFO file creation ##
######################################################


import argparse, glob, os, shutil, subprocess, sys

from Bio import SeqIO
from collections import OrderedDict
from io import StringIO  # python3
from pathlib import Path


class Info:

    def __init__(self, line):

        pos_list, ref_list, alt_list, qual_list = [], [], [], []

        for variation in line.split(','):

            pos_end_idx = variation.index('*')
            pos = int(variation[ : pos_end_idx ])

            ref_end_idx = variation.index('>')
            ref = variation[ pos_end_idx +1 : ref_end_idx ]

            alt_end_idx = variation.index('-')
            alt = variation[ ref_end_idx + 1 : alt_end_idx ]
                        
            qual = variation[ alt_end_idx + 1 : ] 
            
            pos_list.append(pos)
            ref_list.append(ref)
            alt_list.append(alt)
            qual_list.append(qual)

        self.pos_list = pos_list
        self.ref_list = ref_list
        self.alt_list = alt_list
        self.qual_list = qual_list

    def __repr__(self):

        return f'positions: {" ".join(list(map(str, self.pos_list)))}\n' \
                f'refs: {" ".join(self.ref_list)}\n' \
                f'alts: {" ".join(self.alt_list)}\n' \
                f'quals: {" ".join(list(map(str, self.qual_list)))}\n'


class Paf:

	def __init__( self, paf_line ):

		fields = paf_line.strip('\n').split('\t')

		self.qname = fields[0]
		self.qlen = int(fields[1])
		self.qstart = int(fields[2])
		self.qend = int(fields[3])
		self.strand = fields[4]
		self.tname = fields[5]
		self.tlen = int(fields[6])
		self.tstart = int(fields[7])
		self.tend = int(fields[8])
		self.nmatch = int(fields[9])
		self.alen = int(fields[10])
		self.mapq = int(fields[11])

	def __repr__(self):

		return f"query sequence name: {self.qname}\t" \
			   f"query sequence length: {self.qlen}\t" \
			   f"query start: {self.qstart}\t" \
			   f"query end: {self.qend}\t" \
			   f"strand (+ or -): {self.strand}\t" \
			   f"target sequence name: {self.tname}\t" \
			   f"target sequence length: {self.tlen}\t" \
			   f"target start: {self.tstart}\t" \
			   f"target end: {self.tend}\t" \
			   f"number of matching bases in the mapping: {self.nmatch}\t" \
			   f"number of bases including gaps in the mapping: {self.alen}\t" \
			   f"mapping quality:{self.mapq}\n"


class Vcf:

	def __init__( self, vcf_line ):

		fields = vcf_line.split('\t')

		self.chr = get_cds_name_from_allele_name(fields[0])
		self.pos = fields[1]
		self.id = fields[2]
		self.ref = fields[3]
		self.alt = fields[4]
		self.qual = fields[5]
		self.filter = fields[6]
		self.info = fields[7]
		self.sample_format = fields[8]
		self.sample = fields[9]

	def __repr__(self):

		return f"chr: {self.chr}\t" \
			   f"pos: {self.pos}\t" \
			   f"id: {self.id}\t" \
			   f"ref: {self.ref}\t" \
			   f"alt: {self.alt}\t" \
			   f"qual: {self.qual}\t" \
			   f"filter: {self.filter}\t" \
			   f"info: {self.info}\t" \
			   f"format: {self.sample_format}\t" \
			   f"sample: {self.sample}\n"


def get_cds_name_from_allele_name(allele_name: str) -> str:
	"""
	Get {cds_name}_{allele_id} and return {cds_name}

	Parameters
	----------
	allele_name : {cds_name}_{allele_id}

	Returns
	-------
	cds_name : name of CDS for given allele_name
	"""

	cds_name = allele_name.strip("\n").split("_")[0]

	return cds_name


def get_allele_id_from_allele_name(allele_name: str) -> str:
	"""
	Get {cds_name}_{allele_id} and return {allele_id}

	Parameters
	----------
	allele_name : {cds_name}_{allele_id}

	Returns
	-------
	allele_id : allele ID for given allele_name
	"""

	allele_id = allele_name.strip("\n").split("_")[-1]

	return allele_id


def write_allele_defining_variant_list_to_file( cds_name: str, sv_at_edges: str, allele_id: str, pos_dict: dict ) -> None:
	"""
	Write the variants of alleles that are not equal to the reference CDS

	Parameters
	----------
	cds_name : name of cds of which allele will be written
	allele_id : ID of allele of which variant will be written
	pos_dict : positions of variants for the allele
	"""

	# checkpoint : sv_at_edges unsorted variations

	with open( args.reference_info, 'a' ) as out_file:

		out_file.write(f'{cds_name}_{allele_id}\t')

		pos_list = []

		for pos, alt_dict in pos_dict.items():

			alt_list = []

			for alt, qual in alt_dict.items():

				if type(alt) is list:
					alt = ';'.join(alt)

				alt_list.append("-".join([alt, qual]))

			pos_list.append('>'.join([pos, "/".join(alt_list)]))

		out_file.write(",".join(pos_list))

		out_file.write('\n')

		out_file.close()


def get_ref_alt_qual_of_position_s_variant_dict( vcf_file: str, sv_at_edges: list, cds_name: str, allele_id: str ) -> dict:
	"""
	Read {sample}.vcf to create dictionary that contains positions
	of variants of allele of cds.

	Parameters
	----------
	vcf_file : {sample}.vcf contains {allele_id}'s variants for {cds_name}
	sv_at_edges : Info-formatted variations consisting of uncovered variations
	cds_name : name of CDS of which positions will be taken
	allele_id : ID of allele of CDS of which positions will be taken

	Returns
	-------
	pos_dict : Dictionary of variant positions for allele
	"""

	pos_dict = OrderedDict()

	has_variant = False

	with open( vcf_file, 'r' ) as file:

		for line in file.readlines():

			if not line.startswith("#"):

				vcf_line = Vcf(line)

				has_variant = True
	
				pos_ref = f'{vcf_line.pos}*{vcf_line.ref}'

				if not pos_ref in pos_dict.keys():

					pos_dict[pos_ref] = { vcf_line.alt : vcf_line.qual }

				else:    # if there is a proof for a high quality variant take
							# with the highest

					if pos_dict[pos_ref][vcf_line.alt] <= vcf_line.qual:

						pos_dict[pos_ref][vcf_line.alt] = vcf_line.qual

		file.close()

	for info in sv_at_edges:

		has_variant = True
		
		for i in range(len(info.pos_list)):

			pos_ref = f'{info.pos_list[i]}*{info.ref_list[i]}'

			if not pos_ref in pos_dict.keys():

				pos_dict[pos_ref] = { info.alt_list[i] : info.qual_list[i] }

			else:

				if pos_dict[pos_ref][info.alt_list[i]] <= info.qual_list[i]:

					pos_dict[pos_ref][info.alt_list[i]] = info.qual_list[i]

	sorted_pos_dict = OrderedDict(sorted(pos_dict.items(), key=lambda t: int(t[0].split('*')[0])))

	if has_variant:

		write_allele_defining_variant_list_to_file( cds_name, sv_at_edges, allele_id, sorted_pos_dict )

	return sorted_pos_dict


def merge_intervals(intervals: dict) -> list:
	"""
	Merge intersecting intervals

	Parameter
	---------
	intervals: Dictionary with sequence identifiers as keys and a list of
			   lists as values. Each sublist has a start and stop position
			   in the sequence and a dictionary with the coverage for
			   every position in the sequence interval.

	Return
	------
	merged: Dictionary with the result of merging intervals
			that overlapped (coverage data is updated and
			incremented for positions in common).
	"""

	from copy import deepcopy

	merged = [ deepcopy(intervals[0]) ]

	for current in intervals[1:]:

		previous = merged[-1]

		# current and previous intervals intersect
		if current[0] <= previous[1]:

			# determine top position
			previous[1] = max( previous[1], current[1] )

		# current and previous intervals do not intersect
		else:

			merged.append( deepcopy(current) )

	return merged


def read_paf_file(file_name: str) -> None:
	"""
	Return aligned positions in reference and sample

	Parameter
	---------
	file_name : Name of Paf file
	"""

	paf_list = []
	mapped_reference_regions, mapped_sample_regions = {}, {}
	with open( file_name, 'r' ) as file:

		for i, line in enumerate(file.readlines()):
			
			paf_line = Paf(line)
			paf_list.append(paf_line)

			if i == 0:
				cds_len = paf_line.tlen

			mapped_reference_regions.setdefault(paf_line.tname, []).append([paf_line.tstart, paf_line.tend])
			mapped_sample_regions.setdefault(paf_line.tname, []).append([paf_line.qstart, paf_line.qend])

		file.close()

	merged_reference_intervals = { k: merge_intervals(v) for k, v in mapped_reference_regions.items() }
	cds_merged_reference_interval_list = merged_reference_intervals[list(merged_reference_intervals.keys())[0]]
	
	merged_sample_intervals = { k: merge_intervals(v) for k, v in mapped_sample_regions.items() }
	cds_merged_sample_interval_list = merged_sample_intervals[list(merged_sample_intervals.keys())[0]]

	return paf_list, get_nonintersecting_intervals( cds_merged_reference_interval_list, paf_line.tlen ), get_nonintersecting_intervals( cds_merged_sample_interval_list, paf_line.qlen ), paf_line.mapq


def get_nonintersecting_intervals( interval_list : list, cds_len: int ) -> list:
	"""
	Take the length of CDS and aligned regions
	Extract the unaligned regions

	Parameter
	---------
	interval_list : List of aligned regions
	cds_len : Length of CDS to define range

	Return
	------
	nonintersecting_interval_list : List of unaligned regions
	"""

	nonintersecting_interval_list = list()
	for i, interval in enumerate(interval_list):
		
		if i == 0 and interval[0] != 0:
			nonintersecting_interval_list.append( [0, interval[0] - 1] )

		if i < len(interval_list) - 1:
			nonintersecting_interval_list.append( [interval[1] + 1, interval_list[i+1][0] - 1] )

		if i == len(interval_list) - 1 and interval[1] < cds_len:
			nonintersecting_interval_list.append( [ interval[1] + 1, cds_len] )

	return nonintersecting_interval_list


def remove_common_suffices( var1: str, var2: str ) -> [str, str]:
    """
    Take two variations and remove the common suffices

    Parameters
    ----------
    var1 : variation sequence
    var2 : variation sequence

    Returns
    -------
    var1 : updated variation 1 of which common suffix is deleted
    var2 : updated variation 2 of which common suffix is deleted
    """

    check_len = min(len(var1), len(var2))

    var1 = var1[::-1]
    var2 = var2[::-1]

    i=0
    while i < check_len and var1[i] == var2[i]:

        i+=1

    return var1[i:][::-1], var2[i:][::-1]


def remove_common_prefices( pos: int, var1: str, var2: str ) -> [ int, str, str ]:
    """
    Take two variations and remove the common prefices

    Parameter
    ---------
    pos : position of the variation which might be affected by the change
    var1 : variation sequence
    var2 : variation sequence

    Return
    ------
    pos : position of the variation which might be affected by the change
    var1 : updated variation 1 of which common prefix is deleted
    var2 : updated variation 2 of which common prefix is deleted
    """

    check_len = min(len(var1), len(var2))

    i=0
    while i < check_len and var1[i] == var2[i]:

        i+=1

    return pos+i, var1[i:], var2[i:]


def call_variants_of_allele( reference: str, sample: str ) -> str:
	"""
	Align sample sequences onto the reference and create sample.vcf

	Parameter
	---------
	reference: <reference> in <reference>.fasta
	sample: name of sample to be aligned onto the <reference.fasta>

	Return
	------
	sv_at_edges : info-formatted string to add reference files
	"""

	sv_at_edges = []

	# match : 1, mismatch : 4, gap_opening : 6, gap_extend : 1 (Graph aligner's)
	os.system(f"minimap2 -c --cs=long -t {args.threads} -A 1 -B 4 -O 6 -E 1 {reference}.fasta {sample}.fasta 2>/dev/null > {sample}.paf;")
	os.system(f"sort -k6,6 -k8,8n {sample}.paf | paftools.js call -L0 -l0 -f {reference}.fasta -s {sample} - 2>/dev/null > {sample}.vcf;")

	cds_allele = sample.strip('\n').split('/')[-1].split('_')
	cds, allele_id = cds_allele[0], cds_allele[-1]

	reference_seq = str( list( SeqIO.parse( StringIO( open( f'{reference}.fasta', 'r').read() ), 'fasta') )[0].seq )

	sample_seq = str( list( SeqIO.parse( StringIO( open( f'{sample}.fasta', 'r').read() ), 'fasta') )[0].seq )

	if os.path.getsize(f'{sample}.paf') != 0:

		paf_list, unaligned_reference_pos_list, unaligned_sample_pos_list, mapq = read_paf_file(f'{sample}.paf')

		if unaligned_reference_pos_list != [] or unaligned_sample_pos_list != []:

			if unaligned_sample_pos_list == []:

				if unaligned_reference_pos_list[0][0] == 0:

					unaligned_sample_pos_list.insert( 0, [ 0, 0 ] )

				elif unaligned_reference_pos_list[-1][1] == paf_list[0].tlen:

					unaligned_sample_pos_list.append( [ paf_list[0].qlen, paf_list[0].qlen ] )

			# missing interval in the reference : i.e. [[479, 618]] [[0, 71], [562, 582]]
			if len(unaligned_reference_pos_list) < len(unaligned_sample_pos_list):

				# sample[0,...] : i.e.[[479, 618]] [[0, 71], [562, 582]]
				if unaligned_reference_pos_list[0][0] != 0 and unaligned_sample_pos_list[0][0] == 0:

					unaligned_reference_pos_list.insert( 0, [ 0, 0 ] )

			elif len(unaligned_reference_pos_list) > len(unaligned_sample_pos_list):

				if unaligned_reference_pos_list[0][0] == 0 and unaligned_sample_pos_list[0][0] != 0:

					unaligned_sample_pos_list.insert( 0, [ 0, 0 ] )

			if len(unaligned_reference_pos_list) == len(unaligned_sample_pos_list):

				if unaligned_reference_pos_list[0][1] == paf_list[0].tlen and unaligned_sample_pos_list[0][0] == 0:

					unaligned_reference_pos_list.insert( 0, [ 0, 0 ] )
					unaligned_sample_pos_list.append( [ paf_list[0].qlen, paf_list[0].qlen ] )

			for pos_ref, pos_alt in zip(unaligned_reference_pos_list, unaligned_sample_pos_list):

				# take the different-sized variations
				if pos_ref[1] - pos_ref[0] != pos_alt[1] - pos_alt[0]:

					# both at the beginning
					if pos_ref[0] == 0 and pos_alt[0] == 0:

						# deletion at the beginning
						if pos_ref[1] != 0 and pos_alt[1] == 0:

							pos = 1
							ref = reference_seq[ pos_ref[0] : pos_ref[1] + 2 ]
							alt = sample_seq[ pos_alt[0] : pos_alt[1] + 1 ]

						# insertion at the beginning
						elif pos_ref[1] == 0 and pos_alt[1] != 0:

							pos = 1
							ref = reference_seq[0]
							alt = sample_seq[ pos_alt[0] : pos_alt[1] + 2 ] # 0-based -> 1-based & cover the last index

						# alignment at the beginning
						else:

							pos = 1
							ref = reference_seq[ pos_ref[0] : pos_ref[1] + 2 ]
							alt = sample_seq[ pos_alt[0] : pos_alt[1] + 2 ]

					# both at the end
					elif pos_ref[1] == paf_list[0].tlen and pos_alt[1] == paf_list[0].qlen:

						# insertion at the end
						if pos_ref[0] == paf_list[0].tlen:

							pos = pos_ref[0] - 2
							ref = reference_seq[-1]
							alt = sample_seq[ pos_alt[0] - 2 : pos_alt[1] + 1 ]

						# deletion at the end
						elif pos_alt[0] == paf_list[0].qlen:

							pos = pos_ref[0] - 1
							ref = reference_seq[ pos_ref[0] - 2 : ]
							alt = sample_seq[-1]

						# alignment at the end
						else:

							pos = pos_ref[0] - 1
							ref = reference_seq[ pos_ref[0] - 1 : ]
							alt = sample_seq[ pos_alt[0] - 1 : ]

					# alt is inserted at the end
					elif pos_ref[1] == paf_list[0].tlen and pos_alt[0] == 0:

						pos = pos_ref[1] - 1
						ref = reference_seq[-1]
						alt = sample_seq[ : pos_alt[1] + 2 ]

					else:

						pos = pos_ref[0] - 1
						ref = reference_seq[ pos_ref[0] - 1 : pos_ref[1] + 2 ]
						alt = sample_seq[ pos_alt[0] - 1 : pos_alt[1] + 2 ]

					sv_at_edges.append( Info(f'{pos}*{ref}>{alt}-{mapq}') )

	return sv_at_edges


def sort_zip_and_index_vcf_files(vcf_file: str) -> None:
	"""
	Sort, zip, and index vcf file
	
	Parameter
	---------
	vcf_file : Name of vcf file with its directory
	"""

	# @todo: These might be shortened.
	os.system(f"bcftools sort {vcf_file} -Ov -o {vcf_file}. 2>/dev/null")
	os.system(f"mv {vcf_file}. {vcf_file}")
	os.system(f"bgzip -f {vcf_file} > {vcf_file}.gz; tabix -f -p vcf {vcf_file}.gz")
	os.system(f"bcftools concat -a --rm-dups none {vcf_file}.gz -Ov -o {vcf_file} 2>/dev/null")
	os.system(f"bgzip -f {vcf_file} > {vcf_file}.gz; tabix -f -p vcf {vcf_file}.gz")


def remove_redundant_files(sample: str) -> None:
	"""
	Take the name of sample and remove the related files
	Parameter
	---------
	sample : Name of sample with its directory.
	"""

	for extension in [ "fasta", "sam", "bam", "sorted.bam", "bam.bai", "sorted.bam.bai", "vcf" ]:

		os.system( f"rm {sample}.{extension}")


def convert_info_into_vcf( sv_at_edges: list, cds_name: str, allele_name: str ) -> list:
	"""
	Convert info into vcf format

	Parameter
	---------
	sv_at_edges : long indels that cannot be managed by minimap2,
	              especially at the edges. Info-formatted list.
	cds_name : name of cds
	allele_name : name of allele

	Return
	------
	vcf_list : vcf lines representing SVs at edges.
	"""

	vcf_list = []

	for info in sv_at_edges:

		for i in range(len(info.pos_list)):

			vcf_list.append(f'{cds_name}_1\t{info.pos_list[i]}\t.\t{info.ref_list[i]}\t{info.alt_list[i]}\t{info.qual_list[i]}\t.\tQNAME={allele_name};QSTART={info.pos_list[i]};QSTRAND=+\tGT\t1')

	return vcf_list


def insert_sv_at_edges_to_vcf( sample_vcf_file: str, cds_name: str, allele_name: str, sv_at_edges: list ) -> None:
	"""
	Add long indels to sample's vcf file to represent the allele

	Parameter
	---------
	sample_vcf_file : <sample.vcf> file to add the long indels
	sv_at_edges : list of variations (esp. long indels) to be added
	"""

	with open( sample_vcf_file, 'a' ) as file:

		for vcf_line in convert_info_into_vcf( sv_at_edges, cds_name, allele_name ):

			file.write(f'{vcf_line}\n')

		file.close()


def create_allele_dict_for_a_cds( write_dir: str, allele_name: str, cds_dir: str, cds_name: str ) -> dict:
	"""
	Create allele dictionary for given CDS

	Parameters
	----------
	write_dir : core genome directory
	allele_name : {cds_name}_{allele_id}
	cds_dir : alleles' directory that contains sequences of alleles
	cds_name : CDS name to create its allele dictionary

	Returns
	-------
	allele_dict : position dictionary for allele
	"""

	sample = f"{write_dir}/{allele_name}"
	reference = f"{cds_dir}/references/{cds_name}_1"
	allele_id = get_allele_id_from_allele_name(allele_name)

	# create <sample_allele.vcf>
	sv_at_edges = call_variants_of_allele( reference, sample )

	# {sample}.vcf contains {allele_id}'s variants for {cds_name}
	allele_dict = get_ref_alt_qual_of_position_s_variant_dict( f'{sample}.vcf', sv_at_edges, cds_name, allele_id )

	if len(sv_at_edges) != 0:

		insert_sv_at_edges_to_vcf( f'{sample}.vcf', cds_name, allele_name, sv_at_edges )

	sort_zip_and_index_vcf_files(f'{sample}.vcf')

	# remove_redundant_files()

	return allele_dict


def create_cds_list( cds_dir: str, cds_fasta: str, cds_to_merge_list: list ) -> list:
	"""
	Creates CDS list

	Parameters
	----------
	cds_dir : alleles' directory that contains sequences of alleles
	cds_fasta : FASTA file for CDS
	cds_to_merge_list : all CDSs for reference vcf

	Returns
	-------
	cds_to_merge_list : all CDSs for reference vcf
	"""

	cds_dict = {}

	try:

		for sequence in list( SeqIO.parse( StringIO( open( f"{cds_dir}/{cds_fasta}", 'r' ).read() ), 'fasta') ):

			allele_seq_list = []
		
			for seq_record in SeqIO.parse(f"{args.schema_dir}/{cds}", "fasta"):

				ref_allele_cds, ref_allele_id = str(seq_record.id).split('_')[0], str(seq_record.id).split('_')[-1]
				ref_allele_name = f'{ref_allele_cds}_{ref_allele_id}'

				if ref_allele_id == '1':
			
					out_file = open( f"{args.schema_dir}/references/{ref_allele_name}.fasta", 'w' )

					# allele_id
					out_file.write(f'>{ref_allele_name}\n')

					# allele sequence
					out_file.write(f'{str(seq_record.seq)}\n')

					out_file.close()

					Path(f"{cds_dir}/alleles/{ref_allele_cds}").mkdir(exist_ok=True)

			cds_name = get_cds_name_from_allele_name(sequence.id)
			allele_id = get_allele_id_from_allele_name(sequence.id)

			allele_name = f"{cds_name}_{allele_id}"

			if allele_id != '1':

				write_dir = f"{cds_dir}/alleles/{cds_name}"

				with open( f"{write_dir}/{allele_name}.fasta", "w" ) as out_file:

					out_file.write(f">{allele_name}\n{str(sequence.seq)}\n")
					out_file.close()

				cds_dict[cds_name] = create_allele_dict_for_a_cds( write_dir, allele_name, cds_dir, cds_name )

	except FileNotFoundError:

		pass # {cds_fasta} file does not exist.

	try:

		cds_name = cds_fasta

		if cds_fasta.endswith('.fasta'):

			cds_name = cds_fasta[:-6]

			if cds_name.endswith('_short'):
				cds_name = cds_name[:-6]

		elif cds_fasta.endswith('.fa'):

			cds_name = cds_fasta[:-3]

			if cds_name.endswith('_short'):
				cds_name = cds_name[:-6]

		wd = f"{cds_dir}/alleles/{cds_name}" # directory of allele's vcf files

		command_list = []

		no_alleles_of_cds = len(glob.glob(f'{wd}/*_*.vcf.gz'))

		if no_alleles_of_cds == 0: # skip CDS with 1 alleles to avoid redundance
			pass

		elif no_alleles_of_cds == 1:

			cds_to_merge_list.append(cds_name)

			# remove .vcf.gz get the allele_id ..._'allele_id'
			unzip_allele_id = glob.glob(f'{wd}/*_*.vcf.gz')[0][:-7].split('_')[-1]
			command_list.append(f"gunzip {wd}/{cds_name}_{unzip_allele_id}.vcf.gz")
			command_list.append(f"mv {wd}/{cds_name}_{unzip_allele_id}.vcf {wd}/{cds_name}.vcf")

		else:

			cds_to_merge_list.append(cds_name)

			command_list.append(f"bcftools merge {' '.join(glob.glob(f'{wd}/*_*.vcf.gz'))} -O v -o {wd}/{cds_name}.vcf")

		for command in command_list:
			subprocess.call(command, shell=True, stdout=subprocess.DEVNULL)

	except FileNotFoundError:

		pass # There is no vcf file to merge.

	return cds_to_merge_list


def create_reference_vcf_fasta( wd: str, cds_to_merge_list: list ) -> None:
	"""
	Creates FASTA file for reference

	Parameters
	----------
	wd : working directory
	cds_to_merge_list : all CDSs for reference vcf
	"""

	reference_file = open( f"{args.reference_vcf}.temp", 'w' )

	contig_info_set = set()

	for cds in cds_to_merge_list:

		cds_vcf_file_name = f"{wd}/alleles/{cds}/{cds}.vcf"

		with open( cds_vcf_file_name, 'r' ) as cds_file:

			for line in cds_file.readlines():

				line = line.strip('\n')

				if line.startswith("##contig"):

					contig_info_set.add(line)

				elif line.startswith("##INFO"):

					contig_info_set.add(line)

				elif not line.startswith("#"):

					fields = (line.split('\t'))[:8]
					fields.append("GT")
					fields.append("1")

					line = '\t'.join(fields)

					reference_file.write(f"{line}\n")

			cds_file.close()

	reference_file.close()

	header_file = open( f"{wd}/alleles/header.txt", 'w' )

	header = [  '##fileformat=VCFv4.2',
				'##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">',
				'##FORMAT=<ID=PL,Number=G,Type=Float,Description="Phred-scaled Genotype Likelihoods">',
				'##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Approximate read depth (reads with MQ=255 or with bad mates are filtered)">',
				'##FORMAT=<ID=AD,Number=R,Type=Integer,Description="Allelic depths for the ref and alt alleles in the order listed">',
				'##FORMAT=<ID=RO,Number=1,Type=Integer,Description="Reference allele observation count, with partial observations recorded fractionally">',
				'##FORMAT=<ID=GL,Number=3,Type=Float,Description="Likelihoods for RR,RA,AA genotypes (R=ref,A=alt)">',
				'##FORMAT=<ID=QR,Number=1,Type=Integer,Description="Reference allele quality sum in phred">',
				'##FORMAT=<ID=QA,Number=A,Type=Integer,Description="Alternate allele quality sum in phred">',
				'##FORMAT=<ID=AO,Number=A,Type=Integer,Description="Alternate allele observation count">',
				'##INFO=<ID=AC,Number=A,Type=Integer,Description="Allele Counts">',
				'##INFO=<ID=AO,Number=A,Type=Integer,Description="Count of full observations of this alternate haplotype.">',
				'##INFO=<ID=AN,Number=1,Type=Integer,Description="Total number of alleles">'
			 ]

	header.extend(list(contig_info_set))

	header.append("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t"
				  "REFERENCE\n")

	header_file.write("\n".join(header))

	header_file.close()

	os.system(f"cat {wd}/alleles/header.txt {args.reference_vcf}.temp > "
			  f"{args.reference_vcf}")

	os.system(f"rm {args.reference_vcf}.temp; rm {wd}/alleles/header.txt;")

	# merge all CDS fasta files to create reference FASTA

	with open( args.reference_fasta, 'a' ) as f:

		for file in glob.glob(f'{wd}/references/*.fasta'):

			with open( file, 'r' ) as infile:

				f.write(infile.read().strip('\n')+'\n')

		f.close()
	

def get_cds_list() -> list:
	"""
	Returns names of CDSs as list

	Returns
	-------
	cds_list : List of CDS names
	"""

	cds_list = []

	for file_name in os.listdir(args.schema_dir):

		if file_name.endswith('.fasta'):
			cds_list.append(file_name)

	return cds_list


def create_dirs_to_split_sequences_to_call_variants() -> None:
	"""
	delete references and alleles directory including its content,
	which is created by milestone.py create_schema
	"""
	
	try:

		# create reference FASTA directory for each core gene
		Path(f"{args.schema_dir}/references").mkdir(exist_ok=True)

		# create FASTA directory for each core gene's alleles
		Path(f"{args.schema_dir}/alleles").mkdir(exist_ok=True)

	except OSError:

		print(f"Creation of the directories {args.schema_dir}/references "
			  f"and {args.schema_dir}/alleles are failed.")


def delete_sequences_created_by_milestone() -> None:
	"""
	delete references and alleles directory including its content,
	which is created by milestone.py create_schema
	"""
	
	try:

		ref_dir = f"{args.schema_dir}/references"
		shutil.rmtree(ref_dir)

	except OSError as e:

		print(f"Error: {ref_dir} : {e.strerror}")

	try:

		allele_dir = f"{args.schema_dir}/alleles"
		shutil.rmtree(allele_dir)

	except OSError as e:

		print(f"Error: {allele_dir} : {e.strerror}")


if __name__ == "__main__":

	parser = argparse.ArgumentParser(add_help = True)

	parser.add_argument('-s', '--schema_dir', 
						type = str,
						required = True,
						help = 'Directory of allele sequences.')

	parser.add_argument('-v', '--reference_vcf',
						type = str,
						required = True,
						help = 'VCF file of reference genome.')
	
	parser.add_argument('-f', '--reference_fasta',
						type = str,
						required=True,
						help = 'FASTA file of reference genome.')
	
	parser.add_argument('-i', '--reference_info',
						type = str,
						required=True,
						help='Info file of reference genome')
	
	parser.add_argument('-t', '--threads',
						type = str,
						required = False,
						default = '1',
						help='Number of threads [ default: 1 ].')

	args = parser.parse_args()

	cds_list = get_cds_list()

	create_dirs_to_split_sequences_to_call_variants()

	cds_to_merge_list = [] # to merge all CDSs for reference vcf

	for cds in cds_list:

		create_cds_list(args.schema_dir, cds, cds_to_merge_list)

	create_reference_vcf_fasta(args.schema_dir, cds_to_merge_list)

	# it creates reference_info.txt file in case that
	# cds sequences are provided as only references
	# which have no alleles. This step is for the further analysis.
	reference_info_txt_file = Path(args.reference_info)
	reference_info_txt_file.touch(exist_ok=True)

	delete_sequences_created_by_milestone()