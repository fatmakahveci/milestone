#!/usr/bin/env python3

######################################################
## Author: @fatmakhv                                ##
## The latest update: 27/10/2021                    ##
## Aim: Reference FASTA, VCF and INFO file creation ##
######################################################


import argparse, glob, os, shutil, subprocess, sys

from Bio import SeqIO
from collections import OrderedDict
from pathlib import Path

SCRIPT_DIR = Path(__file__).resolve().parent
if str(SCRIPT_DIR) not in sys.path:
	sys.path.insert(0, str(SCRIPT_DIR))

from wgmlst_utils import (
	get_allele_id_from_allele_name,
	get_cds_name_from_allele_name,
	select_reference_record,
)
from reference_pipeline import (
	call_variants_of_allele as pipeline_call_variants_of_allele,
	sort_zip_and_index_vcf_files as pipeline_sort_zip_and_index_vcf_files,
)
from reference_formats import Paf, Vcf
from script_utils import ParsedVariationInfo, run_checked_command


Info = ParsedVariationInfo

def run_command(command: list[str], stdout=None, stderr=None) -> None:
	run_checked_command(command, stdout=stdout, stderr=stderr)


def write_allele_defining_variant_list_to_file( cds_name: str, allele_id: str, pos_dict: dict ) -> None:
	"""
	Write the variants of alleles that are not equal to the reference CDS

	Parameter
	---------
	cds_name : name of cds of which allele will be written
	allele_id : ID of allele of which variant will be written
	pos_dict : positions of variants for the allele
	"""

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


def get_ref_alt_qual_of_position_s_variant_dict( vcf_file: str, cds_name: str, allele_id: str ) -> dict:
	"""
	Read {sample}.vcf to create dictionary that contains positions
	of variants of allele of cds.

	Parameter
	---------
	vcf_file : {sample}.vcf contains {allele_id}'s variants for {cds_name}
	cds_name : name of CDS of which positions will be taken
	allele_id : ID of allele of CDS of which positions will be taken

	Return
	------
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

	if has_variant:

		write_allele_defining_variant_list_to_file( cds_name=cds_name, allele_id=allele_id, pos_dict=pos_dict)

	return pos_dict


def call_variants_of_allele( reference: str, sample: str ) ->  str:
	"""
	Align sample sequences onto the reference and create sample.vcf

	Parameter
	---------
	reference: <reference> in <reference>.fasta
	sample: name of sample to be aligned onto the <reference.fasta>
	"""

	pipeline_call_variants_of_allele(reference, sample, int(args.threads), run_command)


def sort_zip_and_index_vcf_files(vcf_file: str) -> None:
	"""
	Sort, zip, and index vcf file
	
	Parameter
	---------
	vcf_file : Name of vcf file with its directory
	"""

	pipeline_sort_zip_and_index_vcf_files(vcf_file, run_command)


def remove_redundant_files(sample: str) -> None:
	"""
	Take the name of sample and remove the related files

	Parameter
	---------
	sample : Name of sample with its directory.
	"""

	for extension in [ "fasta", "sam", "bam", "sorted.bam", "bam.bai", "sorted.bam.bai", "vcf" ]:
		Path(f"{sample}.{extension}").unlink(missing_ok=True)


def create_allele_dict_for_a_cds(
	write_dir: str,
	allele_name: str,
	cds_dir: str,
	cds_name: str,
	reference_allele_name: str,
) -> dict:
	"""
	Create allele dictionary for given CDS

	Parameter
	---------
	write_dir : core genome directory
	allele_name : {cds_name}_{allele_id}
	cds_dir : alleles' directory that contains sequences of alleles
	cds_name : CDS name to create its allele dictionary
	"""

	sample = f"{write_dir}/{allele_name}"
	reference = f"{cds_dir}/references/{reference_allele_name}"
	allele_id = get_allele_id_from_allele_name( allele_name=allele_name )

	# create <sample_allele.vcf>
	call_variants_of_allele( reference, sample )

	# {sample}.vcf contains {allele_id}'s variants for {cds_name}
	get_ref_alt_qual_of_position_s_variant_dict( f'{sample}.vcf', cds_name, allele_id )

	sort_zip_and_index_vcf_files( vcf_file=f'{sample}.vcf' )

	# remove_redundant_files()


def create_cds_list( cds_dir: str, cds_fasta: str, cds_to_merge_list: list ) -> list:
	"""
	Creates CDS list

	Parameter
	---------
	cds_dir : alleles' directory that contains sequences of alleles
	cds_fasta : FASTA file for CDS
	cds_to_merge_list : all CDSs for reference vcf

	Return
	------
	cds_to_merge_list : all CDSs for reference vcf
	"""

	try:
		with open(f"{cds_dir}/{cds_fasta}", "r", encoding="ascii") as handle:
			locus_records = list(SeqIO.parse(handle, 'fasta'))
		if not locus_records:
			return cds_to_merge_list

		reference_record = select_reference_record(locus_records)
		reference_allele_name = reference_record.id
		reference_cds = get_cds_name_from_allele_name(reference_allele_name)

		with open(f"{args.schema_dir}/references/{reference_allele_name}.fasta", 'w') as out_file:
			out_file.write(f'>{reference_allele_name}\n')
			out_file.write(f'{str(reference_record.seq)}\n')

		Path(f"{cds_dir}/alleles/{reference_cds}").mkdir(exist_ok=True)

		for sequence in locus_records:
			cds_name = get_cds_name_from_allele_name( allele_name=sequence.id )
			allele_id = get_allele_id_from_allele_name( allele_name=sequence.id )
			allele_name = f"{cds_name}_{allele_id}"

			if allele_name != reference_allele_name:
				write_dir = f"{cds_dir}/alleles/{cds_name}"

				with open( f"{write_dir}/{allele_name}.fasta", "w" ) as out_file:

					out_file.write(f">{allele_name}\n{str(sequence.seq)}\n")

				create_allele_dict_for_a_cds(
					write_dir,
					allele_name,
					cds_dir,
					cds_name,
					reference_allele_name,
				)

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

		no_alleles_of_cds = len( glob.glob( f'{wd}/*_*.vcf.gz' ) )

		if no_alleles_of_cds == 0: # skip CDS with 1 alleles to avoid redundance
			pass

		elif no_alleles_of_cds == 1:

			cds_to_merge_list.append(cds_name)

			# remove .vcf.gz get the allele_id ..._'allele_id'
			vcf_stem = Path(glob.glob(f"{wd}/*_*.vcf.gz")[0]).name[:-7]
			unzip_allele_id = get_allele_id_from_allele_name(vcf_stem)
			command_list.append(["gunzip", "-f", f"{wd}/{cds_name}_{unzip_allele_id}.vcf.gz"])
			command_list.append(["mv", f"{wd}/{cds_name}_{unzip_allele_id}.vcf", f"{wd}/{cds_name}.vcf"])

		else:

			cds_to_merge_list.append(cds_name)

			command_list.append(
				["bcftools", "merge", *glob.glob(f"{wd}/*_*.vcf.gz"), "-O", "v", "-o", f"{wd}/{cds_name}.vcf"]
			)

		for command in command_list:
			run_command(command, stderr=subprocess.DEVNULL)

	except FileNotFoundError:

		pass # There is no vcf file to merge.

	return cds_to_merge_list


def create_reference_vcf_fasta( wd: str, cds_to_merge_list: list ) -> None:
	"""
	Creates FASTA file for reference

	Parameter
	---------
	wd : working directory
	cds_to_merge_list : all CDSs for reference vcf
	"""

	contig_info_set = set()
	with open(f"{args.reference_vcf}.temp", "w") as reference_file:
		for cds in cds_to_merge_list:
			cds_vcf_file_name = f"{wd}/alleles/{cds}/{cds}.vcf"
			with open(cds_vcf_file_name, "r") as cds_file:
				for line in cds_file.readlines():
					line = line.strip('\n')
					if line.startswith("##contig"):
						contig_info_set.add(line)
					elif line.startswith("##INFO"):
						contig_info_set.add(line)
					elif not line.startswith("#"):
						fields = ( line.split('\t') )[:8]
						fields.append("GT")
						fields.append("1")
						reference_file.write("\t".join(fields) + "\n")

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

	header.extend( list(contig_info_set) )
	header.append("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tREFERENCE\n")

	with open(f"{wd}/alleles/header.txt", "w") as header_file:
		header_file.write("\n".join(header))

	with open(args.reference_vcf, "w") as output_handle:
		for input_path in [f"{wd}/alleles/header.txt", f"{args.reference_vcf}.temp"]:
			with open(input_path, "r") as input_handle:
				output_handle.write(input_handle.read())

	Path(f"{args.reference_vcf}.temp").unlink(missing_ok=True)
	Path(f"{wd}/alleles/header.txt").unlink(missing_ok=True)

	# merge all CDS fasta files to create reference FASTA

	with open( args.reference_fasta, 'a' ) as f:

		for file in glob.glob( f'{wd}/references/*.fasta' ):

			with open( file, 'r' ) as infile:

				f.write(infile.read().strip('\n')+'\n')
	

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

		create_cds_list( cds_dir=args.schema_dir, cds_fasta=cds, cds_to_merge_list=cds_to_merge_list )

	create_reference_vcf_fasta( wd=args.schema_dir, cds_to_merge_list=cds_to_merge_list )

	# it creates reference_info.txt file in case that
	# cds sequences are provided as only references
	# which have no alleles. This step is for the further analysis.
	reference_info_txt_file = Path(args.reference_info)
	reference_info_txt_file.touch(exist_ok=True)

	delete_sequences_created_by_milestone()
