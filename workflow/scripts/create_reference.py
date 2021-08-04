#!/usr/bin/env python3

"""
-----------------------------------------------
Aim: Reference FASTA VCF and INFO file creation
-----------------------------------------------
Authors: @fatmakhv
The latest update: June 14, 2021
-----------------------------------------------
"""

import argparse, glob, os, shutil, subprocess
from Bio import SeqIO
from io import StringIO  # python3
from pathlib import Path


class vcf:

	def __init__(self, vcf_line):

		fields = vcf_line.split('\t')

		self.chr = get_cds_name_from_allele_name(fields[0])
		self.pos = fields[1]
		self.id = fields[2]
		self.ref = fields[3]
		self.alt = fields[4]
		self.qual = fields[5]
		self.filter = fields[6]
		self.info = fields[7]
		self.format = fields[8]
		self.sample = fields[9]

	def __repr__(self):

		return f"chr: {self.chr}\tpos: {self.pos}\tid: {self.id}\tref: " \
				f"{self.ref}\talt: {self.alt}\tqual: {self.qual}\tfilter: " \
				f"{self.filter}\tinfo: {self.info}\tformat: {self.format}\t" \
				f"sample: {self.sample}"


def write_allele_defining_variant_list_to_file(cds_name: str, allele_id: str,
														pos_dict: dict) -> None:
	"""
	Writes the variants of alleles that are not equal to the reference CDS

	Parameters
	----------
	cds_name : name of cds of which allele will be written
	allele_id : ID of allele of which variant will be written
	pos_dict : positions of variants for the allele
	"""

	with open(args.reference_info, 'a') as out_file:
		
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


def get_ref_alt_qual_of_position_s_variant_dict(vcf_file: str, cds_name: str,
														allele_id: str) -> dict:
	"""
	Reads {sample}.vcf to create dictionary that contains positions
	of variants of allele of cds.

	Parameters
	----------
	vcf_file : {sample}.vcf contains {allele_id}'s variants for {cds_name}
	cds_name : name of CDS of which positions will be taken
	allele_id : ID of allele of CDS of which positions will be taken

	Returns
	-------
	pos_dict : Dictionary of variant positions for allele
	"""

	pos_dict = {}

	has_variant = False

	with open(vcf_file, 'r') as file:

		for line in file.readlines():

			if not line.startswith("#"):

				has_variant = True

				vcf_line = vcf(line)

				pos_ref = f'{vcf_line.pos}*{vcf_line.ref}'

				if not pos_ref in pos_dict.keys():

					pos_dict[pos_ref] = {vcf_line.alt:vcf_line.qual}

				else:	# if there is a proof for a high quality variant take
						# with the highest

					if pos_dict[pos_ref][vcf_line.alt] <= vcf_line.qual:
				
						pos_dict[pos_ref][vcf_line.alt] = vcf_line.qual

		file.close()

	if has_variant:
		write_allele_defining_variant_list_to_file(cds_name, allele_id, pos_dict)

	return pos_dict


def get_allele_id_from_allele_name(allele_name: str) -> str:
	"""
	Gets {cds_name}_{allele_id} and returns {allele_id}

	Parameters
	----------
	allele_name : {cds_name}_{allele_id}

	Returns
	-------
	allele_id : allele ID for given allele_name
	"""

	allele_id = allele_name.strip("\n").split("_")[-1]

	return allele_id


def get_cds_name_from_allele_name(allele_name: str) -> str:
	"""
	Gets {cds_name}_{allele_id} and returns {cds_name}

	Parameters
	----------
	allele_name : {cds_name}_{allele_id}

	Returns
	-------
	cds_name : name of CDS for given allele_name
	"""

	cds_name = allele_name.strip("\n").split("_")[0]

	return cds_name


def create_allele_dict_for_a_cds(write_dir: str, allele_name: str, cg_dir: str,
								 		cds_name: str, threads: str) -> dict:
	"""
	Creates allele dictionary for given CDS

	Parameters
	----------
	write_dir : core genome directory
	allele_name : {cds_name}_{allele_id}
	cg_dir : schema_seed directory that contains sequences of alleles
	cds_name : CDS name to create its allele dictionary
	threads : number of threads to run minimap2

	Returns
	-------
	allele_dict : position dictionary for allele
	"""

	sample = f"{write_dir}/{allele_name}"
	reference = f"{cg_dir}/references/{cds_name}_1"
	allele_id = get_allele_id_from_allele_name(allele_name)

	# create <sample_allele.vcf>
	command_list = []

	command_list.append(f"minimap2 -ax asm5 {reference}.fasta {sample}.fasta -t"
						f" {threads} --cs=long -o {sample}.sam 2>&1")
	command_list.append(f"samtools view -Sb {sample}.sam > {sample}.bam")
	command_list.append(f"samtools index {sample}.bam")
	command_list.append(f"samtools sort {sample}.bam -o {sample}.sorted.bam")
	command_list.append(f"samtools index {sample}.sorted.bam")
	command_list.append(f"bcftools mpileup -O u -f {reference}.fasta {sample}."
						f"sorted.bam | bcftools call --ploidy 1 -Ov -c -v -o "
						f"{sample}.vcf")

	for command in command_list:
		subprocess.call(command, shell=True, stdout=subprocess.DEVNULL)

	# {sample}.vcf contains {allele_id}'s variants for {cds_name}
	allele_dict = get_ref_alt_qual_of_position_s_variant_dict(f'{sample}.vcf',
															cds_name, allele_id)
	
	command_list = []

	command_list.append(f"bgzip -c {sample}.vcf > {sample}.vcf.gz")
	command_list.append(f"tabix -p vcf {sample}.vcf.gz")
	command_list.append(f"rm {sample}.fasta; rm {sample}.sam; rm {sample}.bam;"
						f" rm {sample}.sorted.bam; rm {sample}.bam.bai; rm "
						f"{sample}.sorted.bam.bai; rm {sample}.vcf;")

	for command in command_list:
		subprocess.call(command, shell=True, stdout=subprocess.DEVNULL)

	return allele_dict


def create_cds_list(cg_dir: str, cds_fasta: str, cds_to_merge_list: list,
														threads: str) -> list:
	"""
	Creates CDS list

	Parameters
	----------
	cg_dir : schema_seed directory that contains sequences of alleles
	cds_fasta : FASTA file for CDS
	cds_to_merge_list : all CDSs for reference vcf
	threads : number of threads to run minimap2

	Returns
	-------
	cds_to_merge_list : all CDSs for reference vcf
	"""

	write_dir = f"{cg_dir}/references"

	cds_dict = {}

	try:

		for sequence in list(SeqIO.parse(StringIO(open(f"{cg_dir}/{cds_fasta}",
														'r').read()), 'fasta')):

			cds_name = get_cds_name_from_allele_name(sequence.id)
			allele_id = get_allele_id_from_allele_name(sequence.id)

			allele_name = f"{cds_name}_{allele_id}"

			if allele_id == '1':

				Path(f"{cg_dir}/alleles/{cds_name}").mkdir(exist_ok=True)

			else:
				write_dir = f"{cg_dir}/alleles/{cds_name}"

			with open(f"{write_dir}/{allele_name}.fasta", "w") as out_file:
				out_file.write(f">{allele_name}\n{str(sequence.seq)}\n")
				out_file.close()

			if allele_id != '1':

				cds_dict[cds_name] = create_allele_dict_for_a_cds(write_dir,
										allele_name, cg_dir, cds_name, threads)

	except FileNotFoundError:
		pass # {cds_fasta} file does not exist.

	try:

		cds_name = str(cds_fasta.strip('.fasta')) # remove fasta extension

		wd = f"{cg_dir}/alleles/{cds_name}" # directory of allele's vcf files

		command_list = []

		no_alleles_of_cds = len(glob.glob(f'{wd}/*_*.vcf.gz'))

		if no_alleles_of_cds == 0: # skip CDS with 1 alleles to avoid redundance
			pass

		elif no_alleles_of_cds == 1:

			cds_to_merge_list.append(cds_name)
			command_list.append(f"gunzip {wd}/{cds_name}_2.vcf.gz")
			command_list.append(f"mv {wd}/{cds_name}_2.vcf {wd}/{cds_name}.vcf")

		else:

			cds_to_merge_list.append(cds_name)
			command_list.append(f"bcftools merge "
								f"{' '.join(glob.glob(f'{wd}/*_*.vcf.gz'))} -O"
								f" v -o {wd}/{cds_name}.vcf")

		for command in command_list:
			subprocess.call(command, shell=True, stdout=subprocess.DEVNULL)

	except FileNotFoundError:
		pass # There is no vcf file to merge.

	return cds_to_merge_list


def create_reference_vcf_fasta(wd: str, cds_to_merge_list: list) -> None:
	"""
	Creates FASTA file for reference 

	Parameters
	----------
	wd : working directory
	cds_to_merge_list : all CDSs for reference vcf
	"""

	reference_file = open(f"{args.reference_vcf}.temp", 'w')

	contig_info_set = set()

	for cds in cds_to_merge_list:

		cds_vcf_file_name = f"{wd}/alleles/{cds}/{cds}.vcf"

		with open(cds_vcf_file_name, 'r') as cds_file:

			for line in cds_file.readlines():

				line = line.strip('\n')

				if line.startswith("##contig"):
					contig_info_set.add(line)

				elif line.startswith("##INFO"):
					contig_info_set.add(line)

				elif not line.startswith("#"):
					fields = (line.split('\t'))[:9]
					fields.append("1:60,0")
					line = '\t'.join(fields)
					reference_file.write(f"{line}\n")

			cds_file.close()

	reference_file.close()

	header_file = open(f"{wd}/alleles/header.txt", 'w')

	header = ['##fileformat=VCFv4.2',
			'##FILTER=<ID=PASS,Description="All filters passed">")',
			'##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">',
			'##FORMAT=<ID=PL,Number=G,Type=Integer,Description="Normalized, Phred-scaled likelihoods for genotypes as defined in the VCF specification">',
    		'##INFO=<ID=AC,Number=A,Type=Integer,Description="Allele Counts">',
    		'##INFO=<ID=AN,Number=1,Type=Integer,Description="Total number of '
			'alleles">']

	header.extend(list(contig_info_set))

	header.append("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t"
				  "REFERENCE\n")

	header_file.write("\n".join(header))

	header_file.close()

	os.system(f"cat {wd}/alleles/header.txt {args.reference_vcf}.temp > "
			  f"{args.reference_vcf}")

	os.system(f"rm {args.reference_vcf}.temp; rm {wd}/alleles/header.txt;")

	# merge all CDS fasta files to create reference FASTA
	os.system(f"cat {wd}/references/*.fasta > {args.reference_fasta}")


def get_cg_list(cg_schema_file: str) -> list:
	"""
	Returns names of core genes as list
	
	Parameters
	----------
	cg_schema_file : file that is created by chewBBACA and contains core genes

	Returns
	-------
	cg_list : List of core genes
	"""

	cg_list = []

	with open(cg_schema_file, 'r') as file:

		for line in file.readlines():

			cds = line.strip()
			cg_list.append(f"{cds}")

		file.close()

	return cg_list


if __name__ == "__main__":

	parser = argparse.ArgumentParser(add_help = True)

	parser.add_argument('--cgmlst_dir', required=True,
						help='cgMLST directory')
	parser.add_argument('--schema_seed_dir', required=True,
						help='schema seed directory')
	parser.add_argument('--reference_vcf', required=True,
						help='VCF file of reference genome')
	parser.add_argument('--reference_fasta', required=True,
						help='FASTA file of reference genome')
	parser.add_argument('--reference_info', required=True,
						help='Info file of reference genome')
	parser.add_argument('--threads', required=True,
						help='Number of threads')

	args = parser.parse_args()

	cg_list = get_cg_list(f'{args.cgmlst_dir}/cgMLSTschema.txt')

	try:

		# create reference FASTA directory for each core gene
		Path(f"{args.schema_seed_dir}/references").mkdir(exist_ok=True)

		# create FASTA directory for each core gene's alleles
		Path(f"{args.schema_seed_dir}/alleles").mkdir(exist_ok=True)

	except OSError:

		print(f"Creation of the directories {args.schema_seed_dir}/references "
							f"and {args.schema_seed_dir}/alleles are failed.")

	cds_to_merge_list = [] # to merge all CDSs for reference vcf

	# create reference vcf
	for cds in cg_list:

		create_cds_list(f"{args.schema_seed_dir}", cds, cds_to_merge_list,
																{args.threads})

	create_reference_vcf_fasta(f"{args.schema_seed_dir}", cds_to_merge_list)

	# @todo redundant files will be removed.
	# os.system(f"rm -rf {args.schema_seed_dir}")