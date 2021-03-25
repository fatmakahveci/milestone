#!/usr/bin/env python3

import argparse, glob, os, shutil, subprocess
from Bio import SeqIO
from io import StringIO  # python3
from pathlib import Path


def create_cds_vcf(cg_dir, cds_fasta, cds_to_merge_list):

	write_dir = f"{cg_dir}/references"
	try:

		for sequence in list(SeqIO.parse(StringIO(open(f"{cg_dir}/{cds_fasta}", 'r').read()), 'fasta')):

			allele_name = str(sequence.id).strip("\n").split("_")
			cds_name, allele_id = str(allele_name[0]), str(allele_name[-1])
			allele_name = f"{cds_name}_{allele_id}"

			if allele_id == '1':

				Path(f"{cg_dir}/alleles/{cds_name}").mkdir(exist_ok=True)

			else:
				write_dir = f"{cg_dir}/alleles/{cds_name}"

			with open(f"{write_dir}/{allele_name}.fasta", "w") as out_file:
				out_file.write(f">{allele_name}\n{str(sequence.seq)}\n")
				out_file.close()

			if allele_id != '1':

				command_list = []  # Command List
				sample = f"{write_dir}/{allele_name}"
				reference = f"{cg_dir}/references/{cds_name}_1"
				command_list.append(f"minimap2 -ax asm5 {reference}.fasta {sample}.fasta --cs=long -o {sample}.sam 2>&1")
				command_list.append(f"samtools view -Sb {sample}.sam > {sample}.bam")
				command_list.append(f"samtools index {sample}.bam")
				command_list.append(f"samtools sort {sample}.bam -o {sample}.sorted.bam")
				command_list.append(f"samtools index {sample}.sorted.bam")
				command_list.append(f"bcftools mpileup -O u -f {reference}.fasta {sample}.sorted.bam | bcftools call --ploidy 1 -Ov -c -v -o {sample}.vcf")
				command_list.append(f"bgzip -c {sample}.vcf > {sample}.vcf.gz")
				command_list.append(f"tabix -p vcf {sample}.vcf.gz")
				command_list.append(f"rm {sample}.fasta; rm {sample}.sam; rm {sample}.bam; rm {sample}.sorted.bam; rm {sample}.bam.bai; rm {sample}.sorted.bam.bai; rm {sample}.vcf;")

				for command in command_list:
					subprocess.call(command, shell=True, stdout=subprocess.DEVNULL)

	except FileNotFoundError:
		print(f'Allele fasta file cannot be found')

	try:

		cds_name = str(cds_fasta.strip('.fasta')) # remove fasta extension

		wd = f"{cg_dir}/alleles/{cds_name}" # directory of allele's vcf files

		command_list = []

		no_alleles_of_cds = len(glob.glob(f'{wd}/*_*.vcf.gz'))

		if no_alleles_of_cds == 0: # skip CDS with one alleles to avoid redundance
			pass

		elif no_alleles_of_cds == 1:

			cds_to_merge_list.append(cds_name)
			command_list.append(f"gunzip {wd}/{cds_name}_2.vcf.gz")
			command_list.append(f"mv {wd}/{cds_name}_2.vcf {wd}/{cds_name}.vcf")

		else:

			cds_to_merge_list.append(cds_name)
			command_list.append(f"bcftools merge {' '.join(glob.glob(f'{wd}/*_*.vcf.gz'))} -O v -o {wd}/{cds_name}.vcf")

		for command in command_list:
			subprocess.call(command, shell=True, stdout=subprocess.DEVNULL)

	except FileNotFoundError:
		pass

	return cds_to_merge_list

def create_reference_vcf_fasta(wd, cds_to_merge_list):

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
    	     '##INFO=<ID=AC,Number=A,Type=Integer,Description="Allele Counts">',
    	     '##INFO=<ID=AN,Number=1,Type=Integer,Description="Total number of alleles">']
	header.extend(list(contig_info_set))
	header.append("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tREFERENCE\n")

	header_file.write("\n".join(header))

	header_file.close()

	os.system(f"cat {wd}/alleles/header.txt {args.reference_vcf}.temp > {args.reference_vcf}")

	os.system(f"rm {args.reference_vcf}.temp; rm {wd}/alleles/header.txt;")

	# merge all CDS fasta files to create reference FASTA
	os.system(f"cat {wd}/references/*.fasta > {args.reference_fasta}")

def get_cg_list(cg_schema_file):

	cg_list = []

	with open(cg_schema_file, 'r') as file:

		for line in file.readlines():

			cds = line.strip()
			cg_list.append(f"{cds}")

		file.close()

	return cg_list


if __name__ == "__main__":

	parser = argparse.ArgumentParser()

	parser.add_argument('--cgmlst_dir')
	parser.add_argument('--schema_seed_dir')
	parser.add_argument('--reference_vcf')
	parser.add_argument('--reference_fasta')

	args = parser.parse_args()

	cg_list = get_cg_list(f'{args.cgmlst_dir}/cgMLSTschema.txt')

	try:
		Path(f"{args.schema_seed_dir}/references").mkdir(exist_ok=True)  # create reference FASTA directory for each core gene
		Path(f"{args.schema_seed_dir}/alleles").mkdir(exist_ok=True)  # create FASTA directory for each core gene's alleles

	except OSError:
		print(f"Creation of the directories {args.schema_seed_dir}/references and {args.schema_seed_dir}/alleles are failed.")

	cds_to_merge_list = [] # to merge all CDSs for reference vcf

	# create reference vcf
	for cds in cg_list:
		create_cds_vcf(f"{args.schema_seed_dir}", cds, cds_to_merge_list)

	create_reference_vcf_fasta(f"{args.schema_seed_dir}", cds_to_merge_list)

	# os.system(f"rm -rf {args.schema_seed_dir}")

