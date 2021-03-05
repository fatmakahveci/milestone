#!/usr/bin/env python3

import argparse, glob, os, shutil, subprocess
from Bio import SeqIO
from io import StringIO  # python3
from pathlib import Path


def create_cds_vcf(cg_dir, cds_fasta, cds_to_merge_list):

	try:
		for sequence in list(SeqIO.parse(StringIO(open(f"{os.getcwd()}/{cg_dir}/{cds_fasta}", 'r').read()), 'fasta')):

			allele_name = str(sequence.id).strip("\n").split("_")
			cds_name, allele_id = str(allele_name[0]), str(allele_name[-1])
			allele_name = f"{cds_name}_{allele_id}"

			if allele_id == '1':
				write_dir = f"{cg_dir}/references"
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
		print(f'{write_dir}/{allele_name}.fasta file cannot be found')

	try:
		cds_name = str(cds_fasta.strip('.fasta')) # remove fasta extension
		wd = f"{os.getcwd()}/{cg_dir}/alleles/{cds_name}" # directory of allele's vcf files

		command_list = []

		no_alleles_of_cds = len(glob.glob(f'{wd}/*_*.vcf.gz')) + len(glob.glob(f'{wd}/*.vcf')) # the second is temporary.

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

def create_reference_vcf_fasta(wd, cds_to_merge_list, reference_vcf, reference_fasta):

	reference_file = open(f"{reference_vcf}.temp", 'w')

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

	header_file = open(f"{reference_vcf}header.txt", 'w')

	header = ['##fileformat=VCFv4.2',
    	     '##FILTER=<ID=PASS,Description="All filters passed">")',
    	     '##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">',
    	     '##INFO=<ID=AC,Number=A,Type=Integer,Description="Allele Counts">',
    	     '##INFO=<ID=AN,Number=1,Type=Integer,Description="Total number of alleles">']
	header.extend(list(contig_info_set))
	header.append("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tREFERENCE\n")

	header_file.write("\n".join(header))

	os.system(f"cat {reference_vcf}header.txt {reference_vcf}.temp > {reference_vcf}")
	os.system(f"rm {reference_vcf}.temp; rm {reference_vcf}header.txt;")
	

	# merge all CDS fasta files to create reference FASTA
	os.system(f"cat {wd}/references/*.fasta > {reference_fasta}")


def clean_directory(wd, reference_vcf, reference_fasta):

	for root, dirs, files in os.walk(wd):
		for file in files:
			if not(file == {reference_vcf} or file == {reference_fasta}):
				os.unlink(os.path.join(root, file))
		for dir in dirs:
			shutil.rmtree(os.path.join(root, dir))
	os.system(f"rmdir {wd}/alleles 2&>1")
