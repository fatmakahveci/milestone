#!/usr/bin/env python3

import argparse, glob, os, shutil, subprocess
from Bio import SeqIO
from io import StringIO  # python3
from pathlib import Path


def create_vcf(cg_dir, cds_fasta, cds_to_merge_list):

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

		no_alleles_of_cds = len(glob.glob(f'{wd}/*.vcf.gz'))

		if no_alleles_of_cds == 0: # skip CDS with one alleles to avoid redundance
			pass

		elif no_alleles_of_cds == 1:
			cds_to_merge_list.append(cds_name)
			command_list.append(f"mv {wd}/{cds_name}_2.vcf.gz {wd}/{cds_name}.vcf.gz")

		else:
			cds_to_merge_list.append(cds_name)
			command_list.append(f"bcftools merge {' '.join(glob.glob(f'{wd}/*.vcf.gz'))} -O v -o {wd}/{cds_name}.vcf 2>&1")
			command_list.append(f"bgzip -c {wd}/{cds_name}.vcf > {wd}/{cds_name}.vcf.gz")
			command_list.append(f"tabix {wd}/{cds_name}.vcf.gz")
		for command in command_list:
			subprocess.call(command, shell=True, stdout=subprocess.DEVNULL)
			
	except FileNotFoundError:
		pass

	return cds_to_merge_list