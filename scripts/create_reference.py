#!/usr/bin/env python3

import argparse, os, shutil, subprocess
from Bio import SeqIO
from io import StringIO  # python3
from pathlib import Path


def create_vcf(cg_dir, cds_fasta):

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
				out_file.write(f">{allele_name}\n{str(sequence.seq)}")
				out_file.close()

			if allele_id != '1':
				command_list = []  # Command List
				sample = f"{write_dir}/{allele_name}"
				reference = f"{cg_dir}/references/{cds_name}_1"
				command_list.append(f"minimap2 -ax asm5 {reference}.fasta {sample}.fasta --cs=long > {sample}.sam")
				command_list.append(f"samtools view -Sb {sample}.sam > {sample}.bam")
				command_list.append(f"samtools index {sample}.bam")
				command_list.append(f"samtools sort {sample}.bam -o {sample}.sorted.bam")
				command_list.append(f"samtools index {sample}.sorted.bam")
				command_list.append(f"samtools mpileup -uf {reference}.fasta {sample}.sorted.bam | bcftools call -Ov -c -v > {sample}.vcf")
				command_list.append(f"rm {sample}.fasta; rm {sample}.sam; rm {sample}.bam; rm {sample}.sorted.bam; rm {sample}.bam.bai; rm {sample}.sorted.bam.bai;")

				for command in command_list:
					subprocess.call(command, shell=True, stdout=subprocess.DEVNULL)


	except FileNotFoundError:
		print(f'{write_dir}/{allele_name}.fasta file cannot be found')
