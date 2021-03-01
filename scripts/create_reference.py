#!/usr/bin/env python3

import argparse, os, shutil
from Bio import SeqIO
from io import StringIO  # python3
from pathlib import Path


def split_cg_sequences(cg_dir, cds_fasta):

	try:
		wd = f"{os.getcwd()}/{cg_dir}"

		for sequence in list(SeqIO.parse(StringIO(open(f"{wd}/{cds_fasta}", 'r').read()), 'fasta')):
			raw_seq_id_fields = str(sequence.id).strip("\n").split("_")
			cds_name = raw_seq_id_fields[0]
			allele_id = raw_seq_id_fields[-1]
			allele_name = f"{raw_seq_id_fields[0]}_{allele_id}"

			wr_dir = "alleles"
			if allele_id == '1':
				wr_dir = "references"
				seq_name = raw_seq_id_fields[0]
				allele_dir = f"{wd}/alleles/{cds_name}"
				Path(allele_dir).mkdir(exist_ok=True)
			else:
				wr_dir = f"alleles/{cds_name}"
				seq_name = f"{allele_name}/{cds_name}"
			
			out_file_name = f"{wd}/{wr_dir}/{allele_name}.fasta"
			
			with open(out_file_name, "w") as out_file:
				out_file.write(f">{seq_name}\n{str(sequence.seq)}")
				out_file.close()

	except FileNotFoundError:
		print(f'{cg_dir}/{cds_fasta} file cannot be found')