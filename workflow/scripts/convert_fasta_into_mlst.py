#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
author: @rfm
"""

import argparse, csv, os, re, subprocess

from Bio import SeqIO
from copy import deepcopy
from collections import Counter
from itertools import groupby, chain


def run_minimap2(reference, fastq1, fastq2, output_file, output_type):
	""" Executes minimap2 to map short sequences
		to a reference genome.

		Parameters
		----------
		reference : str
			Path to the FASTA file with the reference.
		map_fasta : str
			Path to FASTA file with the short sequences
			to map against the reference.
		output_file : str
			Path to the output file with mapping results.

		Returns
		-------
		List with following elements:
			stdout : list
				List with stdout from minimap2 in bytes.
			stderr : list
				List with the stderr from minimpa2 in bytes.
	"""

	# -I parameter to control number of target bases loaded into memory
	if output_type == 'paf':
		minimap_args = ['minimap2 -I 100M --cs -cx sr {0} {1} {2} > '
						'{3}'.format(reference, fastq1, fastq2, output_file)]
	elif output_type == 'sam':
		minimap_args = ['minimap2 -I 100M --cs -ax sr {0} {1} {2} > '
						'{3}'.format(reference, fastq1, fastq2, output_file)]

	minimap_proc = subprocess.Popen(minimap_args,
									shell=True,
									stdout=subprocess.PIPE,
									stderr=subprocess.PIPE)

	stderr = minimap_proc.stderr.readlines()
	stdout = minimap_proc.stdout.readlines()

	return [stdout, stderr]


def assign_alleleid(milestone_seqs, chewie_schema):
	""" Determines if the allele predicted by Milestone
		is in the schema and attributes an allele
		identifier.

		Parameters
		----------
		milestone_seqs : dict
			Dictionary with loci identifiers as keys and
			DNA sequences predicted by Milestone as values.
		chewie_schema : str
			Path to the schema's directory.

		Returns
		-------
		milestone_allele_ids : dict
			Dictionary with loci identifiers as keys and
			allele identifiers as values. If an allele
			predicted by Milestone is not in the schema,
			the value will be '?'.
	"""

	milestone_allele_ids = {}
	for locid, seq in milestone_seqs.items():
		locus_file = os.path.join(chewie_schema, '{0}.fasta'.format(locid))
		locus_alleles = {str(rec.seq): (rec.id).split('_')[-1]
						 for rec in SeqIO.parse(locus_file, 'fasta')}

		# get allele ID for allele predicted by Milestone
		# assign '?' if allele is not in Chewie's schema
		milestone_id = locus_alleles.get(seq, '?')
		milestone_allele_ids[locid] = milestone_id

	return milestone_allele_ids


def get_alleles_seqs(loci_alleles, schema_path):
	""" Gets the DNA sequence of one allele for each
		locus identifier that is a key in the input
		dictionary.

		Parameters
		---------
		loci_alleles : dict
			Dictionary with loci ids as keys and alleles
			identifiers as values.
		schema_path : str
			Path to the schema's directory.

		Returns
		-------
		loci_seqs : dict
			Dictionary with loci identifiers as keys and
			DNA sequences as values.
	"""

	loci_seqs = {}
	for locus, alleleid in loci_alleles.items():
		locus_file = os.path.join(schema_path, locus+'.fasta')
		allele_seq = [str(rec.seq)
					  for rec in SeqIO.parse(locus_file, 'fasta')
					  if (rec.id).split('_')[-1] == alleleid]
		if len(allele_seq) > 0:
			loci_seqs[locus] = allele_seq[0]

	return loci_seqs


def main(chewie_matrix, chewie_schema, milestone_results, strain_id,
		 output_mlst, fastq1, fastq2):

	# process chewie's results
	# read AlleleCall matrix
	with open(chewie_matrix, 'r') as infile:
		reader = csv.reader(infile, delimiter='\t')
		chewie_lines = [line for line in reader]

	# get line with loci IDs
	loci = chewie_lines[0]
	# extract line for strain_id
	strain_profile = [l for l in chewie_lines if l[0].split('.fasta')[0] == strain_id][0]

	# create dictionary with chewie classification per locus
	chewie_allele_ids = {loci[i].split('.fasta')[0]: strain_profile[i].replace('INF-', '')
						 for i in range(1, len(loci))}

	# get chewie alleles
	chewie_seqs = get_alleles_seqs(chewie_allele_ids, chewie_schema)

	# read allele sequences predicted by Milestone
	milestone_seqs = {(rec.id).split('_')[0]: str(rec.seq)
					  for rec in SeqIO.parse(milestone_results, 'fasta')}

	# determine allele IDs for alleles predicted by Milestone
	milestone_allele_ids = assign_alleleid(milestone_seqs, chewie_schema)

	# create dictionary with loci IDs as keys and chewie and
	# milestone classifications as values
	merged_classifications = {k: [v, milestone_allele_ids.get(k, '-')]
							  for k, v in chewie_allele_ids.items()}

	# write results to output file
	with open(output_mlst, 'w') as out:
		out.write('locus\tmilestone\n')
		out.write('\n'.join(['{0}\t{1}'.format(k, v1) for k, [v0, v1] in merged_classifications.items()]))


def parse_arguments():

	parser = argparse.ArgumentParser(description=__doc__,
									 formatter_class=argparse.RawDescriptionHelpFormatter)

	parser.add_argument('-cm', type=str, required=True,
						dest='chewie_matrix',
						help='Path to the TSV file with allelic '
							 'profiles determined by chewBBACA.')

	parser.add_argument('-cs', type=str, required=True,
						dest='chewie_schema',
						help='Path to chewies\'s schema.')

	parser.add_argument('-mr', type=str, required=True,
						dest='milestone_results',
						help='Path to a FASTA file with the '
							 'alleles predicted by Milestone.')

	parser.add_argument('-sid', type=str, required=True,
						dest='strain_id',
						help='Identifier of the strain.')

	parser.add_argument('-o', type=str, required=True,
						dest='output_mlst',
						help='Output sample mlst file.')

	parser.add_argument('-fq1', type=str, required=True,
						dest='fastq1',
						help='FASTQ R1 file.')

	parser.add_argument('-fq2', type=str, required=True,
						dest='fastq2',
						help='FASTQ R2 file.')
	args = parser.parse_args()

	return args


if __name__ == '__main__':

	args = parse_arguments()
	main(**vars(args))
