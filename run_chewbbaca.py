#!/usr/bin/env python3

"""
create reference files: python cg_to_vcf.py --cglist ./create_baumannii_cgMLST/cgMLSTschema.txt --fastadir ./baumannii_schema_seed/

Usage:
    run_chewbbaca.py [--taxid=<taxid>] [--training=<file>] [--fasta=<species_ref_fasta>] [--cpu=<no_cpu>] [--genome_samples_dir=<genome_samples_fasta_dir>] [--samples_dir=<samples_fasta_dir>]
    run_chewbbaca.py --help

Options:
    --help                                              Show this help message and exit.
    --taxid=<taxid>                                     Download species complete genome from NCBI if it doesn't provided by user [default: ].
    --training=<file>                                   Prodigal output training file [default: ].
    --fasta=<species_ref_fasta>                         Species gold standard reference FASTA file for training file creation [default: ].
    --cpu=<no_cpu>                                      Number of CPUs for chewBBACA running [default: 4].
    --genome_samples_dir=<genome_samples_fasta_dir>     Samples' FASTA files [default: ].
    --samples_dir=<samples_fasta_dir>                   Samples' FASTA files to be analyzed [default: ].
"""

from docopt import docopt
from pathlib import Path

import datetime as dt, glob, os, shutil, subprocess


# Search taxonomy id: https://www.ncbi.nlm.nih.gov/Taxonomy/Browser/wwwtax.cgi
# taxonomy id -> A.baumannii: 470, S.agalactiae: 1311
def download_species_data_from_ncbi():

	command = "pip install ncbi-genome-download"
	subprocess.call(command, shell = True, stdout=subprocess.DEVNULL)

	command = "ncbi-genome-download -s refseq -l complete -T %s -F fasta bacteria" % taxid
	subprocess.call(command, shell = True, stdout=subprocess.DEVNULL)

	Path("./complete_genome/").mkdir(exist_ok=True)

	os.chdir("./refseq/bacteria/")

	for genome_gz_dir in glob.glob('*'):
	
		os.chdir(genome_gz_dir)
		
		if "gz" in glob.glob("*"):
			os.system("gunzip *gz")

		if "fna" in glob.glob("*") or "fasta" in glob.glob("*") or "fa" in glob.glob("*"):
			os.system("mv *a ../../../complete_genome")

		os.chdir("..")

	os.chdir("../..")

	try:
		shutil.rmtree("./refseq")
	except:
		print('Error while deleting directory')

# create training file
def create_training_file(training_file):

	command = "prodigal -i %s -t %s -p single" % (args["--fasta"], training_file)
	os.system(command)

# create schema
def create_schema():

	command = "chewBBACA.py CreateSchema -i %s --cpu %s -o ./schema_seed --ptf %s" % (genome_samples_dir, args["--cpu"], training_file)
	os.system(command)

# allele calling
def call_allele():

	command = "chewBBACA.py AlleleCall -i %s -g ./schema_seed/ -o ./allele_call --cpu %s --ptf %s" % (genome_samples_dir, args["--cpu"], training_file)
	os.system(command)

# create core genome MLST schema
def create_cgMLST():
	
	results_alleles_tsv = (glob.glob("./allele_call/result*/results_alleles.tsv"))[0]
	repeated_loci = (glob.glob("./allele_call/result*/RepeatedLoci.txt"))[0]
	
	command = "chewBBACA.py ExtractCgMLST -i %s -r %s -p 0.95 -o ./create_cgMLST" % (results_alleles_tsv, repeated_loci)
	os.system(command)

def move_cg_fasta_files():

	Path("./cg_schema_seed/").mkdir(exist_ok=True)

	with open("./create_cgMLST/cgMLSTschema.txt", 'r') as file:

		for line in file.readlines():

			command = "mv %s %s" % ("./schema_seed/"+line.strip(), "./cg_schema_seed/")
			os.system(command)

		file.close()

# create core genome MLST schema for given samples
def create_cgMLST_for_given_samples():

	command = "chewBBACA.py PrepExternalSchema -i ./cg_schema_seed/ -o ./prep_cg_schema_seed/ --cpu 4"
	os.system(command)

	command = "chewBBACA.py AlleleCall -i %s -g ./prep_cg_schema_seed/ -o ./samples_allele_call/ --cpu %s --ptf %s" % (samples_dir, args["--cpu"], training_file)
	os.system(command)

	results_alleles_tsv = (glob.glob("./samples_allele_call/result*/results_alleles.tsv"))[0]
	repeated_loci = (glob.glob("./samples_allele_call/result*/RepeatedLoci.txt"))[0]

	command = "chewBBACA.py ExtractCgMLST -i %s -r %s -p 0.95 -o ./samples_create_cgMLST" % (results_alleles_tsv, repeated_loci)
	os.system(command)

def run_c():

	start_date = dt.datetime.now()
	start_date_str = dt.datetime.strftime(start_date, '%H:%M:%S-%d/%m/%Y')
	print('\nStarted at: {0}\n'.format(start_date_str))

	if not taxid == '':
		print("Complete genome for the species of interest is being downloaded from NCBI database...")
		download_species_data_from_ncbi()

	if args["--training"] == '':
		training_file = 'data.trn'
		print("Prodigal training data is being created...")
		create_training_file(training_file)

	print("chewBBACA schema creation...")
	create_schema()
	print("chewBBACA allele calling...")
	call_allele()
	print("chewBBACA core genome MLST schema creation...")
	create_cgMLST()
	move_cg_fasta_files()
	print("chewBBACA core genome MLST schema creation for given samples...")
	create_cgMLST_for_given_samples()
	
	end_date = dt.datetime.now()
	end_date_str = dt.datetime.strftime(end_date, '%H:%M:%S-%d/%m/%Y')
	
	delta = end_date - start_date
	minutes, seconds = divmod(delta.total_seconds(), 60)


if __name__ == "__main__":
	
	# handling with parameters
	args = docopt(__doc__, version="0.6.2")

	taxid = args["--taxid"]
	training_file = args["--training"]
	samples_dir = args["--samples_dir"]
	genome_samples_dir = args["--genome_samples_dir"]

	run_c()
