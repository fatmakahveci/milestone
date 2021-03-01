#!/usr/bin/env python3

import argparse
import glob
import os
import subprocess
import sys


def create_config(configfile, args):

	output_file = open(configfile, "w")

	output_file.write(f"# milestoneSnake.py created {configfile}.\n\n")
	
	# config files	
	output_file.write(f'workdir:  "{args.working_directory}"\n')
	
	## create wgMLST
	output_file.write(f'training_file: "data/{args.training_file}"\n')
	output_file.write(f'genome_dir: "data/{args.genome_dir}"\n')
	output_file.write(f'schema_seed_dir: "data/schema_seed"\n')

	## allele call
	output_file.write(f'allele_call_dir: "data/allele_call"\n')

	## create cgMLST
	output_file.write(f'cgmlst_dir: "data/create_cgMLST"\n')
	
	# config parameters
	output_file.write('\n')
	output_file.write(f'parameters: \n')
	output_file.write(f'  threads: {args.threads}\n')

	output_file.close()

def run_snakemake(configfile, args, snakefile):

	dryrun = ""
	if args.dryrun:
		dryrun = "--dryrun"

	forceall = ""
	if args.forceall:
		forceall = "--forceall"

	printshellcmds = ""
	if args.printshellcmds:
		printshellcmds = "--printshellcmds"

	ri = ""
	if args.ri:
		ri = "--rerun-incomplete"

	unlock = ""
	if args.unlock:
		unlock="--unlock"

	cmd = (f'snakemake -p --conda-prefix {args.condaprefix} --use-conda --configfile "{configfile}" --cores {args.threads} --snakefile {args.snakefile} {forceall} {dryrun} {printshellcmds} {ri} {unlock}')
	os.system(cmd)

def main():

	parser = argparse.ArgumentParser()
	
	# SNAKEMAKE
	parser.add_argument('-d', '--working_directory', help='Working directory to keep results.', required=True, type=os.path.abspath)
	parser.add_argument('-n', '--dryrun', help='Do not execute anything, and display what would be done.', default=False, action='store_true', required=False)
	parser.add_argument('-p', '--printshellcmds', help='Print out the shell commands that will be executed.', default=False, action='store_true', required=False)
	parser.add_argument('-s', '--snakefile', help='Snake file name', default='Snakefile', required=False)
	parser.add_argument('-t', '--threads', help='Number of threads to use', required=False, type=int)
	parser.add_argument('-r', '--report', help='Create html report', default=False, action='store_true', required=False)
	parser.add_argument('-c', '--condaprefix', help='Path of default conda environment, enables recycling built environments, default is in folder conda_env in repository directory.', default=os.path.join(os.path.dirname(os.path.realpath(__file__)), "conda_env"), required=False)
	parser.add_argument('-f', '--forceall', help='Snakemake force recalculation of all steps', default=False, action='store_true', required=False)
	parser.add_argument('--ri', '--rerun-incomplete', help='Snakemake rerun-incomplete', default=False, action='store_true', required=False)
	parser.add_argument('--unlock', help='Snakemake unlock', default=False, action='store_true', required=False)

	# CHEWIE
	parser.add_argument('-a', '--training_file', help='Prodigal training file name', required=True)
	parser.add_argument('-g', '--genome_dir', help='Assembled genome directory name to create species MLST schema', required=True)
	
	# REFERENCE GENOME CREATION
	

	# GRAPH ALIGNMENT
	# parser.add_argument('--aligner', help='Aligner option, sbg or vg (default: sbg)', default='sbg', required=False)
	
	# CREATION OF MLST SCHEMA OF SAMPLE

	args = parser.parse_args()

	if not os.path.exists(args.working_directory):
		os.makedirs(args.working_directory)

	configfile = os.path.join(args.working_directory, "config.yaml")

	create_config(configfile, args)

	snakefile = os.path.join(os.path.dirname(os.path.realpath(__file__)), "Snakefile")

	if os.path.exists(configfile):
		run_snakemake(configfile, args, snakefile)
	else:
		print("Path to configfile does not exist: " + configfile)


if __name__ == "__main__":
	main()