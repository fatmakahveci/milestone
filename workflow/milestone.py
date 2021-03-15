#!/usr/bin/env python3

import argparse, glob, os, subprocess, sys


def create_config():

	output_file = open(configfile, "w")

	output_file.write(f"# milestone.py created {configfile}.\n\n")

	# config files
	output_file.write(f'workdir:  "{args.directory}"\n')
	output_file.write(f'logs: "{args.directory}/logs/"\n')
	output_file.write(f'envs: "{args.directory}/envs/"\n')

	# config parameters
	output_file.write('\n')
	output_file.write(f'parameters: \n')
	output_file.write(f' threads: {args.threads}\n')
	output_file.write(f'reference: {args.reference}\n')

	## chewBBACA run
	if args.command == 'chewbbaca':

		## create wgMLST
		output_file.write(f'training_file: "data/{args.training_file}"\n')
		output_file.write(f'genome_dir: "data/{args.genome_dir}"\n')
		output_file.write(f'schema_seed_dir: "data/schema_seed"\n')

		## allele call
		output_file.write(f'allele_call_dir: "data/allele_call"\n')

		## create cgMLST
		output_file.write(f'cgmlst_dir: "data/create_cgMLST"\n')

	## alignment run
	elif args.command == 'alignment':
		## graph alignment
		output_file.write(f'samples:\n sample1: "{args.read1}"\n sample2: "{args.read2}"\n')
		# reference file name of VCF and FASTA
		output_file.write(f'aligner: "{args.aligner}"\n')
	
	output_file.close()

def run_snakemake():

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

def create_snakefile():

	output_file = open(snakefile, 'w')
	
	output_file.write('import glob, os, sys\n\n')
	output_file.write('configfile: "config.yaml"\n\n')

	if args.command == 'chewbbaca':
		output_file.write('include: "rules/chewie.smk"\n\n')
		output_file.write('rule all:\n\tinput:\n')
		output_file.write(f'\t\treference_vcf = "data/{args.reference}.vcf",')
	
	elif args.command == 'alignment':
		output_file.write('include: "rules/milestone.smk"\n\n')
		output_file.write('rule all:\n\tinput:\n')
		read_name = ''.join(args.read1).split('_1')[0]
		output_file.write(f'\t\tsample_fasta = "data/{args.aligner}/{read_name}.fasta"')

	output_file.close()


if __name__ == "__main__":

	parent_parser = argparse.ArgumentParser(add_help=False)

	# SNAKEMAKE
	parent_parser.add_argument('-d', '--directory', help='Snakemake - Specify working directory (relative paths in the snakefile will use this as their origin).', default="", required=True, type=os.path.abspath)
	parent_parser.add_argument('-n', '--dryrun', '--dry-run', help="Snakemake - Do not execute anything, and display what would be done. If you have a very large workflow, use --dry-run --quiet to just print a summary of the DAG of jobs.", default=False, action='store_true', required=False)
	parent_parser.add_argument('-p', '--printshellcmds', help='Snakemake - Print out the shell commands that will be executed.', default=False, action='store_true', required=False)
	parent_parser.add_argument('-s', '--snakefile', help="Snakemake - The workflow definition in form of a snakefile.Usually, you should not need to specify this. By default, Snakemake will search for 'Snakefile','snakefile', 'workflow/Snakefile', 'workflow/snakefile' beneath the current working directory, in this order. Only if you definitely want a different layout, you need to use this parameter.", default='Snakefile', required=False)
	parent_parser.add_argument('-t', '--threads', '--set-threads', help='Snakemake - Overwrite thread usage of rules. This allows to fine-tune workflow parallelization. In particular, this is helpful to target certain cluster nodes by e.g. shifting a rule to use more, or less threads than defined in the workflow. Thereby, THREADS has to be a positive integer, and RULE has to be the name of the rule.', required=False, type=int)
	parent_parser.add_argument('--report', help='Snakemake - Create an HTML report with results and statistics. This can be either a .html file or a .zip file. In the former case, all results are embedded into the .html (this only works for small data). In the latter case, results are stored along with a file report.html in the zip archive. If no filename is given, an embedded report.html is the default.', action='store_true', required=False)
	parent_parser.add_argument('-c', '--condaprefix', help='Path of default conda environment, enables recycling built environments, default is in folder conda_env in repository directory.', default=os.path.join(os.path.dirname(os.path.realpath(__file__)), "envs"), required=False)
	parent_parser.add_argument('-F', '--forceall', help='Snakemake - Force the execution of the selected (or the first) rule and all rules it is dependent on regardless of already created output.', default=False, action='store_true', required=False)
	parent_parser.add_argument('--ri', '--rerun-incomplete', help='Snakemake - Re-run all jobs the output of which is recognized as incomplete. ', default=False, action='store_true', required=False)
	parent_parser.add_argument('--unlock', help='Snakemake - Remove a lock on the working directory.', default=False, action='store_true', required=False)
	parent_parser.add_argument('-q', '--quiet', help='Snakemake - Do not output any progress or rule information.', default=False, action='store_true', required=False)

	# REFERENCE
	parent_parser.add_argument('-r', '--reference', help='Name of reference file without extension. (Both VCF and FASTA file name of the reference.)', required=True)

	parser = argparse.ArgumentParser(add_help=False)
	subparsers = parser.add_subparsers(title='commands', dest='command')

	# CHEWIE
	chewbbaca_parser = subparsers.add_parser("chewbbaca", parents=[parent_parser], description='chewBBACA', help='ChewBBACA - Run chewBBACA workflow to create FASTA and VCF files for reference genome.')
	chewbbaca_parser.add_argument('-a', '--training_file', help='Prodigal training file name', required=True)
	chewbbaca_parser.add_argument('-g', '--genome_dir', help='Assembled genome directory name to create species MLST schema', required=True)
	
	# GRAPH ALIGNMENT
	alignment_parser = subparsers.add_parser("alignment", parents=[parent_parser], description='Graph Alignment', help='Graph aligner - Choose VG or SBG GRAF aligners to align reads onto the reference genome.')
	alignment_parser.add_argument('--aligner', help='Aligner option, sbg or vg', default='vg', required=False)
	alignment_parser.add_argument('-e', '--read1', help='Sample first read', required=True)
	alignment_parser.add_argument('-E', '--read2', help='Sample second read', required=True)

	# CREATION OF MLST SCHEMA OF SAMPLE
	# Rafael's code...

	# PARSE ARGUMENTS

	args = parser.parse_args()

	if not os.path.exists(args.directory):
		os.makedirs(args.directory)

	configfile = os.path.join(args.directory, "config.yaml")
	create_config()

	snakefile = os.path.join(os.path.dirname(os.path.realpath(__file__)), args.snakefile)
	create_snakefile()

	if os.path.exists(configfile):
		run_snakemake()
	else:
		print("Path to configfile does not exist: " + configfile)
