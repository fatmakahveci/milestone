#!/usr/bin/env python3

import argparse, os


def create_config(configfile, args):

	output_file = open(configfile, "w")

	output_file.write(f"# milestoneSnake.py created {configfile}.\n\n")

	# config files
	output_file.write(f'workdir:  "{args.directory}"\n')

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
	output_file.write(f' threads: {args.threads}\n')

	output_file.write(f'reference: {args.reference}\n')

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

	cmd = (f'snakemake -p --conda-prefix {args.condaprefix} --use-conda --configfile "{configfile}" --cores {args.threads} {forceall} {dryrun} {printshellcmds} {ri} {unlock}')
	os.system(cmd)


def main():

	parser = argparse.ArgumentParser()

	# SNAKEMAKE
	parser.add_argument('-d', '--directory', help='Snakemake - Specify working directory (relative paths in the snakefile will use this as their origin).', default="", required=True, type=os.path.abspath)
	parser.add_argument('-n', '--dryrun', '--dry-run', help="Snakemake - Do not execute anything, and display what would be done. If you have a very large workflow, use --dry-run --quiet to just print a summary of the DAG of jobs.", default=False, action='store_true', required=False)
	parser.add_argument('-p', '--printshellcmds', help='Snakemake - Print out the shell commands that will be executed.', default=False, action='store_true', required=False)
	parser.add_argument('-t', '--threads', '--set-threads', help='Snakemake - Overwrite thread usage of rules. This allows to fine-tune workflow parallelization. In particular, this is helpful to target certain cluster nodes by e.g. shifting a rule to use more, or less threads than defined in the workflow. Thereby, THREADS has to be a positive integer, and RULE has to be the name of the rule.', required=False, type=int)
	parser.add_argument('--report', help='Snakemake - Create an HTML report with results and statistics. This can be either a .html file or a .zip file. In the former case, all results are embedded into the .html (this only works for small data). In the latter case, results are stored along with a file report.html in the zip archive. If no filename is given, an embedded report.html is the default.', action='store_true', required=False)
	parser.add_argument('-c', '--condaprefix', help='Path of default conda environment, enables recycling built environments, default is in folder conda_env in repository directory.', default=os.path.join(os.path.dirname(os.path.realpath(__file__)), "conda_env"), required=False)
	parser.add_argument('-F', '--forceall', help='Snakemake - Force the execution of the selected (or the first) rule and all rules it is dependent on regardless of already created output.', default=False, action='store_true', required=False)
	parser.add_argument('--ri', '--rerun-incomplete', help='Snakemake - Re-run all jobs the output of which is recognized as incomplete. ', default=False, action='store_true', required=False)
	parser.add_argument('--unlock', help='Snakemake - Remove a lock on the working directory.', default=False, action='store_true', required=False)
	parser.add_argument('-q', '--quiet', help='Snakemake - Do not output any progress or rule information.', default=False, action='store_true', required=False)

	# CHEWIE
	parser.add_argument('-a', '--training_file', help='Prodigal training file name', required=True)
	parser.add_argument('-g', '--genome_dir', help='Assembled genome directory name to create species MLST schema', required=True)
	parser.add_argument('-r', '--reference', help='Name of reference file without extension. (Both VCF and FASTA file name of the reference.)', required=True)

	# CREATION OF MLST SCHEMA OF SAMPLE

	args = parser.parse_args()

	if not os.path.exists(args.directory):
		os.makedirs(args.directory)

	configfile = os.path.join(args.directory, "chewie_config.yaml")

	create_config(configfile, args)

	snakefile = os.path.join(os.path.dirname(os.path.realpath(__file__)), "Snakefile")

	if os.path.exists(configfile):
		run_snakemake(configfile, args, snakefile)
	else:
		print("Path to configfile does not exist: " + configfile)


if __name__ == "__main__":
	main()
