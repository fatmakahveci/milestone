#!/usr/bin/env python3

'''
----------------------------------------------
Aim: Creation of config file and milestone run
----------------------------------------------
Authors: @fatmakhv
The latest update: June 14, 2021
----------------------------------------------
'''

import argparse, glob, os, subprocess, sys


def create_config():

    output_file = open(configfile, "w")

    output_file.write(f"# milestone.py created {configfile}.\n\n")

    # config files
    output_file.write(f'logs: "{log_dir}"\n')

    # wd for code
    output_file.write(f'working_dir: "{code_dir}"\n')

    # config parameters
    output_file.write('\n')
    output_file.write(f'parameters: \n')
    output_file.write(f' threads: {args.threads}\n')
    output_file.write(f'reference: {os.path.join(args.output, f"{args.reference}")}\n')
    output_file.write(f'reference_vcf: {os.path.join(args.output, f"{args.reference}.vcf")}\n')
    output_file.write(f'reference_info_txt: {os.path.join(args.output, f"{args.reference}_info.txt")}\n')
    output_file.write(f'reference_fasta: {os.path.join(args.output, f"{args.reference}.fasta")}\n')

    ## chewBBACA run
    if args.command == 'chewbbaca':

        ## create wgMLST
        output_file.write(f'genome_dir: "{args.genome_dir}"\n')

        ## create cgMLST
        output_file.write(f'cgmlst_dir: "{os.path.join(args.output, "create_cgMLST")}"\n')

        ## allele call
        output_file.write(f'allele_call_dir: "{os.path.join(args.output, "allele_call")}"\n')

        ## alleles
        output_file.write(f'alleles_dir: "{os.path.join(args.output, "alleles")}"\n')

        ## schema seed
        output_file.write(f'schema_seed_dir: "{os.path.join(args.output, "alleles/schema_seed")}"\n')

        ## chewbbaca log file
        output_file.write(f'chewbbaca_log_file: "{chewbbaca_log_file}"\n')

    ## Run to create <sample.mlst> or
    ## [<reference.updated.fasta> and <reference.updated.vcf>]
    elif args.command == 'mlst':

        output_file.write(f'sample: "{args.read1.split("/")[-1].split("_1")[0]}"\n')

        ## graph alignment
        output_file.write(f'samples:\n sample1: "{args.read1}"\n sample2: "{args.read2}"\n')

        ## name of aligner vg or sbg
        output_file.write(f'aligner: "{args.aligner}"\n')
        output_file.write(f'output_dir: {args.output}\n')
        output_file.write(f'aligner_reference: {os.path.join(args.output, f"{args.aligner}/{args.reference}")}\n')

        ## update reference
        output_file.write(f'update_reference: "{os.path.join(args.output, str(args.update_reference))}"\n')

        ## chewbbaca log file
        output_file.write(f'mlst_log_file: "{mlst_log_file}"\n')

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

    cmd = (f'snakemake -p --use-conda --configfile "{configfile}" --cores {args.threads} --snakefile {snakefile} {forceall} {dryrun} {printshellcmds} {ri} {unlock}')

    os.system(cmd)


def create_snakefile():

    output_file = open(snakefile, 'w')

    output_file.write('import glob, os, sys\n\n')
    output_file.write(f'configfile: "{configfile}"\n\n')

    if args.command == 'chewbbaca':

        # create a log file for chewbbaca
        os.system(f'touch {chewbbaca_log_file}')

        output_file.write(f'include: "{rules_dir}/chewie.smk"\n\n')
        output_file.write('rule all:\n\tinput:\n')
        # <reference.fasta> and <reference_info.txt> can be added, but not obligatory.
        output_file.write(f'\t\treference_vcf = "{args.output}/{args.reference}.vcf"')

    elif args.command == 'mlst':

        # create a log file for chewbbaca
        os.system(f'touch {mlst_log_file}')

        output_file.write(f'include: "{rules_dir}/milestone.smk"\n\n')
        output_file.write('rule all:\n\tinput:\n')

        sample_name = args.read1.split("/")[-1].split("_1")[0]
        output_file.write(f'\t\tsample_mlst = "{args.output}/{args.aligner}/{sample_name}_mlst.tsv"\n')

    output_file.close()


def parse_arguments():

    parent_parser = argparse.ArgumentParser(add_help=False)

    ########################################
    # SNAKEMAKE PARAMETERS
    parent_parser.add_argument('-n', '--dryrun', '--dry-run',
        help="Snakemake - Do not execute anything, and display what would be\
             done. If you have a very large workflow, use --dry-run --quiet to\
             just print a summary of the DAG of jobs.",
        default=False, action='store_true', required=False)

    parent_parser.add_argument('-p', '--printshellcmds',
        help='Snakemake - Print out the shell commands that will be executed.',
        default=False, action='store_true', required=False)

    parent_parser.add_argument('-s', '--snakefile',
        help="Snakemake - The workflow definition in form of a snakefile.\
             Usually, you should not need to specify this. By default,\
             Snakemake will search for 'Snakefile','snakefile',\
             'workflow/Snakefile', 'workflow/snakefile' beneath the current\
             working directory, in this order. Only if you definitely want a\
             different layout, you need to use this parameter.",
        default='Snakefile', required=False)

    parent_parser.add_argument('-t', '--threads', '--set-threads',
        help='Snakemake - Overwrite thread usage of rules. This allows to\
             fine-tune workflow parallelization. In particular, this is\
             helpful to target certain cluster nodes by e.g. shifting a rule\
             to use more, or less threads than defined in the workflow.\
             Thereby, THREADS has to be a positive integer, and RULE has to be\
             the name of the rule.',
        required=False, type=int, default=1)

    parent_parser.add_argument('-F', '--forceall',
        help='Snakemake - Force the execution of the selected (or the first)\
             rule and all rules it is dependent on regardless of already\
             created output.',
        default=False, action='store_true', required=False)

    parent_parser.add_argument('--ri', '--rerun-incomplete',
        help='Snakemake - Re-run all jobs the output of which is recognized\
             as incomplete. ',
        default=False, action='store_true', required=False)

    parent_parser.add_argument('--unlock',
        help='Snakemake - Remove a lock on the working directory.',
        default=False, action='store_true', required=False)

    parent_parser.add_argument('-q', '--quiet',
        help='Snakemake - Do not output any progress or rule information.',
        default=False, action='store_true', required=False)

    # REFERENCE FILES PARAMETER
    parent_parser.add_argument('-r', '--reference',
        help='Name of reference file to be given without extension and directory.\
             (Both VCF and FASTA file name of the reference.)',
        required=True)

    # OUTPUT DIRECTORY
    parent_parser.add_argument('-o', '--output', 
        help='Directory to be created for the output files', required=True)
    ########################################

    parser = argparse.ArgumentParser(add_help=True)

    subparsers = parser.add_subparsers(title='commands', dest='command')

    ########################################
    # milestone.py chewbbaca mode PARAMETERS
    chewbbaca_parser = subparsers.add_parser("chewbbaca",
        parents=[parent_parser], description='chewBBACA',
        help='ChewBBACA - Run chewBBACA workflow to create FASTA and VCF files\
             for reference genome.')

    chewbbaca_parser.add_argument('-g', '--genome_dir',
        help='ChewBBACA - Assembled genome directory name to create species\
             MLST schema',
        required=True)
    ########################################

    ########################################
    # milestone.py mlst mode PARAMETERS

    mlst_parser = subparsers.add_parser("mlst", parents=[parent_parser],
        description='Graph Alignment',
        help='Graph Aligner - Choose VG or SBG GRAF aligners to align reads\
             onto the reference genome.')

    mlst_parser.add_argument('--aligner',
        help='Graph Aligner option, sbg or vg',
        default='vg', required=False)

    mlst_parser.add_argument('-e', '--read1', type=str,
        help='Graph Aligner - Sample first read including its directory',
        required=True)

    mlst_parser.add_argument('-E', '--read2', type=str,
        help='Graph Aligner - Sample second read including its directory',
        required=True)

    mlst_parser.add_argument('--ur', '--update_reference',
        help='Reference genome - Update <reference_info.txt> and <reference.vcf> after\
             the alignment of the given sample.',
        dest='update_reference', default=False, action='store_true', required=False)
    ########################################

    args = parser.parse_args()

    return args


if __name__ == "__main__":

    args = parse_arguments()

    code_dir = os.path.dirname(os.path.realpath(__file__))
    rules_dir = os.path.join(code_dir, 'rules')

    # Create directory of the log files
    log_dir = os.path.join(args.output, 'logs')
    chewbbaca_log_file = f'{log_dir}/chewbbaca.log'
    mlst_log_file = f'{log_dir}/mlst.log'
    if not os.path.exists(f'{log_dir}'):
        os.makedirs(log_dir)

    # Create config file in YAML format for snakemake
    configfile = os.path.join(code_dir, 'config.yaml')
    create_config()

    # Create snakefile for snakemake
    snakefile = os.path.join(code_dir, args.snakefile)
    create_snakefile()

    # Check the existence of config file and run snakemake
    if os.path.exists(configfile):
        run_snakemake()

    else:
        print("Path to configfile does not exist: " + configfile)

