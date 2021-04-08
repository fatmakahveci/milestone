#!/usr/bin/env python3

import argparse, glob, os, subprocess, sys


def create_config():

    output_file = open(configfile, "w")

    output_file.write(f"# milestone.py created {configfile}.\n\n")

    # config files
    output_file.write(f'workdir:  "{args.directory}"\n')
    output_file.write(f'logs: "{args.directory}/logs"\n')
    output_file.write(f'data_dir: "data"\n')

    # config parameters
    output_file.write('\n')
    output_file.write(f'parameters: \n')
    output_file.write(f' threads: {args.threads}\n')
    output_file.write(f'reference: {args.reference}\n')
    output_file.write(f'sample: {"".join({args.read1}).split("_1")[0]}\n')

    ## chewBBACA run
    if args.command == 'chewbbaca':

        ## create wgMLST
        output_file.write(f'genome_dir: "data/{args.genome_dir}"\n')

        ## create cgMLST
        output_file.write(f'cgmlst_dir: "data/create_cgMLST"\n')

    ## Run to create <sample.mlst> or
    ## [<reference.updated.fasta> and <reference.updated.vcf>]
    elif args.command == 'mlst':

        ## graph alignment
        output_file.write(f'samples:\n sample1: "{args.read1}"\n sample2: "{args.read2}"\n')

        ## name of aligner vg or sbg
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

    cmd = (f'snakemake -p --conda-prefix milestone --use-conda --configfile "{configfile}" --cores {args.threads} --snakefile {args.snakefile} {forceall} {dryrun} {printshellcmds} {ri} {unlock}')
    os.system(cmd)


def create_snakefile():

    output_file = open(snakefile, 'w')

    output_file.write('import glob, os, sys\n\n')
    output_file.write('configfile: "config.yaml"\n\n')

    if args.command == 'chewbbaca':

        # create a log file for chewbbaca
        os.system(f'touch {args.directory}/logs/chewbbaca.log')

        output_file.write('include: "rules/chewie.smk"\n\n')
        output_file.write('rule all:\n\tinput:\n')
        output_file.write(f'\t\treference_vcf = "data/{args.reference}.vcf"')

    elif args.command == 'mlst':

        # create a log file for chewbbaca
        os.system(f'touch {args.directory}/logs/mlst.log')

        output_file.write('include: "rules/milestone.smk"\n\n')
        output_file.write('rule all:\n\tinput:\n')

        sample_name = ''.join(args.read1).split('_1')[0]

        if args.sample_mlst:
            output_file.write(f'\t\tsample_mlst = "data/{args.aligner}/{sample_name}.tsv"')
        elif args.update_reference:
            output_file.write(f'\t\treference_vcf = "data/{args.reference}_updated.vcf"')
            output_file.write(f'\t\treference_vcf = "data/{args.reference}_updated.fasta"')
        else:
            output_file.write(f'\t\tsample_fasta = "data/{args.aligner}/{sample_name}.fasta"')

    output_file.close()

def parse_arguments():

    parent_parser = argparse.ArgumentParser(add_help=False)

    # SNAKEMAKE PARAMETERS
    parent_parser.add_argument('-d', '--directory',
        help='Snakemake - Specify working directory (relative paths in the\
             snakefile will use this as their origin).',
        default="", required=True, type=os.path.abspath)
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
        required=False, type=int)
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
        help='Name of reference file without extension.\
             (Both VCF and FASTA file name of the reference.)',
        required=True)

    parser = argparse.ArgumentParser(add_help=True)

    subparsers = parser.add_subparsers(title='commands', dest='command')

    # milestone.py chewbbaca mode PARAMETERS
    chewbbaca_parser = subparsers.add_parser("chewbbaca",
        parents=[parent_parser], description='chewBBACA',
        help='ChewBBACA - Run chewBBACA workflow to create FASTA and VCF files\
             for reference genome.')

    chewbbaca_parser.add_argument('-g', '--genome_dir',
        help='ChewBBACA - Assembled genome directory name to create species\
             MLST schema',
        required=True)

    # milestone.py mlst mode PARAMETERS
    mlst_parser = subparsers.add_parser("mlst", parents=[parent_parser],
        description='Graph Alignment',
        help='Graph Alignment - Choose VG or SBG GRAF aligners to align reads\
             onto the reference genome.')

    mlst_parser.add_argument('--aligner',
        help='Graph Alignment - Aligner option, sbg or vg',
        default='vg', required=False)
    mlst_parser.add_argument('-e', '--read1',
        help='Graph Alignment - Sample first read',
        required=True)
    mlst_parser.add_argument('-E', '--read2',
        help='Graph Alignment - Sample second read',
        required=True)
    mlst_parser.add_argument('--ur', '--update_reference',
        help='Reference - Update <reference.fasta> and <reference.vcf> after\
             the alignment of the given sample.',
        dest='update_reference', default=False, required=False)
    mlst_parser.add_argument('--sm', '--sample_mlst',
        help='Sample MLST - Create MLST profile of the given sample.',
        dest='sample_mlst', default=False, action='store_true', required=False)

    args = parser.parse_args()

    return args


if __name__ == "__main__":

    args = parse_arguments()

    # Create directory of the log files
    if not os.path.exists(f'{args.directory}/logs'):
        os.makedirs(f'{args.directory}/logs')

    # Create config file in YAML format for snakemake
    configfile = os.path.join(args.directory, "config.yaml")
    create_config()

    # Create snakefile for snakemake
    snakefile = os.path.join(os.path.dirname(os.path.realpath(__file__)),
        args.snakefile)
    create_snakefile()

    # Check the existence of config file and run snakemake
    if os.path.exists(configfile):
        run_snakemake()

    else:
        print("Path to configfile does not exist: " + configfile)

