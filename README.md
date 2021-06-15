<div align="left"> <h1> <img src="images/milestone.png" alt="milestone_logo"> MILESTONE </h1> </div>

## Table of Contents

<!-- MarkdownTOC -->

- Preprocess
    - Downloading NCBI data to create the reference genome for given species
- Install conda environment
- Install VG
- Creation of reference genome
- Creation of MLST schema of the given sample and update reference genome
- Creates sample's mlst without updating reference files
- Creates sample's mlst and updates reference files

<!-- /MarkdownTOC -->

---

Milestone is an end-to-end sample-based cgMLST profile creation workflow for given bacterial species. It only uses available genome assemblies of the species provided by the user or NCBI's public database, and raw reads of the given sample.

---

## 1. Preprocess

### Downloading NCBI data to create the reference genome for given species

- `$ >> bash download_species_reference_fasta.sh -s <species_name>`

---

## Install conda environment

```bash
$ >> conda config --add channels defaults
$ >> conda config --add channels bioconda
$ >> conda config --add channels conda-forge
$ >> conda create --name milestone chewbbaca=2.7 freebayes=1.3 minimap2=2.17 snakemake=5.32 pysam=0.16 bcftools=1.12
```

- Activate the created environment: `$ >> source activate milestone`
- Deactivate the created environment: `$ >> source deactivate milestone`

---

## Install VG

@todo

---

## 2. Milestone

![milestone workflow](images/milestone_workflow_github.png)

---

### STEP 1 and STEP 2

### *milestone.py chewbbaca*: Creation of reference genome

`$ >> python milestone.py chewbbaca --help` for more detailed information.

usage: milestone.py chewbbaca [-h] -d DIRECTORY [-n] [-p] [-s SNAKEFILE] [-t THREADS] [-F] [--ri] [--unlock] [-q] -r REFERENCE -g GENOME_DIR

chewBBACA

optional arguments:
  -h, --help            show this help message and exit
  -d DIRECTORY, --directory DIRECTORY
                        Snakemake - Specify working directory (relative paths in the snakefile will use this as their origin).
  -n, --dryrun, --dry-run
                        Snakemake - Do not execute anything, and display what would be done. If you have a very large workflow, use --dry-run --quiet to just print a summary of the DAG of jobs.
  -p, --printshellcmds  Snakemake - Print out the shell commands that will be executed.
  -s SNAKEFILE, --snakefile SNAKEFILE
                        Snakemake - The workflow definition in form of a snakefile. Usually, you should not need to specify this. By default, Snakemake will search for 'Snakefile','snakefile', 'workflow/Snakefile', 'workflow/snakefile' beneath the current working directory, in this order. Only if you
                        definitely want a different layout, you need to use this parameter.
  -t THREADS, --threads THREADS, --set-threads THREADS
                        Snakemake - Overwrite thread usage of rules. This allows to fine-tune workflow parallelization. In particular, this is helpful to target certain cluster nodes by e.g. shifting a rule to use more, or less threads than defined in the workflow. Thereby, THREADS has to be a positive
                        integer, and RULE has to be the name of the rule.
  -F, --forceall        Snakemake - Force the execution of the selected (or the first) rule and all rules it is dependent on regardless of already created output.
  --ri, --rerun-incomplete
                        Snakemake - Re-run all jobs the output of which is recognized as incomplete.
  --unlock              Snakemake - Remove a lock on the working directory.
  -q, --quiet           Snakemake - Do not output any progress or rule information.
  -r REFERENCE, --reference REFERENCE
                        Name of reference file without extension. (Both VCF and FASTA file name of the reference.)
  -g GENOME_DIR, --genome_dir GENOME_DIR
                        ChewBBACA - Assembled genome directory name to create species MLST schema

`$ >> python milestone.py chewbbaca -d <input_data_directory> -t <number_of_threads> -g <reference_genome_assemblies_directory> -p -r <to_be_created_reference_file_name_without_extension> --snakefile Snakefile -F`

---

### STEP 3 and STEP 4

### *milestone.py mlst*: Creation of MLST schema of the given sample and update reference genome

`$ >> python milestone.py mlst --help` for more detailed information.

usage: milestone.py mlst [-h] -d DIRECTORY [-n] [-p] [-s SNAKEFILE] [-t THREADS] [-F] [--ri] [--unlock] [-q] -r REFERENCE [--aligner ALIGNER] -e READ1 -E READ2 [--ur]

Graph Alignment

optional arguments:
  -h, --help            show this help message and exit
  -d DIRECTORY, --directory DIRECTORY
                        Snakemake - Specify working directory (relative paths in the snakefile will use this as their origin).
  -n, --dryrun, --dry-run
                        Snakemake - Do not execute anything, and display what would be done. If you have a very large workflow, use --dry-run --quiet to just print a summary of the DAG of jobs.
  -p, --printshellcmds  Snakemake - Print out the shell commands that will be executed.
  -s SNAKEFILE, --snakefile SNAKEFILE
                        Snakemake - The workflow definition in form of a snakefile. Usually, you should not need to specify this. By default, Snakemake will search for 'Snakefile','snakefile', 'workflow/Snakefile', 'workflow/snakefile' beneath the current working directory, in this order. Only if you
                        definitely want a different layout, you need to use this parameter.
  -t THREADS, --threads THREADS, --set-threads THREADS
                        Snakemake - Overwrite thread usage of rules. This allows to fine-tune workflow parallelization. In particular, this is helpful to target certain cluster nodes by e.g. shifting a rule to use more, or less threads than defined in the workflow. Thereby, THREADS has to be a positive
                        integer, and RULE has to be the name of the rule.
  -F, --forceall        Snakemake - Force the execution of the selected (or the first) rule and all rules it is dependent on regardless of already created output.
  --ri, --rerun-incomplete
                        Snakemake - Re-run all jobs the output of which is recognized as incomplete.
  --unlock              Snakemake - Remove a lock on the working directory.
  -q, --quiet           Snakemake - Do not output any progress or rule information.
  -r REFERENCE, --reference REFERENCE
                        Name of reference file without extension. (Both VCF and FASTA file name of the reference.)
  --aligner ALIGNER     Graph Alignment - Aligner option, sbg or vg
  -e READ1, --read1 READ1
                        Graph Alignment - Sample first read
  -E READ2, --read2 READ2
                        Graph Alignment - Sample second read
  --ur, --update_reference
                        Reference - Update <reference.fasta> and <reference.vcf> after the alignment of the given sample.

#### 1. Creates sample's mlst without updating reference files

`$ >> python milestone.py mlst -d <input_data_directory> -t <number_of_threads> -p -r <to_be_created_reference_file_name_without_extension> -e <sample_1.fastq> -E <sample_2.fastq> --aligner <vg/sbg> --snakefile Snakefile -F`

#### 2. Creates sample's mlst and updates reference files

`$ >> python milestone.py mlst -d <input_data_directory> -t <number_of_threads> -p -r <to_be_created_reference_file_name_without_extension> -e <sample_1.fastq> -E <sample_2.fastq> --aligner <vg/sbg> --snakefile Snakefile --update_reference -F`

---



### How to convert allele sequence into MLST

![allele to vcf](images/allele_to_vcf_github.png)
