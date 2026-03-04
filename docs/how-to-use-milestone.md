# How To Use Milestone

This page gives the shortest path to using Milestone without reading the whole repository.

## When To Use It

Use Milestone when you want to:

- work with a bacterial `wgMLST` schema locally
- call alleles from your own samples
- compare two or more `wgMLST` profiles
- check whether a schema looks usable before analysis
- keep schema, comparison, and validation steps in one workflow

## The Main Workflow

Milestone has three common paths.

### 1. Build or import a schema

Choose this path when you are starting with loci or a public schema.

Use:

- `schema_creation` if you have your own CDS FASTA files
- `import_pubmlst_scheme.py` if you want to import a public schema

Output:

- locus FASTA files
- reference files
- schema metadata and manifest files

### 2. Run allele calling

Choose this path when you already have a schema and want a `wgMLST` profile for a sample.

Use:

- `allele_calling`

Output:

- `*_wgmlst.tsv`
- novel allele exports if new alleles are detected
- optional QC and reporting files

### 3. Compare profiles

Choose this path when you already have finished `wgMLST` profiles.

Use:

- `profile_compare` for two profiles
- `profile_compare_batch` for many pairwise comparisons
- `profile_matrix` for a full matrix
- `profile_benchmark` if you want to compare predicted profiles with a truth set

Output:

- comparison tables
- summary JSON files
- matrix and report files

## Fast Start

### Run the web app

```bash
python -m venv .venv
source .venv/bin/activate
pip install -r webapp/requirements.txt
streamlit run webapp/app.py
```

### Run the CLI

Compare two finished profiles:

```bash
python workflow/milestone.py profile_compare \
  --profile-a output/run_a/sample_a_wgmlst.tsv \
  --profile-b output/run_b/sample_b_wgmlst.tsv \
  --output-dir output/compare/sample_a_vs_b
```

Build a distance matrix from multiple profiles:

```bash
python workflow/milestone.py profile_matrix \
  --profiles output/run_a/sample_a_wgmlst.tsv output/run_b/sample_b_wgmlst.tsv output/run_c/sample_c_wgmlst.tsv \
  --output output/compare/matrix_abc
```

## What To Expect From The Results

Milestone gives you schema-based allele calls and `wgMLST` profile comparisons.

That means it can help you:

- see whether profiles differ across comparable loci
- summarize unresolved or differing loci
- organize validation and reporting outputs

It does not automatically give:

- exact biological strain identity
- phylogeny
- outbreak confirmation on its own

## Which Interface To Use

Use the CLI when:

- you want reproducible scripted runs
- you are processing many samples
- you want full control

Use the web app when:

- you want to demo the workflow
- you want a quicker manual interface
- you want to import public schemas or inspect outputs visually

## Good Follow-Up Reads

- [What Milestone Is](/Users/fatmakhv/Desktop/milestone/docs/scientific-positioning.md)
- [Milestone Compared With Other Tools](/Users/fatmakhv/Desktop/milestone/docs/competitive-analysis.md)
- [Public Benchmark Datasets](/Users/fatmakhv/Desktop/milestone/docs/public-benchmark-datasets.md)
