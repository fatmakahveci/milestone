# What Milestone Is

Milestone is a local `wgMLST` workflow for bacterial allele calling, profile comparison, schema quality control, and validation support. It is built for users who want to run and inspect a full schema-based typing workflow in one place instead of stitching together many separate scripts and tools.

Milestone includes:

- schema-based allele calling
- `wgMLST` profile generation
- pairwise and matrix-based profile comparison
- schema quality control
- public schema and benchmark-pack import
- validation corpus construction
- provenance and reporting outputs

In short, Milestone is a workflow for producing, checking, comparing, and documenting `wgMLST` results.

## What Milestone Is Good For

Milestone is useful when you want to:

- build or import a `wgMLST` schema
- call alleles from bacterial samples
- generate locus-by-locus `wgMLST` profiles
- compare two or more profiles
- check whether a schema is usable before analysis
- package results and metadata in a reproducible way

It is especially useful in local research settings, custom analysis pipelines, and projects where reproducibility matters.

## What Milestone Measures

Milestone works at the allele and profile level. It assigns alleles to loci and creates `wgMLST` profiles that can then be compared across comparable loci.

This allows Milestone to:

- report locus-level allele calls
- show where two profiles are the same or different
- summarize unresolved, comparable, and differing loci
- benchmark predicted profiles against truth-set profiles

These outputs are helpful for bacterial typing workflows and for comparing samples under a clearly defined schema.

## What Milestone Does Not Do by Itself

Milestone should not be treated as doing all of the following on its own:

- giving a universal or exact strain identity
- replacing phylogenetic analysis
- proving outbreak linkage
- measuring evolutionary distance directly
- confirming novel alleles biologically without further validation

Those questions need more context, such as species-specific validation, external databases, epidemiological information, phylogenetic analysis, or lab confirmation.

## Main Strengths

### Everything in One Workflow

Milestone brings together several steps that are often spread across different tools: schema handling, allele calling, profile comparison, quality control, benchmarking, and reporting.

### Reproducibility

The workflow includes provenance capture, benchmark packaging, readiness checks, and reporting outputs. That makes it easier to rerun an analysis and understand how a result was produced.

### Local and Flexible Use

Milestone can be run locally and adapted to custom schemas, custom thresholds, and custom validation assets. This makes it practical for exploratory and project-specific work.

### Validation Support

Milestone includes benchmark packs, validation corpus generation, and publication-support outputs. Even when validation is not complete, the workflow already supports organizing and documenting it.

## Current Limits

Milestone still has some important limits.

### Validation Is Not Complete for All Species

The project supports frozen benchmark packs and species-level validation packaging, but it does not yet include a broad, fully curated validation set for every relevant bacterial species.

### External Tool Comparison Is Not Finished

The repository includes paths for comparing Milestone with tools such as `chewBBACA` and `MentaLiST`, but a full same-dataset comparison across all tools is not yet complete in this environment.

### `wgMLST` Results Still Need Careful Interpretation

Milestone can show whether two profiles are distinguishable across comparable loci, but that should not automatically be read as phylogenetic distance or epidemiological relatedness.

## Simple Description

If you need one short way to describe the project, this is accurate:

“Milestone is a local `wgMLST` workflow for bacterial allele calling, profile comparison, validation support, and reproducible reporting.”

These descriptions are too strong and should be avoided:

- “Milestone identifies the exact strain.”
- “Milestone infers phylogeny.”
- “Milestone proves relatedness on its own.”

## Bottom Line

Milestone is useful as a practical, reproducible workflow for schema-based bacterial `wgMLST` analysis. Its main value is that it combines allele calling, profile comparison, schema validation, and reporting support in a single local framework. It is best used as a bacterial typing workflow, not as a standalone biological truth engine.
