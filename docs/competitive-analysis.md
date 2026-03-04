# Milestone Compared With Other Tools

This page explains where Milestone fits next to other `wgMLST`, `cgMLST`, and gene-by-gene typing tools. The goal is simple: help users understand when Milestone is a good choice and when another tool may be a better fit.

## What Milestone Covers

Milestone currently includes:

- local schema-based allele calling
- `wgMLST` profile generation
- pairwise and matrix-based profile comparison
- schema quality control
- public schema import from sources such as PubMLST/BIGSdb
- benchmark-pack and validation support
- a lightweight web interface for demos and operational use

Milestone is strongest as a local, modifiable workflow that keeps multiple analysis steps in one repository.

## Where Other Tools Are Stronger

### chewBBACA

`chewBBACA` is usually the stronger choice when you need a more mature gene-by-gene ecosystem.

Use `chewBBACA` when:

- you want a widely used schema-creation and allele-calling environment
- schema lifecycle management is a priority
- interoperability with established community workflows matters

Milestone differs by keeping schema import, schema QC, comparison, and web access together in one place. It is easier to modify when the goal is method development or custom workflow design.

### MentaLiST

`MentaLiST` is stronger when speed and scale matter most.

Use `MentaLiST` when:

- you need fast calling on very large `cgMLST` or `wgMLST` schemes
- you are processing many isolates
- lower runtime and memory use matter more than workflow flexibility

Milestone differs by being easier to inspect and adapt when you want to change calling logic, QC behavior, or validation steps.

### ARIBA

`ARIBA` is stronger for targeted locus analysis and read-level inspection.

Use `ARIBA` when:

- the main task is targeted locus reconstruction
- read-level QC is central to your workflow
- you are doing focused gene analysis rather than broad schema-based profiling

Milestone differs by focusing on a wider schema-based allele-profile workflow rather than a smaller targeted panel workflow.

### EnteroBase

`EnteroBase` is stronger when you need a hosted public comparison ecosystem.

Use `EnteroBase` when:

- you want access to curated schemes and large historical isolate collections
- surveillance context matters more than local pipeline control
- you need a broader comparison universe beyond your own local dataset

Milestone differs by running fully locally and supporting custom workflows without depending on a hosted platform.

### Ridom SeqSphere+ / MBioSeq Typer

Ridom-style platforms are stronger for routine lab and institutional use.

Use them when:

- you need a stable GUI-first laboratory workflow
- permissions, audit trails, and controlled reporting are important
- the software will be used by a team, lab, or surveillance unit

Milestone differs by being open, scriptable, and easier to rework for experimental or custom research use.

### refMLST

`refMLST` is stronger when there is no strong curated schema for the species you care about.

Use `refMLST` when:

- schema availability is limited
- you want to reduce dependence on a traditional curated allele schema

Milestone differs by supporting custom and imported schemas, but it still works within a schema-based allele-calling model.

### pyMLST

`pyMLST` is stronger when you want a smaller and simpler local setup.

Use `pyMLST` when:

- you want lighter local database handling
- you need a simpler environment for clonal or typing work

Milestone differs by offering more layers in the same repository, including schema QC, profile comparison, validation support, and web execution.

## Where Milestone Is Strong

Milestone is a good fit when:

- you want a fully local workflow
- you need to work with custom loci or custom schemas
- you want schema import, schema QC, allele calling, and profile comparison together
- you want a workflow that is easier to inspect and modify than a larger platform
- you are doing method development, exploratory analysis, or project-specific bacterial typing work

## Where Milestone Is Still Behind

Milestone is not yet as strong as the most established tools in a few areas:

- broad species-by-species validation
- public hosted nomenclature and surveillance context
- enterprise-style multi-user operation
- large-scale performance characterization
- long-term curated schema maintenance

These are important gaps, especially for users who need production surveillance or institution-scale deployment.

## Simple Takeaway

Choose Milestone when you want:

- local control
- transparency
- customizability
- one workflow that covers multiple `wgMLST` tasks

Choose another tool when you need:

- a mature public ecosystem
- very high throughput
- heavy multi-user lab operations
- long-established schema curation infrastructure

## Reference Links

- chewBBACA documentation: https://chewbbaca.readthedocs.io/
- chewBBACA publication: https://www.microbiologyresearch.org/content/journal/mgen/10.1099/mgen.0.000166
- MentaLiST publication: https://www.microbiologyresearch.org/content/journal/mgen/10.1099/mgen.0.000146
- ARIBA repository and documentation: https://github.com/sanger-pathogens/ariba
- EnteroBase documentation: https://enterobase.readthedocs.io/en/latest/
- EnteroBase scheme browser example: https://enterobase.warwick.ac.uk/schemes/Escherichia.wgMLST/
- Ridom SeqSphere+ overview: https://www.ridom.de/seqsphere/
- refMLST publication: https://bmcbioinformatics.biomedcentral.com/articles/10.1186/s12859-024-05913-4
- pyMLST repository: https://github.com/bvalot/pyMLST
- pyMLST publication: https://pmc.ncbi.nlm.nih.gov/articles/PMC10711306/
