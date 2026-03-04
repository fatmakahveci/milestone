# Public Benchmark Dataset Sources

This note lists realistic public sources that can be used to build larger benchmark packs for Milestone.

These are **not** all ready-made, publication-grade truth sets by themselves. In most cases, they are public data sources from which a defensible truth-set benchmark can be assembled by selecting isolates with stable public profiles and associated assemblies or genomes.

## Recommended Public Sources

### Neisseria spp. via PubMLST

Why it is useful:

- large public typing database
- large isolate collection
- large genome collection
- dedicated genome libraries

Relevant official sources:

- Organism page: https://pubmlst.org/organisms/neisseria-spp
- Genome libraries: https://pubmlst.org/organisms/neisseria-spp/genome-libraries
- BIGSdb REST API documentation: https://bigsdb.readthedocs.io/en/latest/rest.html

Why it fits Milestone well:

- stable allele and profile nomenclature are available through the BIGSdb/PubMLST API
- genome-backed isolate records can be selected for benchmark construction

### Escherichia spp. via PubMLST

Why it is useful:

- large typing database
- substantial isolate and genome collections

Relevant official source:

- https://pubmlst.org/organisms/escherichia-spp

Important note:

- for some Escherichia workflows, broader hosted context may also be available through EnteroBase

### Acinetobacter baumannii via PubMLST

Why it is useful:

- sizeable genome collection
- clinically important organism with active public nomenclature

Relevant official source:

- https://pubmlst.org/organisms/acinetobacter-baumannii

### Campylobacter jejuni/coli via PubMLST

Why it is useful:

- very large isolate and genome collections
- useful for stress-testing benchmark-pack collection workflows

Relevant official source:

- https://pubmlst.org/organisms/campylobacter-jejunicoli

### EnteroBase scheme downloads

Why it is useful:

- official support for downloading scheme profiles and allele FASTAs
- useful where the benchmark should align with EnteroBase-hosted nomenclature and schemes

Relevant official sources:

- Scheme download docs: https://enterobase.readthedocs.io/en/latest/features/user-download-schemes.html
- API download docs: https://enterobase.readthedocs.io/en/latest/api/api-download-schemes-assemblies.html

## What Counts As A Better Benchmark Pack

A stronger public benchmark pack should contain:

- a provenance manifest with source URLs and retrieval date
- a clearly defined species and schema
- predicted profiles produced by Milestone
- truth profiles from the public source
- ideally, linked assemblies or genomes used to derive the Milestone predictions

## Repository Importers

Milestone now includes importer scripts for two practical workflows:

- [import_pubmlst_benchmark_pack.py](/Users/fatmakhv/Desktop/milestone/workflow/scripts/import_pubmlst_benchmark_pack.py): build a benchmark-pack layout directly from PubMLST/BIGSdb isolate IDs and scheme allele calls
- [import_enterobase_scheme.py](/Users/fatmakhv/Desktop/milestone/workflow/scripts/import_enterobase_scheme.py): download EnteroBase scheme profile archives and metadata

For PubMLST benchmark-pack creation, use:

- a `seqdef` database for scheme definitions
- an `isolates` database for isolate-backed truth profiles

Milestone now includes two freezing layers:

- [freeze_public_benchmark_pack.py](/Users/fatmakhv/Desktop/milestone/workflow/scripts/freeze_public_benchmark_pack.py): download a PubMLST-backed pack and archive it as a frozen ZIP bundle
- [freeze_benchmark_pack_collection.py](/Users/fatmakhv/Desktop/milestone/workflow/scripts/freeze_benchmark_pack_collection.py): archive an entire local benchmark-pack collection into reproducible ZIP files

## Honest Limitation

The sources above are public and suitable for assembling benchmark packs, but they are not automatically “gold-standard truth sets” without curation. A publication-grade validation dataset still requires:

- explicit isolate selection criteria
- frozen retrieval dates or snapshots
- documentation of the exact public profiles used as truth
- confirmation that the same schema definition was used on both sides of the comparison
