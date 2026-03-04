# Milestone Publication Checklist

Use this checklist before submitting a manuscript, thesis chapter, or technical report based on Milestone.

## Required Provenance

- record the Milestone commit or release identifier used for the analysis
- archive the exact `schema_manifest.json`
- archive the exact `schema_qc_summary.json`
- archive the exact benchmark report (`wgmlst_benchmark_summary.json`)
- archive the exact calling parameters, including `min_locus_coverage`, `translation_table`, and any custom start codons

## Required Scientific Reporting

- state clearly whether the schema is custom, imported, or externally curated
- report the species and schema type (`wgmlst` or `cgmlst`)
- report how unresolved loci (`LNF`, `ASM`, `ALM`, and related states) were handled
- state that profile comparison reflects wgMLST discrimination, not direct phylogenetic distance
- state whether novel allele calls were externally validated or only computationally inferred

## Required Validation Reporting

- report the number of matched truth-set samples used in benchmarking
- report mean concordance and mean non-comparable rate
- report any schema QC warnings that remain unresolved
- avoid presenting exact allele-ID agreement as universal biological equivalence across different schemas or tools

## Claims To Avoid

- "Milestone proves phylogenetic relatedness"
- "Milestone universally identifies strains without external validation"
- "Novel allele calls are biologically confirmed" unless independently verified
