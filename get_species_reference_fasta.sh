#!/bin/bash

# author: fatmakhv

while [[ "$#" -gt 0 ]]; do
    case $1 in
        -s|--species) species=$2; shift ;;
    *) echo "Usage: bash get_species_reference_fasta.sh -s <species_name>"
    exit 1;;
    esac
    shift
done

wget https://ftp.ncbi.nlm.nih.gov/genomes/refseq/bacteria/assembly_summary.txt;

species_link=`awk -F'\t' -v v="$species" '{ if ($12 == "Complete Genome" && $11 == "latest" && $5 ~ /representative/ && $8 ~ v) print $20 }' assembly_summary.txt`;

rm assembly_summary.txt;

sample=`echo ${species_link} | awk -F/ '{ print $NF }'`;

wget "${species_link}/${sample}_genomic.fna.gz";
