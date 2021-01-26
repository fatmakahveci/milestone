#!/usr/bin/env python3

"""
Usage:
    create_reference_files.py (--cglist=<core_genome_list>) (--fastadir=<allele_fasta_files_dir>) [--cg_mlst_tsv=<tsv_file>] [--strain_cluster=<strain_cluster_file>]
    create_reference_files.py --help

Options:
    --help                                 Show this help message and exit.
    --cglist=<core_genome_list>            Core genome list text file.
    --fastadir=<allele_fasta_files_dir>    Allele FASTA files.
    --cg_mlst_tsv=<tsv_file>               Sample set's core genome MLST .tsv file.
    --strain_cluster=<strain_cluster_file> A txt file contained strain and their related samples given as a separate line. i.e. STRAIN1:SAMPLE1,...[default:strain_cluster.txt].
"""

# python create_reference_files.py --cglist ./sagalactiae_cgMLST_files/cgMLSTschema.txt --cg_mlst_tsv ./sagalactiae_cgMLST_files/sagalactiae_cgMLST.tsv --fastadir ./sagalactiae_schema/
from Bio import SeqIO
from docopt import docopt
from io import StringIO  # python3
from pathlib import Path


import os, subprocess


# to remove duplicate alleles and squeezed form
def get_cg_pos_ref_alt_dict(cg, no_alleles):

    allele_locus_ref_alt_dict = dict()

    for allele_id in range(1, no_alleles):

        vcf_file = cg_allele_fasta_dir + cg + "/" + cg + "_" + str(allele_id) + ".vcf"

        with open(vcf_file, 'r') as file:

            for line in file.readlines():

                if not line.startswith("#"):

                    fields = line.strip().split('\t')

                    pos = fields[1]
                    ref = fields[3]
                    alt = fields[4]

                    if not pos in allele_locus_ref_alt_dict.keys():

                        allele_locus_ref_alt_dict[pos] = list()
                        allele_locus_ref_alt_dict[pos].append(ref)
                        allele_locus_ref_alt_dict[pos].append(alt)

                    else:

                        if alt not in allele_locus_ref_alt_dict[pos]:
                            allele_locus_ref_alt_dict[pos].append(alt)

            file.close()

    return dict(sorted(allele_locus_ref_alt_dict.items(), key=lambda x: int(x[0])))

# write graph vcf file for an allele
def write_graph_vcf_file_for_an_allele(cg_name, no_alleles, cg_pos_ref_alt_dict):

    output_file = open(cg_allele_fasta_dir + "/" + cg_name + "_gt.vcf", 'w')

    graph_vcf_line_dict = dict()

    is_vcf_header_written = False

    id_list = list(range(1, no_alleles))

    for allele_id in id_list:

        vcf_file = cg_allele_fasta_dir + cg_name + "/" + cg_name + "_" + str(allele_id) + ".vcf"

        with open(vcf_file, 'r') as file:

            for line in file.readlines():

                if line.startswith("#"):

                    if not is_vcf_header_written:

                        if line.startswith("##"):

                            output_file.write(line.strip() + '\n')

                        else:

                            allele_id_list = [cg_name + '_' + str(id) for id in id_list]
                            fields = line.strip().split('\t')[:-1] + allele_id_list
                            output_file.write('\t'.join(fields) + '\n')

                else:

                    fields = list(line.strip().split('\t')[:])

                    pos = fields[1]
                    alt = fields[4]

                    if pos not in graph_vcf_line_dict.keys():
                        graph_vcf_line_dict[pos] = ['0'] * len(id_list)

                    graph_vcf_line_dict[pos][allele_id-1] = str(cg_pos_ref_alt_dict[pos].index(alt))

            is_vcf_header_written = True

            file.close()

    for key, value in cg_pos_ref_alt_dict.items():
        output_file.write(str(cg_name) + '\t' + str(key) + '\t.\t' + str(value[0]) + '\t' + ','.join(value[1:]) + '\t30\t.\t.\tGT\t' + '\t'.join(graph_vcf_line_dict[key]) + '\n')

    output_file.close()

    # os.system("rm -r %s" % cg_allele_fasta_dir + "/" + cg_name)


# return all sequences with details of each core gene
def get_cg_sequences(cg_fasta_file_name):

    sequence_dict = {}
    try:
        for sequence in list(SeqIO.parse(StringIO(open(cg_fasta_file_name, 'r').read()), 'fasta')):
            raw_seq_id_fields = str(sequence.id).strip("\n").split("_")

            sequence_dict[raw_seq_id_fields[0] + "_" + raw_seq_id_fields[-1]] = str(sequence.seq)
    except FileNotFoundError:
        print('')
    return sorted(sequence_dict.items(), key=lambda item: int(item[0].split("_")[-1]))


# create vcf file for a core gene
def create_vcf_for_cg(cg):

    cg_name, allele_id = cg[0], cg[1]

    ref = cg_ref_fasta_dir + cg_name
    allele = cg_allele_fasta_dir + cg_name + "/" + cg_name + "_" + str(allele_id)

    # REFERENCE SEQUENCE
    ref_fasta = add_ext(ref, "fasta")

    # ALLELE SEQUENCE
    fasta = add_ext(allele, "fasta")
    sam = add_ext(allele, "sam")
    bam = add_ext(allele, "bam")
    sorted_bam = add_ext(allele, "sorted.bam")
    vcf = add_ext(allele, "vcf")
    bam_bai = add_ext(allele, "bam.bai")
    sorted_bam_bai = add_ext(allele, "sorted.bam.bai")

    command_list = []  # Command List
    command_list.append("minimap2 -ax asm5 %s %s --cs=long > %s" % (ref_fasta, fasta, sam))
    command_list.append("samtools view -Sb %s > %s" % (sam, bam))
    command_list.append("samtools index %s" % bam)
    command_list.append("samtools sort %s -o %s" % (bam, sorted_bam))
    command_list.append("samtools index %s" % sorted_bam)
    command_list.append("samtools mpileup -uf %s %s | bcftools call -Ov -c -v > %s" % (ref_fasta, sorted_bam, vcf))
    command_list.append("rm %s; rm %s; rm %s; rm %s; rm %s; rm %s" % (fasta, sam, bam, sorted_bam, bam_bai, sorted_bam_bai))

    for command in command_list:
        subprocess.call(command, shell=True, stdout=subprocess.DEVNULL)

# get core gene names from file
def get_cg_list(cg_list_file_name):

    with open(cg_list_file_name, 'r') as file:

        cg_list = list()
        for line in file.readlines():
            cg_list.append(line.strip())

        file.close()
        return cg_list

# write sequences to reference and allele fasta directories
def split_and_write_reference_and_other_alleles(sequence_dict):

    for allele_name, sequence in sequence_dict:

        if len(sequence_dict) != 1:  # only get genes with more than one allele to be discriminative

            cg_name, allele_id = allele_name.split("_")[0], int(allele_name.split("_")[1])

            cg_fasta_header = cg_name + "_" + str(allele_id)

            seq = ">" + cg_fasta_header + "\n" + sequence + "\n"

            if allele_id == 1:

                out_file_name = cg_ref_fasta_dir + cg_name + ".fasta"
                Path(cg_allele_fasta_dir + cg_name).mkdir(exist_ok=True)  # create separate directory for each gene
                out_ref_file = open(out_file_name, "w")
                out_ref_file.write(seq)
                out_ref_file.close()

            out_file_name = cg_allele_fasta_dir + cg_name + "/" + cg_fasta_header + ".fasta"
            out_ref_file = open(out_file_name, "w")
            out_ref_file.write(seq)
            out_ref_file.close()

            create_vcf_for_cg([cg_name, allele_id])


# write separated reference and remaining alleles to the given directory
def create_ref_and_other_allele_directories():

    try:
        Path(cg_ref_fasta_dir).mkdir(exist_ok=True)  # create reference FASTA directory for each core gene
        Path(cg_allele_fasta_dir).mkdir(exist_ok=True)  # create FASTA directory for each core gene's alleles

    except OSError:
        print("Creation of the directories %s and %s are failed.", cg_ref_fasta_dir, cg_allele_fasta_dir)

    for cg in get_cg_list(cg_list_file_name):

        cg_name = cg.strip(".fasta")

        sequence_dict = get_cg_sequences(fasta_dir + cg)

        no_alleles = len(sequence_dict) + 1
        split_and_write_reference_and_other_alleles(sequence_dict)

        write_graph_vcf_file_for_an_allele(cg_name, no_alleles, get_cg_pos_ref_alt_dict(cg_name, no_alleles))

# add extension at the end of the file name
def add_ext(file_name, ext):
    return file_name + "." + ext


# dict = { ERR100:[1,2,2,...,n], ERR200:[1,2,2,...,n], ... }
def get_strain_allele_profile_dict():
    strain_allele_profile_dict = dict()

    with open(args["--cg_mlst_tsv"], 'r') as file:
        for line in file.readlines():
            fields = line.strip().split('\t')

            strain_allele_profile_dict[fields[0]] = fields[1:]

        file.close()

    return strain_allele_profile_dict


# dict = { ACIBA00015_pos_ref_alt:[id_ACIBA00002_1, id_ACIBA00002_2, ..., id_ACIBA00002_n], ACIBA00035_pos_ref_alt:[id_ACIBA00035_1, id_ACIBA00035_2, ..., id_ACIBA00035_m], ... }
def get_cds_allele_gt_dict():  #

    cds_length_dict = dict()
    cds_allele_gt_dict = dict()

    for allele_id in strain_allele_profile_dict["FILE"]:

        cds_gt_vcf_file = cg_allele_fasta_dir + allele_id[:-6] + "_gt.vcf"  # [:-6] removes .fasta extension from header

        if os.path.isfile(cds_gt_vcf_file):

            with open(cds_gt_vcf_file, 'r') as file:

                for line in file.readlines():

                    if "contig" in line:
                        cds = line[line.index("ID=") + 3:line.index("_1")]
                        length = line[line.rindex("=") + 1:-1]

                        cds_length_dict[cds] = length

                    if not line.startswith("#"):
                        fields = line.strip().split('\t')

                        cds_key = fields[0] + "_" + fields[1] + "_" + fields[3] + "_" + fields[4]

                        cds_allele_gt_dict[cds_key] = fields[9:]

                file.close()

    return [cds_allele_gt_dict, cds_length_dict]


# dict = {}
def get_strain_allele_gt_dict():

    strain_allele_gt_dict = dict()

    for strain, strain_allele_id_list in strain_allele_profile_dict.items():  # { strain1:[cds1_allele_id, cds2_allele_id, ...], strain2:...}

        if strain != "FILE":

            for strain_cds_allele_id in range(len(strain_allele_id_list)):  # strain_cds_allele_id th cds

                cds_allele = strain_allele_profile_dict["FILE"][strain_cds_allele_id].strip(".fasta")

                for key, values in cds_allele_gt_dict.items():

                    if cds_allele == key.split("_")[0]:

                        if key not in strain_allele_gt_dict.keys():

                            key_fields = key.split('_')

                            # sample: key_fields[0] pos: key_fields[1] ref: key_fields[2] alt: key_fields[3]
                            strain_allele_gt_dict[key] = [key_fields[0], key_fields[1], '.', key_fields[2],
                                                          key_fields[3], '.', '.', 'AC=100;AN=100;', 'GT',
                                                          values[int(strain_allele_id_list[strain_cds_allele_id])]]

                        else:
                            strain_allele_gt_dict[key].append(
                                values[int(strain_allele_id_list[strain_cds_allele_id]) - 1])

    return strain_allele_gt_dict


# TODO: change here
# TODO: check positions from 0 or separate genes
def add_header_to_vcf_file():
    header = """##fileformat=VCFv4.2
##FILTER=<ID=PASS,Description="All filters passed">
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
##INFO=<ID=AC,Number=A,Type=Integer,Description="Allele Counts">
##INFO=<ID=AN,Number=1,Type=Integer,Description="Total number of alleles">
"""
    return header


# write results to reference.vcf file
def create_strain_graph_vcf_file():
    output_file = open(cg_allele_fasta_dir + "reference.vcf", 'w')

    add_header_to_vcf_file()
    for key, values in strain_allele_gt_dict.items():
        output_file.write('\t'.join(str(value) for value in values) + '\n')

    output_file.close()


def get_strain_cluster_dict():
    strain_cluster_dict = {}

    with open(args["--strain_cluster"], "r") as file:

        for line in file.readlines():

            strain, samples = line.strip().split(":")

            strain_cluster_dict[strain] = []

            for sample in samples.split(","):
                strain_cluster_dict[strain].append(sample)

        file.close()

    return strain_cluster_dict


def create_graph_reference_files():
    out_vcf_file = open(cg_allele_fasta_dir + "strain_graph.vcf", 'w')

    remaining_cds_set = set()

    with open(cg_allele_fasta_dir + "reference.vcf", 'r') as in_file:

        out_vcf_file.write(add_header_to_vcf_file())

        for strain, samples in strain_cluster_dict.items():
            out_vcf_file.write("##SAMPLE=<ID=%s,GENOMES=%s>\n" % (strain, ";".join(samples)))

        all_samples = []
        for samples in strain_cluster_dict.values():
            all_samples.extend(samples)

        out_vcf_file.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t%s\n" % "\t".join(all_samples))

        for line in in_file.readlines():

            is_line_written = False

            if line != "":

                gt_fields = line.strip("\n").split("\t")[9:]  # no of info field = 9

                gt_to_compare = gt_fields[0]  # ref

                for gt in gt_fields:

                    if gt != gt_to_compare:
                        is_line_written = True

            if is_line_written:

                current_cds = line.split("\t")[0]

                if current_cds not in remaining_cds_set:
                    remaining_cds_set.add(current_cds)

                out_vcf_file.write(line.strip("\n") + '\n')

        in_file.close()

    out_vcf_file.close()

    # add contig line without seeking
    ##contig=<ID=CDS_ID,length=CDS_LENGTH>
    with open(cg_allele_fasta_dir + "strain_graph.vcf", 'r') as in_file:

        with open(cg_allele_fasta_dir + "reference.vcf", 'w') as out_file:

            for line_idx, line in enumerate(in_file.readlines()):

                if line_idx == 2:

                    for cds in remaining_cds_set:
                        out_file.write("##contig=<ID=%s,length=%s>\n" % (cds, cds_length_dict[cds]))

                else:
                    out_file.write(line)

            out_file.close()

        in_file.close()

    os.remove(cg_allele_fasta_dir + "strain_graph.vcf")

    # create fasta file for the reference genome
    with open(cg_allele_fasta_dir + "reference.fasta", 'w') as out_file:

        for cds in remaining_cds_set:

            with open(cg_ref_fasta_dir + cds + ".fasta") as in_file:

                for line in in_file.readlines():
                    if line.startswith(">"):
                        out_file.write(line[:-3] + "\n")

                    else:
                        out_file.write(line)

                in_file.close()

        out_file.close()


if __name__ == "__main__":

    args = docopt(__doc__, version="0.6.2")

    cg_list_file_name, fasta_dir = args["--cglist"], args["--fastadir"]

    wd = ("/").join(fasta_dir.split("/")[:-2])
    cg_ref_fasta_dir = wd + "/cg_reference_fasta/"
    cg_allele_fasta_dir = wd + "/cg_allele_fasta/"

    # create_ref_and_other_allele_directories()
    strain_allele_profile_dict = get_strain_allele_profile_dict()
    cds_allele_gt_dict, cds_length_dict = get_cds_allele_gt_dict()
    strain_allele_gt_dict = get_strain_allele_gt_dict()
    create_strain_graph_vcf_file()
    # strain_cluster_dict = get_strain_cluster_dict()
    # create_graph_reference_files()
