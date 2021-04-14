#!/usr/bin/env python3
# -*- coding: utf-8 -*-

'''
author: @rfm, @fatmakhv
the latest update: April 13, 2021
'''

import argparse, csv, os, subprocess

from Bio import SeqIO


class variant:

    def __init__(self, id, pos, ref, alt):
        self.id = id
        self.pos = pos
        self.ref = ref
        self.alt = alt

    def __repr__(self):
        return f'{self.id}\t{self.pos}\t{self.ref}\t{self.alt}\n'


def get_sequence_dict_from_fasta_file(file_name: str, remove_allele_id: bool) -> dict:
    """ reads FASTA-formatted file and returns the sequences inside
        as a list.

    Parameters
    ----------
    file_name : FASTA-formatted file name to be read
    remove_allele_id : checks whether if allele IDs will be removed or not

    Returns
    -------
    seq_dict : dict of sequences inside <file_name>
    """

    seq_dict = {}

    for rec in SeqIO.parse(file_name, 'fasta'):

        # ERR3464559-protein120_extra-info_2 -> ERR3464559-protein120_2
        cds_id, allele_id = str(rec.id).split('_')[0], str(rec.id).split('_')[-1]

        if not remove_allele_id:
            cds_id = f'{cds_id}_{allele_id}'

        seq_dict[cds_id] = str(rec.seq)

    return seq_dict


def get_chewbbaca_sequences_for_milestone_cds(milestone_cds: str) -> list:
    """ reads allele sequences of given CDS from schema_seed.

    Parameters
    ----------
    milestone_cds : name of CDS

    Returns
    -------
    chewbbaca_sequence_list_of_given_cds : list of sequences created by chewBBACA for given CDS
    """

    chewbbaca_cds_name = f'{args.chewbbaca_path}/{milestone_cds}.fasta'

    chewbbaca_sequence_list_of_given_cds = get_sequence_dict_from_fasta_file(chewbbaca_cds_name, False)

    return chewbbaca_sequence_list_of_given_cds


def write_seq_to_file(seq_id: str, seq: str, seq_file: str) -> None:
    """ gets <seq> FASTA and <seq_id>
        to write FASTA-formatted <seq_file>

    Parameters
    ----------
    seq_id : to write to header
    seq : FASTA sequence to be written to <seq_file>
    seq_file : name of FASTA file to write sequence
    """

    file = open(seq_file, 'w')
    file.write(f'>{seq_id}\n{seq}')
    file.close()


def read_tabular(input_file: str) -> list:
    """ reads tabular file and returns as list.

        Parameters
        ----------
        input_file : path to tab-separated file.

        Returns
        -------
        lines : a list with a sublist per line in the input file.
    """

    with open(input_file, 'r') as infile:

        reader = csv.reader(infile, delimiter='\t')
        lines = [line for line in reader]

    return lines


def call_variants(cds_name: str, sequence: str) -> list:
    """ Call variants from PAF-formatted file
        (insertions, deletions, and substitutions)

    Parameters
    ----------
    cds_name : name of CDS of which variants are called
    sequence : sequence to be processed

    Returns
    -------
    variant_list : called variants will be returned as 'variant' list
    """

    variant_list = []

    idx, pos = 0, 0
    while idx < len(sequence):

        base = sequence[idx]

        if base.upper() in ['A', 'C', 'G', 'T']: # identical
            pos += 1 # skip matching bases

        elif base == '=': # equality ~ placeholder
            pass

        elif base == '*': # substitution
            start_idx = idx + 2
            ref_base = sequence[start_idx-1].upper()
            alt_base = sequence[start_idx].upper()
            variant_list.append(variant(cds_name, pos+1, ref_base, alt_base))
            pos += 2 # reference base is inside of *[ref_base][alt_base]=

        elif base == '+': # insertion
            start_idx = idx + 1
            end_idx = start_idx + sequence[start_idx:].index('=')
            ref_base = sequence[idx-1]
            insertion = sequence[start_idx:end_idx].upper()
            alt_base = ref_base + insertion
            variant_list.append(variant(cds_name, pos, ref_base, alt_base))

        elif base == '-': # deletion
            start_idx = idx + 1
            end_idx = start_idx + sequence[start_idx:].index('=')
            base_before_deletion = sequence[start_idx-2]
            ref_base = base_before_deletion + sequence[start_idx:end_idx].upper()
            alt_base = base_before_deletion
            variant_list.append(variant(cds_name, pos+1, ref_base, alt_base))

        idx += 1

    return variant_list


def write_variants_to_file(variant_list: list) -> None:
    """ gets variants from PAF-formatted alignment and
        appends variants at the end of the <reference.vcf> file.

    Parameters
    ----------
    cds_name : name of CDS to be processed
    """

    with open(args.reference_vcf, 'a') as vcf_file:

        for variant in variant_list:
            vcf_file.write(f'{variant.id}_1\t{variant.pos}\t.\t{variant.ref}\t{variant.alt}\t{30}\t.\tDP=1;SGB=-0.379885;MQ0F=0;AF1=1;AC1=1;DP4=0,0,1,0;MQ=60;FQ=-999   GT:PL\t1:60,0\n')

        vcf_file.close()


def find_variants(milestone_cds: str, chewbbaca_seq: str, milestone_seq: str, milestone_allele_id: int) -> None:
    """ compares milestone_cds sequence to chewbbaca_cds
        alleles minimap2 to find variants

    Parameters
    ----------
    milestone_cds : CDS ID to work on
    chewbbaca_seq : reference sequence of CDS created by chewBBACA
    milestone_seq : sequence of CDS created by milestone
    milestone_allele_id : (novel) allele ID to be written into <sample.mlst.tsv>
    """

    chewbbaca_seq_fasta = f'{args.output_dir}/{milestone_cds}_chewbbaca.fasta'
    write_seq_to_file(milestone_cds, chewbbaca_seq, chewbbaca_seq_fasta)

    milestone_seq_fasta = f'{args.output_dir}/{milestone_cds}_milestone.fasta'
    write_seq_to_file(milestone_cds, milestone_seq, milestone_seq_fasta)

    cds_paf_file_name = f'{args.output_dir}/{milestone_cds}.paf'

    cmd = (f'minimap2 -cx asm5 {chewbbaca_seq_fasta} {milestone_seq_fasta} -t {args.threads} --cs=long -o {cds_paf_file_name} 2>&1')
    subprocess.call(cmd, shell=True, stdout=subprocess.DEVNULL)

    os.remove(chewbbaca_seq_fasta) if os.path.exists(chewbbaca_seq_fasta) else None
    os.remove(milestone_seq_fasta) if os.path.exists(milestone_seq_fasta) else None

    paf_line = read_tabular(cds_paf_file_name)

    # case: len(paf_lines) = 0
    if len(paf_line) != 0: # if the alignment info exists.

        paf_line = paf_line[0] # list of list -> list
        variant_list = call_variants(cds_name=paf_line[0], sequence=paf_line[-1][5:]) # '[5:] skip cs:z:'

        if args.reference_vcf != "": # if --update_reference is used to run this script

            write_variants_to_file(variant_list)

            for _ in variant_list: # write novel allele sequences among the chewbbaca's schema_seed
                # milestone_sequence_id = novel sequence id
                milestone_sequence_id = f'{milestone_cds}_{milestone_allele_id}'
                write_seq_to_file(milestone_sequence_id, milestone_seq, f'{args.chewbbaca_path}/{milestone_sequence_id}.fasta')

    os.remove(cds_paf_file_name) if os.path.exists(cds_paf_file_name) else None

def compare_milestone_seq_to_chewbbaca_allele_seqs(milestone_seq_dict: dict) -> None:
    """ compares each CDS sequence created by milestone to the allele sequences
        of the CDS created by chewBBACA. If it exist among them, it writes the same
        allele ID. Otherwise, it writes allele ID as (#alleles+1) to show its novelty.

    Parameters
    ----------
    milestone_seq_dict : { CDS1: <seq>, CDS2: <seq>, ... }
    """

    sample_mlst_file = open(f'{args.output_dir}/{args.sample_id}.mlst.tsv', 'w')

    # header of the <sample.mlst.tsv>
    sample_mlst_file.write(f'CDS\tAllele ID\n')

    for milestone_cds, milestone_seq in milestone_seq_dict.items():

        chewbbaca_sequence_list_of_given_cds = get_chewbbaca_sequences_for_milestone_cds(milestone_cds)

        if milestone_seq in chewbbaca_sequence_list_of_given_cds.values():

            for chewbbaca_allele_name, chewbbaca_allele_seq in chewbbaca_sequence_list_of_given_cds.items():

                if milestone_seq == chewbbaca_allele_seq:

                    # chewbbaca_allele_name = <chewbbaca_CDS_name>_[<chewbbaca_allele_id>==(milestone_allele_id)]
                    milestone_allele_id = chewbbaca_allele_name.split('_')[-1]

        else:

            chewbbaca_cds_reference = chewbbaca_sequence_list_of_given_cds[f'{milestone_cds}_1']

            # change allele ID for the novel allele to add new allele
            # novel allele id = milestone_allele_id
            milestone_allele_id = len(chewbbaca_sequence_list_of_given_cds) + 1

            find_variants(milestone_cds, chewbbaca_cds_reference, milestone_seq, milestone_allele_id)

        sample_mlst_file.write(f'{milestone_cds}\t{milestone_allele_id}\n')

    sample_mlst_file.close()


def parse_arguments():

    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)

    parser.add_argument('--chewbbaca_path', type=str, required=True,
                        help='Path to directory containing FASTA sequences of '
                             'alleles of CDSs created by chewBBACA.')

    parser.add_argument('--milestone_path', type=str, required=True,
                        help='Path to directory containing FASTA sequences of '
                             'alleles of CDSs created by milestone.')

    parser.add_argument('--sample_id', type=str, required=True,
                        help='Identifier of the sample (without extensions).')

    parser.add_argument('--output_dir', type=str, required=True,
                        help='Output directory for <sample.mlst.tsv>')

    parser.add_argument('--threads', type=int, required=True,
                        help='Number of threads to run minimap2')

    # after_completion - if --update_reference is given to run this script
    parser.add_argument('--reference_vcf', type=str, required=False, default="",
                        help='Graph genome reference VCF file to update for called variants.')

    args = parser.parse_args()

    return args


if __name__ == '__main__':

    args = parse_arguments()

    milestone_seq_dict = get_sequence_dict_from_fasta_file(args.milestone_path, True)

    compare_milestone_seq_to_chewbbaca_allele_seqs(milestone_seq_dict)
