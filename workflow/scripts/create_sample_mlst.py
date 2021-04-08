#!/usr/bin/env python3
# -*- coding: utf-8 -*-

'''
author: @rfm, @fatmakhv
'''

import argparse

from Bio import SeqIO


def get_sequence_dict_from_fasta_file(file_name: str, remove_allele_id: bool) -> dict:
    """ It reads FASTA-formatted file and returns the sequences inside
        as a list.

    Parameters
    ----------
    file_name : FASTA-formatted file name to be read

    Returns
    -------
    seq_dict : dict of sequences inside <file_name>

    """

    seq_dict = {}

    for rec in SeqIO.parse(file_name, 'fasta'):
        # ERR3464559-protein120_extra-info_2 -> ERR3464559-protein120_2
        cds_id, allele_id = str(rec.id).split('_')[0],str(rec.id).split('_')[-1]
        if not remove_allele_id:
            cds_id = f'{cds_id}_{allele_id}'
        seq_dict[cds_id] = str(rec.seq)

    return seq_dict


def get_chewbbaca_sequences_for_milestone_cds(milestone_cds: str) -> list:
    """ Read allele sequences of given CDS from schema_seed

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


def compare_milestone_seq_to_chewbbaca_allele_seqs(milestone_seq_dict: dict) -> None:

    for milestone_cds, milestone_seq in milestone_seq_dict.items():

        chewbbaca_sequence_list_of_given_cds = get_chewbbaca_sequences_for_milestone_cds(milestone_cds)
        
        # @todo create sample.mlst.tsv
        if milestone_seq not in chewbbaca_sequence_list_of_given_cds.values():
            pass
            # @todo compare milestone_cds sequence to chewbbaca_cds alleles minimap2 to find variants
            # @todo milestone_cds sequence add to schema_seed/cds alleles
            # @todo milestone_cds.vcf add to reference.vcf
            # @todo add milestone_allele_id to sample.mlst.tsv allele_id = {{number of alleles of chewbbaca CDS}+1}

        # else:
            # @todo add milestone_allele_id = '1' to sample.mlst.tsv


def parse_arguments():

    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)

    parser.add_argument('--chewbbaca_path', type=str, required=True,
                        help='Path to directory containing FASTA sequences of '
                             'alleles of CDSs created by chewBBACA.')

    parser.add_argument('--milestone_path', type=str, required=True,
                        help='Path to directory containing FASTA sequences of '
                             'alleles of CDSs created by milestone.')

    parser.add_argument('--strain_id', type=str, required=True,
                        help='Identifier of the strain.')

    parser.add_argument('--sample_mlst',type=str, required=True,
                        help='Output file name representing <sample.tsv>')

    args = parser.parse_args()

    return args


if __name__ == '__main__':

    args = parse_arguments()

    milestone_seq_dict = get_sequence_dict_from_fasta_file(args.milestone_path, True)

    compare_milestone_seq_to_chewbbaca_allele_seqs(milestone_seq_dict)
