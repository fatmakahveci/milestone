#!/usr/bin/env python3
# -*- coding: utf-8 -*-

'''
author: @rfm, @fatmakhv
the latest update: April 09, 2021
'''

import argparse, csv, subprocess

from Bio import SeqIO


def get_sequence_dict_from_fasta_file(file_name: str, remove_allele_id: bool) -> dict:
    """ reads FASTA-formatted file and returns the sequences inside
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


def write_seq_to_file(seq_id, seq, seq_file):
    file = open(seq_file, 'w')
    file.write(f'>{seq_id}\n{seq}')
    file.close()


def read_tabular(input_file, delimiter='\t'):
    """ Read tabular file.
        Parameters
        ----------
        input_file : str
            Path to a tabular file.
        delimiter : str
            Delimiter used to separate file fields.
        Returns
        -------
        lines : list
            A list with a sublist per line in the input file.
            Each sublist has the fields that were separated by
            the defined delimiter.
    """

    with open(input_file, 'r') as infile:
        reader = csv.reader(infile, delimiter=delimiter)
        lines = [line for line in reader]

    return lines


def find_variants(milestone_cds: str, chewbbaca_seq: str, milestone_seq: str, novel_allele_id: int, threads: int, output_dir: str) -> None:
    """ compares milestone_cds sequence to chewbbaca_cds alleles minimap2 to find variants

    Parameters
    ----------
    milestone_cds : CDS ID to work on
    chewbbaca_seq : reference sequence of CDS created by chewBBACA
    milestone_seq : sequence of CDS created by milestone
    novel_allele_id : allele ID to be written into <sample.mlst.tsv>
    threads : number of threads to run minimap2
    """

    chewbbaca_seq_fasta = f'{output_dir}/{milestone_cds}_chewbbaca.fasta'
    write_seq_to_file(milestone_cds, chewbbaca_seq, chewbbaca_seq_fasta)

    milestone_seq_fasta = f'{output_dir}/{milestone_cds}_milestone.fasta'
    write_seq_to_file(milestone_cds, milestone_seq, milestone_seq_fasta)

    cds_paf_file_name = f'{output_dir}/{milestone_cds}.paf'

    cmd = (f'minimap2 -cx asm5 {chewbbaca_seq_fasta} {milestone_seq_fasta} -t {threads} --cs=long -o {cds_paf_file_name} 2>&1')
    subprocess.call(cmd, shell=True, stdout=subprocess.DEVNULL)

    paf_lines = read_tabular(cds_paf_file_name)
    print(paf_lines)


def compare_milestone_seq_to_chewbbaca_allele_seqs(milestone_seq_dict: dict, sample_id: str, output_dir: str, threads: int) -> None:
    """ compares each CDS sequence created by milestone to the allele sequences
        of the CDS created by chewBBACA. If it exist among them, it writes the same
        allele ID. Otherwise, it writes allele ID as (#alleles+1) to show its novelty.

    Parameters
    ----------
    milestone_seq_dict : { CDS1: <seq>, CDS2: <seq>, ... }
    sample_id : Sample ID given as parameter to use in output file name.
    output_dir : Output directory given as parameter to use for <sample.mlst.tsv>.
    """

    sample_mlst_file = open(f'{output_dir}/{sample_id}.mlst.tsv', 'w')

    # header of the <sample.mlst.tsv>
    sample_mlst_file.write(f'CDS\tAllele ID\n')

    for milestone_cds, milestone_seq in milestone_seq_dict.items():

        chewbbaca_sequence_list_of_given_cds = get_chewbbaca_sequences_for_milestone_cds(milestone_cds)

        if milestone_seq in chewbbaca_sequence_list_of_given_cds.values():

            for chewbbaca_allele_name, chewbbaca_allele_seq in chewbbaca_sequence_list_of_given_cds.items():

                if milestone_seq == chewbbaca_allele_seq:

                    # chewbbaca_allele_name = <chewbbaca_CDS_name>_<chewbbaca_allele_id>
                    chewbbaca_allele_id = chewbbaca_allele_name.split('_')[-1]

        else:

            chewbbaca_cds_reference = chewbbaca_sequence_list_of_given_cds[f'{milestone_cds}_1']

            # change allele ID for the novel allele to add new allele
            novel_allele_id = len(chewbbaca_sequence_list_of_given_cds) + 1

            find_variants(milestone_cds, chewbbaca_cds_reference, milestone_seq, novel_allele_id, threads, output_dir)

            # @todo milestone_cds sequence add to schema_seed/cds alleles
            # @todo milestone_cds.vcf add to reference.vcf
            # @todo add milestone_allele_id to sample.mlst.tsv allele_id = {{number of alleles of chewbbaca CDS}+1}

        sample_mlst_file.write(f'{milestone_cds}\t{chewbbaca_allele_id}\n')

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

    args = parser.parse_args()

    return args


if __name__ == '__main__':

    args = parse_arguments()

    milestone_seq_dict = get_sequence_dict_from_fasta_file(args.milestone_path, True)

    compare_milestone_seq_to_chewbbaca_allele_seqs(milestone_seq_dict, args.sample_id, args.output_dir, args.threads)
