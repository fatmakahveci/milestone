#!/usr/bin/env python3

"""
-------------------------------------------------------
Aim: High quality variant selection based on base ratio
-------------------------------------------------------
Authors: @fatmakhv
The latest update: June 01, 2021
-------------------------------------------------------
"""

import argparse


class Sam:

    def __init__(self, sam_line):
        fields = sam_line.strip('\n').split('\t')

        self.query_name = fields[0]
        self.read_info = fields[1]
        self.target_sequence_name = fields[2]
        self.position = int(fields[3])
        self.mapping_quality = fields[4]
        self.cigar = fields[5]
        self.rnext = fields[6]
        self.pnext = fields[7]
        self.tlen = fields[8]
        self.sequence = fields[9]
        self.read_quality = fields[10]

    def __repr__(self):
        return f'QNAME: {self.query_name}\t' \
               f'FLAG: {self.read_info}\t' \
               f'RNAME: {self.target_sequence_name}\t' \
               f'POS: {self.position}\t' \
               f'MAPQ: {self.mapping_quality}\t' \
               f'CIGAR: {self.cigar}\t' \
               f'RNEXT: {self.rnext}\t' \
               f'PNEXT: {self.pnext}\t' \
               f'TLEN: {self.tlen}\t' \
               f'SEQ: {self.sequence}\t' \
               f'QUAL: {self.read_quality}'


def reference_cds_seq() -> str:
    """
    Reads <reference.fasta> and returns sequence for CDS

    Returns
    -------
    cds_ref_seq: Sequence of CDS reference
    """

    from Bio import SeqIO

    for sequence in list(SeqIO.parse(args.reference_fasta, 'fasta')):

        if sequence.id == args.cds:

            cds_ref_seq = str(sequence.seq)

            return cds_ref_seq


def mapping_quality_check(sam) -> bool:
    """
    Compares Sam's mapping quality to user-specified parameter

    Returns
    -------
    is_passed : True if its mapping quality is higher than parameter
    """

    is_passed = True if int(sam.mapping_quality) >= args.mapping_quality else \
        False

    return is_passed


def matching_sequence_region(sequence: str, start: int, end: int) -> str:
    """
    @todo

    Parameters
    ----------
    sequence : @todo
    start : @todo
    end : @todo

    Returns
    -------
    matching_part : @todo
    """

    matching_part = sequence[start-1:end-1]

    return matching_part


def cigar(sam: Sam, base_dict: dict) -> None:
    """
    @todo

    Parameters
    ----------
    sam : @todo
    """

    import re

    match_pos_list = []
    current_position = 0
    pos_in_ref = args.pos - 1
    sam_start_pos = sam.position - 1

    # match/mismatch loci in sam sequence
    for case, number_of_case in zip(re.findall("[A-Z]", sam.cigar),
                                    map(int, re.findall("[0-9]+", sam.cigar))):
        if case == 'M':
            match_start = current_position
            current_position += number_of_case
            match_end = current_position - 1
            match_pos_list.append([match_start, match_end])

        else:
            current_position += number_of_case

    for start, end in match_pos_list:

        pos_in_sam = pos_in_ref - sam_start_pos + start

        if start <= pos_in_sam <= end:
            base_dict[sam.sequence[pos_in_sam]] += 1

    return base_dict

def base_ratio_check() -> None:
    """

    Returns
    -------

    """

    # avoid redundant operations after cds is found
    is_cds_found = False

    base_dict = { 'A': 0, 'C': 0, 'G': 0, 'T': 0 }

    with open(args.sam, 'r') as file:

        for line in file.readlines():

            # skip header
            if not line.startswith('@'):

                sam = Sam(line)

                # get sample related sam lines
                if sam.target_sequence_name == args.cds:

                    is_cds_found = True

                    if mapping_quality_check(sam):
                        cigar(sam, base_dict)

                # avoid redundant operations after cds is found
                if sam.target_sequence_name != args.cds and is_cds_found:
                    break

        file.close()
    print(base_dict)

if __name__ == "__main__":

    parser = argparse.ArgumentParser()

    # @todo check: sorted sam is required
    parser.add_argument('--sam', type = str, required =  True,
                        help = 'Sample\'s Sam file to be analyzed.')
    parser.add_argument('--mapping_quality', type = int, required = False,
                        default = 30,
                        help = 'Minimum quality for read mapping.')
    parser.add_argument('--cds', type = str, required = True,
                        help = 'Name of CDS to be analyzed.' )
    parser.add_argument('--pos', type = int, required = True,
                        help = 'Position to be analyzed.')
    parser.add_argument('--reference_fasta', type = str, required = False,
                        help = 'Reference Fasta file for quality checks.')

    args = parser.parse_args()

    base_ratio_check()