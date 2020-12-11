#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""

@author: rfm
"""


import os
import csv
import argparse

from Bio import SeqIO


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
        infile.close()

    return lines


def join_list(lst, link):
    """ Joins all elements in a list into a single string.

        Parameters
        ----------
        lst : list
            List with elements to be joined.
        link : str
            Character used to join list elements.

        Returns
        -------
        joined_list : str
            A single string with all elements in the input
            list joined by the character chosen as link.
    """

    joined_list = link.join(lst)

    return joined_list


def write_to_file(text, output_file, write_mode, end_char):
    """ Writes a single string to a file.

        Parameters
        ----------
        text : str
            A single string to write to the output file.
        output_file : str
            Path to the output file.
        write_mode : str
            Write mode can be 'w', writes text and overwrites
            any text in file, or 'a', appends text to text
            already in file.
        end_char : str
            Character added to the end of the file.
    """

    with open(output_file, write_mode) as outfile:
        outfile.write(text+end_char)
        outfile.close()


def write_lines(lines, output_file):
    """ Writes a list of strings to a file. The strings
        are joined with newlines before being written to
        file.

        Parameters
        ----------
        lines : list
            List with the lines/strings to write to the
            output file.
        output_file : str
            Path to the output file.
    """

    joined_lines = join_list(lines, '\n')

    write_to_file(joined_lines, output_file, 'a', '\n')


def main(chewie_matrix, chewie_schema, milestone_results, strain_id,
         output_file):

    # process chewie's results
    chewie_lines = read_tabular(chewie_matrix)

    loci = chewie_lines[0]
    strain_profile = [l for l in chewie_lines
                      if l[0].split('.fasta')[0] == strain_id][0]

    # create dictionary with data from chewie's results
    chewie_allele_ids = {loci[i].split('.fasta')[0]: strain_profile[i].replace('INF-', '')
                         for i in range(1, len(loci))}

    # process milestone's results
    milestone_seqs = {(rec.id).split('_')[0]: str(rec.seq)
                      for rec in SeqIO.parse(milestone_results, 'fasta')}

    # get alleles in schema and compare with milestone's alleles
    milestone_allele_ids = {}
    for locid, seq in milestone_seqs.items():
        locus_file = os.path.join(chewie_schema, '{0}.fasta'.format(locid))
        locus_alleles = {str(rec.seq): (rec.id).split('_')[-1]
                         for rec in SeqIO.parse(locus_file, 'fasta')}

        milestone_id = locus_alleles.get(seq, '?')
        milestone_allele_ids[locid] = milestone_id

    # merge results
    merged_results = {k: [v, milestone_allele_ids.get(k, '-')]
                      for k, v in chewie_allele_ids.items()}

    # write results to output file
    table_header = 'locus\tchewie\tmilestone'
    table_lines = ['{0}\t{1}\t{2}'.format(k, v[0], v[1])
                   for k, v in merged_results.items()]
    table_lines = [table_header] + table_lines
    write_lines(table_lines, output_file)

    # get results that differ
    diffs = {k: v for k, v in merged_results.items() if v[0] != v[1]}

    # print diffs before completing process
    diffs_header = '{:<25}  {:^6}  {:^9}'.format('locus', 'chewie', 'milestone')
    diffs_lines = ['{:<25}  {:^6}  {:^9}'.format(k, v[0], v[1])
                   for k, v in diffs.items()]
    diffs_lines = [diffs_header] + diffs_lines
    diffs_text = '\n'.join(diffs_lines)
    print(diffs_text)


def parse_arguments():

    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)

    parser.add_argument('-cm', type=str, required=True,
                        dest='chewie_matrix',
                        help='Path to the TSV file with allelic '
                             'profiles determined by chewBBACA.')

    parser.add_argument('-cs', type=str, required=True,
                        dest='chewie_schema',
                        help='Path to chewies\'s schema.')

    parser.add_argument('-mr', type=str, required=True,
                        dest='milestone_results',
                        help='Path to a FASTA file with the '
                             'alleles predicted by Milestone.')

    parser.add_argument('-sid', type=str, required=True,
                        dest='strain_id',
                        help='Identifier of the strain.')

    parser.add_argument('-o', type=str, required=True,
                        dest='output_file',
                        help='Path to the output file. The '
                             'output file has three columns: '
                             'locus identifier, chewie allele, '
                             'milestone allele.')

    args = parser.parse_args()

    return [args.chewie_matrix, args.chewie_schema,
            args.milestone_results, args.strain_id,
            args.output_file]


if __name__ == '__main__':

    args = parse_arguments()
    main(*args)
