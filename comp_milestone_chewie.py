#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""

@author: rfm
"""


import os
import re
import csv
import argparse
import subprocess
from copy import deepcopy
from itertools import groupby
from collections import Counter

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

    with open(output_file, write_mode) as out:
        out.write(text+end_char)


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


def run_minimap2(reference, fastq1, fastq2, output_file):
    """ Executes minimap2 to map short sequences
        to a reference genome.

        Parameters
        ----------
        reference : str
            Path to the FASTA file with the reference.
        map_fasta : str
            Path to FASTA file with the short sequences
            to map against the reference.
        output_file : str
            Path to the output file with mapping results.

        Returns
        -------
        List with following elements:
            stdout : list
                List with stdout from minimap2 in bytes.
            stderr : list
                List with the stderr from minimpa2 in bytes.
    """

    # -I parameter to control number of target bases loaded into memory
    minimap_args = ['minimap2 -I 100M --cs -cx sr {0} {1} {2} > '
                    '{3}'.format(reference, fastq1, fastq2, output_file)]

    minimap_proc = subprocess.Popen(minimap_args,
                                    shell=True,
                                    stdout=subprocess.PIPE,
                                    stderr=subprocess.PIPE)

    stderr = minimap_proc.stderr.readlines()
    stdout = minimap_proc.stdout.readlines()

    return [stdout, stderr]


def regex_matcher(string, pattern):
    """ Finds substrings that match a regex pattern.

        Parameters
        ----------
        string : str
            Input string.
        pattern : str
            Pattern to match. Patterns should
            include 'r' before string to match
            or characters after backslashes will
            be escaped.

        Returns
        -------
        matches : list
            List with substrings that matched the pattern.
    """

    matches = re.findall(pattern, string)

    return matches


def single_position_coverage(coverage_info, start):
    """ Determine if positions in a subsequence are
        covered based on information in the cs field
        in a PAF file created by minimpa2.

        Parameters
        ----------
        coverage_info : list
            List with subsequent operations extracted
            from the cd field in a PAF file created by
            minimap2.
        start : int
            Subsequence start position in the complete
            sequence.

        Returns
        -------
        coverage : dict
            Dictionary with sequence positions as keys
            and coverage for each position as values.
    """

    coverage = {}
    mismatches = []
    for m in coverage_info:
        # subsequence part with exact matches
        if m[0] == ':':
            # create dctionary entries with coverage = 1
            new_cov = {i: 1 for i in range(start, start+int(m[1:]))}
            coverage = {**coverage, **new_cov}
            # increment start position
            start = start + int(m[1:])
        # position with substitution
        elif m[0] == '*':
            coverage[start] = 0
            mismatches.extend([(start, m)])
            start += 1
        # position with deletion
        elif m[0] == '-':
            # coverage 0 for missing bases
            new_cov = {i: 0 for i in range(start, start+len(m[1:]))}
            coverage = {**coverage, **new_cov}
            mismatches.extend([(start, m)])
            start = start + len(m[1:])
        # insertion
        elif m[0] == '+':
            # do not add coverage values for positions because
            # insertion does not exist in reference
            mismatches.extend([(start, m)])

    return [coverage, mismatches]


def merge_intervals(intervals):
    """ Merges intersecting intervals.

        Parameters
        ----------
        intervals : dict
            Dictionary with sequence identifiers as keys
            and a list of lists as values. Each sublist has
            a start and stop position in the sequence and
            a dictionary with the coverage for every position
            in the sequence interval.

        Returns
        -------
        merged : list
            Dictionary with the result of merging intervals
            that overlapped (coverage data is updated and
            incremented for positions in common).
    """

    merged = [deepcopy(intervals[0])]
    for current in intervals[1:]:
        previous = merged[-1]
        # current and previous intervals intersect
        if current[0] <= previous[1]:
            # determine top position
            previous[1] = max(previous[1], current[1])
            # merge coverage dictionaries
            previous_cov = previous[2]
            current_cov = current[2]
            for k, v in current_cov.items():
                if k not in previous_cov:
                    previous_cov[k] = v
                else:
                    previous_cov[k] += v
            previous[2] = previous_cov
        # current and previous intervals do not intersect
        else:
            merged.append(deepcopy(current))

    return merged


def determine_breadth_coverage(intervals, total_bases):
    """ Determines the percentage and total number of covered
        bases according to a set of coverage intervals.

        Parameters
        ----------
        intervals : dict
            Dictionary with sequence identifiers as keys
            and a list of lists as values. Each sublist has
            a start and stop position in the sequence and
            a dictionary with the coverage for every position
            in the sequence interval.
        total_bases : int
            Total number of bases in the reference.

        Returns
        -------
        List with following elements:
            breadth_coverage : float
                Percentage of covered bases.
            covered_bases : int
                Total number of covered bases.
    """

    # determine breadth of coverage
    covered_bases = 0
    for k, v in intervals.items():
        for e in v:
            covered_bases += sum([1 for p, c in e[2].items() if c > 0])

    breadth_coverage = covered_bases / total_bases

    return [breadth_coverage, covered_bases]


def determine_depth_coverage(intervals, total_len):
    """ Determine depth of coverage for a sequence.

        Parameters
        ----------
        intervals : dict
            Dictionary with sequence identifiers as keys
            and a list of lists as values. Each sublist has
            a start and stop position in the sequence and
            a dictionary with the coverage for every position
            in the sequence interval.
        total_len : int
            Total length of the sequence.

        Returns
        -------
        List with following elements:
            positions_depth : dict
                Dictonary with sequence positions and keys
                and coverage for each position as values.
            counts : dict
                Dictionary with coverage values as keys and
                the total number of positions with that coverage
                value as values.
    """

    # create dictionary to add coverage for all positions
    positions = list(range(0, total_len))
    positions_depth = {p: 0 for p in positions}
    # increment coverage values based on intervals
    for i in intervals:
        for p, c in i[2].items():
            positions_depth[p] += c

    # determine coverage distribution
    counts = sorted(Counter(positions_depth.values()).most_common(),
                    key=lambda x: x[0])

    return [positions_depth, counts]


def determine_missing_intervals(intervals, identifier, total_len):
    """ Determines sequence intervals that are not covered by any
        probes.

        Parameters
        ----------
        intervals : dict
            Dictionary with sequence identifiers as keys
            and a list of lists as values. Each sublist has
            a start and stop position in the sequence and
            a dictionary with the coverage for every position
            in the sequence interval.
        identifier : str
            Sequence identifier.
        total_len : int
            Total length of the sequence.

        Returns
        -------
        List with following elements:
            missing_regions : dict
                Dictionary with sequence identifiers as keys
                a list of lists as values. Each sublist has
                the start and stop positions for a sequence
                interval that is not covered by probes.
            not_covered : int
                Total number of bases not covered by probes.
    """

    start = 0
    not_covered = 0
    missing_regions = {identifier: []}
    for i in intervals:
        diff = i[0] - start
        if diff > 0:
            missing_regions[identifier].append([start, start+diff])
            not_covered += diff
            start += diff

        # create groups of equal values
        values_groups = [list(j) for i, j in groupby(i[2].values())]
        for g in values_groups:
            if g[0] == 0:
                missing_regions[identifier].append([start, start+len(g)])
                not_covered += len(g)
                start += len(g)
            else:
                start += len(g)

    # add terminal region
    if start != total_len:
        missing_regions[identifier].append([start, total_len])
        not_covered += total_len - start

    return [missing_regions, not_covered]


def diffs_info(paf_lines, seqs, identity, miss_perc):
    """
    """

    # filter out based on percentage of read that has aligned
    # keep reads that align well at allele extremities?
    # only keep the best alignment for each read?

    valid = deepcopy(paf_lines)

    # compute alignment identity
    for i in range(len(valid)):
        valid[i].append(int(valid[i][9]) / int(valid[i][10]))

    # filter out alignments below defined identity
    invalid = [line for line in valid if line[-1] < identity]
    valid = [line for line in valid if line[-1] >= identity]

    # match alignment string with regex
    pattern = r':[0-9]+|\*[a-z][a-z]|\+[a-z]+|-[a-z]+'
    for i in range(len(valid)):
        current = valid[i][-2]
        valid[i].append(regex_matcher(current, pattern))

    # split results per locus
    loci_results = {}
    for l in valid:
        loci_results.setdefault(l[5], []).append(l)

    # for the results of each locus
    # get information about positions that match to determine coverage
    for i in range(len(valid)):
        current = valid[i][-1]
        start = int(valid[i][7])
        valid[i].append(single_position_coverage(current, start))

    # get indels and SNPs for each locus
    loci_miss = {}
    for l in valid:
        loci_miss.setdefault(l[5], []).extend(l[-1][1])

    # indels and SNPs counts
    loci_miss = {k: [v, Counter(v)] for k, v in loci_miss.items()}

    # identify subsequences that are well covered by baits
    covered_intervals = {}
    for l in valid:
        covered_intervals.setdefault(l[5], []).append([int(l[7]), int(l[8]), l[-1][0]])

    # sort covered intervals
    covered_intervals_sorted = {k: sorted(v, key=lambda x: x[0])
                                for k, v in covered_intervals.items()}

    # merge overlapping intervals
    # deepcopy to avoid altering original intervals
    merged_intervals = {k: merge_intervals(v)
                        for k, v in covered_intervals_sorted.items()}

    # select indels or SNPs based on frequency
    selected_miss = {}
    for k, v in loci_miss.items():
        counts = v[1]
        selected_counts = {}
        for p, c in counts.items():
            coverage_values = merged_intervals[k]
            for region in coverage_values:
                if p[0] >= region[0] and p[0] <= region[1]:
                    pos_cov = region[2][p[0]]
                    if pos_cov == 0 or c / pos_cov >= miss_perc:
                        selected_counts[p] = (c, pos_cov)
        selected_miss[k] = selected_counts

    for k, v in merged_intervals.items():
        coverage = determine_breadth_coverage({k: v}, len(seqs[k]))
        # determine subsequences that are not covered
        missing = determine_missing_intervals(v, k, len(seqs[k]))
        depth_info = determine_depth_coverage(v, len(seqs[k]))
        # determine mean depth of coverage
        depth_counts = depth_info[1]
        depth_sum = sum([c[0]*c[1] for c in depth_counts])
        mean_depth = round(depth_sum/len(seqs[k]), 4)

        merged_intervals[k].append([coverage, missing, mean_depth])

    return [merged_intervals, invalid, selected_miss]


# final table should have breatdh of coverage per allele, mean depth of coverage,
# number of positions with coverage below certain value, number of positions with
# 0 coverage, different variants supported by reads.
# change minimizer size used in minimap2
def diffs_coverage(diffs, chewie_schema, milestone_seqs, strain_id, output_dir,
                   fastq1, fastq2, identity, miss_perc):
    """
    """

    # get alleles predicted by milestone
    milestone_diffs = {k: milestone_seqs[k] for k, v in diffs.items() if v[1] != '0'}
    # write milstone alleles to file
    milestone_recs = ['>{0}\n{1}'.format(k, v) for k, v in milestone_diffs.items()]
    milestone_file = os.path.join(output_dir, 'milestone_diffs.fasta')
    write_lines(milestone_recs, milestone_file)

    # map reads against milestone alleles
    milestone_paf = os.path.join(output_dir, 'milestone.paf')
    run_minimap2(milestone_file, fastq1, fastq2, milestone_paf)

    # read
    milestone_paf_lines = read_tabular(milestone_paf)

    # get coverage info
    milestone_covinfo = diffs_info(milestone_paf_lines, milestone_seqs, identity, miss_perc)

    # get chewie alleles
    chewie_seqs = {}
    for k, v in diffs.items():
        allele_id = v[0]
        if allele_id != '0':
            locus_file = os.path.join(chewie_schema, k+'.fasta')
            # get allele sequence
            allele_seq = [str(rec.seq)
                          for rec in SeqIO.parse(locus_file, 'fasta')
                          if (rec.id).split('_')[-1] == allele_id]
            chewie_seqs[k] = allele_seq[0]

    chewie_recs = ['>{0}\n{1}'.format(k, v) for k, v in chewie_seqs.items()]
    chewie_file = os.path.join(output_dir, 'chewie_diffs.fasta')
    write_lines(chewie_recs, chewie_file)

    # map reads against chewie alleles
    chewie_paf = os.path.join(output_dir, 'chewie.paf')
    run_minimap2(chewie_file, fastq1, fastq2, chewie_paf)

    # read
    chewie_paf_lines = read_tabular(chewie_paf)

    # get coverage info
    chewie_covinfo = diffs_info(chewie_paf_lines, chewie_seqs, identity, miss_perc)

    return [chewie_covinfo, milestone_covinfo]


def main(chewie_matrix, chewie_schema, milestone_results, strain_id,
         output_dir, fastq1, fastq2, identity, missing_percentage):

    if os.path.isdir(output_dir) is not True:
        os.mkdir(output_dir)

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
    output_file = os.path.join(output_dir, strain_id+'.tsv')
    table_header = 'locus\tchewie\tmilestone'
    table_lines = ['{0}\t{1}\t{2}'.format(k, v[0], v[1])
                   for k, v in merged_results.items()]
    table_lines = [table_header] + table_lines
    write_lines(table_lines, output_file)

    # get results that differ
    diffs = {k: v for k, v in merged_results.items() if v[0] != v[1]}

    # determine coverage for differences
    final_info = diffs_coverage(diffs, chewie_schema, milestone_seqs, strain_id,
                                output_dir, fastq1, fastq2, identity, missing_percentage)

    # print diffs before completing process
    line_template = '{:<25}\t{:^6}\t{:^6}\t{:^6}\t{:^6}\t{:^6}\t{:<}'
    diffs_header = line_template.format('locus', 'tool', 'allele', 'length',
                                        'breadth_cov', 'depth_cov', 'miss')
    diffs_lines = {}
    tools = ['chewie', 'milestone']
    for i, t in enumerate(final_info):
        for k, v in t[0].items():
            allele_length = v[-1][0][1]
            if len(v[-1][1]) > 0:
                allele_length += v[-1][1][-1]
            breadth_of_coverage = v[-1][0][0]
            depth_of_coverage = v[-1][-1]
            miss_cases = t[2][k]
            miss_line = ':'.join(['({0},{1})-({2}, {3})'.format(p[0], p[1], c[0], c[1]) for p, c in miss_cases.items()])

            diffs_lines.setdefault(k, []).append(line_template.format(k, tools[i], diffs[k][i], allele_length, breadth_of_coverage, depth_of_coverage, miss_line))

    diffs_lines_ordered = []
    for k, v in diffs_lines.items():
        diffs_lines_ordered += v

    diffs_file = output_file.replace('.tsv', '_diffs.tsv')
    write_lines([diffs_header]+diffs_lines_ordered, diffs_file)


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
                        dest='output_dir',
                        help='Path to the output directory '
                             'that will be created to store output '
                             'files.')

    parser.add_argument('-fq1', type=str, required=True,
                        dest='fastq1',
                        help='FASTQ R1 file.')

    parser.add_argument('-fq2', type=str, required=True,
                        dest='fastq2',
                        help='FASTQ R2 file.')

    parser.add_argument('--idt', type=float, required=False,
                        dest='identity',
                        help='Minimum identity percentage for mapped '
                             'reads.')

    parser.add_argument('--mp', type=float, required=False,
                        default=0.10,
                        dest='missing_percentage',
                        help='Minimum frequency of indels or SNPs to be '
                             'reported.')

    args = parser.parse_args()

    return args


if __name__ == '__main__':

    args = parse_arguments()
    main(**vars(args))
