#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
"""


import os
import re
import csv
import argparse
import subprocess
from copy import deepcopy
from collections import Counter
from itertools import groupby, chain
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
def run_minimap2(reference, fastq1, fastq2, output_file, output_type):
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
    if output_type == 'paf':
        minimap_args = ['minimap2 -I 100M --cs -cx sr {0} {1} {2} > '
                        '{3}'.format(reference, fastq1, fastq2, output_file)]
    elif output_type == 'sam':
        minimap_args = ['minimap2 -I 100M --cs -ax sr {0} {1} {2} > '
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
                Dictonary with sequence positions as keys
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
    counts = {c[0]: c[1] for c in counts}
    return [positions_depth, counts]
def determine_missing_intervals(intervals, identifier, total_len):
    """ Determines sequence intervals that are not covered by any
        reads.
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
                interval that is not covered by reads.
            not_covered : int
                Total number of bases not covered by reads.
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
def select_variants(variants, intervals, minimum_frequency, variant_frequency):
    """ Identifies variants with high frequency.
        Parameters
        ----------
        variants : dict
            A dictionary with loci identifiers as keys and a list
            as value with identified variants and the frequency for
            each variant.
        intervals : dict
            A dictionary with loci identifiers as keys and a
            list as value. Each list has a set of intervals with
            start positions, stop position and a dictionary with
            the depth of coverage per position.
        minimum_frequency : float
            Minimum coverage/frequency of a SNP or indel to be
            selected.
        variant_frequency :
            Variants with frequency equal or greater than this value
            are selected as probable variants.
        Returns
        -------
        selected_variants: dict
            A dictionary with loci identifiers as keys and a
            dictionary with position:variant tuples as keys and
            with a position_coverage:variant_coverage tuple as
            value. Stores variants with frequency equal or
            greater than value passed to minimum_frequency.
        high_frequency : dict
            Same structure as previous returned variable but
            only stores variants with frequency value equal
            or greater than value passed to variant_frequency.
    """
    # select indels or SNPs based on frequency
    high_frequency = {}
    selected_variants = {}
    # for each locus and detected indels/SNPs
    for k, v in variants.items():
        counts = v[1]
        selected = {}
        probable = {}
        # for each detected indel/SNP
        for p, c in counts.items():
            # get coverage values for locus positions
            coverage_values = intervals[k]
            # determine which region the indel/SNP is in
            for region in coverage_values:
                if p[0] >= region[0] and p[0] <= region[1]:
                    pos_cov = region[2][p[0]]
                    # store indel/SNP info if coverage value for the position
                    # is 0 or the indel/SNP has high frequency
                    frequency = (c/(c+pos_cov))
                    if pos_cov == 0 or frequency >= minimum_frequency:
                        # key has position and variant (ref to reads)
                        # value has position coverage and indel/SNP frequency
                        selected[p] = (pos_cov, c)
                    if pos_cov == 0 or frequency >= variant_frequency:
                        probable[p] = (pos_cov, c)
        selected_variants[k] = selected
        high_frequency[k] = probable
    return [selected_variants, high_frequency]
def process_paf(paf_lines, seqs, identity, minimum_frequency, low_coverage,
                variant_frequency):
    """ Processes results in a PAF file in order to determine
        coverage statistics and identify variants.
        Parameters
        ----------
        paf_lines : list
            A list with all lines in the PAF file generated
            by minimap2.
        seqs : dict
            Dictionary with loci identifiers as keys and DNA
            sequences of the alleles predicted by Chewie or
            Milestone as values.
        identity : float
            Minimum sequence identity between reference sequences
            and mapped reads. Mappings below this value are excluded.
        minimum_frequency : float
            Minimum coverage/frequency of a SNP or indel to be
            selected.
        low_coverage : int
            Positions with coverage value below this value are
            identified as low coverage regions.
        variant_frequency : float
            Variants with frequency equal or greater than this value
            are sleected as probable variants.
        Returns
        -------
        A list with the following elements:
            - A dictionary with loci identifiers as keys and a
            list as value. Each list has a set of intervals with
            start positions, stop position and a dictionary with
            the depth of coverage per position.
            - A list with the mappings that were below the identity
            threshold.
            - A dictionary with loci identifiers as keys and a
            dictionary with position:variant tuples as keys and
            with a position_coverage:variant_coverage tuple as
            value.
            - A dictionary with loci identifiers as keys and a list
            as value. Each list has the following elements: a list
            with the breadth of coverage value and the number of
            covered bases. A list with a dictionary with the locus
            identifier as key and a list with missing positions as
            value and the number of missing positions. Mean depth
            of coverage. Number of positions with 0 coverage.
            Number of positions below minimum coverage threshold.
    """
    mapped_reads = deepcopy(paf_lines)
    # compute alignment identity
    for i in range(len(mapped_reads)):
        mapped_reads[i].append(int(mapped_reads[i][9]) / int(mapped_reads[i][10]))
    # filter out alignments below defined identity
    # this kepps reads that aligned well at sequence extremities
    invalid_mappings = [line for line in mapped_reads if line[-1] < identity]
    mapped_reads = [line for line in mapped_reads if line[-1] >= identity]
    # match alignment string with regex
    pattern = r':[0-9]+|\*[a-z][a-z]|\+[a-z]+|-[a-z]+'
    for i in range(len(mapped_reads)):
        current = mapped_reads[i][-2]
        mapped_reads[i].append(regex_matcher(current, pattern))
    # get information about positions that match to determine coverage
    for i in range(len(mapped_reads)):
        current = mapped_reads[i][-1]
        start = int(mapped_reads[i][7])
        mapped_reads[i].append(single_position_coverage(current, start))
    # get indels and SNPs for each locus
    loci_miss = {}
    for l in mapped_reads:
        loci_miss.setdefault(l[5], []).extend(l[-1][1])
    # indels and SNPs counts
    loci_miss = {k: [v, Counter(v)] for k, v in loci_miss.items()}
    # identify subsequences that are well covered by reads
    covered_intervals = {}
    for l in mapped_reads:
        covered_intervals.setdefault(l[5], []).append([int(l[7]), int(l[8]), l[-1][0]])
    # sort covered intervals
    covered_intervals_sorted = {k: sorted(v, key=lambda x: x[0])
                                for k, v in covered_intervals.items()}
    # merge overlapping intervals
    # deepcopy to avoid altering original intervals
    merged_intervals = {k: merge_intervals(v)
                        for k, v in covered_intervals_sorted.items()}
    # identify variants and highly probable variants
    selected_miss, probable_variants = select_variants(loci_miss,
                                                       merged_intervals,
                                                       minimum_frequency,
                                                       variant_frequency)
    coverage_stats = {k: [] for k in merged_intervals}
    for k, v in merged_intervals.items():
        # breadth of coverage and number of covered bases
        coverage = determine_breadth_coverage({k: v}, len(seqs[k]))
        # determine subsequences that are not covered
        missing = determine_missing_intervals(v, k, len(seqs[k]))
        depth_info = determine_depth_coverage(v, len(seqs[k]))
        # determine mean depth of coverage
        depth_counts = depth_info[1]
        # get number of positions with 0 coverage
        zero_coverage = depth_counts.get(0, 0)
        # get number of positions below low coverage threshold
        below_coverage = sum([v for k, v in depth_counts.items() if k < low_coverage])
        # get mean depth of coverage
        depth_sum = sum([k*v for k, v in depth_counts.items()])
        mean_depth = round(depth_sum/len(seqs[k]), 4)
        coverage_stats[k].extend([coverage, missing, mean_depth,
                                  zero_coverage, below_coverage])
    return [merged_intervals, invalid_mappings, selected_miss,
            coverage_stats, probable_variants]
def map_reads(reference_fasta, fastq1, fastq2, output_dir, sam_output,
              output_prefix):
    """ Executes minimap2 to map reads against sequences in a
        Fasta file.
        Parameters
        ----------
        reference_fasta : str
            Reference Fasta file the reads will be mapped
            against.
        fastq1 : str
            Path to the R1 FASTQ file.
        fastq2 : str
            Path to the R2 FASTQ file.
        output_dir : str
            Path to the output directory where files are
            created.
        sam_output : bool
            True if mapping results should also be generated
            in SAM format.
        output_prefix : str
            Prefix to include in the name of output files.
        Returns
        -------
        paf_lines : list
            A list with all lines in the PAF file generated
            by minimap2.
    """
    # map reads and create output in PAF format
    paf_outfile = os.path.join(output_dir, output_prefix+'.paf')
    paf_std = run_minimap2(reference_fasta, fastq1, fastq2, paf_outfile, 'paf')
    # read data from mapped reads
    paf_lines = read_tabular(paf_outfile)
    if sam_output is True:
        # create SAM file to view in IGV
        sam_outfile = os.path.join(output_dir, output_prefix+'.sam')
        sam_std = run_minimap2(reference_fasta, fastq1, fastq2,
                               sam_outfile, 'sam')
    return paf_lines
def differences_coverage(differences, chewie_schema, milestone_seqs,
                         chewie_seqs, strain_id, output_dir, fastq1,
                         fastq2, identity, minimum_frequency, low_coverage,
                         variant_frequency, sam_output):
    """ Determines coverage statistics and detected variants for
        loci with predictions that differ between Chewie and Milestone
        results.
        Parameters
        ----------
        differences : dict
            Dictionary with loci identifiers as keys and a list
            as value. The list has the allele identifier attributed
            by Chewie (0 if missing data) and the allele identifier
            of the allele predicted by Milestone (? if allele is
            not in Chewie's schema).
        chewie_schema : str
            Path to the schema's directory.
        milestone_seqs : dict
            Dictionary with loci identifiers as keys and DNA
            sequences of the alleles predicted by Milestone.
        strain_id : str
            Strain identifier.
        output_dir : str
            Path to the output directory where files will
            be created.
        fastq1 : str
            Path to the R1 FASTQ file.
        fastq2 : str
            Path to the R2 FASTQ file.
        identity : float
            Minimum sequence identity between reference sequences
            and mapped reads. Mappings below this value are excluded.
        minimum_frequency : float
            Minimum coverage/frequency of a SNP or indel to be
            selected.
        low_coverage : int
            Positions with coverage value below this value are
            identified as low coverage regions.
        variant_frequency : float
            Variants with frequency equal or greater than this value
            are sleected as probable variants.
        Returns
        -------
        A list with the variables returned by the `process_paf`
        function for Chewie's results and Milestone's results.
    """
    # get alleles predicted by milestone
    milestone_diffs = {k: milestone_seqs[k] for k, v in differences.items() if v[1] != '0'}
    # write milestone alleles to file
    milestone_recs = ['>{0}\n{1}'.format(k, v) for k, v in milestone_diffs.items()]
    milestone_file = os.path.join(output_dir, 'milestone_diffs.fasta')
    write_lines(milestone_recs, milestone_file)
    # map reads against milestone alleles
    milestone_lines = map_reads(milestone_file, fastq1, fastq2,
                                output_dir, sam_output, 'milestone')
    # get coverage info
    milestone_covinfo = process_paf(milestone_lines, milestone_seqs,
                                    identity, minimum_frequency,
                                    low_coverage, variant_frequency)
    chewie_diffs = {k: chewie_seqs[k] for k, v in differences.items() if v[0] != '0'}
    # write chewie alleles to file
    chewie_recs = ['>{0}\n{1}'.format(k, v) for k, v in chewie_diffs.items()]
    chewie_file = os.path.join(output_dir, 'chewie_diffs.fasta')
    write_lines(chewie_recs, chewie_file)
    # map reads against chewie alleles
    chewie_lines = map_reads(chewie_file, fastq1, fastq2,
                             output_dir, sam_output, 'chewie')
    # get coverage info
    chewie_covinfo = process_paf(chewie_lines, chewie_seqs,
                                 identity, minimum_frequency,
                                 low_coverage, variant_frequency)
    # each element in list has coverage info for one tool
    # coverage per position
    # invalid mappings based on identity threshold
    # detected variants that have frequency above minimum_frequency
    # breadth of coverage, covered bases, regions not covered, depth
    # of coverage
    # positions with 0 coverage and positions with coverage below low_coverage
    # variants with frequency above variant_frequency
    return [chewie_covinfo, milestone_covinfo]
def create_diff_lines(coverage_info, differences):
    """ Creates lines for the output TSV file with
        the list of differences and coverage statsitics.
        Parameters
        ----------
        coverage_info : list
            A list with the variables returned by the
            `differences_coverage` function.
        differences : dict
            Dictionary with loci identifiers as keys and a list
            as value. The list has the allele identifier attributed
            by Chewie (0 if missing data) and the allele identifier
            of the allele predicted by Milestone (? if allele is
            not in Chewie's schema).
        Returns
        -------
        differences_lines : list
            List with two elements: the string/header of
            the TSV file and a dictionary with loci
            identifiers as keys and strings/lines for the
            classifications that differe between Chewie
            and Milestone.
    """
    line_template = '{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}'
    columns_titles = ['locus', 'tool', 'allele', 'length', 'breadth_cov',
                      'depth_cov', 'zero_cov', 'below_cov',
                      'high_variants', 'variants']
    header = line_template.format(*columns_titles)
    differences_lines = [header, {}]
    tools = ['chewie', 'milestone']
    for i, t in enumerate(coverage_info):
        tool = tools[i]
        for k, v in t[3].items():
            allele_length = v[0][1] + v[1][-1]
            breadth_of_coverage = v[0][0]
            depth_of_coverage = v[2]
            zero_cov = v[3]
            below_cov = v[4]
            high_cov_variants = t[4][k]
            high_cov_variants_line = ':'.join(['({0},{1})({2},{3})'.format(*p, *c)
                                              for p, c in high_cov_variants.items()])
            variants = t[2][k]
            variants_line = ':'.join(['({0},{1})({2},{3})'.format(*p, *c)
                                      for p, c in variants.items()])
            line = line_template.format(k, tool, differences[k][i],
                                        allele_length, breadth_of_coverage,
                                        depth_of_coverage, zero_cov, below_cov,
                                        high_cov_variants_line, variants_line)
            differences_lines[1].setdefault(k, []).append(line)
    return differences_lines
def assign_alleleid(milestone_seqs, chewie_schema):
    """ Determines if the allele predicted by Milestone
        is in the schema and attributes an allele
        identifier.
        Parameters
        ----------
        milestone_seqs : dict
            Dictionary with loci identifiers as keys and
            DNA sequences predicted by Milestone as values.
        chewie_schema : str
            Path to the schema's directory.
        Returns
        -------
        milestone_allele_ids : dict
            Dictionary with loci identifiers as keys and
            allele identifiers as values. If an allele
            predicted by Milestone is not in the schema,
            the value will be '?'.
    """
    milestone_allele_ids = {}
    for locid, seq in milestone_seqs.items():
        locus_file = os.path.join(chewie_schema, '{0}.fasta'.format(locid))
        locus_alleles = {str(rec.seq): (rec.id).split('_')[-1]
                         for rec in SeqIO.parse(locus_file, 'fasta')}
        # get allele ID for allele predicted by Milestone
        # assign '?' if allele is not in Chewie's schema
        milestone_id = locus_alleles.get(seq, '?')
        milestone_allele_ids[locid] = milestone_id
    return milestone_allele_ids
def get_alleles_seqs(loci_alleles, schema_path):
    """ Gets the DNA sequence of one allele for each
        locus identifier that is a key in the input
        dictionary.
        Parameters
        ---------
        loci_alleles : dict
            Dictionary with loci ids as keys and alleles
            identifiers as values.
        schema_path : str
            Path to the schema's directory.
        Returns
        -------
        loci_seqs : dict
            Dictionary with loci identifiers as keys and
            DNA sequences as values.
    """
    loci_seqs = {}
    for locus, alleleid in loci_alleles.items():
        locus_file = os.path.join(schema_path, locus+'.fasta')
        allele_seq = [str(rec.seq)
                      for rec in SeqIO.parse(locus_file, 'fasta')
                      if (rec.id).split('_')[-1] == alleleid]
        if len(allele_seq) > 0:
            loci_seqs[locus] = allele_seq[0]
    return loci_seqs
def main(chewie_matrix, chewie_schema, milestone_results, strain_id,
         output_dir, fastq1, fastq2, identity, minimum_frequency,
         low_coverage, variant_frequency, sam_output):
    if os.path.isdir(output_dir) is not True:
        os.mkdir(output_dir)
    # process chewie's results
    # read AlleleCall matrix
    chewie_lines = read_tabular(chewie_matrix)
    # get line with loci IDs
    loci = chewie_lines[0]
    # extract line for strain_id
    strain_profile = [l for l in chewie_lines
                      if l[0].split('.fasta')[0] == strain_id][0]
    # create dictionary with chewie classification per locus
    chewie_allele_ids = {loci[i].split('.fasta')[0]: strain_profile[i].replace('INF-', '')
                         for i in range(1, len(loci))}
    # get chewie alleles
    chewie_seqs = get_alleles_seqs(chewie_allele_ids, chewie_schema)
    # read allele sequences predicted by Milestone
    milestone_seqs = {(rec.id).split('_')[0]: str(rec.seq)
                      for rec in SeqIO.parse(milestone_results, 'fasta')}
    # determine allele IDs for alleles predicted by Milestone
    milestone_allele_ids = assign_alleleid(milestone_seqs, chewie_schema)
    # create dictionary with loci IDs as keys and chewie and
    # milestone classifications as values
    merged_classifications = {k: [v, milestone_allele_ids.get(k, '-')]
                              for k, v in chewie_allele_ids.items()}
    # write results to output file
    classifications_file = os.path.join(output_dir, strain_id+'.tsv')
    classifications_header = 'locus\tchewie\tmilestone'
    classifications_lines = ['{0}\t{1}\t{2}'.format(k, v[0], v[1])
                             for k, v in merged_classifications.items()]
    classifications_lines = [classifications_header] + classifications_lines
    write_lines(classifications_lines, classifications_file)
    # identify loci for which Chewie and Milestone
    # attributed different classifications
    differences = {k: v
                   for k, v in merged_classifications.items()
                   if v[0] != v[1]}
    # get more information about loci predictions that chewie
    # and milestone do not agree on
    # determine coverage for each allele and SNPs or indels that
    # have a frequency higher than defined threshold
    coverage_info = differences_coverage(differences, chewie_schema,
                                         milestone_seqs, chewie_seqs,
                                         strain_id, output_dir, fastq1,
                                         fastq2, identity, minimum_frequency,
                                         low_coverage, variant_frequency,
                                         sam_output)
    # save info about differences
    differences_lines = create_diff_lines(coverage_info, differences)
    # group
    ordered_lines = list(chain(*[v for k, v in differences_lines[1].items()]))
    ordered_lines = [differences_lines[0]] + ordered_lines
    differences_file = classifications_file.replace('.tsv', '_diffs.tsv')
    write_lines(ordered_lines, differences_file)
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
    parser.add_argument('--mf', type=float, required=False,
                        default=0.10,
                        dest='minimum_frequency',
                        help='Minimum relative frequency of indels or '
                             'SNPs that are reported.')
    parser.add_argument('--lc', type=int, required=False,
                        default=10,
                        dest='low_coverage',
                        help='Minimum coverage value. Positions with '
                             'smaller values are reported as low coverage '
                             'positions.')
    parser.add_argument('--vf', type=float, required=False,
                        default=0.30,
                        dest='variant_frequency',
                        help='Indels or SNPs with a relative frequency '
                             'that is equal or greater that this value are '
                             'reported as highly probable.')
    parser.add_argument('--sam', required=False, action='store_true',
                        dest='sam_output',
                        help='Create SAM file with mapping results.')
    args = parser.parse_args()
    return args
if __name__ == '__main__':
    args = parse_arguments()
    main(**vars(args))
