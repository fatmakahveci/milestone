#!/usr/bin/env python3

"""
author: fatmakhv
date: January 18, 2021

usage (given CDS): python coverage_stats.py -b <sample.sorted.bam> -mq <min_mapping_quality_score> --cds <cds_name> --print_pos_depth_cds --print_mean_depth_cds
usage (all CDSs): python coverage_stats.py -b <sample.sorted.bam> -mq <min_mapping_quality_score> --print_all_cds_mean_depth
"""

import argparse
import numpy as np
import os
import pandas as pd
import pysam
import sys


def check_file_exists(file_name: str) -> None:
    """ If the file doesn't exist, it exits. """

    try:
        open(file_name, 'r').close()

    except FileNotFoundError:
        print("\"{}\" file does not exist.".format(file_name))
        sys.exit(1) # 1 := There was an issue.


def check_and_create_bam_index_file(bam_file_name: str) -> None:
    """ If bam file is not indexed, it creates an index file. """

    bam_index_file = os.path.abspath(bam_file_name) + '.bai'

    if not os.path.isfile(bam_index_file) and not os.path.isfile(os.path.abspath(bam_file_name)[:-4] + '.bai'):
        print("{} file is not indexed.".format(bam_file_name))
        pysam.index(bam_file_name, bam_index_file)
        print("{} file is created.".format(bam_file_name+".bai"))


def bam_file_preprocess(file_name: str) -> None:
    """ Check existence of the bam file and whether the bam file is indexed """

    check_file_exists(file_name)
    check_and_create_bam_index_file(file_name)
    create_bam_file_with_reads_more_than_mq(file_name, args.mq)


def sort_bam(bam_file_name: str) -> None:
    """ It sorts bam file and returns none """

    pysam.sort("-o", "output.bam", bam_file_name)


def change_file_extension(in_file_name: str, file_ext: str) -> str:
    """ Replace the file's extension with the given extension """

    return in_file_name[:in_file_name.rindex('.')+1]+file_ext


def create_bam_file_with_reads_more_than_mq(bam_file_name: str, mq: int) -> str:
    """ Create an output bam file file containing reads of which
        mapping quality is more than the threshold """

    output_bam_file_name = change_file_extension(bam_file_name, f"q{str(mq)}.bam")

    bam_in = pysam.AlignmentFile(bam_file_name, "rb")
    bam_out = pysam.Samfile(output_bam_file_name, 'wb', template=bam_in)

    for read in bam_in.fetch():

        if read.flag in [99, 147, 83, 163] and read.mapping_quality >= 25: # correctly mapped reads, taken from: https://www.samformat.info/sam-format-flag 
            bam_out.write(read)
            
    bam_in.close()
    bam_out.close()

    return output_bam_file_name

def find_mean_read_depth_for_each_cds(depth_df: pd.DataFrame) -> pd.DataFrame:
    """ Find read depth for each cds """
    
    return pd.pivot_table(depth_df, values='depth', index='cds', aggfunc=np.mean, fill_value='0')
    

def find_read_depth_for_each_position_of_cds(bam_file_name: str, cds_name: str) -> pd.DataFrame:
    """ Find read depth for each position """
    
    depth_df = pd.DataFrame([line.split('\t') for line in pysam.depth(bam_file_name).split('\n') if not line.strip() == "" ], columns = ['cds', 'pos', 'depth'])
    depth_df = depth_df.astype({'cds': str, 'pos': int, 'depth': float})

    if cds_name == "":
        return depth_df
    
    return depth_df[depth_df['cds'] == cds_name]
    

def analyze_bam_file(unfiltered_file_name: str) -> None:
    """ Read sam file for statistics """

    file_name = create_bam_file_with_reads_more_than_mq(unfiltered_file_name, args.mq)
    
    if args.print_all_cds_mean_depth:

        print(find_mean_read_depth_for_each_cds(find_read_depth_for_each_position_of_cds(file_name, "")).to_string())

    if args.print_pos_depth_cds or args.print_mean_depth_cds:

        cds_pos_depth_df = find_read_depth_for_each_position_of_cds(file_name, args.cds)

        if args.print_pos_depth_cds:
            print(cds_pos_depth_df.to_string(index=False))

        if args.print_mean_depth_cds:
            print(find_mean_read_depth_for_each_cds(cds_pos_depth_df).to_string(index=False))
    

def parse_arguments() -> None:
    """ Argument parser """

    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)

    parser.add_argument('-b', '--bam', type=str,
                        required=True, dest='bam_file_name',
                        help='Name of bam file to be analyzed.')

    parser.add_argument('-c', '--cds', type=str,
                        required=False, dest='cds',
                        help='Statistics info of for a CDS'
                             'given by its ID')

    parser.add_argument('-mq', '--mapping_quality', type=int,
                        required=False, dest='mq', default=25,
                        help='Filter the reads mapping quality'
                             'is lower than the threshold.')

    parser.add_argument('-ppdc', '--print_pos_depth_cds',
                        action='store_true', dest='print_pos_depth_cds',
                        help='Create depth file and print read '
                        'depth for each position of given cds')

    parser.add_argument('-pmdc', '--print_mean_depth_cds',
                        action='store_true', dest='print_mean_depth_cds',
                        help='Create depth file and print read '
                        'depth of given cds')

    parser.add_argument('-pacmd', '--print_all_cds_mean_depth',
                        action='store_true', dest='print_all_cds_mean_depth',
                        help='Get mean read depth for all CDSs')

    args = parser.parse_args()

    return args


def main():

    bam_file_preprocess(args.bam_file_name)
    analyze_bam_file(args.bam_file_name)


if __name__ == '__main__':

    args = parse_arguments()

    sys.exit(main())
