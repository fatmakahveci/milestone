#!/usr/bin/env python3

# Author: fatmakhv

## >> pip3 install pysam pysamstats matplotlib docopt
## >> python coverage_stats --bam ERR3464558.sorted.bam --mq 25 --cds ERR3464558-protein1_1 --print_aln_mq --plot

"""
------------------
Usage:
  coverage_stats.py --bam=<bam_file> --cds=<cds_id> [--mq=<mapping_quality_phred>] [--print_aln_mq] [--plot]

Options:
  --help                        Show this message.
  --mq=<mapping_quality_phred>  It filters the reads mapping quality is lower than the threshold [default: 30].
  --cds=<cds_id>                It gets the id of cds of which statistics info is searched for [default: ""].
  --print_aln_mq                To print the lines of which mapping quality is more than given [default: False].
  --plot                        Show the plot for given cds [default: False].
------------------
"""

from docopt import docopt
import os.path, pysam, pysamstats, sys, matplotlib.pyplot as plt


# to check whether whether the file is found.
def check_file_exists(file_name):

	try:
		open(file_name, 'r').close()
	except FileNotFoundError:
		print("\"{}\" file does not exist.".format(file_name))
		sys.exit()


# if bam file is not indexed, it creates an index file
def generate_index_if_needed(bam_file_name):

	bam_index_file = os.path.abspath(bam_file_name) + '.bai'

	if not os.path.isfile(bam_index_file) and not os.path.isfile(os.path.abspath(bam_file_name)[:-4] + '.bai'):
		print("{} file is not indexed.".format(bam_file_name))
		pysam.index(bam_file_name, bam_index_file)
		print("{} file is created.".format(bam_file_name+".bai"))

	return True


# bam existence check and indexing
def bam_file_preprocess():

	check_file_exists(bam_file_name)
	generate_index_if_needed(bam_file_name)


# taken from: https://www.samformat.info/sam-format-flag
def get_flag_stats(flag_info, flag_id):

	if flag_id in [73, 133, 89, 121, 165, 181, 101, 117, 153, 185, 69, 137]:
		flag_info[1] += 1
	elif flag_id in [77, 141]:
		flag_info[2] += 1
	elif flag_id in [99, 147, 83, 163]:
		flag_info[3] += 1
	elif flag_id in [67, 131, 115, 179]:
		flag_info[4] += 1
	elif flag_id in [81, 161, 97, 145, 65, 129, 113, 177]:
		flag_info[5] += 1
	else:
		flag_info[6] += 1

# get detailed information for flags in lines
def get_flag_info(flag_stats_id):
	if flag_stats_id == 1:
		return "Number of lines in sam file where one of the reads is unmapped:"
	elif flag_stats_id == 2:
		return "Number of lines in sam file where both reads are unmapped:"
	elif flag_stats_id == 3:
		return "Number of lines in sam file where mapped within the insert size and in correct orientation:"
	elif flag_stats_id == 4:
		return "Number of lines in sam file where mapped within the insert size but in wrong orientation:"
	elif flag_stats_id == 5:
		return "Number of lines in sam file where mapped uniquely, but with wrong insert size:"
	else:
		return "Others, not a common flag:"


# read sam file for statistics
def read_sam_file():

	sam_file = pysam.AlignmentFile(bam_file_name, "rb")

	flag_info = { 1:0, 2:0, 3:0, 4:0, 5:0, 6:0 }

	flag_info_detailed_dict = {}

	cds = args["--cds"]
	for read in sam_file.fetch(cds):
		if read.mapping_quality >= int(args["--mq"]): # Illumina

			get_flag_stats(flag_info, read.flag)
			if args["--print_aln_mq"]:
				print(cds, read.pos, read.flag, read.seq)

	for key, value in flag_info.items():
		print(get_flag_info(key), value)

	sam_file.close()

	if args["--plot"]:
		bam_file = pysam.AlignmentFile(bam_file_name)
		a = pysamstats.load_coverage(bam_file, chrom=cds)
		plt.plot(a.pos, a.reads_pp)
		plt.show()
		bam_file.close()

# main method
if __name__ == "__main__":

	args = docopt(__doc__)

	bam_file_name = args["--bam"]

	bam_file_preprocess()

	read_sam_file()
