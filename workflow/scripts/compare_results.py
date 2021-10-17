#!/usr/bin/env python3

import argparse


def get_result_dict() -> dict:

	result_dict = dict()

	with open(args.file, 'r') as file:

		for line in file.readlines():

			fields = line.rstrip('\n').replace(' ','\t').split('\t')

			result_dict[fields[0]] = { 'chew' : fields[1], args.aligner : fields[2] }

		file.close()

	return result_dict


def compare_results() -> list:

	common_lnf, diff_lnf_in_chew, diff_lnf_in_aligner, common_asm, diff_asm_in_chew, diff_asm_in_aligner, common_ref, diff_ref_in_chew, diff_ref_in_aligner, common_niph, diff_niph_in_chew, diff_niph_in_aligner, common_exc, diff_exc, common_plot, diff_plot_in_chew, diff_plot_in_aligner = 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0

	alleles_with_different_ids = []

	for cds, alleles in get_result_dict().items():

		if alleles['chew'] == 'LNF':

			if alleles[args.aligner] == 'LNF':

				common_lnf += 1

			elif alleles[args.aligner] != 'LNF':

				diff_lnf_in_chew += 1

		elif alleles['chew'] == 'ASM':

			if alleles[args.aligner] == 'ASM':

				common_asm += 1

			elif alleles[args.aligner] != 'ASM':

				diff_asm_in_chew += 1

		elif alleles['chew'] == '1':

			if alleles[args.aligner] == '1':

				common_ref += 1

			elif alleles[args.aligner] != '1':

				diff_ref_in_chew += 1

		elif alleles['chew'] in [ 'NIPH', 'NIPHEM' ]:

			if alleles[args.aligner] in [ 'NIPH', 'NIPHEM' ]:

				common_niph += 1

			elif alleles[args.aligner] not in [ 'NIPH', 'NIPHEM' ]:

				diff_niph_in_chew += 1

		elif 'PLOT' in alleles['chew']:

			if 'PLOT' in alleles[args.aligner]:

				common_plot += 1

			elif 'PLOT' not in alleles[args.aligner]:

				diff_plot_in_chew += 1

		else:

			if alleles['chew'] == alleles[args.aligner]:

				common_exc += 1

			else:

				if alleles[args.aligner] == 'LNF':

					diff_lnf_in_aligner += 1

				elif alleles[args.aligner] == 'ASM':

					diff_asm_in_aligner += 1

				elif alleles[args.aligner] in [ 'NIPH', 'NIPHEM' ]:

					diff_niph_in_aligner += 1

				elif 'PLOT' in alleles[args.aligner]:

					diff_plot_in_aligner += 1

				elif alleles[args.aligner] == '1':

					diff_ref_in_aligner += 1

				else:
					alleles_with_different_ids.append([cds, alleles['chew'], alleles[args.aligner]])
					diff_exc += 1

	return [ common_lnf, diff_lnf_in_chew, diff_lnf_in_aligner, common_asm, diff_asm_in_chew, diff_asm_in_aligner, common_ref, diff_ref_in_chew, diff_ref_in_aligner, common_niph, diff_niph_in_chew, diff_niph_in_aligner, common_exc, diff_exc, common_plot, diff_plot_in_chew, diff_plot_in_aligner, alleles_with_different_ids ]


def print_results():

	common_lnf, diff_lnf_in_chew, diff_lnf_in_aligner, common_asm, diff_asm_in_chew, diff_asm_in_aligner, common_ref, diff_ref_in_chew, diff_ref_in_aligner, common_niph, diff_niph_in_chew, diff_niph_in_aligner, common_exc, diff_exc, common_plot, diff_plot_in_chew, diff_plot_in_aligner, alleles_with_different_ids = compare_results()

	print('===========================================================')
	print(f'Total number of CDSs: {sum(compare_results()[:-1])}')
	print('===========================================================')
	print(f'Number of LNFs (common in both {args.aligner} and chewBBACA): {common_lnf}')
	print(f'Number of LNFs (LNF in chewBBACA ({args.aligner} is not LNF)):    {diff_lnf_in_chew}')
	print(f'Number of LNFs (LNF in {args.aligner} (chewBBACA is not LNF)):    {diff_lnf_in_aligner}')
	print('-----------------------------------------------------------')
	print(f'Number of ASMs (common in both {args.aligner} and chewBBACA): {common_asm}')
	print(f'Number of ASMs (ASM in chewBBACA ({args.aligner} is not ASM)):    {diff_asm_in_chew}')
	print(f'Number of ASMs (ASM in {args.aligner} (chewBBACA is not ASM)):    {diff_asm_in_aligner}')
	print('-----------------------------------------------------------')
	print(f'Number of Allele IDs Equal to the Reference Allele ID (common in both {args.aligner} and chewBBACA): {common_ref}')
	print(f'Number of Allele IDs Equal to the Reference Allele ID (Reference Allele ID in chewBBACA, different in {args.aligner}): {diff_ref_in_chew}')
	print(f'Number of Allele IDs Equal to the Reference Allele ID (Reference Allele ID in {args.aligner}, different in chewBBACA: {diff_ref_in_aligner}')
	print('-----------------------------------------------------------')
	print(f'Number of NIPHs (common in both {args.aligner} and chewBBACA): {common_niph}')
	print(f'Number of NIPHs (NIPH in chewBBACA ({args.aligner} is not NIPH)): {diff_niph_in_chew}')
	print(f'Number of NIPHs (NIPH in {args.aligner} (chewBBACA is not NIPH)):{diff_niph_in_aligner}')
	print('-----------------------------------------------------------')
	print(f'Number of PLOTs (common in both {args.aligner} and chewBBACA): {common_plot}')
	print(f'Number of PLOTs: {diff_plot_in_chew}')
	print(f'Number of PLOTs: {diff_plot_in_aligner}')
	print('-----------------------------------------------------------')
	print(f'Number of allele IDs(common in both {args.aligner} and chewBBACA): {common_exc}')
	print(f'Number of different allele IDs: {diff_exc}')

	print('===========================================================')
	print(f'Different alleles for both chewBBACA and {args.aligner}:')
	print('-----------------------------------------------------------')
	print(f'CDS\tchewBBACA\t{args.aligner}')
	print('-----------------------------------------------------------')
	for alleles in alleles_with_different_ids:
		print("\t".join(alleles))
	print('===========================================================')

if __name__ == "__main__":

	parser = argparse.ArgumentParser()

	parser.add_argument('-f', '--file', type=str, required=True, help='Name of the file. field1=<cds_name>, field2=<first_tsv>, field3=<second_tsv>')
	parser.add_argument('-a', '--aligner', type=str, required=True, help='The aligner to create the outputs. one of them [ sbg, args.aligner ]')
	args = parser.parse_args()

	print_results()

#######################
# Commands to compare #
#######################

# awk '{ print $1"\t"$2 }' ERR1624739.tsv > ERR1624739.tsv.;
# sort ERR1624739.tsv. > ERR1624739_chewbbaca.tsv

# sort ERR1624739_mlst.tsv > ERR1624739_sbg.tsv.;
# mv ERR1624739_sbg.tsv. ERR1624739_sbg.tsv;

# join -1 1 -2 1 ERR1624739_chewbbaca.tsv ERR1624739_sbg.tsv > ERR1624739_chewbbaca_sbg.tsv;

# python compare_results.py -f ERR1624739_chewbbaca_sbg.tsv -a sbg;