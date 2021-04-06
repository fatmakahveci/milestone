#!/usr/bin/env python3
# -*- coding: utf-8 -*-

'''
author: @rfm, @fatmakhv
'''


import argparse, csv, os, subprocess

from Bio import SeqIO


def read_tabular(input_file: str) -> list:
    """ Read tabular file.

        Parameters
        ----------
        input_file :
            Path to a tabular file.

        Returns
        -------
        lines :
            A list with a sublist per line in the input file.
            Each sublist has the fields that were separated by
            the defined delimiter.
    """

    with open(input_file, 'r') as infile:
        reader = csv.reader(infile, delimiter='\t')
        lines = [line for line in reader]

    return lines


def compare_novel_allele(cds: str, milestone_seq: str, chewbbaca_seq: str) -> None:
    """ Compares the 

        Parameters
        ----------
        cds :
            name of CDS to be compared to the reference sequence created by
            chewBBACA
        milestone_seq :
            FASTA sequence of CDS created by milestone
        chewbbaca_seq :
            reference FASTA sequence of CDS created by chewbbaca
    """

    # {source: sequence_from_the_source}
    seq_dict = {cds: milestone_seq, 'ref': chewbbaca_seq}
    
    for source in seq_dict.keys():

        file_name = f'{source}.fasta'
        file = open(file_name, 'w')
        file.write(f'>{source}\n{seq_dict[source]}')
        file.close()

    # HERE CORRECTION IS NEEDED
    command_list = []
    command_list.append(f"minimap2 -ax asm5 {cds}.fasta ref.fasta --cs=long -o {cds}.sam 2>&1")

    for command in command_list:
        subprocess.call(command, shell=True, stdout=subprocess.DEVNULL)

    for source in seq_dict.keys():
        os.remove(f'{source}.fasta')

def assign_alleleid(milestone_seqs: dict, chewbbaca_schema: str) -> dict:
    """ Determines if the allele predicted by Milestone
        is in the schema and attributes an allele
        identifier.

        Parameters
        ----------
        milestone_seqs :
            Dictionary with loci identifiers as keys and
            DNA sequences predicted by Milestone as values.
        chewbbaca_schema :
            Path to the schema's directory.

        Returns
        -------
        milestone_allele_ids :
            Dictionary with loci identifiers as keys and
            allele identifiers as values. If an allele
            predicted by Milestone is not in the schema,
            the value will be '?'.
    """

    milestone_allele_ids = {}
    for milestone_locid, milestone_seq in milestone_seqs.items():

        locus_file = os.path.join(chewbbaca_schema, '{0}.fasta'.format(milestone_locid))
        locus_alleles = {str(rec.seq): (rec.id).split('_')[-1]
                         for rec in SeqIO.parse(locus_file, 'fasta')}

        # get allele ID for allele predicted by Milestone
        # assign '?' if allele is not in Chewie's schema
        milestone_id = locus_alleles.get(milestone_seq, '?')
        milestone_allele_ids[milestone_locid] = milestone_id

        is_allele_in_schema = all(map(str.isdigit, milestone_id))

        # a chewBBACA allele
        # a milestone allele
        if is_allele_in_schema:
            milestone_allele_ids[milestone_locid] = milestone_id

        # a chewBBACA allele
        # "not" a milestone allele
        elif milestone_id == '-':
            milestone_allele_ids[milestone_locid] = 'LNF'

        # "not" a chewBBACA allele
        # a milestone allele
        else:
            # chewBBACA schema_seed, reference of given CDS := CDS_x_1
            ref_seq = list(locus_alleles.keys())[list(locus_alleles.values()).index('1')]
            compare_novel_allele(milestone_locid, milestone_seq, ref_seq)

    return milestone_allele_ids


def get_alleles_seqs(loci_alleles: dict, schema_path: str) -> dict:
    """ Gets the DNA sequence of one allele for each
        locus identifier that is a key in the input
        dictionary.

        Parameters
        ----------
        loci_alleles :
            Dictionary with loci ids as keys and alleles
            identifiers as values.
        schema_path :
            Path to the schema's directory.

        Returns
        -------
        loci_seqs :
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


def write_milestone_allele_ids_to_file(output_file: str, milestone_allele_ids: dict) -> None:
    """ Definition ...

        Parameters
        ----------
        output_file :
            ....
        milestone_allele_ids :
            ...

    """
    with open(output_file, 'w') as file:

        file.write(f'cds\tallele id\n')

        for cds, allele_id in milestone_allele_ids.items():

            file.write(f'{cds}\t{allele_id}\n')

        file.close()


def main(chewbbaca_matrix: str, chewbbaca_schema: str, milestone_results: str,
         strain_id: str, output_dir: str) -> None:
    """ Definition ...

        Parameters
        ----------
        chewbbaca_matrix :
            ...
        chewbbaca_schema :
            ...
        milestone_results:
            ...
        strain_id :
            ...
        output_dir :
            ...

    """
    if os.path.isdir(output_dir) is not True:
        os.mkdir(output_dir)

    # process chewbbaca's results
    # read AlleleCall matrix
    chewbbaca_lines = read_tabular(chewbbaca_matrix)

    # get line with loci IDs
    loci = chewbbaca_lines[0]
    # extract line for strain_id
    strain_profile = [l for l in chewbbaca_lines
                      if l[0].split('.fasta')[0] == strain_id][0]

    # create dictionary with chewbbaca classification per locus
    chewbbaca_allele_ids = {loci[i].split('.fasta')[0]: strain_profile[i].replace('INF-', '')
                         for i in range(1, len(loci))}

    # # get chewbbaca alleles
    # chewbbaca_seqs = get_alleles_seqs(chewbbaca_allele_ids, chewbbaca_schema)

    # read allele sequences predicted by Milestone
    milestone_seqs = {(rec.id).split('_')[0]: str(rec.seq)
                      for rec in SeqIO.parse(milestone_results, 'fasta')}

    # determine allele IDs for alleles predicted by Milestone
    milestone_allele_ids = assign_alleleid(milestone_seqs, chewbbaca_schema)

    write_milestone_allele_ids_to_file(os.path.join(output_dir, 'milestone_mlst.tsv'), milestone_allele_ids)


def parse_arguments():

    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)

    parser.add_argument('-cm', type=str, required=True,
                        dest='chewbbaca_matrix',
                        help='Path to the TSV file with allelic '
                             'profiles determined by chewBBACA.')

    parser.add_argument('-cs', type=str, required=True,
                        dest='chewbbaca_schema',
                        help='Path to chewbbacas\'s schema.')

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

    args = parser.parse_args()

    return args


if __name__ == '__main__':

    args = parse_arguments()
    main(**vars(args))
