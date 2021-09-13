###########################################
## author: @fatmakhv                     ##
## date: 23/08/2021                      ##
## aim: create info file for the sample  ##
###########################################

import argparse

from Bio import Align
from Bio.Seq import Seq
from collections import Counter


class Coverage:

    def __init__(self, coverage_line):

        fields = coverage_line.strip('\n').split('\t')

        self.cds = fields[0]
        self.length = int(fields[2])
        self.number_of_reads = int(fields[3])
        self.covered_bases = int(fields[4])
        self.coverage = float(fields[5])
        self.mean_depth = float(fields[6])
        self.mean_baseq = float(fields[7])
        self.mean_mapq = float(fields[8])

    def __repr__(self):
        return f"cds: {self.cds}\nlength: {self.length}\nnumber of reads: {self.number_of_reads}\ncovered bases: {self.covered_bases}\ncoverage: {self.coverage}\nmean depth: {self.mean_depth}\nmean baseq: {self.mean_baseq}\nmean_mapq: {self.mean_mapq}\n"


class Info:

    def __init__(self, line):

        pos_list, ref_list, alt_list, qual_list = [], [], [], []

        for variation in line.split(','):

            pos_end_idx = variation.index('*')
            pos = int(variation[ : pos_end_idx ])

            ref_end_idx = variation.index('>')
            ref = variation[ pos_end_idx +1 : ref_end_idx ]

            alt_end_idx = variation.index('-')
            alt = variation[ ref_end_idx + 1 : alt_end_idx ]
                        
            qual = variation[ alt_end_idx + 1 : ] 
            
            if len(ref) > 1 or len(alt) > 1:

                reduced_pos_list, reduced_ref_list, reduced_alt_list, reduced_qual_list = reduce_redundant_bases(pos, ref, alt, qual)

                pos_list.extend(reduced_pos_list)
                ref_list.extend(reduced_ref_list)
                alt_list.extend(reduced_alt_list)
                qual_list.extend(reduced_qual_list)

            else:

                pos_list.append(pos)
                ref_list.append(ref)
                alt_list.append(alt)
                qual_list.append(qual)

        self.pos_list = pos_list
        self.ref_list = ref_list
        self.alt_list = alt_list
        self.qual_list = qual_list

    def __repr__(self):

        return f'positions: {" ".join(list(map(str, self.pos_list)))}\n' \
                f'refs: {" ".join(self.ref_list)}\n' \
                f'alts: {" ".join(self.alt_list)}\n' \
                f'quals: {" ".join(list(map(str, self.qual_list)))}\n'


class Vcf:

    def __init__(self, vcf_line):

        fields = vcf_line.strip('\n').split('\t')

        self.chr = get_cds_name_from_allele_name(fields[0])
        self.pos = int(fields[1])
        self.id = fields[2]
        self.ref = fields[3]
        self.alt = fields[4].split(',')[0]  # in case that there is more than one alt
        self.qual = float(fields[5])
        self.filter = fields[6]
        self.info = fields[7]
        self.format = fields[8]
        self.sample = fields[9]

    def __str__(self):
        return f"{self.chr}\t{self.pos}\t{self.id}\t{self.ref}\t{self.alt}\t{self.qual}\t{self.filter}\t{self.info}\t{self.format}\t{self.sample}\n"

    def __repr__(self):
        return f"chr: {self.chr}\npos: {str(self.pos)}\nid: {self.id}\nref: {self.ref}\nalt: {self.alt}\nqual: {str(self.qual)}\nfilter: {self.filter}\ninfo: {self.info}\nformat: {self.format}\nsample: {self.sample}"


def compare_different_sized_variations(seq1: str, seq2: str) -> list:

    aligner = Align.PairwiseAligner()
    alignments = aligner.align(Seq(seq1), Seq(seq2))

    alignment = alignments[len(alignments)-1] # -1 is not working.
    seq1_aln, matching, seq2_aln = str(alignment).strip('\n').split('\n')

    idx_ref_alt_list = list()

    for idx in range(len(seq1_aln)):

        if matching[idx] != '|':

            if matching[idx] == '-':

                if seq2_aln[idx] == '-':

                    ref = seq1_aln[idx]
                    alt = '.'

                if seq1_aln[idx] == '-':

                    ref = '.'
                    alt = seq2_aln[idx]

                idx_ref_alt_list.append([idx, ref, alt])

    return idx_ref_alt_list


def reduce_redundant_bases(pos: int, ref: str, alt: str, qual: float) -> [ list, list, list, list ]:
    """
    Reduces redundant bases

    i.e.
    Input: 2333 AGAAAAAAACGAAAAAAACGAAAAAAACT AGAAAAAAACGAAAAAAACT
    Output: [2352] ['G'] ['T']

    Parameter
    ---------
    pos : int
    ref : str
    alt : str
    qual : float
    
    Return
    ------
    pos_list : list[int]
    ref_list : list[str]
    alt_list : list[str]
    qual_list : list[float]
    """

    pos_list, ref_list, alt_list, qual_list = [], [], [], []

    if len(ref) == len(alt):

        for idx in range( min(len(ref), len(alt)) ):

            if ref[idx] == alt[idx]:

                continue

            else:
                pos_list.append(pos+idx)
                ref_list.append(ref[idx])
                alt_list.append(alt[idx])
                qual_list.append(qual)

    else:

        for idx_ref_alt in compare_different_sized_variations(ref, alt):

            pos_list.append(pos+idx_ref_alt[0])
            ref_list.append(idx_ref_alt[1])
            alt_list.append(idx_ref_alt[2])
            qual_list.append(qual)

    return [ pos_list, ref_list, alt_list, qual_list ]


def get_cds_name_from_allele_name(allele_name: str) -> str:
    """
    Get <cds-name> from <cds-name_allele-id>

    Parameters
    ----------
    allele_name : <cds-name_allele-id>

    Returns
    -------
    cds_name : <cds-name>
    """

    cds_name = allele_name.strip("\n").split("_")[0]

    return cds_name


def get_allele_id_from_allele_name(allele_name: str) -> str:
    """
    Get <allele-id> from <cds-name_allele-id>

    Parameters
    ----------
    allele_name : <cds-name_allele-id>

    Returns
    -------
    allele_id : <allele-id>
    """

    allele_id = allele_name.strip("\n").split("_")[-1]

    return allele_id


def create_sample_variant_dict() -> dict:
    """
    @todo

    Return
    ------
    sample_variant_dict : @todo
    """

    sample_variant_dict = {}

    with open(args.sample_vcf, 'r') as file:

        for line in file.readlines():

            if not line.startswith("#"):

                vcf_line = Vcf(line)

                info_line = f'{vcf_line.pos}*{vcf_line.ref}>{vcf_line.alt}-{vcf_line.qual}'

                if vcf_line.chr not in sample_variant_dict.keys():

                    # allele_field # s := sample
                    sample_variant_dict[vcf_line.chr] = Info(info_line)

                else:

                    if len(vcf_line.ref) > 1 or len(vcf_line.alt) > 1:

                        reduced_pos_list, reduced_ref_list, reduced_alt_list, reduced_qual_list = reduce_redundant_bases(vcf_line.pos, vcf_line.ref, vcf_line.alt, vcf_line.qual)

                        sample_variant_dict[vcf_line.chr].pos_list.extend(reduced_pos_list)
                        sample_variant_dict[vcf_line.chr].ref_list.extend(reduced_ref_list)
                        sample_variant_dict[vcf_line.chr].alt_list.extend(reduced_alt_list)
                        sample_variant_dict[vcf_line.chr].qual_list.extend(reduced_qual_list)

                    else:

                        sample_variant_dict[vcf_line.chr].pos_list.append(vcf_line.pos)
                        sample_variant_dict[vcf_line.chr].ref_list.append(vcf_line.ref)
                        sample_variant_dict[vcf_line.chr].alt_list.append(vcf_line.alt)
                        sample_variant_dict[vcf_line.chr].qual_list.append(vcf_line.qual)

    return sample_variant_dict


def read_reference_info_txt(info_file: str) -> dict:
    """
    @todo

    Parameter
    ---------
    info_file : @todo

    Return
    ------
    reference_allele_variant_dict : @todo
    """

    reference_allele_variant_dict = {}

    with open(info_file, 'r') as file:

        for line in file.readlines():

            if "*" in line: # to check if it has any variation

                fields = line.strip('\n').split('\t')

                cds = get_cds_name_from_allele_name(fields[0])
                allele_id = get_allele_id_from_allele_name(fields[0])
            
                if cds not in reference_allele_variant_dict.keys():

                    reference_allele_variant_dict[cds] = {}
            
                reference_allele_variant_dict[cds][allele_id] = Info(fields[1])

        file.close()

    return reference_allele_variant_dict


def get_cds_coverage_info() -> dict:
    """
    Gets the output of `samtools depth` and creates cds dictionary
    for its results

    Return
    ------
    cds_depth_dict : {cds_1: coverage_info_for_cds_1, ...}
    """

    cds_depth_dict = {}

    # @todo sum all the mapped reads to each cds region
    with open(args.sample_depth, 'r') as file:

        for line in file.readlines():

            if not line.startswith('#'):

                coverage = Coverage(line)
                cds_depth_dict[coverage.cds] = coverage

        file.close()

    return cds_depth_dict


def get_reference_cds_seq_dict() -> dict:
    """
    Reads <reference.fasta> and returns cds_seq_dict

    Returns
    -------
    cds_seq_dict : { cds1: seq1, cds2: seq2, ... }
    """
    from Bio import SeqIO

    cds_seq_dict = {}

    for sequence in list(SeqIO.parse(args.reference_fasta, "fasta")):

        cds_seq_dict[sequence.id] = str(sequence.seq)

    return cds_seq_dict


def insert_variants_into_sequence(cds_reference: str, pos_list: list, ref_list: list, alt_list: list) -> str:
    """
    Takes reference sequence of CDS and inserts variants
    to create sequence with variants for CDS

    Parameters
    ----------
    cds_reference : reference sequence for given CDS
    pos_list : list of variant positions for given CDS
    ref_list : reference bases of variants in pos_list for given CDS
    alt_list : alternate bases of variants in pos_list for given CDS

    Returns
    -------
    cds_reference : CDS sequence with variants
    """
    
    # reverse list to avoid their effect on each other
    for pos, ref, alt in zip(pos_list[::-1], ref_list[::-1], alt_list[::-1]):

        if type(alt) == list:
            alt = alt[0]

        ref_start = int(pos) - 1
        ref_end = ref_start + len(ref)

        cds_reference = cds_reference[:ref_start]+alt+cds_reference[ref_end:]

    return cds_reference


def quality_check(seq: str, ref_seq: str) -> bool:
    """
    Checks its length is 3n, the first three base is for start codon,
    and the last three base is for stop codon

    Parameter
    ---------
    seq : FASTA sequence
    ref_seq : FASTA sequence for the reference of seq

    Return
    ------
    is_passed : Returns True if all is valid
    """

    # It couldn't pass quality checks.
    allele_id = "NQ"

    # check if cds length is multiple of 3
    # stop codons: seq[-3:] - 49% TAG (likely for high GC),
    #                                  32% TAA (likely for low GC), 19% TGA
    # start codons: seq[:3] - 90% MET (ATG)
    if ( len(seq) % 3 == 0) and \
            (seq[-3:] in ['TAG', 'TAA', 'TGA']) and \
            (seq[:3] in ['ATG', 'CTG', 'GTG', 'TTG'] ):

        allele_id = "Q" # It passed quality checks.
    
    elif len(seq) < ( len(ref_seq) * 0.8 ):

        allele_id = "ASM"

    elif len(seq) > ( len(ref_seq) * 1.2 ):

        allele_id = "ALM"

    
    return allele_id


def compare_ref_to_sample_variations(cds: str, cds_seq_dict: dict, reference_info : Info, sample_cds_info : Info) -> int:
    """
    @todo

    Parameter
    ---------
    cds : @todo
    cds_seq_dict : @todo
    reference_info : @todo
    sample_cds_info : @todo

    Return
    ------
    allele_id : @todo
    is_novel : @todo
    """

    is_novel = True
    allele_id = 'not_known'

    # both reference and sample alleles contain only snps
    # so direct comparison is possible
    for cds_id, ref_info in reference_info.items():

        if Counter(sample_cds_info.pos_list) == Counter(ref_info.pos_list):

            # get allele ID from the chewbbaca results
            if Counter(sample_cds_info.ref_list) == Counter(ref_info.ref_list) and Counter(sample_cds_info.alt_list) == Counter(ref_info.alt_list):

                allele_id = cds_id
                is_novel = False # allele is found among the reference alleles

                break

            # detailed search
            else:

                is_equal = False

                for idx, [ref_pos, sample_pos] in enumerate(zip(ref_info.pos_list, sample_cds_info.pos_list)):

                    # remove CC from ref_ref: ACC ref_alt: TCC and create ref_ref: A and ref_alt: T 
                    if sample_cds_info.alt_list[idx][1:] == sample_cds_info.ref_list[idx][1:] and len(sample_cds_info.alt_list[idx]) > 1 and len(sample_cds_info.ref_list[idx]) > 1:

                        sample_cds_info.alt_list[idx] = sample_cds_info.alt_list[idx][0]
                        sample_cds_info.ref_list[idx] = sample_cds_info.ref_list[idx][0]

                    # skip for the equal bases on the same position
                    if sample_cds_info.ref_list[idx] == ref_info.ref_list[idx] and sample_cds_info.alt_list[idx] == ref_info.alt_list[idx]:

                        is_equal = True # base is found on this allele

                    # detailed search for the equal variations that look like different
                    else:

                        ref_ref = ref_info.ref_list[idx]
                        ref_alt = ref_info.alt_list[idx]
                        # ref_ref <-> ref_alt to compare using the longest sequence

                        if len(ref_ref) < len(ref_alt):

                            ref_ref, ref_alt = ref_alt, ref_ref

                        sample_ref = sample_cds_info.ref_list[idx]
                        sample_alt = sample_cds_info.alt_list[idx]

                        # sample_ref <-> sample_alt to compare using the longest sequence
                        if len(sample_ref) < len(sample_alt):

                            sample_ref, sample_alt = sample_alt, sample_ref

                        ref = ref_ref[len(ref_alt):]

                        if sample_ref[:sample_ref.index(ref)]+sample_ref[sample_ref.index(ref)+len(ref):] == sample_alt:

                            is_equal = True # base is found on this allele
    
                    if not is_equal: # detailed search

                        cds_reference = cds_seq_dict[f'{cds}_1']

                        # sequence_with_variations = insert...
                        # allele_id = ALM ASM or NQ-Q quality
                        allele_id = quality_check(seq = insert_variants_into_sequence(cds_reference = cds_reference, pos_list = sample_cds_info.pos_list, ref_list = sample_cds_info.ref_list, alt_list = sample_cds_info.alt_list), ref_seq = cds_reference)

                        if allele_id == 'Q':

                            # novel allele
                            allele_id = str( max( list( map(int, reference_info.keys()) ) ) + 1 )
                            is_novel = True # this is a new allele that is not found on the reference

                    if is_equal: # all bases are equal to this reference allele's

                        allele_id = cds_id
                        is_novel = False # allele is found among the reference alleles

                        break       

    return [ allele_id, is_novel ]


def take_allele_id_for_sample_from_chewbbaca_alleles() -> dict:
    """
    @todo

    Return
    ------
    sample_allele_dict : { cds_1 : allele_ID_1, ... }
    novel_allele_variant_dict : @todo
    """
    
    sample_variant_dict = create_sample_variant_dict()

    reference_allele_variant_dict = read_reference_info_txt(info_file=args.reference_info)

    cds_seq_dict = get_reference_cds_seq_dict()

    novel_allele_variant_dict = {}

    sample_allele_dict = {}

    cds_coverage_dict = get_cds_coverage_info()

    for cds, coverage in cds_coverage_dict.items():

        sample_cds = cds.split('_')[0]

        if coverage.coverage <= 0.98:

            sample_allele_dict[cds] = 'LNF' # The CDS is not covered by the reads.

        else: # sample reads are mapped to the reference

            if sample_cds not in sample_variant_dict.keys():

                sample_allele_dict[cds] = '1' # The sample does not have any variation in the CDS so it equals to the reference

            else:

                if sample_cds not in reference_allele_variant_dict.keys(): # a novel allele for the reference

                    if len(sample_variant_dict[sample_cds].pos_list) != 0:

                        sample_allele_dict[cds] = '2' # The reference doesn't have any variation for the CDS. It is also a novel allele for the reference. The reference has only the reference sequence for the CDS.

                    else:

                        sample_allele_dict[cds] = '1'

                    if args.update_reference == 'True':

                        if len(sample_variant_dict[sample_cds].pos_list) != 0:

                            write_variations_to_reference_vcf_file(sample_cds, sample_allele_dict[cds], sample_variant_dict[sample_cds])
                            write_variations_to_reference_info_file(sample_cds, sample_allele_dict[cds], sample_variant_dict[sample_cds])

                else: # both sample and reference has the variations of this allele

                    sample_allele_dict[cds], is_novel = compare_ref_to_sample_variations(cds = sample_cds, cds_seq_dict = cds_seq_dict, reference_info = reference_allele_variant_dict[sample_cds], sample_cds_info = sample_variant_dict[sample_cds])

                    if is_novel == True and args.update_reference == 'True':

                        if len(sample_variant_dict[sample_cds].pos_list) == 0:

                            sample_allele_dict[cds] = '1' # sample doesn't have variations for the CDS so it equals to the reference CDS

                        else: # sample has variations for the CDS

                            sample_allele_dict[cds] = str( max( list( map(int, reference_allele_variant_dict[sample_cds].keys()) ) ) + 1 )
                            
                            write_variations_to_reference_vcf_file(sample_cds, sample_allele_dict[cds], sample_variant_dict[sample_cds]) 
                            write_variations_to_reference_info_file(sample_cds, sample_allele_dict[cds], sample_variant_dict[sample_cds])

    return sample_allele_dict


def write_variations_to_reference_info_file(cds: str, allele_id: str, cds_variant: Info) -> None:
    """
    @todo

    Parameter
    ---------
    cds : @todo
    allele_id : @todo
    cds_variant : @todo
    """

    with open(args.reference_info, 'a') as file:

        line = []
        for pos, ref, alt, qual in zip(cds_variant.pos_list,
                                        cds_variant.ref_list,
                                        cds_variant.alt_list,
                                        cds_variant.qual_list):

            if type(alt) is list:
                alt = ";".join(alt)

            line.append(f'{pos}*{ref[0]}>{alt}-{qual}')

        file.write(f'{cds}_{allele_id}\t{",".join(line)}\n')

        file.close()

def write_variations_to_reference_vcf_file(cds: str, allele_id: str, cds_variant: Info) -> None:
    """
    Take the variations from sample variation list
    Write from sample.vcf to reference.vcf

    Parameter
    ---------
    cds : @todo
    cds_variant : @todo
    """

    cds_vcf_line_list = []

    reference_vcf_file = open(args.reference_vcf, 'a')

    info_line_list = []

    with open(args.sample_vcf, 'r') as file:

        for line in file.readlines():

            if f'{cds}_1' in line and not line.startswith('#'):

                vcf_line = Vcf(line)
                
                vcf_line.chr = cds + '_1'
                vcf_line.info = "."

                if len(vcf_line.ref) > 1 or len(vcf_line.alt) > 1:

                    reduced_pos_list, reduced_ref_list, reduced_alt_list, reduced_qual_list = reduce_redundant_bases(vcf_line.pos, vcf_line.ref, vcf_line.alt, vcf_line.qual)

                    for i in range(len(reduced_ref_list)):

                        new_vcf_line = vcf_line

                        new_vcf_line.info = "."

                        new_vcf_line.pos = reduced_pos_list[i]
                        new_vcf_line.ref = reduced_ref_list[i]
                        new_vcf_line.alt = reduced_alt_list[i]
                        new_vcf_line.qual = reduced_qual_list[i]

                        reference_vcf_file.write(str(new_vcf_line))

                else:

                    reference_vcf_file.write(str(vcf_line))
                    
        file.close()

    reference_vcf_file.close()


if __name__ == "__main__":

    parser = argparse.ArgumentParser(add_help=True)

    parser.add_argument('--reference_fasta',  type=str, required=True,  help='Reference\'s fasta file name with its directory')
    parser.add_argument('--reference_info',   type=str, required=True,  help='Reference\'s info file name with its directory')
    parser.add_argument('--reference_vcf',    type=str, required=True,  help='Reference\'s vcf file name with its directory')
    parser.add_argument('--sample_depth',     type=str, required=True,  help='Sample\'s depth file name with its directory')
    parser.add_argument('--sample_vcf',       type=str, required=True,  help='Sample\'s vcf file name with its directory')
    parser.add_argument('--update_reference', type=str, required=False, help='Update reference\'s vcf and info file for the further analysis.')

    args = parser.parse_args()

    with open(f'{args.sample_vcf[:-4]}_mlst.tsv', 'w') as file:

        for sample_cds, allele_id in take_allele_id_for_sample_from_chewbbaca_alleles().items():

            file.write(f'{sample_cds}\t{allele_id}\n')

        file.close()
