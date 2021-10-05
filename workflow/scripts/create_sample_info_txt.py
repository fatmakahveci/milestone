#################################################################################
## author: @fatmakhv                                                           ##
## date: 04/10/2021                                                            ##
## aim: create info file for the sample and write novel alleles to schema seed ##
#################################################################################

import argparse, glob, os, pysam


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


def create_sample_variation_dict(ref_allele_id: str) -> dict:
    """
    Creates variation dictionary

    Return
    ------
    sample_variation_dict : variations of sample dictionary for alleles
    """

    sample_variation_dict = {}

    with open(args.sample_vcf, 'r') as file:

        for line in file.readlines():

            if not line.startswith("#"):

                vcf_line = Vcf(line)

                # dominant_base := base_ratio_check(vcf_line.pos, f'{vcf_line.chr}_{ref_allele_id}')
                if vcf_line.ref[0] == base_ratio_check(vcf_line.pos, f'{vcf_line.chr}_{ref_allele_id}'):
                    continue

                info_line = f'{vcf_line.pos}*{vcf_line.ref}>{vcf_line.alt}-{vcf_line.qual}'

                if vcf_line.chr not in sample_variation_dict.keys():

                    # allele_field # s := sample
                    sample_variation_dict[vcf_line.chr] = Info(info_line)

                else:

                    sample_variation_dict[vcf_line.chr].pos_list.append(vcf_line.pos)
                    sample_variation_dict[vcf_line.chr].ref_list.append(vcf_line.ref)
                    sample_variation_dict[vcf_line.chr].alt_list.append(vcf_line.alt)
                    sample_variation_dict[vcf_line.chr].qual_list.append(vcf_line.qual)

    return sample_variation_dict


def read_reference_info_txt(info_file: str) -> dict:
    """
    Read reference_info.txt
    Create dictionary for allele of variations in CDS

    Parameter
    ---------
    info_file : Name of reference_info.txt to write

    Return
    ------
    reference_allele_variation_dict : allele of variations in CDS
    """

    reference_allele_variation_dict = {}

    with open(info_file, 'r') as file:

        for line in file.readlines():

            if "*" in line: # to check if it has any variation

                fields = line.strip('\n').split('\t')

                cds = get_cds_name_from_allele_name(fields[0])
                allele_id = get_allele_id_from_allele_name(fields[0])
            
                if cds not in reference_allele_variation_dict.keys():

                    reference_allele_variation_dict[cds] = {}
            
                reference_allele_variation_dict[cds][allele_id] = Info(fields[1])

        file.close()

    return reference_allele_variation_dict


def get_cds_coverage_info() -> dict:
    """
    Gets the output of `samtools depth` and creates cds dictionary
    for its results

    Return
    ------
    cds_depth_dict : {cds_1: coverage_info_for_cds_1, ...}
    """

    cds_depth_dict = {}


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


def insert_variations_into_sequence(cds_reference: str, pos_list: list, ref_list: list, alt_list: list) -> str:
    """
    Takes reference sequence of CDS and inserts variations
    to create sequence with variations for CDS

    Parameters
    ----------
    cds_reference : reference sequence for given CDS
    pos_list : list of variation positions for given CDS
    ref_list : reference bases of variations in pos_list for given CDS
    alt_list : alternate bases of variations in pos_list for given CDS

    Returns
    -------
    cds_reference : CDS sequence with variations
    """
    
    # reverse list to avoid their effect on each other
    for pos, ref, alt in zip(pos_list[::-1], ref_list[::-1], alt_list[::-1]):

        pos -= 1

        if type(alt) == list:

            alt = alt[0]

        cds_reference = cds_reference[:pos]+alt.replace('.','')+cds_reference[pos+len(ref):]

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

    if seq == ref_seq:

        return "EQ"

    if len(seq) < len(ref_seq) * 0.8:

        return "ASM"

    if len(seq) > len(ref_seq) * 1.2:

        return "ALM"

    for i in range(0, len(seq)-3, 3):

        if seq[i:i+3] == 'TAG' or seq[i:i+3] == 'TAA' or seq[i:i+3] == 'TGA':

            return "LNF" # in-frame stop codon

    # check if cds length is multiple of 3
    # stop codons: seq[-3:] - 49% TAG (likely for high GC),
    #                                  32% TAA (likely for low GC), 19% TGA
    # start codons: seq[:3] - 90% MET (ATG)
    if seq[-3:] in ['TAG', 'TAA', 'TGA'] and seq[:3] in ['ATG', 'CTG', 'GTG', 'TTG']:

        return "Q" # It passed quality checks.

    return "LNF"


def compare_ref_to_sample_variations(cds: str, ref_allele_id: str, cds_seq_dict: dict, reference_info : Info, sample_cds_info : Info) -> int:
    """
    Compare reference variations to sample variations

    Parameter
    ---------
    cds : Name of CDS
    cds_seq_dict : Reference sequence dict of CDSs
    reference_info : reference variations
    sample_cds_info : sample_cds variations

    Return
    ------
    allele_id : allele ID after comparison
    is_novel : check if allele is novel or not
    """

    from collections import Counter

    # cds length check
    diff_len = 0

    for ref, alt in zip(sample_cds_info.ref_list, sample_cds_info.alt_list):

        if alt == '.':

            diff_len += len(ref)

        else:

            diff_len += len(ref) - len(alt)

    if diff_len % 3 != 0:

        is_novel = False
        allele_id = 'LNF' # incorrect length

    else:

        is_novel = True

        cds_reference = cds_seq_dict[f'{cds}_{ref_allele_id}']

        allele_id = quality_check(insert_variations_into_sequence(cds_reference, sample_cds_info.pos_list, sample_cds_info.ref_list, sample_cds_info.alt_list), cds_reference)
        
        if allele_id == "EQ":

            is_novel = False
            allele_id = ref_allele_id

        elif allele_id == 'Q':

            # both reference and sample alleles contain only snps
            # so direct comparison is possible
            for cds_id, ref_info in reference_info.items():

                # if positions are not equal there is nothing to compare
                if Counter(sample_cds_info.pos_list) == Counter(ref_info.pos_list):

                    # get allele ID from the chewbbaca results
                    if Counter(sample_cds_info.ref_list) == Counter(ref_info.ref_list) and Counter(sample_cds_info.alt_list) == Counter(ref_info.alt_list):

                        is_novel = False # allele is found among the reference alleles
                        allele_id = cds_id
                    
                        # allele is found on the reference
                        break
                
                    else:

                        is_novel = True
                        allele_id = str(novel_allele_id_of_cds_dict[cds])

                else:

                    is_novel = True
                    allele_id = str(novel_allele_id_of_cds_dict[cds])

        else:

            is_novel = False

    return [ allele_id, is_novel ]


def take_allele_id_for_sample_from_chewbbaca_alleles() -> dict:
    """
    Returns allele ID for samples coding sequences

    Return
    ------
    sample_allele_dict : { cds_1 : allele_ID_1, ... }
    """
    
    ref_allele_id = '1' # cds.split('_')[1]
    
    sample_variation_dict = create_sample_variation_dict(ref_allele_id)

    reference_allele_variation_dict = read_reference_info_txt(info_file=args.reference_info)

    cds_seq_dict = get_reference_cds_seq_dict()

    sample_allele_dict = {}

    cds_coverage_dict = get_cds_coverage_info()

    for cds, coverage in cds_coverage_dict.items():

        sample_cds = cds.split('_')[0]

        if coverage.coverage <= 0.98:

            # CDS is not covered by the reads.
            sample_allele_dict[cds] = 'LNF'

        # Sample reads are mapped to the reference.
        else:

            if sample_cds not in sample_variation_dict.keys():

                # Sample does not have any variation in CDS so it equals to the reference.
                sample_allele_dict[cds] = ref_allele_id

            else:

                is_novel = False

                # It is a novel allele for the reference.
                if sample_cds not in reference_allele_variation_dict.keys():
                    
                    # Sample has variations for this CDS.
                    if len( sample_variation_dict[sample_cds].pos_list ) != 0:

                        # cds length check
                        diff_len = 0

                        for ref, alt in zip( sample_variation_dict[sample_cds].ref_list, sample_variation_dict[sample_cds].alt_list ):

                            if alt == '.':

                                diff_len += len(ref)

                            else:

                                diff_len += len(ref) - len(alt)

                        if diff_len % 3 != 0:

                            is_novel = False
                            sample_allele_dict[cds] = 'LNF' # incorrect length

                        else:

                            cds_reference = cds_seq_dict[f'{sample_cds}_{ref_allele_id}']

                            sample_allele_dict[cds] = quality_check(insert_variations_into_sequence(cds_reference, sample_variation_dict[sample_cds].pos_list, sample_variation_dict[sample_cds].ref_list, sample_variation_dict[sample_cds].alt_list), cds_reference)

                            if sample_allele_dict[cds] == "EQ": # sample equals to the reference

                                is_novel = False
                                sample_allele_dict[cds] = ref_allele_id

                            if sample_allele_dict[cds] == 'Q': # it passed 3n, start, stop, in-frame stop

                                is_novel = True

                                # The reference doesn't have any variation for the CDS.
                                # It is also a novel allele for the reference.
                                # The reference has only the reference sequence for the CDS.
                                sample_allele_dict[cds] = '2'

                            if sample_allele_dict[cds] in [ 'ASM', 'ALM', 'LNF', 'EQ']:

                                is_novel = False
                    else:

                        # Reference allele variation set is empty so sample's CDS equals to the reference.
                        sample_allele_dict[cds] = ref_allele_id

                    if args.update_reference == 'True' and is_novel == True:

                        if len(sample_variation_dict[sample_cds].pos_list) != 0:

                            write_allele_sequence_to_schema_seed(sample_cds, ref_allele_id, cds_seq_dict[f'{sample_cds}_{ref_allele_id}'], sample_variation_dict[sample_cds])
                            write_variations_to_reference_vcf_file(sample_cds, ref_allele_id, sample_variation_dict[sample_cds])
                            write_variations_to_reference_info_file(sample_cds, sample_allele_dict[cds], sample_variation_dict[sample_cds])

                else: # both sample and reference has the variations of this allele

                    sample_allele_dict[cds], is_novel = compare_ref_to_sample_variations(cds = sample_cds, ref_allele_id = ref_allele_id, cds_seq_dict = cds_seq_dict, reference_info = reference_allele_variation_dict[sample_cds], sample_cds_info = sample_variation_dict[sample_cds])

                    if is_novel == True and args.update_reference == 'True':

                        if len(sample_variation_dict[sample_cds].pos_list) == 0:

                            sample_allele_dict[cds] = ref_allele_id # sample doesn't have variations for the CDS so it equals to the reference CDS

                        else: # sample has variations for the CDS

                            sample_allele_dict[cds] = novel_allele_id_of_cds_dict[sample_cds]

                            write_allele_sequence_to_schema_seed(sample_cds, sample_allele_dict[cds], cds_seq_dict[f'{sample_cds}_{ref_allele_id}'], sample_variation_dict[sample_cds])
                            write_variations_to_reference_vcf_file(sample_cds, ref_allele_id, sample_variation_dict[sample_cds]) 
                            write_variations_to_reference_info_file(sample_cds, sample_allele_dict[cds], sample_variation_dict[sample_cds])

    return sample_allele_dict


def write_allele_sequence_to_schema_seed(sample_cds: str, cds_allele_id: str, sample_ref_seq: str, sample_cds_variation: Info) -> None:
    """
    Write variations to the reference_info.txt file.

    Parameter
    ---------
    sample_cds : Name of CDS
    cds_allele_id : Allele ID of CDS
    sample_ref_seq : Sequence of CDS reference
    sample_cds_variation : Variant list of allele ID from CDS dict
    """

    cds_allele_seq_with_variation = insert_variations_into_sequence(sample_ref_seq, sample_cds_variation.pos_list, sample_cds_variation.ref_list, sample_cds_variation.alt_list)

    with open(os.path.join(args.schema_dir, sample_cds+'.fasta'), 'a') as file:

        file.write(f'>{sample_cds}_{cds_allele_id}\n')
        file.write(f'{cds_allele_seq_with_variation}\n')

        file.close()


def write_variations_to_reference_info_file(cds: str, allele_id: str, cds_variation: Info) -> None:
    """
    Write variations to the reference_info.txt file.

    Parameter
    ---------
    cds : Name of CDS
    allele_id : Allele ID of CDS
    cds_variation : Variant list of allele ID from CDS dict
    """
    
    with open(args.reference_info, 'a') as file:

        line = []
        for pos, ref, alt, qual in zip(cds_variation.pos_list,
                                        cds_variation.ref_list,
                                        cds_variation.alt_list,
                                        cds_variation.qual_list):

            if type(alt) is list:

                alt = ";".join(alt)

            line.append(f'{pos}*{ref}>{alt}-{qual}')

        file.write(f'{cds}_{allele_id}\t{",".join(line)}\n')

        file.close()


def write_variations_to_reference_vcf_file(cds: str, ref_allele_id: str, cds_variation: Info) -> None:
    """
    Take the variations from sample variation list
    Write from sample.vcf to reference.vcf

    Parameter
    ---------
    cds : CDS name
    ref_allele_id : Allele ID of reference sequence of CDS
    cds_variation : list of variations in allele of CDS
    """

    cds_vcf_line_list = []

    reference_vcf_file = open(args.reference_vcf, 'a')

    info_line_list = []

    with open(args.sample_vcf, 'r') as file:

        for line in file.readlines():

            if f'{cds}_1' in line and not line.startswith('#'):

                vcf_line = Vcf(line)

                if vcf_line.pos in cds_variation.pos_list:
                
                    vcf_line.chr = cds + '_' + ref_allele_id
                    vcf_line.info = "."
                    reference_vcf_file.write(str(vcf_line))
                    
        file.close()

    reference_vcf_file.close()


def get_allele_ids_of_cds_in_reference_info_txt() -> dict:
    """
    Reads reference_info.txt and returns allele ID dictionary
    for each CDS in reference_info.txt

    Returns
    -------
    novel_allele_id_of_cds_dict: {CDS: {allele_id: allele_info} ...}
    """

    novel_allele_id_of_cds_dict = {}

    for cds, alleles in read_reference_info_txt(args.reference_info).items():

        novel_allele_id_of_cds_dict[cds] = max(map(int, alleles.keys())) + 1

    return novel_allele_id_of_cds_dict


def base_dict_for_sam( base_dict: dict, pos_seq_cigar: list, variation_pos: int) -> dict:
    """
    Returns number of bases aligned to given position in sam file's line

    Parameters
    ----------
    base_dict : dictionary of bases to count
    pos_seq_cigar : position, sequence, and cigar sequence in sample sam file's line
    variation_pos: variation position on CDS sequence

    Returns
    -------
    base_dict : dictionary of bases to count 
    """

    import re

    pos, seq, cigar = pos_seq_cigar[0], pos_seq_cigar[1], pos_seq_cigar[2]

    match_pos_list = []
    current_position = 0
    pos_in_ref = variation_pos - 1
    sam_start_pos = pos - 1

    # match/mismatch loci in sam sequence
    for case, number_of_case in zip(re.findall("[A-Z]", cigar),
                                    map(int, re.findall("[0-9]+", cigar))):

        if case == 'M':

            match_start = current_position
            current_position += number_of_case
            match_end = current_position - 1
            match_pos_list.append([match_start, match_end])

        elif case == 'D':

            current_position -= number_of_case

        else:

            current_position += number_of_case

    for start, end in match_pos_list:
        
        pos_in_sam = pos_in_ref - sam_start_pos + start

        if start >= 0 and start <= pos_in_sam < end:

            base_dict[seq[pos_in_sam]] += 1

    return base_dict


def base_ratio_check(variation_pos: int, cds_name: str) -> str:
    """
    Get the dominant base for the locus on CDS
    
    Parameters
    ----------
    variation_pos: variation position on CDS sequence
    cds_name: CDS sequence name to be investigated

    Returns
    -------
    dominant_base : base name with max count
    """

    # avoid redundant operations after cds is found
    is_cds_found = False

    base_dict = { 'A': 0, 'C': 0, 'G': 0, 'T': 0 }

    for pos_seq_cigar in sample_sam_dict[cds_name]:

        base_dict_for_sam( base_dict, pos_seq_cigar, variation_pos )

    dominant_base = max( base_dict, key=base_dict.get )

    return dominant_base


def get_sample_sam_dict() -> dict:
    """
    Read Sample's SAM file and returns sequence, cigar, and position info

    Returns
    -------
    sample_sam_dict : { CDS: [pos, seq, cigar], [pos, seq, cigar], ..., CDS: ...}
    """

    sample_sam_dict = dict()

    with open(args.sample_sam, 'r') as file:

        for line in file.readlines():

            # skip header
            if not line.startswith('@'):

                sam = Sam(line)

                if sam.target_sequence_name not in sample_sam_dict.keys():

                    sample_sam_dict[sam.target_sequence_name] = []
                    sample_sam_dict[sam.target_sequence_name].append([ sam.position, sam.sequence, sam.cigar ])

                else:

                    sample_sam_dict[sam.target_sequence_name].append([ sam.position, sam.sequence, sam.cigar ])

        file.close()

    return sample_sam_dict


if __name__ == "__main__":

    parser = argparse.ArgumentParser(add_help=True)

    parser.add_argument('--reference_fasta',  type=str, required=True,  help='Reference\'s fasta file name with its directory')
    parser.add_argument('--reference_info',   type=str, required=True,  help='Reference\'s info file name with its directory')
    parser.add_argument('--reference_vcf',    type=str, required=True,  help='Reference\'s vcf file name with its directory')
    parser.add_argument('--sample_depth',     type=str, required=True,  help='Sample\'s depth file name with its directory')
    parser.add_argument('--sample_vcf',       type=str, required=True,  help='Sample\'s vcf file name with its directory')
    parser.add_argument('--update_reference', type=str, required=False, help='Update reference\'s vcf and info file for the further analysis.')
    parser.add_argument('--schema_dir',       type=str, required=True,  help='Directory of schema\'s to write novel alleles')
    parser.add_argument('--sample_sam',       type=str, required=True,  help='Sample\'s sam file name with its directory')
    
    args = parser.parse_args()

    novel_allele_id_of_cds_dict = get_allele_ids_of_cds_in_reference_info_txt()
    
    sample_sam_dict = get_sample_sam_dict()

    with open(f'{args.sample_vcf[:-4]}_mlst.tsv', 'w') as file:

        for sample_cds, allele_id in take_allele_id_for_sample_from_chewbbaca_alleles().items():

            file.write(f'{sample_cds.split("_")[0]}\t{allele_id}\n')

        file.close()