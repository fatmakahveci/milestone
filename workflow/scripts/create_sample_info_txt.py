#################################################################################
## author: @fatmakhv                                                           ##
## date: 14/10/2021                                                            ##
## aim: create info file for the sample and write novel alleles to schema seed ##
#################################################################################

import argparse, os


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


def remove_common_suffices(var1: str, var2: str) -> [str, str]:
    """
    Take two variations and remove the common suffices

    Parameters
    ----------
    var1 : variation sequence
    var2 : variation sequence

    Returns
    -------
    var1 : updated variation 1 of which common suffix is deleted
    var2 : updated variation 2 of which common suffix is deleted
    """

    check_len = min(len(var1), len(var2))

    var1 = var1[::-1]
    var2 = var2[::-1]

    i=0
    while i < check_len and var1[i] == var2[i]:

        i+=1

    return var1[i:][::-1], var2[i:][::-1]


def remove_common_prefices(pos: int, var1: str, var2: str) -> [str, str]:
    """
    Take two variations and remove the common prefices

    Parameters
    ----------
    pos : position of the variation which might be affected by the change
    var1 : variation sequence
    var2 : variation sequence

    Returns
    -------
    var1 : updated variation 1 of which common prefix is deleted
    var2 : updated variation 2 of which common prefix is deleted
    """

    check_len = min(len(var1), len(var2))

    i=0
    while i < check_len and var1[i] == var2[i]:

        i+=1

    return pos-i, var1[i:], var2[i:]


def remove_common_mid(pos: int, var1: str, var2: str, qual: int) -> [list, int]:
    """
    Remove the common substrings between ref and alt i.e. TAAG GAAC -> T G - G C

    Parameters
    ----------
    pos : position of the variation which might be affected by the change
    var1 : variation sequence
    var2 : variation sequence
    qual : quality score

    Returns
    -------
    pos : updated positions after the removal
    var1 : updated variation 1 of which common prefix is deleted
    var2 : updated variation 2 of which common prefix is deleted
    qual : quality score
    number_of_common_mids : number of splits for the further investigation
    """

    split_idx_list = []
    
    pos_list = []
    ref_list = []
    alt_list = []
    qual_list = []

    i = 0 
    while i < len(var1):

        if var1[i] != var2[i]:

            split_idx_list.append(i)

        i += 1

    for i in split_idx_list:

        pos_list.append(pos+i)
        ref_list.append(var1[i])
        alt_list.append(var2[i])
        qual_list.append(qual)

    return pos_list, ref_list, alt_list, qual_list, len(split_idx_list)


def remove_redundance(variations: Info) -> Info:
    """
    Remove common prefices and suffices from both reference and alternate

    Parameters
    ----------
    variations : Info

    Returns
    -------
    variations : updated variations with removed suffices and prefices
    """

    variations_from_splits = []

    number_of_splits = 0

    variations_to_replace = []

    for i in range(len(variations.pos_list)):

        # ACG TG -> AC T
        if len(variations.ref_list[i]) > 1 and len(variations.alt_list[i]) > 1:

            variations.ref_list[i], variations.alt_list[i] = remove_common_suffices(variations.ref_list[i], variations.alt_list[i])

        # reducing suffices might have effect on the length
        # ACG AT -> CG T
        if len(variations.ref_list[i]) > 1 and len(variations.alt_list[i]) > 1:

            variations.pos_list[i], variations.ref_list[i], variations.alt_list[i] = remove_common_prefices(variations.pos_list[i], variations.ref_list[i], variations.alt_list[i])

        # AGGT -> A  T
        if len(variations.ref_list[i]) == len(variations.alt_list[i]) and len(variations.ref_list[i]) > 2:

            variations_to_replace.append([i, remove_common_mid(variations.pos_list[i], variations.ref_list[i], variations.alt_list[i], variations.qual_list[i])])

    # if there is a split in variations of cds
    if len(variations_to_replace) > 0:

        for variation in variations_to_replace[::-1]:

            # remove the current redundant variation in that index
            del variations.pos_list[variation[0]]
            del variations.ref_list[variation[0]]
            del variations.alt_list[variation[0]]
            del variations.qual_list[variation[0]]

            idx_to_put = variation[0]
            var = variation[1]

            for i in range(len(var[0]))[::-1]:

                variations.pos_list.insert(idx_to_put, var[0][i])
                variations.ref_list.insert(idx_to_put, var[1][i])
                variations.alt_list.insert(idx_to_put, var[2][i])
                variations.qual_list.insert(idx_to_put, var[3][i])

    return variations


def merge_variations(variations: Info) -> Info:
    """
    Take the variations for reference_info.txt file
    Merge the variations of which positions are intersected.

    Parameters
    ----------
    variations : variations in reference_info.txt file

    Returns
    -------
    variations : merged variations for reference_info.txt file
    """

    # no need to merge single variation
    if len(variations.pos_list) == 1:

        return variations

    merged_variations_list = []

    pos_list, ref_list, alt_list, qual_list = variations.pos_list, variations.ref_list, variations.alt_list, variations.qual_list

    number_of_variations = len(pos_list)

    is_the_last_variation_added = False

    for i in range(number_of_variations):

        # if i < number_of_variations - 1: # 0 1 2

            # if pos_list[i] + len(ref_list[i]) == pos_list[i+1]:

            #     merged_variations_list.append(f"{pos_list[i]}*{ref_list[i]+ref_list[i+1]}>{alt_list[i]+alt_list[i+1]}-{min(qual_list[i],qual_list[i+1])}")

            #     if i == number_of_variations - 2: # added item in the last position to the last-1^th position

            #         is_the_last_variation_added = True

            # else:

    #             merged_variations_list.append(f"{pos_list[i]}*{ref_list[i]}>{alt_list[i]}-{qual_list[i]}")

    # if not is_the_last_variation_added:

        merged_variations_list.append(f"{pos_list[i]}*{ref_list[i]}>{alt_list[i]}-{qual_list[i]}")

    merged_list = []

    var_list = remove_redundance(Info(",".join(merged_variations_list)))

    for i in range(len(var_list.pos_list)):

        merged_list.append(f"{var_list.pos_list[i]}*{var_list.ref_list[i]}>{var_list.alt_list[i]}-{var_list.qual_list[i]}")

    return Info(",".join(merged_list))


def get_from_info(info: str, subinfo: str) -> str:
    """
    Return ...;<subinfo_name>="<subinfo>";... from INFO field in VCF file

    Parameters
    ----------
    info : INFO field in VCF file
    
    Returns
    -------
    <subinfo> : subinfo in INFO field
    """

    try:

        start = info.index(subinfo) + len(subinfo)
        end = info.index( ";", start )

        return info[start:end]

    except ValueError:
        return ""

def get_cigar(info: str) -> str:
    """
    Return ...;<CIGAR>="<cigar>";... from INFO field in VCF file

    Parameters
    ----------
    info : INFO field in VCF file
    
    Returns
    -------
    CIGAR : cigar in INFO field
    """

    # in case that multiple cigar take the first
    cigar = get_from_info(info, "CIGAR=").split(',')[0]
    
    for sep in [ 'M', 'D', 'I', 'S', 'H', '=', 'X']:
        cigar = cigar.replace(sep, '.')

    # return cigar, len(cigar)
    return get_from_info(info, "CIGAR="), sum(list(map(int, cigar.rstrip('.').split('.'))))


def get_var_type(info: str) -> str:
    """
    Return ...;<TYPE=>="<type>";... from TYPE field in VCF file

    Parameters
    ----------
    info : INFO field in VCF file
    
    Returns
    -------
    <type> : TYPE in INFO field
    """

    return get_from_info(info, "TYPE=")


def resolve_cigar(vcf_line: str, cigar: str) -> None:
    """
    Resolves the cigar and returns the corrected variations

    Parameters
    ----------
    vcf_line : vcf line containing complex cigar
    cigar : cigar sequence in vcf_line 
    """

    pos_list = []
    ref_list = []
    alt_list = []
    qual_list = []

    cigar_list = list(zip([int(s) for s in cigar if s.isdigit()], [s for s in cigar if not s.isdigit()]))

    i, j = 0, 0 # i := index in ref and  j := index in alt

    for k, [case_count, case] in enumerate(cigar_list):

        if case == 'D':

            pos_list.append(vcf_line.pos+i-1)
            ref_list.append(vcf_line.ref[i-1:i+case_count])
            alt_list.append(vcf_line.alt[j-1])
            qual_list.append(vcf_line.qual)

            i += case_count

        elif case == 'I':

            pos_list.append(vcf_line.pos+i-1)
            ref_list.append(vcf_line.ref[i-1])
            alt_list.append(vcf_line.alt[j-1:j+case_count])
            qual_list.append(vcf_line.qual)

            j += case_count

        elif case == 'X':

            # X is the last case
            if k == len(cigar_list) - 1:

                for n in range(case_count):

                    pos_list.append(vcf_line.pos+i+n)
                    ref_list.append(vcf_line.ref[i+n])
                    alt_list.append(vcf_line.alt[j+n])
                    qual_list.append(vcf_line.qual)

            # X is in-between
            else:

                # ...3X1M...
                if cigar_list[k+1][1] == 'M':

                    # ...nX or ...1X1M...
                    for n in range(case_count):

                        pos_list.append(vcf_line.pos+i+n)
                        ref_list.append(vcf_line.ref[i+n])
                        alt_list.append(vcf_line.alt[j+n])
                        qual_list.append(vcf_line.qual)

                # ...3X1I... or ...3X1D...
                elif cigar_list[k+1][1] == 'I' or cigar_list[k+1][1] == 'D':

                    # ...nX or ...1X1M...
                    for n in range(case_count):

                        # do not add the last item it will be used in the next D or I
                        pos_list.append(vcf_line.pos+i-1)
                        ref_list.append(vcf_line.ref[i+n-1])
                        alt_list.append(vcf_line.alt[j+n-1])
                        qual_list.append(vcf_line.qual)

            i += case_count
            j += case_count

        else:

            i += case_count
            j += case_count

    return pos_list, ref_list, alt_list, qual_list


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

            is_variation_to_be_added = False

            if not line.startswith("#"):

                vcf_line = Vcf(line)

                if get_var_type(vcf_line.info) == 'complex':

                    cigar, cigar_len = get_cigar(vcf_line.info)

                    if len(vcf_line.ref) != cigar_len:

                        is_variation_to_be_added = True

                        variation_pos_list, variation_ref_list, variation_alt_list, variation_qual_list = resolve_cigar(vcf_line, cigar)

                    # its type was expected as mnp but it is classified as complex 
                    else:

                        if len(vcf_line.ref) == len(vcf_line.alt) and len(vcf_line.ref) == 2:

                            variation_pos_list.append(vcf_line.pos)
                            variation_ref_list.append(vcf_line.ref[0])
                            variation_alt_list.append(vcf_line.alt[0])
                            variation_qual_list.append(vcf_line.qual)

                            variation_pos_list.append(vcf_line.pos+1)
                            variation_ref_list.append(vcf_line.ref[1])
                            variation_alt_list.append(vcf_line.alt[1])
                            variation_qual_list.append(vcf_line.qual)

                elif get_var_type(vcf_line.info) == 'mnp':

                    is_variation_to_be_added = True

                    cigar, cigar_len = get_cigar(vcf_line.info)

                    variation_pos_list, variation_ref_list, variation_alt_list, variation_qual_list = resolve_cigar(vcf_line, cigar)

                elif vcf_line.sample.startswith('0'):

                    continue

                elif vcf_line.sample.startswith('1'):

                    info_line = f'{vcf_line.pos}*{vcf_line.ref}>{vcf_line.alt}-{vcf_line.qual}'

                else:

                    info_line = f'{vcf_line.pos}*{vcf_line.ref}>{vcf_line.alt.split(",")[vcf_line.sample[0]]}-{vcf_line.qual}'

                if vcf_line.chr not in sample_variation_dict.keys():

                    # allele_field # s := sample
                    sample_variation_dict[vcf_line.chr] = Info(info_line)

                else:

                    if is_variation_to_be_added:

                        for i in range(len(variation_pos_list)):

                            sample_variation_dict[vcf_line.chr].pos_list.append(variation_pos_list[i])
                            sample_variation_dict[vcf_line.chr].ref_list.append(variation_ref_list[i])
                            sample_variation_dict[vcf_line.chr].alt_list.append(variation_alt_list[i])
                            sample_variation_dict[vcf_line.chr].qual_list.append(variation_qual_list[i])

                    else:

                        sample_variation_dict[vcf_line.chr].pos_list.append(vcf_line.pos)
                        sample_variation_dict[vcf_line.chr].ref_list.append(vcf_line.ref)
                        sample_variation_dict[vcf_line.chr].alt_list.append(vcf_line.alt)
                        sample_variation_dict[vcf_line.chr].qual_list.append(vcf_line.qual)

    for sample_cds, variations in sample_variation_dict.items():

        sample_variation_dict[sample_cds] = merge_variations(sample_variation_dict[sample_cds])

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

    if len(seq) <= len(ref_seq) * 0.8:

        return "ASM"

    if len(seq) >= len(ref_seq) * 1.2:

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

    import glob, shutil

    temp_sample_vcf_dir = f'{args.sample_vcf}_dir'

    if not os.path.isdir(temp_sample_vcf_dir):

        os.mkdir(temp_sample_vcf_dir)

    ref_allele_id = '1' # cds.split('_')[1]
    
    sample_variation_dict = create_sample_variation_dict(ref_allele_id)

    reference_allele_variation_dict = read_reference_info_txt(info_file=args.reference_info)

    # { CDS1_ref: seq1, ... }
    cds_seq_dict = get_reference_cds_seq_dict()

    sample_allele_dict = {}

    for cds, coverage in get_cds_coverage_info().items():

        sample_cds = cds.split('_')[0]

        if coverage.coverage <= 60:

            # CDS is not covered by the reads.
            sample_allele_dict[cds] = 'LNF'
            is_novel = False

        elif 60 < coverage.coverage <= 80:

            # CDS is not covered by the reads.
            sample_allele_dict[cds] = 'ASM'
            is_novel = False

        # Sample reads are mapped to the reference.
        else:

            if sample_cds not in sample_variation_dict.keys():

                # Sample does not have any variation in CDS so it equals to the reference.
                sample_allele_dict[cds] = ref_allele_id
                is_novel = False

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

                    #         write_allele_sequence_to_schema_seed(sample_cds, ref_allele_id, cds_seq_dict[f'{sample_cds}_{ref_allele_id}'], sample_variation_dict[sample_cds])
                            write_variations_to_reference_vcf_file(sample_cds, temp_sample_vcf_dir, ref_allele_id, sample_variation_dict[sample_cds])
                            write_variations_to_reference_info_file(sample_cds, sample_allele_dict[cds], sample_variation_dict[sample_cds])

                else: # both sample and reference has the variations of this allele

                    sample_allele_dict[cds], is_novel = compare_ref_to_sample_variations(cds = sample_cds, ref_allele_id = ref_allele_id, cds_seq_dict = cds_seq_dict, reference_info = reference_allele_variation_dict[sample_cds], sample_cds_info = sample_variation_dict[sample_cds])

                    if is_novel == True and args.update_reference == 'True':

                        if len(sample_variation_dict[sample_cds].pos_list) == 0:

                            sample_allele_dict[cds] = ref_allele_id # sample doesn't have variations for the CDS so it equals to the reference CDS

                        else: # sample has variations for the CDS

                            sample_allele_dict[cds] = novel_allele_id_of_cds_dict[sample_cds]

                            # write_allele_sequence_to_schema_seed(sample_cds, sample_allele_dict[cds], cds_seq_dict[f'{sample_cds}_{ref_allele_id}'], sample_variation_dict[sample_cds])
                            write_variations_to_reference_vcf_file(sample_cds, temp_sample_vcf_dir, ref_allele_id, sample_variation_dict[sample_cds]) 
                            write_variations_to_reference_info_file(sample_cds, sample_allele_dict[cds], sample_variation_dict[sample_cds])

    # if args.update_reference == 'True':

    #     os.system(f"bcftools concat {' '.join(glob.glob(f'{temp_sample_vcf_dir}/*.vcf.gz'))} --threads {args.threads} -Oz -o {args.sample_vcf}.gz")
    #     os.system(f"tabix -f -p vcf {args.sample_vcf}.gz")
    #     os.system(f"bcftools concat -a --threads {args.threads} {args.reference_vcf}.gz {args.sample_vcf}.gz -Ov -o {args.reference_vcf}")
    #     os.system(f"bcftools sort {args.reference_vcf} -Oz -o {args.reference_vcf}.gz")
    #     os.system(f"bcftools norm {args.reference_vcf}.gz -m +any -Ov -o {args.reference_vcf}")
    #     os.system(f"bgzip -f {args.reference_vcf} && tabix -f -p vcf {args.reference_vcf}.gz")

    try:

        shutil.rmtree(temp_sample_vcf_dir)

    except OSError as e:

        print("Error: %s - %s." % (e.filename, e.strerror))

    return sample_allele_dict


def write_allele_sequence_to_schema_seed(sample_cds: str, cds_allele_id: str, sample_ref_seq: str, sample_cds_variation: Info) -> None:
    """
    Write variations to the schema file.

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


def write_variations_to_reference_vcf_file(cds: str, temp_sample_vcf_dir: str, ref_allele_id: str, cds_variation: Info) -> None:
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

    temp_sample_vcf_file_name = os.path.join(temp_sample_vcf_dir, f'{cds}.vcf')

    temp_sample_vcf_file = open(temp_sample_vcf_file_name , 'w')

    info_line_list = []

    with open(args.sample_vcf, 'r') as file:

        for line in file.readlines():

            if line.startswith('#'):

                if line.startswith('#CHROM'):

                    temp_sample_vcf_file.write("\t".join(line.rstrip('\n').split('\t')[:-1])+"\tREFERENCE\n")

                else:

                    temp_sample_vcf_file.write(line)

            else:

                if f'{cds}_1' in line:

                    vcf_line = Vcf(line)

                    if vcf_line.pos in cds_variation.pos_list:
                
                        vcf_line.chr = cds + '_' + ref_allele_id
                        vcf_line.info = "."
                        vcf_line.format = "GT"
                        vcf_line.sample = '1'

                        temp_sample_vcf_file.write(str(vcf_line))
                    
        file.close()

    temp_sample_vcf_file.close()

    os.system(f'bgzip -f {temp_sample_vcf_file_name} && tabix -p vcf {temp_sample_vcf_file_name}.gz')


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
    parser.add_argument('--threads',          type=str, required=True,  help='Number of threads')

    args = parser.parse_args()

    novel_allele_id_of_cds_dict = get_allele_ids_of_cds_in_reference_info_txt()
    
    with open(f'{args.sample_vcf[:-4]}_mlst.tsv', 'w') as file:

        for sample_cds, allele_id in take_allele_id_for_sample_from_chewbbaca_alleles().items():

            file.write(f'{sample_cds.split("_")[0]}\t{allele_id}\n')

        file.close()
