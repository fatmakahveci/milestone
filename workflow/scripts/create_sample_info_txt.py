#!/usr/bin/env python3

#################################################################################
## author: @fatmakhv                                                           ##
## date: 01/11/2021                                                            ##
## aim: create info file for the sample and write novel alleles to schema seed ##
#################################################################################


import os, subprocess, sys
from collections import Counter
from pathlib import Path
from types import SimpleNamespace

SCRIPT_DIR = Path(__file__).resolve().parent
if str(SCRIPT_DIR) not in sys.path:
    sys.path.insert(0, str(SCRIPT_DIR))

from wgmlst_utils import (
    REFERENCE_ALLELE_ID,
    describe_quality_check,
    get_allele_id_from_allele_name,
    get_cds_name_from_allele_name,
    get_locus_length_mode_from_fasta,
    get_next_novel_allele_id_from_fasta,
    get_reference_alignment_target_name as shared_get_reference_alignment_target_name,
    get_reference_allele_name as shared_get_reference_allele_name,
    quality_check as shared_quality_check,
)
from script_utils import ParsedVariationInfo, run_checked_command
from sample_variation_io import (
    create_sample_variation_dict as shared_create_sample_variation_dict,
    get_cigar as shared_get_cigar,
    get_cigar_info as shared_get_cigar_info,
    get_sample_sam_dict as shared_get_sample_sam_dict,
    get_var_type as shared_get_var_type,
    resolve_cigar as shared_resolve_cigar,
)
from reference_update_io import (
    write_allele_sequence_to_schema_seed as shared_write_allele_sequence_to_schema_seed,
    write_variations_to_reference_info_file as shared_write_variations_to_reference_info_file,
    write_variations_to_reference_vcf_file as shared_write_variations_to_reference_vcf_file,
)
from reference_merge_pipeline import merge_reference_with_novel_vcfs
from sample_call_cli_args import parse_cli_args
from sample_call_orchestrator import run_sample_calling_workflow
from sample_formats import Coverage, Sam, Vcf
from variation_utils import (
    build_variation_signature as shared_build_variation_signature,
    merge_variations as shared_merge_variations,
    remove_common_mid as shared_remove_common_mid,
    remove_common_prefixes as shared_remove_common_prefixes,
    remove_common_suffixes as shared_remove_common_suffixes,
    remove_redundance as shared_remove_redundance,
    write_novel_allele_exports as shared_write_novel_allele_exports,
)


Info = ParsedVariationInfo


DEFAULT_MIN_LOCUS_COVERAGE = 60.0
args = SimpleNamespace(
    schema_dir="",
    reference_info="",
    sample_sam="",
    sample_vcf="",
    sample_depth="",
    reference_fasta="",
    reference_vcf="",
    update_reference="False",
    min_locus_coverage=DEFAULT_MIN_LOCUS_COVERAGE,
    translation_table=11,
    allowed_start_codons=None,
    threads="1",
)
sample_sam_dict: dict = {}
novel_allele_id_of_cds_dict: dict = {}


def _runtime_args():
    return globals().get("args", args)


def _init_runtime_state() -> None:
    global novel_allele_id_of_cds_dict, sample_sam_dict
    novel_allele_id_of_cds_dict = get_allele_ids_of_cds_in_reference_info_txt()
    sample_sam_dict = get_sample_sam_dict()


def parse_allowed_start_codons(raw_value: str | None) -> set[str] | None:
    if not raw_value:
        return None
    codons = {item.strip().upper() for item in raw_value.split(",") if item.strip()}
    return codons or None


def run_command(command: list[str]) -> None:
    run_checked_command(command, stderr=subprocess.DEVNULL)


def base_ratio_check( variation_pos: int, cds_name: str ) -> str:
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

    base_dict = { 'A': 0, 'C': 0, 'G': 0, 'T': 0 }

    for pos_seq_cigar in sample_sam_dict[cds_name]:

        base_dict_for_sam( base_dict=base_dict, pos_seq_cigar=pos_seq_cigar, variation_pos=variation_pos )

    dominant_base = max( base_dict, key = base_dict.get )

    return dominant_base


def base_dict_for_sam( base_dict: dict, pos_seq_cigar: list, variation_pos: int ) -> dict:
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
    for case, number_of_case in zip(
        re.findall("[A-Z]", cigar),
        map(int, re.findall("[0-9]+", cigar)),
        strict=False,
    ):

        if case == 'M':

            match_start = current_position
            current_position += number_of_case
            match_end = current_position - 1
            match_pos_list.append( [ match_start, match_end ] )

        elif case == 'D':

            current_position -= number_of_case

        else:

            current_position += number_of_case

    for start, end in match_pos_list:
        
        pos_in_sam = pos_in_ref - sam_start_pos + start

        if start >= 0 and start <= pos_in_sam < end:

            base_dict[ seq[pos_in_sam] ] += 1

    return base_dict


def get_allele_ids_of_cds_in_reference_info_txt() -> dict:
    """
    Return the next numeric novel allele ID for loci present in the schema.

    The previous implementation derived IDs only from `reference_info.txt`,
    which silently failed for loci whose reference allele had no listed
    variation or used a non-numeric identifier.
    Return
    ------
    novel_allele_id_of_cds_dict: {CDS: {allele_id: allele_info} ...}
    """

    novel_allele_id_of_cds_dict = {}

    runtime = _runtime_args()
    for schema_fasta in Path(runtime.schema_dir).glob("*.fasta"):
        cds = schema_fasta.stem
        novel_allele_id_of_cds_dict[cds] = get_next_novel_allele_id(cds)

    for cds in read_reference_info_txt(info_file=runtime.reference_info).keys():
        novel_allele_id_of_cds_dict.setdefault(cds, get_next_novel_allele_id(cds))

    return novel_allele_id_of_cds_dict


def get_numeric_allele_ids_for_locus(cds: str) -> list[int]:
    runtime = _runtime_args()
    locus_fasta = Path(runtime.schema_dir) / f"{cds}.fasta"
    from wgmlst_utils import get_numeric_allele_ids_from_fasta
    return get_numeric_allele_ids_from_fasta(locus_fasta)


def get_next_novel_allele_id(cds: str) -> int:
    runtime = _runtime_args()
    return get_next_novel_allele_id_from_fasta(Path(runtime.schema_dir) / f"{cds}.fasta")


def get_reference_alignment_target_name(cds: str) -> str:
    return shared_get_reference_alignment_target_name(cds, sample_sam_dict)


def remove_common_suffixes( var1: str, var2: str ) -> [ str, str ]:
    """
    Take two variations and remove the common suffixes
    Parameters
    ----------
    var1 : variation sequence
    var2 : variation sequence
    Returns
    -------
    var1 : updated variation 1 with common suffix removed
    var2 : updated variation 2 with common suffix removed
    """

    return list(shared_remove_common_suffixes(var1, var2))


def remove_common_prefixes( pos: int, var1: str, var2: str ) -> [ int, str, str ]:
    """
    Take two variations and remove the common prefixes
    Parameter
    ---------
    pos : position of the variation which might be affected by the change
    var1 : variation sequence
    var2 : variation sequence
    Return
    ------
    pos : position of the variation which might be affected by the change
    var1 : updated variation 1 with common prefix removed
    var2 : updated variation 2 with common prefix removed
    """

    return list(shared_remove_common_prefixes(pos, var1, var2))


def remove_common_mid( pos: int, var1: str, var2: str, qual: int ) -> [ list, list, list, list, int ]:
    """
    Remove the common substrings between ref and alt i.e. TAAG GAAC -> T G - G C
    Parameter
    ---------
    pos : position of the variation which might be affected by the change
    var1 : variation sequence
    var2 : variation sequence
    qual : quality score
    Return
    ------
    pos : updated positions after the removal
    var1 : updated variation 1 of which common prefix is deleted
    var2 : updated variation 2 of which common prefix is deleted
    qual : quality score
    number_of_common_mids : number of splits for the further investigation
    """

    return list(shared_remove_common_mid(pos, var1, var2, qual))


def remove_redundance(variations: Info) -> Info:
    """
    Remove common prefixes and suffixes from both reference and alternate
    Parameters
    ----------
    variations : Info
    Returns
    -------
    variations : updated variations with removed suffixes and prefixes
    """

    return shared_remove_redundance(variations)


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

    return shared_merge_variations(variations)


def get_var_type(info: str) -> str:
    """
    Return ...;TYPE="<type>";... from INFO field in VCF file
    Parameters
    ----------
    info : INFO field in VCF file
    
    Returns
    -------
    <type> : TYPE in INFO field
    """

    return shared_get_var_type(info)


def get_cigar(info: str) -> str:
    """
    Return ...;CIGAR="<cigar>";... from INFO field in VCF file
    Parameter
    ---------
    info : INFO field in VCF file
    
    Return
    ------
    <cigar> : CIGAR in INFO field
    """

    return shared_get_cigar(info)


def get_cigar_info(info: str) -> [ str, int ]:
    """
    Return ...;<CIGAR>="<cigar>";... from INFO field in VCF file
    Parameter
    ---------
    info : INFO field in VCF file
    
    Return
    ------
    CIGAR : cigar in INFO field
    cigar_len : length of CIGAR
    """

    return shared_get_cigar_info(info)


def resolve_cigar( vcf_line: str, cigar: str ) -> [ list, list, list, list ]:
    """
    Resolve the cigar and returns the corrected variations
    Parameter
    ---------
    vcf_line : vcf line containing complex cigar
    cigar : cigar sequence in vcf_line
    Return
    ------
    pos_list : positions of variations in VCF file
    ref_list : references of variations in VCF file
    alt_list : alternates of variations in VCF file
    qual_list : qualities of variations in VCF file
    """

    return shared_resolve_cigar(
        vcf_line,
        cigar,
        base_ratio_lookup=base_ratio_check,
        reference_target_lookup=get_reference_alignment_target_name,
    )


def get_sample_sam_dict() -> dict:
    """
    Read Sample's SAM file and returns sequence, cigar, and position info
    Returns
    -------
    sample_sam_dict : { CDS: [pos, seq, cigar], [pos, seq, cigar], ..., CDS: ...}
    """

    runtime = _runtime_args()
    return shared_get_sample_sam_dict(runtime.sample_sam, Sam)


def create_sample_variation_dict() -> dict:
    """
    Creates variation dictionary
    Return
    ------
    sample_variation_dict : variations of sample dictionary for alleles
    """

    return shared_create_sample_variation_dict(
        sample_vcf_path=_runtime_args().sample_vcf,
        vcf_class=Vcf,
        variation_factory=Info,
        merge_variations=merge_variations,
        base_ratio_lookup=base_ratio_check,
        reference_target_lookup=get_reference_alignment_target_name,
    )


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

                cds = get_cds_name_from_allele_name( allele_name=fields[0] )
                allele_id = get_allele_id_from_allele_name( allele_name=fields[0] )
            
                if cds not in reference_allele_variation_dict.keys():

                    reference_allele_variation_dict[cds] = {}
            
                reference_allele_variation_dict[cds][allele_id] = Info(fields[1])

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

    with open(_runtime_args().sample_depth, 'r') as file:

        for line in file.readlines():

            if not line.startswith('#'):

                coverage = Coverage(line)
                cds_depth_dict[coverage.cds] = coverage

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

    for sequence in list( SeqIO.parse( _runtime_args().reference_fasta, "fasta" ) ):

        cds_seq_dict[sequence.id] = str(sequence.seq)

    return cds_seq_dict


def get_reference_allele_name(cds: str, cds_seq_dict: dict) -> str:
    """
    Return the schema's designated reference allele name for a locus.

    Prefer allele ID 1 when present, otherwise infer the unique reference
    sequence already embedded in the reference FASTA for that locus.
    """
    try:
        return shared_get_reference_allele_name(cds, cds_seq_dict.keys())
    except ValueError as exc:
        raise ValueError(
            f"Unable to determine a unique reference allele for locus '{cds}'. "
            f"Candidates found in reference FASTA: "
            f"{', '.join(sorted(name for name in cds_seq_dict if get_cds_name_from_allele_name(name) == cds)) or 'none'}."
        ) from exc


def get_reference_allele_id(cds: str, cds_seq_dict: dict) -> str:
    return get_allele_id_from_allele_name(get_reference_allele_name(cds, cds_seq_dict))


def get_locus_length_mode(cds: str, ref_seq: str) -> int:
    """
    Return the modal allele length for a locus based on the schema FASTA.
    Fallback to the reference sequence length if the locus FASTA is missing.
    """

    runtime = _runtime_args()
    return get_locus_length_mode_from_fasta(Path(runtime.schema_dir) / f"{cds}.fasta", ref_seq)


def insert_variations_into_sequence( cds_reference: str, pos_list: list, ref_list: list, alt_list: list ) -> str:
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
    for pos, ref, alt in zip(pos_list[::-1], ref_list[::-1], alt_list[::-1], strict=False):

        pos -= 1

        if type(alt) == list:

            alt = alt[0]

        cds_reference = cds_reference[ : pos ] + alt.replace( '.', '' ) + cds_reference[ pos+len(ref) : ]

    return cds_reference


def quality_check(
    seq: str,
    ref_seq: str,
    locus_mode_length: int,
    translation_table: int | None = None,
    allowed_start_codons: set[str] | None = None,
) -> str:
    """
    Validate whether a CDS-like sequence keeps an intact ORF and whether
    its length is within the tolerated range for the locus.
    Parameter
    ---------
    seq : FASTA sequence
    ref_seq : FASTA sequence for the reference of seq
    locus_mode_length : modal allele length for the locus
    Return
    ------
    quality label : EQ, Q, ASM, ALM, or LNF
    """

    runtime_args = _runtime_args()
    return shared_quality_check(
        seq,
        ref_seq,
        locus_mode_length,
        translation_table=translation_table or runtime_args.translation_table,
        allowed_start_codons=allowed_start_codons or parse_allowed_start_codons(runtime_args.allowed_start_codons),
    )


def quality_check_details(
    seq: str,
    ref_seq: str,
    locus_mode_length: int,
    translation_table: int | None = None,
    allowed_start_codons: set[str] | None = None,
) -> dict[str, object]:
    runtime_args = _runtime_args()
    runtime_translation_table = getattr(runtime_args, "translation_table", 11)
    runtime_start_codons = getattr(runtime_args, "allowed_start_codons", None)
    return describe_quality_check(
        seq,
        ref_seq,
        locus_mode_length,
        translation_table=translation_table or runtime_translation_table,
        allowed_start_codons=allowed_start_codons or parse_allowed_start_codons(runtime_start_codons),
    )


def build_novel_allele_sequence(sample_cds: str, cds_seq_dict: dict, sample_cds_variation: Info) -> str:
    reference_allele_name = get_reference_allele_name(sample_cds, cds_seq_dict)
    cds_reference = cds_seq_dict[reference_allele_name]
    return insert_variations_into_sequence(
        cds_reference=cds_reference,
        pos_list=sample_cds_variation.pos_list,
        ref_list=sample_cds_variation.ref_list,
        alt_list=sample_cds_variation.alt_list,
    )


def summarize_variations(sample_cds_variation: Info) -> str:
    return ",".join(
        f"{pos}:{ref}>{alt}"
        for pos, ref, alt in zip(
            sample_cds_variation.pos_list,
            sample_cds_variation.ref_list,
            sample_cds_variation.alt_list,
            strict=False,
        )
    )


def build_variation_signature(variation_info: Info) -> Counter[tuple[int, str, str]]:
    return shared_build_variation_signature(variation_info)


def write_novel_allele_exports(novel_alleles: list[dict[str, str]]) -> tuple[Path, Path]:
    return shared_write_novel_allele_exports(_runtime_args().sample_vcf, novel_alleles)


def compare_ref_to_sample_variations( cds: str, cds_seq_dict: dict, reference_info : Info, sample_cds_info : Info ) -> int:
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

    # cds length check
    diff_len = 0

    for ref, alt in zip(sample_cds_info.ref_list, sample_cds_info.alt_list, strict=False):

        if alt == '.':

            diff_len += len(ref)

        else:

            diff_len += len(ref) - len(alt)

    if diff_len % 3 != 0:

        is_novel = False
        allele_id = 'LNF' # incorrect length

    else:

        is_novel = True

        reference_allele_name = get_reference_allele_name(cds, cds_seq_dict)
        reference_allele_id = get_allele_id_from_allele_name(reference_allele_name)
        cds_reference = cds_seq_dict[reference_allele_name]
        locus_mode_length = get_locus_length_mode(cds, cds_reference)

        allele_id = quality_check(
            seq=insert_variations_into_sequence(
                cds_reference=cds_reference,
                pos_list=sample_cds_info.pos_list,
                ref_list=sample_cds_info.ref_list,
                alt_list=sample_cds_info.alt_list,
            ),
            ref_seq=cds_reference,
            locus_mode_length=locus_mode_length,
        )
        
        if allele_id == "EQ":

            is_novel = False
            allele_id = reference_allele_id

        elif allele_id == 'Q':

            sample_signature = build_variation_signature(sample_cds_info)
            for cds_id, ref_info in reference_info.items():

                if sample_signature == build_variation_signature(ref_info):
                    is_novel = False
                    allele_id = cds_id
                    break
            else:
                is_novel = True
                allele_id = str(novel_allele_id_of_cds_dict[cds])

        else:

            is_novel = False

    return [ allele_id, is_novel ]


def _variation_length_delta(variation_info: Info) -> int:
    diff_len = 0
    for ref, alt in zip(variation_info.ref_list, variation_info.alt_list, strict=False):
        if alt == ".":
            diff_len += len(ref)
        else:
            diff_len += len(ref) - len(alt)
    return diff_len


def _resolve_novel_reference_locus(
    sample_cds: str,
    cds_seq_dict: dict,
    sample_variation: Info,
) -> tuple[str, bool, dict[str, str] | None]:
    if len(sample_variation.pos_list) == 0:
        return get_reference_allele_id(sample_cds, cds_seq_dict), False, None

    if _variation_length_delta(sample_variation) % 3 != 0:
        return "LNF", False, None

    reference_allele_name = get_reference_allele_name(sample_cds, cds_seq_dict)
    reference_allele_id = get_allele_id_from_allele_name(reference_allele_name)
    cds_reference = cds_seq_dict[reference_allele_name]
    locus_mode_length = get_locus_length_mode(sample_cds, cds_reference)

    candidate_sequence = insert_variations_into_sequence(
        cds_reference=cds_reference,
        pos_list=sample_variation.pos_list,
        ref_list=sample_variation.ref_list,
        alt_list=sample_variation.alt_list,
    )
    qc_result = quality_check_details(
        seq=candidate_sequence,
        ref_seq=cds_reference,
        locus_mode_length=locus_mode_length,
    )
    status = str(qc_result["status"])

    if status == "EQ":
        return reference_allele_id, False, None
    if status != "Q":
        return status, False, None

    allele_id = str(get_next_novel_allele_id(sample_cds))
    record = {
        "locus": sample_cds,
        "allele_id": allele_id,
        "sequence": candidate_sequence,
        "status": "novel_candidate",
        "reference_allele": reference_allele_name,
        "variation_summary": summarize_variations(sample_variation),
        "orf_reason": str(qc_result["reason"]),
    }
    return allele_id, True, record


def _resolve_existing_reference_locus(
    sample_cds: str,
    cds_seq_dict: dict,
    sample_variation: Info,
    reference_allele_variation_dict: dict,
) -> tuple[str, bool, dict[str, str] | None]:
    allele_id, is_novel = compare_ref_to_sample_variations(
        cds=sample_cds,
        cds_seq_dict=cds_seq_dict,
        reference_info=reference_allele_variation_dict[sample_cds],
        sample_cds_info=sample_variation,
    )
    allele_id = str(allele_id)

    if is_novel and len(sample_variation.pos_list) != 0:
        record = {
            "locus": sample_cds,
            "allele_id": allele_id,
            "sequence": build_novel_allele_sequence(sample_cds, cds_seq_dict, sample_variation),
            "status": "novel_candidate",
            "reference_allele": get_reference_allele_name(sample_cds, cds_seq_dict),
            "variation_summary": summarize_variations(sample_variation),
            "orf_reason": "passes_orf_qc",
        }
        return allele_id, is_novel, record

    return allele_id, is_novel, None


def _update_reference_assets_for_novel(
    sample_cds: str,
    sample_variation: Info,
    cds_seq_dict: dict,
    temp_sample_vcf_dir: str,
    allele_id: str,
) -> None:
    if len(sample_variation.pos_list) == 0:
        return

    write_variations_to_reference_vcf_file(
        cds=sample_cds,
        reference_allele_name=get_reference_allele_name(sample_cds, cds_seq_dict),
        temp_sample_vcf_dir=temp_sample_vcf_dir,
        cds_variation=sample_variation,
    )
    write_variations_to_reference_info_file(
        cds=sample_cds,
        allele_id=allele_id,
        cds_variation=sample_variation,
    )


def take_allele_id_for_sample_from_chewbbaca_alleles() -> dict:
    """
    Returns allele ID for samples coding sequences
    Return
    ------
    sample_allele_dict : { cds_1 : allele_ID_1, ... }
    """

    import glob, shutil

    runtime = _runtime_args()
    temp_sample_vcf_dir = f'{runtime.sample_vcf}_dir'

    if not os.path.isdir(temp_sample_vcf_dir):

        os.mkdir(temp_sample_vcf_dir)
    
    sample_variation_dict = create_sample_variation_dict()

    reference_allele_variation_dict = read_reference_info_txt(info_file=runtime.reference_info)

    # { CDS1_ref: seq1, ... }
    cds_seq_dict = get_reference_cds_seq_dict()

    sample_allele_dict = {}
    novel_alleles: list[dict[str, str]] = []

    for cds, coverage in get_cds_coverage_info().items():

        sample_cds = get_cds_name_from_allele_name(cds)

        if coverage.coverage <= runtime.min_locus_coverage:

            # CDS is not covered by the reads.
            sample_allele_dict[cds] = 'LNF'
            is_novel = False

        else:

            if sample_cds not in sample_variation_dict.keys():

                # Sample does not have any variation in CDS so it equals to the reference.
                sample_allele_dict[cds] = get_reference_allele_id(sample_cds, cds_seq_dict)
                is_novel = False

            else:

                is_novel = False

                # It is a novel allele for the reference.
                sample_variation = sample_variation_dict[sample_cds]
                if sample_cds not in reference_allele_variation_dict.keys():
                    allele_id, is_novel, novel_record = _resolve_novel_reference_locus(
                        sample_cds=sample_cds,
                        cds_seq_dict=cds_seq_dict,
                        sample_variation=sample_variation,
                    )
                    sample_allele_dict[cds] = allele_id
                    if novel_record is not None:
                        novel_alleles.append(novel_record)
                else:  # both sample and reference has the variations of this allele
                    allele_id, is_novel, novel_record = _resolve_existing_reference_locus(
                        sample_cds=sample_cds,
                        cds_seq_dict=cds_seq_dict,
                        sample_variation=sample_variation,
                        reference_allele_variation_dict=reference_allele_variation_dict,
                    )
                    sample_allele_dict[cds] = allele_id
                    if novel_record is not None:
                        novel_alleles.append(novel_record)

                if is_novel and runtime.update_reference == 'True':
                    if sample_cds in reference_allele_variation_dict.keys():
                        if len(sample_variation.pos_list) == 0:
                            sample_allele_dict[cds] = get_reference_allele_id(sample_cds, cds_seq_dict)
                        else:
                            sample_allele_dict[cds] = str(novel_allele_id_of_cds_dict[sample_cds])

                    _update_reference_assets_for_novel(
                        sample_cds=sample_cds,
                        sample_variation=sample_variation,
                        cds_seq_dict=cds_seq_dict,
                        temp_sample_vcf_dir=temp_sample_vcf_dir,
                        allele_id=sample_allele_dict[cds],
                    )

    if runtime.update_reference == 'True':
        temp_vcfs = glob.glob(f'{temp_sample_vcf_dir}/*.vcf.gz')
        merge_reference_with_novel_vcfs(
            temp_vcfs=temp_vcfs,
            sample_vcf_path=runtime.sample_vcf,
            reference_vcf_path=runtime.reference_vcf,
            threads=runtime.threads,
            run_command=run_command,
        )

    if novel_alleles:
        write_novel_allele_exports(novel_alleles)

    try:

        shutil.rmtree(temp_sample_vcf_dir)

    except OSError as e:

        print("Error: %s - %s." % (e.filename, e.strerror))

    return sample_allele_dict


def write_allele_sequence_to_schema_seed( sample_cds: str, cds_allele_id: str, sample_ref_seq: str, sample_cds_variation: Info ) -> None:
    """
    Write variations to the schema file.
    Parameter
    ---------
    sample_cds : Name of CDS
    cds_allele_id : Allele ID of CDS
    sample_ref_seq : Sequence of CDS reference
    sample_cds_variation : Variant list of allele ID from CDS dict
    """

    cds_allele_seq_with_variation = insert_variations_into_sequence( cds_reference=sample_ref_seq, pos_list=sample_cds_variation.pos_list, ref_list=sample_cds_variation.ref_list, alt_list=sample_cds_variation.alt_list )

    runtime = _runtime_args()
    shared_write_allele_sequence_to_schema_seed(
        schema_dir=runtime.schema_dir,
        sample_cds=sample_cds,
        cds_allele_id=cds_allele_id,
        sequence_with_variation=cds_allele_seq_with_variation,
    )


def write_variations_to_reference_info_file( cds: str, allele_id: str, cds_variation: Info ) -> None:
    """
    Write variations to the reference_info.txt file.
    Parameter
    ---------
    cds : Name of CDS
    allele_id : Allele ID of CDS
    cds_variation : Variant list of allele ID from CDS dict
    """
    
    shared_write_variations_to_reference_info_file(
        reference_info_path=_runtime_args().reference_info,
        cds=cds,
        allele_id=allele_id,
        cds_variation=cds_variation,
    )


def write_variations_to_reference_vcf_file(
    cds: str,
    reference_allele_name: str,
    temp_sample_vcf_dir: str,
    cds_variation: Info,
) -> None:
    """
    Take the variations from sample variation list
    Write from sample.vcf to reference.vcf
    Parameter
    ---------
    cds : CDS name
    temp_sample_vcf_dir : Temporary directory to put sample vcf files
    cds_variation : list of variations in allele of CDS
    """

    shared_write_variations_to_reference_vcf_file(
        sample_vcf_path=_runtime_args().sample_vcf,
        cds=cds,
        reference_allele_name=reference_allele_name,
        temp_sample_vcf_dir=temp_sample_vcf_dir,
        cds_variation=cds_variation,
        vcf_class=Vcf,
        run_command=run_command,
    )


if __name__ == "__main__":
    parsed_args = parse_cli_args(DEFAULT_MIN_LOCUS_COVERAGE)
    args = parsed_args

    run_sample_calling_workflow(
        sample_vcf_path=parsed_args.sample_vcf,
        init_state=_init_runtime_state,
        collect_allele_calls=take_allele_id_for_sample_from_chewbbaca_alleles,
        get_cds_name_from_allele_name=get_cds_name_from_allele_name,
    )
