from __future__ import annotations

from collections import Counter
from pathlib import Path
import re
from typing import Iterable


COMMON_BACTERIAL_START_CODONS = {"ATG", "CTG", "GTG", "TTG"}
STOP_CODONS = {"TAG", "TAA", "TGA"}
REFERENCE_ALLELE_ID = "1"
LENGTH_MODE_TOLERANCE = 0.20
TRANSLATION_TABLE_DEFAULT = 11
BASE1 = "TTTTTTTTTTTTTTTTCCCCCCCCCCCCCCCCAAAAAAAAAAAAAAAAGGGGGGGGGGGGGGGG"
BASE2 = "TTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGG"
BASE3 = "TCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAG"
TRANSLATION_TABLE_ROWS = {
    1: {
        "aas": "FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG",
        "starts": "---M------**--*----M---------------M----------------------------",
    },
    2: {
        "aas": "FFLLSSSSYY**CCWWLLLLPPPPHHQQRRRRIIMMTTTTNNKKSS**VVVVAAAADDEEGGGG",
        "starts": "----------**--------------------MMMM----------**---M------------",
    },
    3: {
        "aas": "FFLLSSSSYY**CCWWTTTTPPPPHHQQRRRRIIMMTTTTNNKKSSRRVVVVAAAADDEEGGGG",
        "starts": "----------**----------------------MM---------------M------------",
    },
    4: {
        "aas": "FFLLSSSSYY**CCWWLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG",
        "starts": "--MM------**-------M------------MMMM---------------M------------",
    },
    5: {
        "aas": "FFLLSSSSYY**CCWWLLLLPPPPHHQQRRRRIIMMTTTTNNKKSSSSVVVVAAAADDEEGGGG",
        "starts": "---M------**--------------------MMMM---------------M------------",
    },
    6: {
        "aas": "FFLLSSSSYYQQCC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG",
        "starts": "--------------*--------------------M----------------------------",
    },
    9: {
        "aas": "FFLLSSSSYY**CCWWLLLLPPPPHHQQRRRRIIIMTTTTNNNKSSSSVVVVAAAADDEEGGGG",
        "starts": "----------**-----------------------M---------------M------------",
    },
    10: {
        "aas": "FFLLSSSSYY**CCCWLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG",
        "starts": "----------**-----------------------M----------------------------",
    },
    11: {
        "aas": "FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG",
        "starts": "---M------**--*----M------------MMMM---------------M------------",
    },
    12: {
        "aas": "FFLLSSSSYY**CC*WLLLSPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG",
        "starts": "----------**--*----M---------------M----------------------------",
    },
    13: {
        "aas": "FFLLSSSSYY**CCWWLLLLPPPPHHQQRRRRIIMMTTTTNNKKSSGGVVVVAAAADDEEGGGG",
        "starts": "---M------**----------------------MM---------------M------------",
    },
    14: {
        "aas": "FFLLSSSSYYY*CCWWLLLLPPPPHHQQRRRRIIIMTTTTNNNKSSSSVVVVAAAADDEEGGGG",
        "starts": "-----------*-----------------------M----------------------------",
    },
    15: {
        "aas": "FFLLSSSSYY*QCC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG",
        "starts": "----------*---*--------------------M----------------------------",
    },
    16: {
        "aas": "FFLLSSSSYY*LCC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG",
        "starts": "----------*---*--------------------M----------------------------",
    },
    21: {
        "aas": "FFLLSSSSYY**CCWWLLLLPPPPHHQQRRRRIIMMTTTTNNNKSSSSVVVVAAAADDEEGGGG",
        "starts": "----------**-----------------------M---------------M------------",
    },
    22: {
        "aas": "FFLLSS*SYY*LCC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG",
        "starts": "------*---*---*--------------------M----------------------------",
    },
    24: {
        "aas": "FFLLSSSSYY**CCWWLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSSKVVVVAAAADDEEGGGG",
        "starts": "---M------**-------M---------------M---------------M------------",
    },
    25: {
        "aas": "FFLLSSSSYY**CCGWLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG",
        "starts": "---M------**-----------------------M---------------M------------",
    },
}
ALLELE_ID_PATTERN = re.compile(r"^[A-Za-z0-9][A-Za-z0-9._-]*$")


def split_allele_name(allele_name: str) -> tuple[str, str]:
    normalized = allele_name.strip("\n")
    if "_" not in normalized:
        return normalized, normalized
    locus, allele_id = normalized.rsplit("_", 1)
    return locus, allele_id


def get_cds_name_from_allele_name(allele_name: str) -> str:
    return split_allele_name(allele_name)[0]


def get_allele_id_from_allele_name(allele_name: str) -> str:
    return split_allele_name(allele_name)[1]


def _reference_sort_key(record) -> tuple[int, object, str]:
    allele_id = get_allele_id_from_allele_name(record.id)
    if allele_id == REFERENCE_ALLELE_ID:
        return (0, 0, record.id)
    if allele_id.isdigit():
        return (2, int(allele_id), record.id)
    return (1, allele_id.lower(), record.id)


def select_reference_record(records: list):
    if records:
        return min(records, key=_reference_sort_key)
    raise ValueError("Cannot select a reference allele from an empty locus FASTA.")


def get_reference_allele_name(cds: str, allele_names: Iterable[str]) -> str:
    preferred = f"{cds}_{REFERENCE_ALLELE_ID}"
    allele_name_list = list(allele_names)
    if preferred in allele_name_list:
        return preferred

    candidates = [
        allele_name
        for allele_name in allele_name_list
        if get_cds_name_from_allele_name(allele_name) == cds
    ]
    if len(candidates) == 1:
        return candidates[0]

    raise ValueError(
        f"Unable to determine a unique reference allele for locus '{cds}'. "
        f"Candidates found: {', '.join(sorted(candidates)) or 'none'}."
    )


def get_reference_alignment_target_name(cds: str, target_names: Iterable[str]) -> str:
    preferred = f"{cds}_{REFERENCE_ALLELE_ID}"
    target_name_list = list(target_names)
    if preferred in target_name_list:
        return preferred

    candidates = [
        target_name
        for target_name in target_name_list
        if get_cds_name_from_allele_name(target_name) == cds
    ]
    if len(candidates) == 1:
        return candidates[0]

    raise ValueError(
        f"Unable to determine a unique alignment target for locus '{cds}'. "
        f"Targets found: {', '.join(sorted(candidates)) or 'none'}."
    )


def get_locus_length_mode_from_fasta(locus_fasta: str | Path, ref_seq: str) -> int:
    locus_path = Path(locus_fasta)
    if not locus_path.exists():
        return len(ref_seq)

    allele_lengths = []
    current_sequence: list[str] = []
    with locus_path.open("r", encoding="ascii") as handle:
        for line in handle:
            stripped = line.strip()
            if not stripped:
                continue
            if stripped.startswith(">"):
                if current_sequence:
                    allele_lengths.append(len("".join(current_sequence)))
                    current_sequence = []
                continue
            current_sequence.append(stripped)
    if current_sequence:
        allele_lengths.append(len("".join(current_sequence)))

    if not allele_lengths:
        return len(ref_seq)
    return Counter(allele_lengths).most_common(1)[0][0]


def get_translation_table_settings(
    translation_table: int = TRANSLATION_TABLE_DEFAULT,
    allowed_start_codons: set[str] | None = None,
) -> tuple[set[str], set[str]]:
    settings = TRANSLATION_TABLE_ROWS.get(translation_table)
    if settings is None:
        raise ValueError(f"Unsupported translation table: {translation_table}")
    codons = [f"{a}{b}{c}" for a, b, c in zip(BASE1, BASE2, BASE3)]
    stop_codons = {
        codon
        for codon, aa in zip(codons, settings["aas"])
        if aa == "*"
    }
    start_codons = allowed_start_codons or {
        codon
        for codon, marker in zip(codons, settings["starts"])
        if marker == "M"
    }
    return {codon.upper() for codon in start_codons}, stop_codons


def describe_quality_check(
    seq: str,
    ref_seq: str,
    locus_mode_length: int,
    translation_table: int = TRANSLATION_TABLE_DEFAULT,
    allowed_start_codons: set[str] | None = None,
) -> dict[str, object]:
    start_codons, stop_codons = get_translation_table_settings(
        translation_table=translation_table,
        allowed_start_codons=allowed_start_codons,
    )
    for i in range(0, len(seq) - 3, 3):
        if seq[i : i + 3] in stop_codons:
            return {
                "status": "LNF",
                "reason": "internal_stop_codon",
                "reference_match": False,
                "sequence_length": len(seq),
                "length_delta": len(seq) - locus_mode_length,
            }

    if len(seq) % 3 != 0:
        return {
            "status": "LNF",
            "reason": "length_not_divisible_by_three",
            "reference_match": False,
            "sequence_length": len(seq),
            "length_delta": len(seq) - locus_mode_length,
        }

    if seq[:3] not in start_codons:
        return {
            "status": "LNF",
            "reason": "invalid_start_codon",
            "reference_match": False,
            "sequence_length": len(seq),
            "length_delta": len(seq) - locus_mode_length,
        }

    if seq[-3:] not in stop_codons:
        return {
            "status": "LNF",
            "reason": "invalid_stop_codon",
            "reference_match": False,
            "sequence_length": len(seq),
            "length_delta": len(seq) - locus_mode_length,
        }

    if seq == ref_seq:
        return {
            "status": "EQ",
            "reason": "matches_reference",
            "reference_match": True,
            "sequence_length": len(seq),
            "length_delta": len(seq) - locus_mode_length,
        }

    lower_bound = locus_mode_length * (1 - LENGTH_MODE_TOLERANCE)
    upper_bound = locus_mode_length * (1 + LENGTH_MODE_TOLERANCE)
    if len(seq) < lower_bound:
        return {
            "status": "ASM",
            "reason": "below_length_mode_tolerance",
            "reference_match": False,
            "sequence_length": len(seq),
            "length_delta": len(seq) - locus_mode_length,
        }
    if len(seq) > upper_bound:
        return {
            "status": "ALM",
            "reason": "above_length_mode_tolerance",
            "reference_match": False,
            "sequence_length": len(seq),
            "length_delta": len(seq) - locus_mode_length,
        }
    return {
        "status": "Q",
        "reason": "passes_orf_qc",
        "reference_match": False,
        "sequence_length": len(seq),
        "length_delta": len(seq) - locus_mode_length,
    }


def quality_check(
    seq: str,
    ref_seq: str,
    locus_mode_length: int,
    translation_table: int = TRANSLATION_TABLE_DEFAULT,
    allowed_start_codons: set[str] | None = None,
) -> str:
    return describe_quality_check(
        seq,
        ref_seq,
        locus_mode_length,
        translation_table=translation_table,
        allowed_start_codons=allowed_start_codons,
    )["status"]


def get_numeric_allele_ids_from_fasta(locus_fasta: str | Path) -> list[int]:
    locus_path = Path(locus_fasta)
    if not locus_path.exists():
        return []

    numeric_ids = []
    with locus_path.open("r", encoding="ascii") as handle:
        for raw_line in handle:
            line = raw_line.strip()
            if not line.startswith(">"):
                continue
            allele_id = get_allele_id_from_allele_name(line[1:])
            if allele_id.isdigit():
                numeric_ids.append(int(allele_id))
    return numeric_ids


def get_next_novel_allele_id_from_fasta(locus_fasta: str | Path) -> int:
    numeric_ids = get_numeric_allele_ids_from_fasta(locus_fasta)
    if numeric_ids:
        return max(numeric_ids) + 1
    return 2


def is_valid_allele_identifier(allele_id: str) -> bool:
    return bool(ALLELE_ID_PATTERN.match(allele_id))
