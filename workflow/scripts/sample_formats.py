from __future__ import annotations

from wgmlst_utils import get_cds_name_from_allele_name


class Coverage:
    def __init__(self, coverage_line):
        fields = coverage_line.strip("\n").split("\t")

        self.cds = fields[0]
        self.length = int(fields[2])
        self.number_of_reads = int(fields[3])
        self.covered_bases = int(fields[4])
        self.coverage = float(fields[5])
        self.mean_depth = float(fields[6])
        self.mean_baseq = float(fields[7])
        self.mean_mapq = float(fields[8])

    def __repr__(self):
        return f"cds: {self.cds}\n" \
               f"length: {self.length}\n" \
               f"number of reads: {self.number_of_reads}\n" \
               f"covered bases: {self.covered_bases}\n" \
               f"coverage: {self.coverage}\n" \
               f"mean depth: {self.mean_depth}\n" \
               f"mean baseq: {self.mean_baseq}\n" \
               f"mean_mapq: {self.mean_mapq}\n"


class Sam:
    def __init__(self, sam_line):
        fields = sam_line.strip("\n").split("\t")

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
        return f"QNAME: {self.query_name}\t" \
               f"FLAG: {self.read_info}\t" \
               f"RNAME: {self.target_sequence_name}\t" \
               f"POS: {self.position}\t" \
               f"MAPQ: {self.mapping_quality}\t" \
               f"CIGAR: {self.cigar}\t" \
               f"RNEXT: {self.rnext}\t" \
               f"PNEXT: {self.pnext}\t" \
               f"TLEN: {self.tlen}\t" \
               f"SEQ: {self.sequence}\t" \
               f"QUAL: {self.read_quality}"


class Vcf:
    def __init__(self, vcf_line):
        fields = vcf_line.strip("\n").split("\t")

        self.chr = get_cds_name_from_allele_name(allele_name=fields[0])
        self.pos = int(fields[1])
        self.id = fields[2]
        self.ref = fields[3]
        self.alt = fields[4].split(",")[0]
        self.qual = float(fields[5])
        self.filter = fields[6]
        self.info = fields[7]
        self.sample_format = fields[8]
        self.sample = fields[9]

    def __str__(self):
        return f"{self.chr}\t" \
               f"{self.pos}\t" \
               f"{self.id}\t" \
               f"{self.ref}\t" \
               f"{self.alt}\t" \
               f"{self.qual}\t" \
               f"{self.filter}\t" \
               f"{self.info}\t" \
               f"{self.sample_format}\t" \
               f"{self.sample}\n"

    def __repr__(self):
        return f"chr: {self.chr}\n" \
               f"pos: {str(self.pos)}\n" \
               f"id: {self.id}\n" \
               f"ref: {self.ref}\n" \
               f"alt: {self.alt}\n" \
               f"qual: {str(self.qual)}\n" \
               f"filter: {self.filter}\n" \
               f"info: {self.info}\n" \
               f"format: {self.sample_format}\n" \
               f"sample: {self.sample}"
