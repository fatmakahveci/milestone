from __future__ import annotations

from wgmlst_utils import get_cds_name_from_allele_name


class Paf:
    def __init__(self, paf_line):
        fields = paf_line.strip("\n").split("\t")

        self.qname = fields[0]
        self.qlen = int(fields[1])
        self.qstart = int(fields[2])
        self.qend = int(fields[3])
        self.strand = fields[4]
        self.tname = fields[5]
        self.tlen = int(fields[6])
        self.tstart = int(fields[7])
        self.tend = int(fields[8])
        self.nmatch = int(fields[9])
        self.alen = int(fields[10])
        self.mapq = int(fields[11])

    def __repr__(self):
        return f"query sequence name: {self.qname}\t" \
               f"query sequence length: {self.qlen}\t" \
               f"query start: {self.qstart}\t" \
               f"query end: {self.qend}\t" \
               f"strand (+ or -): {self.strand}\t" \
               f"target sequence name: {self.tname}\t" \
               f"target sequence length: {self.tlen}\t" \
               f"target start: {self.tstart}\t" \
               f"target end: {self.tend}\t" \
               f"number of matching bases in the mapping: {self.nmatch}\t" \
               f"number of bases including gaps in the mapping: {self.alen}\t" \
               f"mapping quality:{self.mapq}\n"


class Vcf:
    def __init__(self, vcf_line):
        fields = vcf_line.split("\t")

        self.chr = get_cds_name_from_allele_name(allele_name=fields[0])
        self.pos = fields[1]
        self.id = fields[2]
        self.ref = fields[3]
        self.alt = fields[4]
        self.qual = fields[5]
        self.filter = fields[6]
        self.info = fields[7]
        self.sample_format = fields[8]
        self.sample = fields[9]

    def __repr__(self):
        return f"chr: {self.chr}\t" \
               f"pos: {self.pos}\t" \
               f"id: {self.id}\t" \
               f"ref: {self.ref}\t" \
               f"alt: {self.alt}\t" \
               f"qual: {self.qual}\t" \
               f"filter: {self.filter}\t" \
               f"info: {self.info}\t" \
               f"format: {self.sample_format}\t" \
               f"sample: {self.sample}\n"
