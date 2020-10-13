from typing import List

from varcomb.core import VCF, Info, Location, VCFrow


def _parse_vcf_line(line: str) -> VCFrow:
    elements = line.split('\t')
    loc = Location(chrom=elements[0], pos=int(elements[1]))
    return VCFrow(loc=loc, id=elements[2], ref=elements[3], alt=elements[4],
                  qual=elements[5], filter=elements[6], info=Info(elements[7]),
                  format=elements[8], samples=elements[9:])


def parse_vcf_file(stream: List[str]) -> VCF:
    return VCF(rows=[_parse_vcf_line(row) for row in stream if not row.startswith('#')],
               header=[row for row in stream if row.startswith('#')])
