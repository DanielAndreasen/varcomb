import gzip

from varcomb.exceptions import VCFFileNotSupported


def read_vcf(fname):
    if fname.endswith('.vcf'):
        with open(fname, 'r') as f:
            return [row for row in f.read().split('\n') if row]
    elif fname.endswith('.vcf.gz'):
        with gzip.open(fname, 'rb') as f:
            return [row for row in f.read().decode().split('\n') if row]
    raise VCFFileNotSupported(f'{fname} does not end on ".vcf" or ".vcf.gz"')
