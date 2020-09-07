import logging

import click

from varcomb.parsers import parse_vcf_file
from varcomb.utilities import read_vcf


@click.group()
@click.option('--log_file')
@click.pass_context
def client(cts, log_file, type=click.Path()):
    logging.basicConfig(format='[%(levelname)s] %(asctime)s %(message)s', datefmt='%Y/%m/%d %H:%M:%S', level=logging.INFO, filename=log_file)


@client.command()
@click.option('--vcf_file1', required=True, type=click.Path())
@click.option('--vcf_file2', required=True, type=click.Path())
@click.option('--vcf_out', required=True, type=click.Path())
@click.pass_context
def merge_vcfs(ctx, vcf_file1, vcf_file2, vcf_out):
    logging.info(f'Reading file: {vcf_file1}')
    vcf1 = parse_vcf_file(read_vcf(vcf_file1))
    logging.info(f'Reading file: {vcf_file2}')
    vcf2 = parse_vcf_file(read_vcf(vcf_file2))
    vcf = vcf1 + vcf2
    vcf = vcf.remove_true_duplicates()
    vcf.to_file(vcf_out)


def run():
    client(obj={})
