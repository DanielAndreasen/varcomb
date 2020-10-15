![CI](https://github.com/DanielAndreasen/varcomb/workflows/CI/badge.svg?branch=master)

# varcomb

This is a tool to work with VCF files, and to combine VCF files from different callers.


# Installation
The current recommend way of installing `varcomb` is using `pip`:

```bash
pip install git+https://github.com/DanielAndreasen/varcomb
```

# Test
Read the tests to understand what this tool does. You can run the tests by cloning this repository:

```bash
git clone https://github.com/DanielAndreasen/varcomb
cd varcomb
python setup.py test
```

# Usage
After installation, `varcomb` is available from the command line:

```bash
$ varcomb merge-vcfs \
    --vcf_file1 vcf1.vcf.gz \
    --vcf_file2 vcf2.vcf.gz \
    --vcf_out combined.vcf \
    --ann_vcf1 FIRST \ # Optional
    --ann_vcf2 SECOND  # Optional
```
This will generate a combined VCF called `combined.vcf`. Note that `varcomb` does not
compress the final VCF file. It is recommended the user does this:

```bash
$ bgzip combined.vcf
$ tabix combined.vcf.gz
```
which will also generate the index file `combined.vcf.gz.tbi`.
