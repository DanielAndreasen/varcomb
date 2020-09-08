import os
import unittest

from varcomb.parsers import parse_vcf_file


class TestMergeVCF(unittest.TestCase):
    def setUp(self):
        self.vcf_file1 = [
            '\t'.join(['chr1', '16688', '.', 'G', 'A', '.', '', '', '', '', '']),
            '\t'.join(['chr1', '186478', '.', 'A', 'G', '.', '', '', '', '', '']),
            '\t'.join(['chr1', '16577291', '.', 'T', 'C', '.', '', '', '', '', ''])]

        self.vcf1 = parse_vcf_file(self.vcf_file1)
        self.vcf2 = parse_vcf_file(self.vcf_file1)

    def test_merge_two_vcfs(self):
        vcf = self.vcf1 + self.vcf2
        self.assertEqual(vcf, parse_vcf_file(self.vcf_file1 * 2))

    def test_remove_true_duplicates(self):
        vcf = self.vcf1 + self.vcf2
        vcf_dup = vcf.remove_true_duplicates()
        self.assertEqual(vcf_dup, self.vcf1)

    def test_remove_duplicates_at_same_loc(self):
        vcf_loc_dup = ['\t'.join(['chr1', '186478', '.', 'T', 'G', '.', '', '', '', '', ''])]
        vcf = parse_vcf_file(vcf_loc_dup) + self.vcf1
        vcf_rem_true_dup = vcf.remove_true_duplicates()
        vcf_rem_loc_dup = vcf.remove_loc_dup()

        self.assertEqual(len(vcf), len(self.vcf1) + len(vcf_loc_dup))
        self.assertEqual(len(vcf_rem_true_dup), len(vcf))
        self.assertEqual(vcf_rem_true_dup, vcf)

        self.assertEqual(len(vcf_rem_loc_dup), len(self.vcf1))


class TestMergeVCFToFile(unittest.TestCase):
    def setUp(self):
        self.vcf_file1 = [
            '\t'.join(['chr1', '16688', '.', 'G', 'A', '.', '', '', '', '', '']),
            '\t'.join(['chr1', '186478', '.', 'A', 'G', '.', '', '', '', '', '']),
            '\t'.join(['chr1', '16577291', '.', 'T', 'C', '.', '', '', '', '', ''])]

        self.vcf1 = parse_vcf_file(self.vcf_file1)
        self.vcf2 = parse_vcf_file(self.vcf_file1)

    def tearDown(self):
        fname = self.fname
        if os.path.exists(fname):
            os.remove(fname)

    def test_save_to_file(self):
        vcf = self.vcf1 + self.vcf2
        vcf_merged = vcf.remove_true_duplicates()
        self.fname = 'test.vcf'
        vcf_merged.to_file(self.fname)
        self.assertTrue(os.path.exists(self.fname))
