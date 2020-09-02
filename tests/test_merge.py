import unittest

from varcomb.parsers import parse_vcf_file


class TestMergeVCF(unittest.TestCase):
    def setUp(self):
        self.vcf_file1 = [
            '\t'.join(['chr1', '16688', '.', 'G', 'A', '.', 'map_qual;normal_artifact;strand_bias', 'AS_FilterStatus=map_qual,strand_bias;AS_SB_TABLE=327,30|0,13;DP=380;ECNT=3;GERMQ=93;MBQ=37,34;MFRL=260,332;MMQ=22,22;MPOS=21;NALOD=-1.391e+00;NLOD=43.74;POPAF=6.00;ROQ=80;TLOD=7.68', 'GT:AD:AF:DP:F1R2:F2R1:SB', '0/0:180,4:0.026:184:99,1:80,3:165,15,0,4', '0/1:177,9:0.052:186:98,4:77,5:162,15,0,9']),
            '\t'.join(['chr1', '186478', '.', 'A', 'G', '.', 'map_qual;strand_bias;weak_evidence', 'AS_FilterStatus=weak_evidence,map_qual,strand_bias;AS_SB_TABLE=200,282|6,0;DP=499;ECNT=1;GERMQ=93;MBQ=20,30;MFRL=206,389;MMQ=46,24;MPOS=35;NALOD=1.64;NLOD=50.18;POPAF=6.00;ROQ=61;TLOD=3.43', 'GT:AD:AF:DP:F1R2:F2R1:SB', '0/0:246,1:0.010:247:134,0:109,1:101,145,1,0', '0/1:236,5:0.033:241:109,2:123,3:99,137,5,0']),
            '\t'.join(['chr1', '16577291', '.', 'T', 'C', '.', 'haplotype;normal_artifact', 'AS_FilterStatus=SITE;AS_SB_TABLE=594,487|13,10;DP=1121;ECNT=3;GERMQ=93;MBQ=20,20;MFRL=204,172;MMQ=48,42;MPOS=52;NALOD=-4.662e+01;NLOD=27.62;POPAF=6.00;ROQ=51;TLOD=24.07', 'GT:AD:AF:DP:F1R2:F2R1:PGT:PID:PS:SB', '0|0:433,15:0.035:448:209,12:210,3:0|1:16577281_C_G:16577281:237,196,9,6', '0|1:648,8:0.011:656:323,6:316,2:0|1:16577281_C_G:16577281:357,291,4,4'])]
        self.vcf_file2 = self.vcf_file1

        self.vcf1 = parse_vcf_file(self.vcf_file1)
        self.vcf2 = parse_vcf_file(self.vcf_file2)

    def test_merge_two_vcfs(self):
        vcf = self.vcf1 + self.vcf2
        self.assertEqual(vcf, parse_vcf_file(self.vcf_file1 + self.vcf_file2))

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
