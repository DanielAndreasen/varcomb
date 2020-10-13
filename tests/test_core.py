import unittest

from varcomb.core import VCF, Info, Location, VCFrow
from varcomb.exceptions import LocationShiftError


class TestCoreLocation(unittest.TestCase):

    def test_location(self):
        loc = Location(chrom='chr1', pos=42)
        self.assertEqual(loc.chrom, 'chr1')
        self.assertEqual(loc.pos, 42)

    def test_locations_distance_same_chrom(self):
        loc1 = Location(chrom='chr7', pos=1337)
        loc2 = Location(chrom='chr7', pos=1001)
        distance = loc1 - loc2
        inverse_distance = loc2 - loc1
        self.assertEqual(distance, 1337 - 1001)
        self.assertEqual(inverse_distance, 1001 - 1337)

    def test_locations_distance_different_chrom(self):
        loc1 = Location(chrom='chr4', pos=1337)
        loc2 = Location(chrom='chr7', pos=1001)
        distance = loc1 - loc2
        inverse_distance = loc2 - loc1
        self.assertEqual(distance, None)
        self.assertEqual(inverse_distance, None)

    def test_location_shift(self):
        loc = Location(chrom='chr1', pos=42)
        shift = -2
        shifted_loc1 = loc + shift
        shifted_loc2 = loc.shift(shift)
        self.assertIsInstance(shifted_loc1, Location)
        self.assertEqual(shifted_loc1.pos, loc.pos + shift)
        self.assertEqual(shifted_loc1, shifted_loc2)

    def test_location_shift_error(self):
        loc = Location(chrom='chr1', pos=42)
        shift = 'not_valid'
        with self.assertRaises(LocationShiftError):
            loc + shift

    def test_location_less_than_same_chrom(self):
        loc1 = Location(chrom='chr7', pos=1337)
        loc2 = Location(chrom='chr7', pos=1001)
        self.assertEqual(loc2 < loc1, True)
        self.assertEqual(loc1 > loc2, True)

    def test_location_less_than_different_chrom(self):
        loc1 = Location(chrom='chr7', pos=42)
        loc2 = Location(chrom='chr20', pos=1001)
        self.assertEqual(loc1 < loc2, True)

    def test_location_less_than_sex_chrom(self):
        loc1 = Location(chrom='7', pos=42)
        loc2 = Location(chrom='X', pos=1001)
        self.assertEqual(loc1 < loc2, True)

    def test_location_greater_than_sex_chrom(self):
        loc1 = Location(chrom='Y', pos=42)
        loc2 = Location(chrom='4', pos=1001)
        self.assertEqual(loc2 > loc1, False)


class TestCoreInfo(unittest.TestCase):

    def setUp(self):
        self.loc = Location(chrom='chr1', pos=42)
        self.info = Info('k1=v1;k2=v2;k3=v3')
        self.vcfrow = VCFrow(loc=self.loc, id='1234', ref='A', alt='G', qual='.',
                             filter='PASS', info=self.info, format='format',
                             samples=['sample1', 'sample2'])

    def test_infofield(self):
        self.assertIsInstance(self.vcfrow.info, Info)
        self.assertEqual(len(self.vcfrow.info), 3)
        self.assertIn('k1', self.vcfrow.info.keys())
        self.assertIn('v2', self.vcfrow.info.values())

    def test_infofield_add(self):
        self.vcfrow.info['k4'] = 'v4'
        self.assertEqual(len(self.vcfrow.info), 4)
        self.assertIn('k4', self.vcfrow.info.keys())
        self.assertIn('v4', self.vcfrow.info.values())

    def test_infofield_get(self):
        self.assertEqual(self.vcfrow.info['k1'], 'v1')
        self.assertEqual(self.vcfrow.info['k2'], 'v2')
        self.assertEqual(self.vcfrow.info['k3'], 'v3')

    def test_infofield_to_str(self):
        self.assertIsInstance(str(self.vcfrow.info), str)
        self.assertEqual(str(self.vcfrow.info), 'k1=v1;k2=v2;k3=v3')

    def test_infofield_no_value(self):
        info = Info('k1=v1;k2=v2;k3=v3;field1;field2')
        self.assertEqual(len(info), 5)
        self.assertEqual(info['field1'], True)
        self.assertEqual(info['field2'], True)


class TestCoreVCFrow(unittest.TestCase):

    def setUp(self):
        self.loc = Location(chrom='chr1', pos=42)
        self.vcfrow = VCFrow(loc=self.loc, id='1234', ref='A', alt='G', qual='.',
                             filter='PASS', info='info', format='format',
                             samples=['sample1', 'sample2'])

    def test_vcfrow(self):
        self.assertEqual(self.vcfrow.loc, self.loc)
        self.assertEqual(self.vcfrow.id, '1234')
        self.assertEqual(self.vcfrow.ref, 'A')
        self.assertEqual(self.vcfrow.alt, 'G')
        self.assertEqual(self.vcfrow.qual, '.')
        self.assertEqual(self.vcfrow.filter, 'PASS')
        self.assertEqual(self.vcfrow.info, 'info')
        self.assertEqual(self.vcfrow.format, 'format')
        self.assertEqual(self.vcfrow.samples, ['sample1', 'sample2'])

    def test_format_line(self):
        actual = self.vcfrow._format_row()
        row = [self.vcfrow.loc.chrom, self.vcfrow.loc.pos, self.vcfrow.id, self.vcfrow.ref,
               self.vcfrow.alt, self.vcfrow.qual, self.vcfrow.filter, self.vcfrow.info, self.vcfrow.format]
        expected = '\t'.join(map(str, row + self.vcfrow.samples))
        self.assertEqual(actual, expected)


class TestCoreVCF(unittest.TestCase):

    def setUp(self):
        loc1 = Location(chrom='chr2', pos=21)
        loc2 = Location(chrom='chr16', pos=50)
        self.vcfrow1 = VCFrow(loc=loc1, id='id1', ref='G', alt='A', qual='.',
                              filter='PASS', info='info', format='format',
                              samples=['sample1', 'sample2'])
        self.vcfrow2 = VCFrow(loc=loc2, id='id2', ref='T', alt='ATTGC', qual='.',
                              filter='PASS', info='info', format='format',
                              samples=['sample1', 'sample2'])
        self.vcf = VCF(rows=[self.vcfrow1, self.vcfrow2])

    def test_vcf(self):
        self.assertIsInstance(self.vcf.rows, list)
        self.assertEqual(len(self.vcf), len(self.vcf.rows))
        self.assertEqual(self.vcf[0], self.vcf.rows[0])
        self.assertEqual(self.vcf[0], self.vcfrow1)

    def test_vcf_with_header(self):
        header = '# This is a header'
        vcf = VCF(rows=[self.vcfrow1, self.vcfrow2], header=header)
        self.assertIsNotNone(vcf.header)
        self.assertEqual(len(vcf), 2)

    def test_different_vcfs(self):
        vcfrow = self.vcfrow1
        vcfrow.ref = 'T'
        new_vcf = VCF(rows=[vcfrow])
        self.assertNotEqual(self.vcf, new_vcf)

    def test_get_from_chrom(self):
        vcf = self.vcf + self.vcf
        chrom = 'chr2'
        vcf_chrom2 = vcf.get_from_chrom(chrom)
        self.assertEqual(len(vcf_chrom2), 2)
        self.assertEqual(vcf_chrom2[0].loc.chrom, chrom)

    def test_get_near_location(self):
        chrom = 'chr16'
        pos = 55
        vcf1 = self.vcf.get_near_location(chrom, pos, tol=10)
        vcf2 = self.vcf.get_near_location(chrom, pos, tol=2)
        vcf3 = self.vcf.get_near_location('chr2', pos, tol=10)
        self.assertEqual(len(vcf1), 1)
        self.assertEqual(len(vcf2), 0)
        self.assertEqual(len(vcf3), 0)
