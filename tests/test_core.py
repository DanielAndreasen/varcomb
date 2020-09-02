import unittest

from varcomb.core import VCF, Location, VCFrow
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


class TestCoreVCF(unittest.TestCase):

    def setUp(self):
        loc1 = Location(chrom='chr2', pos=21)
        loc2 = Location(chrom='chr16', pos=543245)
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
