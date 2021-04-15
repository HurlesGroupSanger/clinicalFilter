import unittest

from tests.test_utils import create_test_candidate_vars
from tests.test_utils import create_test_person
from tests.test_utils import create_test_family

from filtering.inheritance_filtering import InheritanceFiltering

class TestAutosomalInheritanceFilter(unittest.TestCase):

    def setUp(self):
        self.maxDiff = None
        self.vardata = {'chrom': '2', 'pos': '100000', 'ref': 'A', 'alt': 'GG',
                        'consequence': 'start_lost', 'ensg': 'ensg',
                        'symbol': 'MECP2', 'feature': 'feature',
                        'canonical': 'YES', 'mane': 'MANE', 'hgnc_id': '1234',
                        'max_af': '0', 'max_af_pops': '.', 'ddd_af': '0',
                        'revel': '.', 'polyphen': '.', 'hgvsc': '.',
                        'hgvsp': '.', 'sex': 'XY', 'denovo_snv': False,
                        'denovo_indel': False, 'gt': '0/1', 'gq': '50', 'pid': '.'}

        self.child = create_test_person('fam', 'child_id', 'dad_id', 'mum_id',
                                        'XY', '2', '/vcf/path')
        self.mum = create_test_person('fam', 'mum_id', '0', '0', 'XX', '1',
                                      '/vcf/path')
        self.mum_aff = create_test_person('fam', 'mum_id', '0', '0', 'XX', '2',
                                      '/vcf/path')
        self.dad = create_test_person('fam', 'dad_id', '0', '0', 'XY', '1',
                                      '/vcf/path')
        self.dad_aff = create_test_person('fam', 'dad_id', '0', '0', 'XY', '2',
                                      '/vcf/path')
        self.family_unaff = create_test_family(self.child, self.mum, self.dad)


    def test_biallelic_heterozygous_parents_filter(self):
        pass

    def test_biallelic_homozygous_parents_filter(self):
        pass

    def test_monoallelic_heterozygous_parents_filter(self):
        pass

    def test_monoallelic_homozygous_parents_filter(self):
        pass

    def test_imprinted_heterozygous_parents_filter(self):
        pass

    def test_imprinted_homozygous_parents_filter(self):
        pass

    def test_mosaic_heterozygous_parents_filter(self):
        pass

    def test_mosaic_homozygous_parents_filter(self):
        pass

    def test_biallelic_heterozygous_no_parents_filter(self):
        pass

    def test_biallelic_homozygous_no_parents_filter(self):
        pass

    def test_monoallelic_heterozygous_no_parents_filter(self):
        pass

    def test_monoallelic_homozygous_no_parents_filter(self):
        pass

    def test_imprinted_heterozygous_no_parents_filter(self):
        pass

    def test_imprinted_homozygous_no_parents_filter(self):
        pass

    def test_mosaic_heterozygous_no_parents_filter(self):
        pass

    def test_mosaic_homozygous_no_parents_filter(self):
        pass

