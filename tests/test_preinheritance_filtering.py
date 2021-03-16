import unittest

from tests.test_utils import create_test_snv
from filtering.preinheritance_filtering import preinheritance_filter

class TestPreInheritanceFilter(unittest.TestCase):

    def setUp(self):
        self.maxDiff = None
        self.vardata = {'chrom': '1', 'pos': '100000', 'ref': 'A', 'alt': 'G',
                        'consequence': 'missense_variant', 'ensg': 'ensg',
                        'symbol': 'KMTD2', 'feature': 'feature',
                        'canonical': 'YES', 'mane': 'MANE', 'hgnc_id': '123',
                        'max_af': '0', 'max_af_pops': '.', 'ddd_af': '0',
                        'revel': '1', 'polyphen': '.', 'hgvsc': '.',
                        'hgvsp': '.', 'gender': 'M', 'denovo_snv': False,
                        'denovo_indel': False, 'gt': '0/1', 'gq': '50'}

    def test_min_gq(self):
        # if GQ < 40 a variant should fail if in an autosome
        testvar = create_test_snv(self.vardata)
        variants = {'child': {'1_100000_A_G': testvar},
                    'mum': {}, 'dad': {}}
        variants_per_gene = preinheritance_filter(variants)
        self.assertEqual(variants_per_gene, {'123': {'1_100000_A_G': {
            'child': variants['child']['1_100000_A_G']}}})

        self.vardata['gq'] = '30'
        testvar = create_test_snv(self.vardata)
        variants = {'child': {'1_100000_A_G': testvar},
                    'mum': {}, 'dad': {}}
        variants_per_gene = preinheritance_filter(variants)
        self.assertEqual(variants_per_gene, {})

        self.vardata['chrom'] = 'X'
        testvar = create_test_snv(self.vardata)
        variants = {'child': {'X_100000_A_G': testvar},
                    'mum': {}, 'dad': {}}
        variants_per_gene = preinheritance_filter(variants)
        self.assertEqual(variants_per_gene, {'123': {'X_100000_A_G': {
            'child': variants['child']['X_100000_A_G']}}})

    def test_max_ddd_af(self):
        # fail if ddd_af > 0.005
        self.vardata['ddd_af'] = '0.5'
        testvar = create_test_snv(self.vardata)
        variants = {'child': {'1_100000_A_G': testvar},
                    'mum': {}, 'dad': {}}
        variants_per_gene = preinheritance_filter(variants)
        self.assertEqual(variants_per_gene, {})

    def test_cq (self):
        # fail if consequence is not functional
        self.vardata['consequence'] = 'synonymous'
        testvar = create_test_snv(self.vardata)
        variants = {'child': {'1_100000_A_G': testvar},
                    'mum': {}, 'dad': {}}
        variants_per_gene = preinheritance_filter(variants)
        self.assertEqual(variants_per_gene, {})

        self.vardata['consequence'] = 'synonymous&missense_variant'
        testvar = create_test_snv(self.vardata)
        variants = {'child': {'1_100000_A_G': testvar},
                    'mum': {}, 'dad': {}}
        variants_per_gene = preinheritance_filter(variants)
        self.assertEqual(variants_per_gene, {'123': {'1_100000_A_G': {
            'child': variants['child']['1_100000_A_G']}}})

    def test_revel(self):
        pass




