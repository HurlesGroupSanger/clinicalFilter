import unittest

from tests.test_utils import create_test_candidate_vars
from tests.test_utils import create_test_person
from tests.test_utils import create_test_family

from filtering.postinheritance_filter import postinheritance_filter

class TestPostInheritanceFilter(unittest.TestCase):

    def setUp(self):
        self.maxDiff = None
        self.vardata = {'chrom': '2', 'pos': '100000', 'ref': 'A', 'alt': 'GG',
                        'consequence': 'start_lost', 'ensg': 'ensg',
                        'symbol': 'MECP2', 'feature': 'feature',
                        'canonical': 'YES', 'mane': 'MANE', 'hgnc_id': '1234',
                        'max_af': '0', 'max_af_pops': '.', 'ddd_af': '0',
                        'revel': '.', 'polyphen': '.', 'hgvsc': '.',
                        'hgvsp': '.', 'gender': 'M', 'denovo_snv': False,
                        'denovo_indel': False, 'gt': '0/1', 'gq': '50', 'pid': '.'}

        self.child = create_test_person('fam', 'child_id', 'dad_id', 'mum_id',
                                        'M', '2', '/vcf/path')
        self.mum = create_test_person('fam', 'mum_id', '0', '0', 'F', '1',
                                      '/vcf/path')
        self.dad = create_test_person('fam', 'dad_id', '0', '0', 'M', '1',
                                      '/vcf/path')


    def test_maf_filter_both_parents(self):
        # MAF threshold for non-biallelic variants where family has both parents is
        # 0.0005, where there are no parents threshold is 0.0001

        family = create_test_family(self.child, self.mum, self.dad)
        test_single_vars = {'2_100000_A_G': {'variant': self.vardata, 'mode':'Monoallelic'}}
        test_compound_hets = {}
        test_candidate_vars = create_test_candidate_vars(test_single_vars, test_compound_hets)
        test_candidate_vars_copy = test_candidate_vars.copy()

        postinheritance_filter(test_candidate_vars, family)

        self.assertEqual(test_candidate_vars, test_candidate_vars_copy)

        # modify test variant so that it fails on ddd allele frequency

        self.vardata['ddd_af'] = '0.1'
        test_single_vars = {
            '2_100000_A_G': {'variant': self.vardata, 'mode': 'Monoallelic'}}
        test_compound_hets = {}
        test_candidate_vars = create_test_candidate_vars(test_single_vars,
                                                         test_compound_hets)
        test_candidate_vars_empty = create_test_candidate_vars({},{})
        postinheritance_filter(test_candidate_vars, family)

        self.assertEqual(test_candidate_vars, test_candidate_vars_empty)

        # modify test variant so that it fails on max allele frequency and not ddd

        self.vardata['ddd_af'] = '0.00008'
        self.vardata['max_af'] = '0.008'
        test_single_vars = {
            '2_100000_A_G': {'variant': self.vardata, 'mode': 'Monoallelic'}}
        test_compound_hets = {}
        test_candidate_vars = create_test_candidate_vars(test_single_vars,
                                                         test_compound_hets)
        postinheritance_filter(test_candidate_vars, family)

        self.assertEqual(test_candidate_vars, test_candidate_vars_empty)

        # modify test variant so that it fails on both max and ddd allele frequency

        self.vardata['ddd_af'] = '0.008'
        self.vardata['max_af'] = '0.008'
        test_single_vars = {
            '2_100000_A_G': {'variant': self.vardata, 'mode': 'Monoallelic'}}
        test_compound_hets = {}
        test_candidate_vars = create_test_candidate_vars(test_single_vars,
                                                         test_compound_hets)
        postinheritance_filter(test_candidate_vars, family)

        self.assertEqual(test_candidate_vars, test_candidate_vars_empty)


    def test_maf_filter_no_parents(self):
        # MAF threshold for non-biallelic variants where family has both parents is
        # 0.0005, where there are no parents threshold is 0.0001
        family = create_test_family(self.child, None, None)

        family = create_test_family(self.child, self.mum, self.dad)
        test_single_vars = {'2_100000_A_G': {'variant': self.vardata, 'mode':'Monoallelic'}}
        test_compound_hets = {}
        test_candidate_vars = create_test_candidate_vars(test_single_vars, test_compound_hets)
        test_candidate_vars_copy = test_candidate_vars.copy()

        postinheritance_filter(test_candidate_vars, family)

        self.assertEqual(test_candidate_vars, test_candidate_vars_copy)

        # modify test variant so that it fails on ddd allele frequency

        self.vardata['ddd_af'] = '0.1'
        test_single_vars = {
            '2_100000_A_G': {'variant': self.vardata, 'mode': 'Monoallelic'}}
        test_compound_hets = {}
        test_candidate_vars = create_test_candidate_vars(test_single_vars,
                                                         test_compound_hets)
        test_candidate_vars_empty = create_test_candidate_vars({},{})
        postinheritance_filter(test_candidate_vars, family)

        self.assertEqual(test_candidate_vars, test_candidate_vars_empty)

        # modify test variant so that it fails on max allele frequency and not ddd

        self.vardata['ddd_af'] = '0.00008'
        self.vardata['max_af'] = '0.008'
        test_single_vars = {
            '2_100000_A_G': {'variant': self.vardata, 'mode': 'Monoallelic'}}
        test_compound_hets = {}
        test_candidate_vars = create_test_candidate_vars(test_single_vars,
                                                         test_compound_hets)
        postinheritance_filter(test_candidate_vars, family)

        self.assertEqual(test_candidate_vars, test_candidate_vars_empty)

        # modify test variant so that it fails on both max and ddd allele frequency

        self.vardata['ddd_af'] = '0.008'
        self.vardata['max_af'] = '0.008'
        test_single_vars = {
            '2_100000_A_G': {'variant': self.vardata, 'mode': 'Monoallelic'}}
        test_compound_hets = {}
        test_candidate_vars = create_test_candidate_vars(test_single_vars,
                                                         test_compound_hets)
        postinheritance_filter(test_candidate_vars, family)

        self.assertEqual(test_candidate_vars, test_candidate_vars_empty)



