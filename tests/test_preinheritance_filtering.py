"""
Copyright (c) 2021 Genome Research Limited
Author: Ruth Eberhardt <re3@sanger.ac.uk>

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in
all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
THE SOFTWARE.
"""

import unittest

from tests.test_utils import create_test_snv
from filtering.preinheritance_filtering import PreInheritanceFiltering


class TestPreInheritanceFilter(unittest.TestCase):

    def setUp(self):
        self.maxDiff = None
        self.vardata = {'chrom': '1', 'pos': '100000', 'ref': 'A', 'alt': 'G',
                        'consequence': 'missense_variant', 'ensg': 'ensg',
                        'symbol': 'KMTD2', 'feature': 'feature',
                        'canonical': 'YES', 'mane': 'MANE', 'hgnc_id': '123',
                        'max_af': '0', 'max_af_pops': '.', 'ddd_af': '0',
                        'revel': '1', 'polyphen': '.', 'hgvsc': '.',
                        'hgvsp': '.', 'sex': 'XY', 'dnm': '.',
                        'gt': '0/1', 'gq': '50', 'ac_XX': '2', 'ac_XY': '5',
                        'nhomalt_XX': '379', 'nhomalt_XY': '4',
                        'ddd_father_af': '.'}
        self.vardataX = self.vardata.copy()
        self.vardataX['chrom'] = 'X'
        self.vardataX['ddd_father_af'] = '0.004'

    def test_min_gq(self):
        # if GQ < 40 a variant should fail if in an autosome
        testvar = create_test_snv(self.vardata)
        variants = {'child': {'1_100000_A_G': testvar},
                    'mum': {}, 'dad': {}}
        preinheritancefilter = PreInheritanceFiltering(variants)
        variants_per_gene = preinheritancefilter.preinheritance_filter()
        self.assertEqual(variants_per_gene, {'123': {'1_100000_A_G': {
            'child': variants['child']['1_100000_A_G']}}})

        self.vardata['gq'] = '30'
        testvar = create_test_snv(self.vardata)
        variants = {'child': {'1_100000_A_G': testvar},
                    'mum': {}, 'dad': {}}
        preinheritancefilter = PreInheritanceFiltering(variants)
        variants_per_gene = preinheritancefilter.preinheritance_filter()
        self.assertEqual(variants_per_gene, {})

        self.vardata['chrom'] = 'X'
        testvar = create_test_snv(self.vardata)
        variants = {'child': {'X_100000_A_G': testvar},
                    'mum': {}, 'dad': {}}
        preinheritancefilter = PreInheritanceFiltering(variants)
        variants_per_gene = preinheritancefilter.preinheritance_filter()
        self.assertEqual(variants_per_gene, {'123': {'X_100000_A_G': {
            'child': variants['child']['X_100000_A_G']}}})

    def test_max_ddd_af(self):
        # fail if ddd_af > 0.005
        self.vardata['ddd_af'] = '0.5'
        testvar = create_test_snv(self.vardata)
        variants = {'child': {'1_100000_A_G': testvar},
                    'mum': {}, 'dad': {}}
        preinheritancefilter = PreInheritanceFiltering(variants)
        variants_per_gene = preinheritancefilter.preinheritance_filter()
        self.assertEqual(variants_per_gene, {})

    def test_cq(self):
        # fail if consequence is not functional
        self.vardata['consequence'] = 'synonymous'
        testvar = create_test_snv(self.vardata)
        variants = {'child': {'1_100000_A_G': testvar},
                    'mum': {}, 'dad': {}}
        preinheritancefilter = PreInheritanceFiltering(variants)
        variants_per_gene = preinheritancefilter.preinheritance_filter()
        self.assertEqual(variants_per_gene, {})

        self.vardata['consequence'] = 'synonymous&missense_variant'
        testvar = create_test_snv(self.vardata)
        variants = {'child': {'1_100000_A_G': testvar},
                    'mum': {}, 'dad': {}}
        preinheritancefilter = PreInheritanceFiltering(variants)
        variants_per_gene = preinheritancefilter.preinheritance_filter()
        self.assertEqual(variants_per_gene, {'123': {'1_100000_A_G': {
            'child': variants['child']['1_100000_A_G']}}})

    def test_revel_filter(self):
        # test revel filter
        self.vardata['revel'] = '0.3'
        testvar = create_test_snv(self.vardataX)
        variants = {'child': {'1_100000_A_G': testvar},
                    'mum': {}, 'dad': {}}
        preinheritancefilter = PreInheritanceFiltering(variants)
        variants_per_gene = preinheritancefilter.preinheritance_filter()
        self.assertEqual(variants_per_gene, {})


if __name__ == '__main__':
    unittest.main()
