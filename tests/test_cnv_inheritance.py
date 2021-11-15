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
import copy

from tests.test_utils import create_test_person
from tests.test_utils import create_test_family
from tests.test_utils import create_test_cnv

from filtering.filter import CNVFiltering

class TestAutosomalInheritanceFilter(unittest.TestCase):

    def setUp(self):
        self.maxDiff = None
        self.candidate_variants = {}
        self.candidate_variants['single_variants'] = {}
        self.candidate_variants['compound_hets'] = {}
        self.child = create_test_person('fam', 'child_id', 'dad_id', 'mum_id',
                                        'XY', '2', '/vcf/path')
        self.child_female = create_test_person('fam', 'child_id', 'dad_id', 'mum_id',
                                        'XX', '2', '/vcf/path')
        self.mum = create_test_person('fam', 'mum_id', '0', '0', 'XX', '1',
                                      '/vcf/path')
        self.mum_aff = create_test_person('fam', 'mum_id', '0', '0', 'XX', '2',
                                          '/vcf/path')
        self.dad = create_test_person('fam', 'dad_id', '0', '0', 'XY', '1',
                                      '/vcf/path')
        self.dad_aff = create_test_person('fam', 'dad_id', '0', '0', 'XY', '2',
                                          '/vcf/path')
        self.family_both_unaff = create_test_family(self.child, self.mum,
                                                    self.dad)
        self.family_both_unaff_female = create_test_family(self.child_female, self.mum,
                                                    self.dad)
        self.family_mum_aff = create_test_family(self.child, self.mum_aff,
                                                 self.dad)
        self.family_dad_aff = create_test_family(self.child, self.mum,
                                                 self.dad_aff)
        self.family_both_aff = create_test_family(self.child, self.mum_aff,
                                                  self.dad_aff)

        self.genes_biallelic = {
            '1234': {'chr': '1', 'start': '10971836', 'end': '10984446',
                     'symbol': 'MECP2', 'status': {'Probable DD gene'},
                     'mode': {'Biallelic'},
                     'mechanism': {'Loss of function'}}
        }
        self.genes_monoallelic = copy.deepcopy(self.genes_biallelic)
        self.genes_monoallelic['1234']['mode'] = {'Monoallelic'}
        self.genes_hemizygous = copy.deepcopy(self.genes_biallelic)
        self.genes_hemizygous['1234']['mode'] = {'Hemizygous'}
        self.genes_hemizygous['1234']['chr'] = {'X'}
        self.genes_X_linked_dominant = copy.deepcopy(self.genes_hemizygous)
        self.genes_X_linked_dominant['1234']['mode'] = {'X-linked dominant'}
        self.genes_increased_dose = copy.deepcopy(self.genes_hemizygous)
        self.genes_increased_dose['1234']['mechanism'] = {'Increased gene dosage'}

        self.cn0vardata_bi = {'chrom': '1', 'pos': '10971936', 'ref': 'T', 'alt': '<DEL>',
             'consequence': 'transcript_ablation', 'ensg': 'ENSG01234',
             'symbol': 'MECP1', 'feature': 'ENST01234', 'canonical': 'YES',
             'mane': 'NM01234', 'hgnc_id': 'HGNC:123', 'cnv_type': 'DEL', 'cnv_filter': 'Pass',
             'hgnc_id_all': 'HGNC:1234|HGNC:2|HGNC:3', 'cnv_end':'11071936',
             'symbol_all': 'MECP1|MECP2|MECP3', 'sex': 'XY', 'cn': '0',
             'cnv_inh':'biparental_inh', 'cnv_length':'5000', 'gt': '.'}
        self.cn0vardata_denovo = copy.deepcopy(self.cn0vardata_bi)
        self.cn0vardata_denovo['cnv_inh']= 'not_inherited'
        self.cn0vardata_pat = copy.deepcopy(self.cn0vardata_bi)
        self.cn0vardata_pat['cnv_inh']= 'paternal_inh'
        self.cn0vardata_mat = copy.deepcopy(self.cn0vardata_bi)
        self.cn0vardata_mat['cnv_inh']= 'maternal_inh'
        self.cn1vardata_bi = copy.deepcopy(self.cn0vardata_bi)
        self.cn1vardata_bi['cn'] = '1'
        self.cn1vardata_mat = copy.deepcopy(self.cn0vardata_mat)
        self.cn1vardata_mat['cn'] = '1'
        self.xvardata_mat = copy.deepcopy(self.cn0vardata_mat)
        self.xvardata_mat['chrom'] = 'X'
        self.cn0vardata_long = copy.deepcopy(self.cn0vardata_bi)
        self.cn0vardata_long['cnv_length'] = '10000001'
        self.dupvardata_surround = copy.deepcopy(self.cn0vardata_bi)
        self.dupvardata_surround['alt'] = '<DUP>'
        self.dupvardata_surround['cn'] = '3'
        self.dupvardata_surround['pos'] = '10970000'
        self.dupvardata_surround['cnv_end'] = '10990000'
        self.dupvardata_surround['cnv_length'] = '20000'
        self.dupvardata_surround['cnv_type'] = 'DUP'
        self.cn3vardata_denovo = copy.deepcopy(self.cn0vardata_denovo)
        self.cn3vardata_denovo['cnv_type'] = 'DUP'
        self.cn3vardata_denovo['alt'] = '<DUP>'
        self.cn3vardata_denovo['cn'] = '3'
        self.cn3vardata_mat = copy.deepcopy(self.cn3vardata_denovo)
        self.cn3vardata_mat['cnv_inh'] = 'maternal_inh'
        self.intragenic_dup = copy.deepcopy(self.dupvardata_surround)
        self.intragenic_dup['cnv_end'] = '10980000'
        self.intragenic_dup['cnv_length'] = '10000'


    def test_inh_matches_parent_aff_status(self):
        #pass if paternal and dad affected
        testcnv1 = create_test_cnv(self.cn0vardata_pat)
        testvars1 = {'child':{'1_10971936_A_DEL':testcnv1}}
        cnvfilter_dadaff = CNVFiltering(testvars1, self.family_dad_aff, self.genes_monoallelic, None,
                                 None, self.candidate_variants)
        cnvfilter_dadaff.cnv_filter()
        self.assertEqual(cnvfilter_dadaff.inhmatch, True)
        #fail if paternal and dad unaffected
        cnvfilter_dadunaff = CNVFiltering(testvars1, self.family_both_unaff,
                                        self.genes_monoallelic, None,
                                        None, self.candidate_variants)
        cnvfilter_dadunaff.cnv_filter()
        self.assertEqual(cnvfilter_dadunaff.inhmatch, False)
        # pass if maternal and mum affected
        testcnv2 = create_test_cnv(self.cn0vardata_mat)
        testvars2 = {'child': {'1_10971936_A_DEL': testcnv2}}
        cnvfilter_mumaff = CNVFiltering(testvars2, self.family_mum_aff,
                                        self.genes_monoallelic, None,
                                        None, self.candidate_variants)
        cnvfilter_mumaff.cnv_filter()
        self.assertEqual(cnvfilter_mumaff.inhmatch, True)
        # fail if maternal and mum unaffected
        cnvfilter_mumunaff = CNVFiltering(testvars2, self.family_both_unaff,
                                          self.genes_monoallelic, None,
                                          None, self.candidate_variants)
        cnvfilter_mumunaff.cnv_filter()
        self.assertEqual(cnvfilter_mumunaff.inhmatch, False)
        #pass if biparental and cn = 0
        testcnv3 = create_test_cnv(self.cn0vardata_bi)
        testvars3 = {'child': {'1_10971936_A_DEL': testcnv3}}
        cnvfilter_bothunaff = CNVFiltering(testvars3, self.family_both_unaff,
                                          self.genes_monoallelic, None,
                                          None, self.candidate_variants)
        cnvfilter_bothunaff.cnv_filter()
        self.assertEqual(cnvfilter_bothunaff.inhmatch, True)
        #pass if biparental and either/both parents affected
        testcnv4 = create_test_cnv(self.cn1vardata_bi)
        testvars4 = {'child': {'1_10971936_A_DEL': testcnv4}}
        cnvfilter_dadaff2 = CNVFiltering(testvars4, self.family_dad_aff,
                                        self.genes_monoallelic, None,
                                        None, self.candidate_variants)
        cnvfilter_mumaff2 = CNVFiltering(testvars4, self.family_mum_aff,
                                         self.genes_monoallelic, None,
                                         None, self.candidate_variants)
        cnvfilter_bothaff2 = CNVFiltering(testvars4, self.family_both_aff,
                                         self.genes_monoallelic, None,
                                         None, self.candidate_variants)
        cnvfilter_dadaff2.cnv_filter()
        cnvfilter_mumaff2.cnv_filter()
        cnvfilter_bothaff2.cnv_filter()
        self.assertEqual(cnvfilter_dadaff2.inhmatch, True)
        self.assertEqual(cnvfilter_mumaff2.inhmatch, True)
        self.assertEqual(cnvfilter_bothaff2.inhmatch, True)
        #fail if biparental, cn = 1 and parents both unaff
        cnvfilter_bothunaff2 = CNVFiltering(testvars4, self.family_both_unaff,
                                          self.genes_monoallelic, None,
                                          None, self.candidate_variants)
        cnvfilter_bothunaff2.cnv_filter()
        self.assertEqual(cnvfilter_bothunaff2.inhmatch, False)
        #pass if male, X, maternal inh and mum unaffected
        testcnv5 = create_test_cnv(self.xvardata_mat)
        testvars5 = {'child': {'1_10971936_A_DEL': testcnv5}}
        cnvfilter_unaff3 = CNVFiltering(testvars5, self.family_both_unaff,
                                         self.genes_monoallelic, None,
                                         None, self.candidate_variants)
        cnvfilter_unaff3.cnv_filter()
        self.assertEqual(cnvfilter_unaff3.inhmatch, True)
        # pass if male, X, maternal inh and mum affected(?)
        cnvfilter_mumaff3 = CNVFiltering(testvars5, self.family_mum_aff,
                                        self.genes_monoallelic, None,
                                        None, self.candidate_variants)
        cnvfilter_mumaff3.cnv_filter()
        self.assertEqual(cnvfilter_mumaff3.inhmatch, True)

    def test_non_ddg2p_filter(self):
        #needs to be tested on a variant/family that firstly passes inheritance_filter
        #cnv < 1M fail
        testcnvshort = create_test_cnv(self.cn0vardata_bi)
        testvarsshort = {'child': {'1_10971936_A_DEL': testcnvshort}}
        cnvfiltershort = CNVFiltering(testvarsshort, self.family_dad_aff,
                                        self.genes_monoallelic, None,
                                        None, self.candidate_variants)
        cnvfiltershort.cnv_filter()
        self.assertEqual(cnvfiltershort.passnonddg2p, False)
        #cnv > 1M pass
        testcnvlong = create_test_cnv(self.cn0vardata_long)
        testvarslong = {'child': {'1_10971936_A_DEL': testcnvlong}}
        cnvfilterlong = CNVFiltering(testvarslong, self.family_dad_aff,
                                      self.genes_monoallelic, None,
                                      None, self.candidate_variants)
        cnvfilterlong.cnv_filter()
        self.assertEqual(cnvfilterlong.passnonddg2p, True)

    def test_ddg2p_filter(self):
        #need to test with CNVs that pass inheritance filter but are <1M
        # - fail duplications which completely surround monoallelic, hemizygous and
        #   x-linked dominant genes with loss of function mechanism
        testcnvdup1 = create_test_cnv(self.dupvardata_surround)
        testvarsdup1 = {'child': {'1_10971936_A_DEL': testcnvdup1}}
        cnvfilterdup1 = CNVFiltering(testvarsdup1, self.family_dad_aff,
                                      self.genes_monoallelic, None,
                                      None, self.candidate_variants)
        cnvfilterdup1.cnv_filter()
        self.assertEqual(cnvfilterdup1.passddg2p, False)

        # - Biallelic gene pass if copy number (CN) = 0 and mechanism in
        #   Uncertain", "Loss of function", "Dominant negative"
        testcnv0 = create_test_cnv(self.cn0vardata_bi)
        testvars0 = {'child': {'1_10971936_A_DEL': testcnv0}}
        cnvfilter0 = CNVFiltering(testvars0, self.family_dad_aff,
                                     self.genes_biallelic, None,
                                     None, self.candidate_variants)
        cnvfilter0.cnv_filter()
        self.assertEqual(cnvfilter0.passddg2p, True)

        # - Monoallelic, X-linked dominant or Hemizygous in male pass if CN=0,
        #   1 or 3 and any mechanism
        testcnv2 = create_test_cnv(self.cn0vardata_denovo)
        testvars2 = {'child': {'1_10971936_A_DEL': testcnv2}}
        cnvfilter2 = CNVFiltering(testvars2, self.family_both_unaff,
                                  self.genes_hemizygous, None,
                                  None, self.candidate_variants)
        cnvfilter2.cnv_filter()
        self.assertEqual(cnvfilter2.passddg2p, True)

        cnvfilter3 = CNVFiltering(testvars2, self.family_both_unaff_female,
                                  self.genes_hemizygous, None,
                                  None, self.candidate_variants)
        cnvfilter3.cnv_filter()
        self.assertEqual(cnvfilter3.passddg2p, False)

        # - Hemizygous in female pass if CN=3 and mechanism = "Increased gene
        #   dosage"
        testcnv4 = create_test_cnv(self.cn3vardata_denovo)
        testvars4 = {'child': {'1_10971936_A_DUP': testcnv4}}
        cnvfilter4 = CNVFiltering(testvars4, self.family_both_unaff_female,
                                  self.genes_hemizygous, None,
                                  None, self.candidate_variants)
        cnvfilter4.cnv_filter()
        self.assertEqual(cnvfilter4.passddg2p, False)

        cnvfilter5 = CNVFiltering(testvars4, self.family_both_unaff_female,
                                  self.genes_increased_dose, None,
                                  None, self.candidate_variants)
        cnvfilter5.cnv_filter()
        self.assertEqual(cnvfilter5.passddg2p, True)
        # - Pass intragenic DUP in monoallelic or X-linked dominant gene with
        #   loss of function mechanism and any part of the gene is outside of the CNV boundary
        testcnv6 = create_test_cnv(self.intragenic_dup)
        testvars6 = {'child': {'1_10971936_A_DUP': testcnv6}}
        cnvfilter6 = CNVFiltering(testvars6, self.family_both_aff,
                                  self.genes_X_linked_dominant, None,
                                  None, self.candidate_variants)
        cnvfilter6.cnv_filter()
        self.assertEqual(cnvfilter6.passddg2p, True)

        testcnv7 = create_test_cnv(self.dupvardata_surround)
        testvars7 = {'child': {'1_10971936_A_DUP': testcnv7}}
        cnvfilter7 = CNVFiltering(testvars7, self.family_both_aff,
                                  self.genes_X_linked_dominant, None,
                                  None, self.candidate_variants)
        cnvfilter7.cnv_filter()
        self.assertEqual(cnvfilter7.passddg2p, True)

    def test_candidate_compound_het_filter(self):
        # add var to candidate compound hets if:
        # cn = 1 or 3 and biallelic and DDG2P gene
        testcnv1 = create_test_cnv(self.cn1vardata_mat)
        testvars1 = {'child': {'1_10971936_A_DEL': testcnv1}}
        cnvfilter1 = CNVFiltering(testvars1, self.family_both_unaff,
                                  self.genes_biallelic, None,
                                  None, self.candidate_variants)
        cnvfilter1.cnv_filter()
        self.assertEqual(cnvfilter1.posscomphet, True)

        testcnv2 = create_test_cnv(self.cn3vardata_mat)
        testvars2 = {'child': {'1_10971936_A_DEL': testcnv2}}
        cnvfilter2 = CNVFiltering(testvars2, self.family_both_unaff,
                                  self.genes_biallelic, None,
                                  None, self.candidate_variants)
        cnvfilter2.cnv_filter()
        self.assertEqual(cnvfilter2.posscomphet, True)

        cnvfilter3 = CNVFiltering(testvars2, self.family_both_unaff,
                                  self.genes_monoallelic, None,
                                  None, self.candidate_variants)
        cnvfilter3.cnv_filter()
        self.assertEqual(cnvfilter3.posscomphet, False)
        # or if cn = 1, male, hemizygous and DDG2P gene
        cnvfilter4 = CNVFiltering(testvars1, self.family_both_unaff,
                                  self.genes_hemizygous, None,
                                  None, self.candidate_variants)
        cnvfilter4.cnv_filter()
        self.assertEqual(cnvfilter4.posscomphet, True)

        cnvfilter5 = CNVFiltering(testvars2, self.family_both_unaff,
                                  self.genes_hemizygous, None,
                                  None, self.candidate_variants)
        cnvfilter5.cnv_filter()
        self.assertEqual(cnvfilter5.posscomphet, False)

        cnvfilter6 = CNVFiltering(testvars1, self.family_both_unaff_female,
                                  self.genes_hemizygous, None,
                                  None, self.candidate_variants)
        cnvfilter6.cnv_filter()
        self.assertEqual(cnvfilter6.posscomphet, False)


if __name__ == '__main__':
    unittest.main()