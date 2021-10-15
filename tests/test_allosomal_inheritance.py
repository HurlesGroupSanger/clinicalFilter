import unittest
import copy

from tests.test_utils import create_test_person
from tests.test_utils import create_test_family
from tests.test_utils import create_test_variants_per_gene

from filtering.inheritance_filtering import InheritanceFiltering
from filtering.inheritance_report import InheritanceReport


class TestAllosomalInheritanceFilter(unittest.TestCase):
    def setUp(self):
        self.maxDiff = None
        self.candidate_variants = {}
        self.candidate_variants['single_variants'] = {}
        self.candidate_variants['compound_hets'] = {}
        self.inhreport = None
        self.inhreport = InheritanceReport()
        self.hetvardata = {'chrom': 'X', 'pos': '1097183', 'ref': 'A',
                           'alt': 'GG',
                           'consequence': 'start_lost', 'ensg': 'ensg',
                           'symbol': 'DDX3X', 'feature': 'feature',
                           'canonical': 'YES', 'mane': 'MANE',
                           'hgnc_id': '1234', 'ddd_father_af': '.',
                           'max_af': '0', 'max_af_pops': '.', 'ddd_af': '0',
                           'revel': '.', 'polyphen': '.', 'hgvsc': '.',
                           'hgvsp': '.', 'sex': 'XY', 'dnm': '.',
                           'gt': '0/1', 'gq': '50',
                           'pid': '.', 'protein_position': '123', 'ad': '4,4',
                           'ac_XX': '0',
                           'ac_XY': '0', 'an_XX': '0', 'an_XY': '0',
                           'nhomalt_XX': '0', 'nhomalt_XY': '0'
                           }
        self.homaltvardata = self.hetvardata.copy()
        self.homaltvardata['gt'] = '1/1'
        self.homrefvardata = self.hetvardata.copy()
        self.homrefvardata['gt'] = '0/0'

        self.variants_100 = {'child': {'X_1097183_A_GG': self.hetvardata},
                             'mum': {'X_1097183_A_GG': self.homrefvardata},
                             'dad': {'X_1097183_A_GG': self.homrefvardata}}
        self.variants_110 = {'child': {'X_1097183_A_GG': self.hetvardata},
                             'mum': {'X_1097183_A_GG': self.hetvardata},
                             'dad': {'X_1097183_A_GG': self.homrefvardata}}
        self.variants_102 = {'child': {'X_1097183_A_GG': self.hetvardata},
                             'mum': {'X_1097183_A_GG': self.homrefvardata},
                             'dad': {'X_1097183_A_GG': self.homaltvardata}}
        self.variants_120 = {'child': {'X_1097183_A_GG': self.hetvardata},
                             'mum': {'X_1097183_A_GG': self.homaltvardata},
                             'dad': {'X_1097183_A_GG': self.homrefvardata}}
        self.variants_112 = {'child': {'X_1097183_A_GG': self.hetvardata},
                             'mum': {'X_1097183_A_GG': self.hetvardata},
                             'dad': {'X_1097183_A_GG': self.homaltvardata}}
        self.variants_122 = {'child': {'X_1097183_A_GG': self.hetvardata},
                             'mum': {'X_1097183_A_GG': self.homaltvardata},
                             'dad': {'X_1097183_A_GG': self.homaltvardata}}
        self.variants_200 = {'child': {'X_1097183_A_GG': self.homaltvardata},
                             'mum': {'X_1097183_A_GG': self.homrefvardata},
                             'dad': {'X_1097183_A_GG': self.homrefvardata}}
        self.variants_210 = {'child': {'X_1097183_A_GG': self.homaltvardata},
                             'mum': {'X_1097183_A_GG': self.hetvardata},
                             'dad': {'X_1097183_A_GG': self.homrefvardata}}
        self.variants_202 = {'child': {'X_1097183_A_GG': self.homaltvardata},
                             'mum': {'X_1097183_A_GG': self.homrefvardata},
                             'dad': {'X_1097183_A_GG': self.homaltvardata}}
        self.variants_220 = {'child': {'X_1097183_A_GG': self.homaltvardata},
                             'mum': {'X_1097183_A_GG': self.homaltvardata},
                             'dad': {'X_1097183_A_GG': self.homrefvardata}}
        self.variants_212 = {'child': {'X_1097183_A_GG': self.homaltvardata},
                             'mum': {'X_1097183_A_GG': self.hetvardata},
                             'dad': {'X_1097183_A_GG': self.homaltvardata}}
        self.variants_222 = {'child': {'X_1097183_A_GG': self.homaltvardata},
                             'mum': {'X_1097183_A_GG': self.homaltvardata},
                             'dad': {'X_1097183_A_GG': self.homaltvardata}}
        self.variants_1 = {'child': {'X_1097183_A_GG': self.hetvardata}}
        self.variants_2 = {'child': {'X_1097183_A_GG': self.homaltvardata}}

        self.hetvardataY = {'chrom': 'Y', 'pos': '1097183', 'ref': 'A',
                           'alt': 'GG',
                           'consequence': 'start_lost', 'ensg': 'ensg',
                           'symbol': 'SRY', 'feature': 'feature',
                           'canonical': 'YES', 'mane': 'MANE',
                           'hgnc_id': '1234', 'ddd_father_af': '.',
                           'max_af': '0', 'max_af_pops': '.', 'ddd_af': '0',
                           'revel': '.', 'polyphen': '.', 'hgvsc': '.',
                           'hgvsp': '.', 'sex': 'XY', 'dnm': '.',
                           'gt': '0/1', 'gq': '50',
                           'pid': '.', 'protein_position': '123', 'ad': '4,4',
                            'ac_XX': '0',
                            'ac_XY': '0', 'an_XX': '0', 'an_XY': '0',
                            'nhomalt_XX': '0', 'nhomalt_XY': '0'
                            }
        self.homaltvardataY = self.hetvardataY.copy()
        self.homaltvardataY['gt'] = '1/1'
        self.homrefvardataY = self.hetvardataY.copy()
        self.homrefvardataY['gt'] = '0/0'

        self.variants_100Y = {'child': {'Y_1097183_A_GG': self.hetvardataY},
                             'mum': {'Y_1097183_A_GG': self.homrefvardataY},
                             'dad': {'Y_1097183_A_GG': self.homrefvardataY}}
        self.variants_110Y = {'child': {'Y_1097183_A_GG': self.hetvardataY},
                             'mum': {'Y_1097183_A_GG': self.hetvardataY},
                             'dad': {'Y_1097183_A_GG': self.homrefvardataY}}
        self.variants_102Y = {'child': {'Y_1097183_A_GG': self.hetvardataY},
                             'mum': {'Y_1097183_A_GG': self.homrefvardataY},
                             'dad': {'Y_1097183_A_GG': self.homaltvardataY}}
        self.variants_120Y = {'child': {'Y_1097183_A_GG': self.hetvardataY},
                             'mum': {'Y_1097183_A_GG': self.homaltvardataY},
                             'dad': {'Y_1097183_A_GG': self.homrefvardataY}}
        self.variants_112Y = {'child': {'Y_1097183_A_GG': self.hetvardataY},
                             'mum': {'Y_1097183_A_GG': self.hetvardataY},
                             'dad': {'Y_1097183_A_GG': self.homaltvardataY}}
        self.variants_122Y = {'child': {'Y_1097183_A_GG': self.hetvardataY},
                             'mum': {'Y_1097183_A_GG': self.homaltvardataY},
                             'dad': {'Y_1097183_A_GG': self.homaltvardataY}}
        self.variants_200Y = {'child': {'Y_1097183_A_GG': self.homaltvardataY},
                             'mum': {'Y_1097183_A_GG': self.homrefvardataY},
                             'dad': {'Y_1097183_A_GG': self.homrefvardataY}}
        self.variants_210Y = {'child': {'Y_1097183_A_GG': self.homaltvardataY},
                             'mum': {'Y_1097183_A_GG': self.hetvardataY},
                             'dad': {'Y_1097183_A_GG': self.homrefvardataY}}
        self.variants_202Y = {'child': {'Y_1097183_A_GG': self.homaltvardataY},
                             'mum': {'Y_1097183_A_GG': self.homrefvardataY},
                             'dad': {'Y_1097183_A_GG': self.homaltvardataY}}
        self.variants_220Y = {'child': {'Y_1097183_A_GG': self.homaltvardataY},
                             'mum': {'Y_1097183_A_GG': self.homaltvardataY},
                             'dad': {'Y_1097183_A_GG': self.homrefvardataY}}
        self.variants_212Y = {'child': {'Y_1097183_A_GG': self.homaltvardataY},
                             'mum': {'Y_1097183_A_GG': self.hetvardataY},
                             'dad': {'Y_1097183_A_GG': self.homaltvardataY}}
        self.variants_222Y = {'child': {'Y_1097183_A_GG': self.homaltvardataY},
                             'mum': {'Y_1097183_A_GG': self.homaltvardataY},
                             'dad': {'Y_1097183_A_GG': self.homaltvardataY}}
        self.variants_1Y = {'child': {'Y_1097183_A_GG': self.hetvardataY}}
        self.variants_2Y = {'child': {'Y_1097183_A_GG': self.homaltvardataY}}

        self.genes_hemizygous = {
        '1234': {'chr': 'X', 'start': '1097083', 'end': '1098083',
                 'symbol': 'DDX3X', 'status': {'Probable DD gene'},
                 'mode': {'Hemizygous'},
                 'mechanism': {'Loss of function'}}
        }
        self.genes_X_linked_dominant = copy.deepcopy(self.genes_hemizygous)
        self.genes_X_linked_dominant['1234']['mode'] = {'X-linked dominant'}
        self.genes_X_linked_over_dominance = copy.deepcopy(self.genes_hemizygous)
        self.genes_X_linked_over_dominance['1234']['mode'] = {'X-linked over-dominance'}
        self.genes_monoallelic_Y_hemizygous = copy.deepcopy(self.genes_hemizygous)
        self.genes_monoallelic_Y_hemizygous['1234']['mode'] = {'monoallelic_Y_hem'}
        self.genes_monoallelic_Y_hemizygous['1234']['chr'] = 'Y'
        self.genes_monoallelic_Y_hemizygous['1234']['symbol'] = 'SRY'

        self.child_XX = create_test_person('fam', 'child_id', 'dad_id', 'mum_id',
                                        'XX', '2', '/vcf/path')
        self.child_XXY = create_test_person('fam', 'child_id', 'dad_id',
                                           'mum_id',
                                           'XXY', '2', '/vcf/path')
        self.child_XY = create_test_person('fam', 'child_id', 'dad_id',
                                           'mum_id',
                                           'XY', '2', '/vcf/path')
        self.child_XX2 = create_test_person('fam', 'child_id', '0',
                                           '0',
                                           'XX', '2', '/vcf/path')
        self.child_XY2 = create_test_person('fam', 'child_id', '0',
                                            '0',
                                            'XY', '2', '/vcf/path')
        self.mum = create_test_person('fam', 'mum_id', '0', '0', 'XX', '1',
                                      '/vcf/path')
        self.mum_aff = create_test_person('fam', 'mum_id', '0', '0', 'XX', '2',
                                          '/vcf/path')
        self.dad = create_test_person('fam', 'dad_id', '0', '0', 'XY', '1',
                                      '/vcf/path')
        self.dad_aff = create_test_person('fam', 'dad_id', '0', '0', 'XY', '2',
                                          '/vcf/path')

        self.family_XX_both_unaff = create_test_family(self.child_XX, self.mum,
                                                    self.dad)
        self.family_XX_mum_aff = create_test_family(self.child_XX, self.mum_aff,
                                                 self.dad)
        self.family_XX_dad_aff = create_test_family(self.child_XX, self.mum,
                                                 self.dad_aff)
        self.family_XX_both_aff = create_test_family(self.child_XX, self.mum_aff,
                                                  self.dad_aff)
        self.family_XX_no_parents = create_test_family(self.child_XX2, None, None)
        self.family_XY_both_unaff = create_test_family(self.child_XY, self.mum,
                                                       self.dad)
        self.family_XY_mum_aff = create_test_family(self.child_XY, self.mum_aff,
                                                    self.dad)
        self.family_XY_dad_aff = create_test_family(self.child_XY, self.mum,
                                                    self.dad_aff)
        self.family_XY_both_aff = create_test_family(self.child_XY, self.mum_aff,
                                                     self.dad_aff)
        self.family_XY_no_parents = create_test_family(self.child_XY2, None,
                                                       None)

    def test_allosomal_no_parents(self):
        # in X all should pass - only testing one for now
        variants_per_gene_1 = create_test_variants_per_gene(self.variants_1,
                                                               self.family_XX_no_parents)
        inheritancefilter_1 = InheritanceFiltering(variants_per_gene_1,
                                                      self.family_XX_no_parents,
                                                      self.genes_hemizygous,
                                                      None, None, self.candidate_variants,
                                                      self.inhreport)
        inheritancefilter_1.inheritance_filter_genes()
        test_candidate_variants_1 = {'single_variants': {
            'X_1097183_A_GG': {
                'mode': {'hemizygous'},
                'variant':
                    variants_per_gene_1[
                        '1234'][
                        'X_1097183_A_GG'][
                        'child'],
                'hgncid': '1234'}}, 'compound_hets': {}}

        self.assertEqual(inheritancefilter_1.candidate_variants,
                         test_candidate_variants_1)

        # in Y should pass if 1/1

        variants_per_gene_2 = create_test_variants_per_gene(self.variants_2Y,
                                                            self.family_XY_no_parents)
        inheritancefilter_2 = InheritanceFiltering(variants_per_gene_2,
                                                   self.family_XY_no_parents,
                                                   self.genes_monoallelic_Y_hemizygous, None,
                                                   None, self.candidate_variants,
                                                      self.inhreport)
        inheritancefilter_2.inheritance_filter_genes()
        test_candidate_variants_2 = {'single_variants': {
            'Y_1097183_A_GG': {
                'mode': {'monoallelic_Y_hemizygous'},
                'variant':
                    variants_per_gene_2[
                        '1234'][
                        'Y_1097183_A_GG'][
                        'child'],
                'hgncid': '1234'}}, 'compound_hets': {}}

        self.assertEqual(inheritancefilter_2.candidate_variants,
                         test_candidate_variants_2)

        # Y fail if 0/1
        variants_per_gene_1 = create_test_variants_per_gene(self.variants_1Y,
                                                            self.family_XY_no_parents)
        inheritancefilter_1 = InheritanceFiltering(variants_per_gene_1,
                                                   self.family_XY_no_parents,
                                                   self.genes_monoallelic_Y_hemizygous,
                                                   None,
                                                   None, self.candidate_variants,
                                                      self.inhreport)
        inheritancefilter_1.inheritance_filter_genes()
        test_candidate_variants_1 = {'single_variants': {}, 'compound_hets': {}}

        self.assertEqual(inheritancefilter_1.candidate_variants,
                         test_candidate_variants_1)

    def test_gn_hemizygous_gt_homozygous_parents_filter(self):
        # both aff, mum 0/0, dad 0/0 fail
        variants_per_gene_200 = create_test_variants_per_gene(self.variants_200,
                                                            self.family_XX_both_aff)
        inheritancefilter_200 = InheritanceFiltering(variants_per_gene_200,
                                                   self.family_XX_both_aff,
                                                   self.genes_hemizygous, None,
                                                   None, self.candidate_variants,
                                                      self.inhreport)
        inheritancefilter_200.inheritance_filter_genes()
        test_candidate_variants_200 = {'single_variants': {}, 'compound_hets': {}}

        self.assertEqual(inheritancefilter_200.candidate_variants,
                         test_candidate_variants_200)
        # both aff, mum 0/1, dad 0/0 fail
        variants_per_gene_210 = create_test_variants_per_gene(self.variants_210,
                                                            self.family_XX_both_aff)
        inheritancefilter_210 = InheritanceFiltering(variants_per_gene_210,
                                                   self.family_XX_both_aff,
                                                   self.genes_hemizygous, None,
                                                   None, self.candidate_variants,
                                                      self.inhreport)
        inheritancefilter_210.inheritance_filter_genes()
        test_candidate_variants_210 = {'single_variants': {}, 'compound_hets': {}}

        self.assertEqual(inheritancefilter_210.candidate_variants,
                         test_candidate_variants_210)
        # both aff, mum 1/1, dad 0/0 fail
        variants_per_gene_220 = create_test_variants_per_gene(self.variants_220,
                                                            self.family_XX_both_aff)
        inheritancefilter_220 = InheritanceFiltering(variants_per_gene_220,
                                                   self.family_XX_both_aff,
                                                   self.genes_hemizygous, None,
                                                   None, self.candidate_variants,
                                                      self.inhreport)
        inheritancefilter_220.inheritance_filter_genes()
        test_candidate_variants_220 = {'single_variants': {}, 'compound_hets': {}}

        self.assertEqual(inheritancefilter_220.candidate_variants,
                         test_candidate_variants_220)
        # both aff, mum 0/0, dad 1/1 fail
        variants_per_gene_202 = create_test_variants_per_gene(self.variants_202,
                                                            self.family_XX_both_aff)
        inheritancefilter_202 = InheritanceFiltering(variants_per_gene_202,
                                                   self.family_XX_both_aff,
                                                   self.genes_hemizygous, None,
                                                   None, self.candidate_variants,
                                                      self.inhreport)
        inheritancefilter_202.inheritance_filter_genes()
        test_candidate_variants_202 = {'single_variants': {}, 'compound_hets': {}}

        self.assertEqual(inheritancefilter_202.candidate_variants,
                         test_candidate_variants_202)
        # both aff, mum 0/1, dad 1/1 fail
        variants_per_gene_212 = create_test_variants_per_gene(self.variants_212,
                                                            self.family_XX_both_aff)
        inheritancefilter_212 = InheritanceFiltering(variants_per_gene_212,
                                                   self.family_XX_both_aff,
                                                   self.genes_hemizygous, None,
                                                   None, self.candidate_variants,
                                                      self.inhreport)
        inheritancefilter_212.inheritance_filter_genes()
        test_candidate_variants_212 = {'single_variants': {}, 'compound_hets': {}}

        self.assertEqual(inheritancefilter_212.candidate_variants,
                         test_candidate_variants_212)
        # both aff, mum 1/1, dad 1/1 fail
        variants_per_gene_222 = create_test_variants_per_gene(self.variants_222,
                                                            self.family_XX_both_aff)
        inheritancefilter_222 = InheritanceFiltering(variants_per_gene_222,
                                                   self.family_XX_both_aff,
                                                   self.genes_hemizygous, None,
                                                   None, self.candidate_variants,
                                                      self.inhreport)
        inheritancefilter_222.inheritance_filter_genes()
        test_candidate_variants_222 = {'single_variants': {}, 'compound_hets': {}}

        self.assertEqual(inheritancefilter_222.candidate_variants,
                         test_candidate_variants_222)

        # mum aff, mum 0/0, dad 0/0 fail
        variants_per_gene_200 = create_test_variants_per_gene(self.variants_200,
                                                              self.family_XX_mum_aff)
        inheritancefilter_200 = InheritanceFiltering(variants_per_gene_200,
                                                     self.family_XX_mum_aff,
                                                     self.genes_hemizygous,
                                                     None,
                                                     None, self.candidate_variants,
                                                      self.inhreport)
        inheritancefilter_200.inheritance_filter_genes()
        test_candidate_variants_200 = {'single_variants': {},
                                       'compound_hets': {}}

        self.assertEqual(inheritancefilter_200.candidate_variants,
                         test_candidate_variants_200)
        # mum aff, mum 0/1, dad 0/0 fail
        variants_per_gene_210 = create_test_variants_per_gene(self.variants_210,
                                                              self.family_XX_mum_aff)
        inheritancefilter_210 = InheritanceFiltering(variants_per_gene_210,
                                                     self.family_XX_mum_aff,
                                                     self.genes_hemizygous,
                                                     None,
                                                     None, self.candidate_variants,
                                                      self.inhreport)
        inheritancefilter_210.inheritance_filter_genes()
        test_candidate_variants_210 = {'single_variants': {},
                                       'compound_hets': {}}

        self.assertEqual(inheritancefilter_210.candidate_variants,
                         test_candidate_variants_210)
        # mum aff, mum 1/1, dad 0/0 fail
        variants_per_gene_220 = create_test_variants_per_gene(self.variants_220,
                                                              self.family_XX_mum_aff)
        inheritancefilter_220 = InheritanceFiltering(variants_per_gene_220,
                                                     self.family_XX_mum_aff,
                                                     self.genes_hemizygous,
                                                     None,
                                                     None, self.candidate_variants,
                                                      self.inhreport)
        inheritancefilter_220.inheritance_filter_genes()
        test_candidate_variants_220 = {'single_variants': {},
                                       'compound_hets': {}}

        self.assertEqual(inheritancefilter_220.candidate_variants,
                         test_candidate_variants_220)
        # mum aff, mum 0/0, dad 1/1 fail
        variants_per_gene_202 = create_test_variants_per_gene(self.variants_202,
                                                              self.family_XX_mum_aff)
        inheritancefilter_202 = InheritanceFiltering(variants_per_gene_202,
                                                     self.family_XX_mum_aff,
                                                     self.genes_hemizygous,
                                                     None,
                                                     None, self.candidate_variants,
                                                      self.inhreport)
        inheritancefilter_202.inheritance_filter_genes()
        test_candidate_variants_202 = {'single_variants': {},
                                       'compound_hets': {}}

        self.assertEqual(inheritancefilter_202.candidate_variants,
                         test_candidate_variants_202)
        # mum aff, mum 0/1, dad 1/1 fail
        variants_per_gene_212 = create_test_variants_per_gene(self.variants_212,
                                                              self.family_XX_mum_aff)
        inheritancefilter_212 = InheritanceFiltering(variants_per_gene_212,
                                                     self.family_XX_mum_aff,
                                                     self.genes_hemizygous,
                                                     None,
                                                     None, self.candidate_variants,
                                                      self.inhreport)
        inheritancefilter_212.inheritance_filter_genes()
        test_candidate_variants_212 = {'single_variants': {},
                                       'compound_hets': {}}

        self.assertEqual(inheritancefilter_212.candidate_variants,
                         test_candidate_variants_212)
        # mum aff, mum 1/1, dad 1/1 fail
        variants_per_gene_222 = create_test_variants_per_gene(self.variants_222,
                                                              self.family_XX_mum_aff)
        inheritancefilter_222 = InheritanceFiltering(variants_per_gene_222,
                                                     self.family_XX_mum_aff,
                                                     self.genes_hemizygous,
                                                     None,
                                                     None, self.candidate_variants,
                                                      self.inhreport)
        inheritancefilter_222.inheritance_filter_genes()
        test_candidate_variants_222 = {'single_variants': {},
                                       'compound_hets': {}}

        self.assertEqual(inheritancefilter_222.candidate_variants,
                         test_candidate_variants_222)

        # dad aff, mum 0/0, dad 0/0 fail
        variants_per_gene_200 = create_test_variants_per_gene(self.variants_200,
                                                              self.family_XX_dad_aff)
        inheritancefilter_200 = InheritanceFiltering(variants_per_gene_200,
                                                     self.family_XX_dad_aff,
                                                     self.genes_hemizygous,
                                                     None,
                                                     None, self.candidate_variants,
                                                      self.inhreport)
        inheritancefilter_200.inheritance_filter_genes()
        test_candidate_variants_200 = {'single_variants': {},
                                       'compound_hets': {}}

        self.assertEqual(inheritancefilter_200.candidate_variants,
                         test_candidate_variants_200)
        # dad aff, mum 0/1, dad 0/0 fail
        variants_per_gene_210 = create_test_variants_per_gene(self.variants_210,
                                                              self.family_XX_dad_aff)
        inheritancefilter_210 = InheritanceFiltering(variants_per_gene_210,
                                                     self.family_XX_dad_aff,
                                                     self.genes_hemizygous,
                                                     None,
                                                     None, self.candidate_variants,
                                                      self.inhreport)
        inheritancefilter_210.inheritance_filter_genes()
        test_candidate_variants_210 = {'single_variants': {},
                                       'compound_hets': {}}

        self.assertEqual(inheritancefilter_210.candidate_variants,
                         test_candidate_variants_210)
        # dad aff, mum 1/1, dad 0/0 fail
        variants_per_gene_220 = create_test_variants_per_gene(self.variants_220,
                                                              self.family_XX_dad_aff)
        inheritancefilter_220 = InheritanceFiltering(variants_per_gene_220,
                                                     self.family_XX_dad_aff,
                                                     self.genes_hemizygous,
                                                     None,
                                                     None, self.candidate_variants,
                                                      self.inhreport)
        inheritancefilter_220.inheritance_filter_genes()
        test_candidate_variants_220 = {'single_variants': {},
                                       'compound_hets': {}}

        self.assertEqual(inheritancefilter_220.candidate_variants,
                         test_candidate_variants_220)
        # dad aff, mum 0/0, dad 1/1 fail
        variants_per_gene_202 = create_test_variants_per_gene(self.variants_202,
                                                              self.family_XX_dad_aff)
        inheritancefilter_202 = InheritanceFiltering(variants_per_gene_202,
                                                     self.family_XX_dad_aff,
                                                     self.genes_hemizygous,
                                                     None,
                                                     None, self.candidate_variants,
                                                      self.inhreport)
        inheritancefilter_202.inheritance_filter_genes()
        test_candidate_variants_202 = {'single_variants': {},
                                       'compound_hets': {}}

        self.assertEqual(inheritancefilter_202.candidate_variants,
                         test_candidate_variants_202)
        # dad aff, mum 0/1, dad 1/1 fail
        variants_per_gene_212 = create_test_variants_per_gene(self.variants_212,
                                                              self.family_XX_dad_aff)
        inheritancefilter_212 = InheritanceFiltering(variants_per_gene_212,
                                                     self.family_XX_dad_aff,
                                                     self.genes_hemizygous,
                                                     None,
                                                     None, self.candidate_variants,
                                                      self.inhreport)
        inheritancefilter_212.inheritance_filter_genes()
        test_candidate_variants_212 = {'single_variants': {},
                                       'compound_hets': {}}

        self.assertEqual(inheritancefilter_212.candidate_variants,
                         test_candidate_variants_212)
        # dad aff, mum 1/1, dad 1/1 fail
        variants_per_gene_222 = create_test_variants_per_gene(self.variants_222,
                                                              self.family_XX_dad_aff)
        inheritancefilter_222 = InheritanceFiltering(variants_per_gene_222,
                                                     self.family_XX_dad_aff,
                                                     self.genes_hemizygous,
                                                     None,
                                                     None, self.candidate_variants,
                                                      self.inhreport)
        inheritancefilter_222.inheritance_filter_genes()
        test_candidate_variants_222 = {'single_variants': {},
                                       'compound_hets': {}}

        self.assertEqual(inheritancefilter_222.candidate_variants,
                         test_candidate_variants_222)

        # both unaff, mum 0/0, dad 0/0 fail
        variants_per_gene_200 = create_test_variants_per_gene(self.variants_200,
                                                              self.family_XX_both_unaff)
        inheritancefilter_200 = InheritanceFiltering(variants_per_gene_200,
                                                     self.family_XX_both_unaff,
                                                     self.genes_hemizygous,
                                                     None,
                                                     None, self.candidate_variants,
                                                      self.inhreport)
        inheritancefilter_200.inheritance_filter_genes()
        test_candidate_variants_200 = {'single_variants': {},
                                       'compound_hets': {}}

        self.assertEqual(inheritancefilter_200.candidate_variants,
                         test_candidate_variants_200)
        # both unaff, mum 0/1, dad 0/0 fail
        variants_per_gene_210 = create_test_variants_per_gene(self.variants_210,
                                                              self.family_XX_both_unaff)
        inheritancefilter_210 = InheritanceFiltering(variants_per_gene_210,
                                                     self.family_XX_both_unaff,
                                                     self.genes_hemizygous,
                                                     None,
                                                     None, self.candidate_variants,
                                                      self.inhreport)
        inheritancefilter_210.inheritance_filter_genes()
        test_candidate_variants_210 = {'single_variants': {},
                                       'compound_hets': {}}

        self.assertEqual(inheritancefilter_210.candidate_variants,
                         test_candidate_variants_210)
        # both unaff, mum 1/1, dad 0/0 fail
        variants_per_gene_220 = create_test_variants_per_gene(self.variants_220,
                                                              self.family_XX_both_unaff)
        inheritancefilter_220 = InheritanceFiltering(variants_per_gene_220,
                                                     self.family_XX_both_unaff,
                                                     self.genes_hemizygous,
                                                     None,
                                                     None, self.candidate_variants,
                                                      self.inhreport)
        inheritancefilter_220.inheritance_filter_genes()
        test_candidate_variants_220 = {'single_variants': {},
                                       'compound_hets': {}}

        self.assertEqual(inheritancefilter_220.candidate_variants,
                         test_candidate_variants_220)
        # both unaff, mum 0/0, dad 1/1 fail
        variants_per_gene_202 = create_test_variants_per_gene(self.variants_202,
                                                              self.family_XX_both_unaff)
        inheritancefilter_202 = InheritanceFiltering(variants_per_gene_202,
                                                     self.family_XX_both_unaff,
                                                     self.genes_hemizygous,
                                                     None,
                                                     None, self.candidate_variants,
                                                      self.inhreport)
        inheritancefilter_202.inheritance_filter_genes()
        test_candidate_variants_202 = {'single_variants': {},
                                       'compound_hets': {}}

        self.assertEqual(inheritancefilter_202.candidate_variants,
                         test_candidate_variants_202)
        # both unaff, mum 0/1, dad 1/1 fail
        variants_per_gene_212 = create_test_variants_per_gene(self.variants_212,
                                                              self.family_XX_both_unaff)
        inheritancefilter_212 = InheritanceFiltering(variants_per_gene_212,
                                                     self.family_XX_both_unaff,
                                                     self.genes_hemizygous,
                                                     None,
                                                     None, self.candidate_variants,
                                                      self.inhreport)
        inheritancefilter_212.inheritance_filter_genes()
        test_candidate_variants_212 = {'single_variants': {},
                                       'compound_hets': {}}

        self.assertEqual(inheritancefilter_212.candidate_variants,
                         test_candidate_variants_212)
        # both unaff, mum 1/1, dad 1/1 fail
        variants_per_gene_222 = create_test_variants_per_gene(self.variants_222,
                                                              self.family_XX_both_unaff)
        inheritancefilter_222 = InheritanceFiltering(variants_per_gene_222,
                                                     self.family_XX_both_unaff,
                                                     self.genes_hemizygous,
                                                     None,
                                                     None, self.candidate_variants,
                                                      self.inhreport)
        inheritancefilter_222.inheritance_filter_genes()
        test_candidate_variants_222 = {'single_variants': {},
                                       'compound_hets': {}}

        self.assertEqual(inheritancefilter_222.candidate_variants,
                         test_candidate_variants_222)

    def test_gn_X_linked_dominant_gt_homozygous_parents_filter(self):
        # both aff, mum 0/0, dad 0/0 fail
        variants_per_gene_200 = create_test_variants_per_gene(self.variants_200,
                                                              self.family_XX_both_aff)
        inheritancefilter_200 = InheritanceFiltering(variants_per_gene_200,
                                                     self.family_XX_both_aff,
                                                     self.genes_X_linked_dominant,
                                                     None,
                                                     None, self.candidate_variants,
                                                      self.inhreport)
        inheritancefilter_200.inheritance_filter_genes()
        test_candidate_variants_200 = {'single_variants': {},
                                       'compound_hets': {}}

        self.assertEqual(inheritancefilter_200.candidate_variants,
                         test_candidate_variants_200)
        # both aff, mum 0/1, dad 0/0 fail
        variants_per_gene_210 = create_test_variants_per_gene(self.variants_210,
                                                              self.family_XX_both_aff)
        inheritancefilter_210 = InheritanceFiltering(variants_per_gene_210,
                                                     self.family_XX_both_aff,
                                                     self.genes_X_linked_dominant,
                                                     None,
                                                     None, self.candidate_variants,
                                                      self.inhreport)
        inheritancefilter_210.inheritance_filter_genes()
        test_candidate_variants_210 = {'single_variants': {},
                                       'compound_hets': {}}

        self.assertEqual(inheritancefilter_210.candidate_variants,
                         test_candidate_variants_210)
        # both aff, mum 1/1, dad 0/0 fail
        variants_per_gene_220 = create_test_variants_per_gene(self.variants_220,
                                                              self.family_XX_both_aff)
        inheritancefilter_220 = InheritanceFiltering(variants_per_gene_220,
                                                     self.family_XX_both_aff,
                                                     self.genes_X_linked_dominant,
                                                     None,
                                                     None, self.candidate_variants,
                                                      self.inhreport)
        inheritancefilter_220.inheritance_filter_genes()
        test_candidate_variants_220 = {'single_variants': {},
                                       'compound_hets': {}}

        self.assertEqual(inheritancefilter_220.candidate_variants,
                         test_candidate_variants_220)
        # both aff, mum 0/0, dad 1/1 fail
        variants_per_gene_202 = create_test_variants_per_gene(self.variants_202,
                                                              self.family_XX_both_aff)
        inheritancefilter_202 = InheritanceFiltering(variants_per_gene_202,
                                                     self.family_XX_both_aff,
                                                     self.genes_X_linked_dominant,
                                                     None,
                                                     None, self.candidate_variants,
                                                      self.inhreport)
        inheritancefilter_202.inheritance_filter_genes()
        test_candidate_variants_202 = {'single_variants': {},
                                       'compound_hets': {}}

        self.assertEqual(inheritancefilter_202.candidate_variants,
                         test_candidate_variants_202)
        # both aff, mum 0/1, dad 1/1 fail
        variants_per_gene_212 = create_test_variants_per_gene(self.variants_212,
                                                              self.family_XX_both_aff)
        inheritancefilter_212 = InheritanceFiltering(variants_per_gene_212,
                                                     self.family_XX_both_aff,
                                                     self.genes_X_linked_dominant,
                                                     None,
                                                     None, self.candidate_variants,
                                                      self.inhreport)
        inheritancefilter_212.inheritance_filter_genes()
        test_candidate_variants_212 = {'single_variants': {},
                                       'compound_hets': {}}

        self.assertEqual(inheritancefilter_212.candidate_variants,
                         test_candidate_variants_212)
        # both aff, mum 1/1, dad 1/1 fail
        variants_per_gene_222 = create_test_variants_per_gene(self.variants_222,
                                                              self.family_XX_both_aff)
        inheritancefilter_222 = InheritanceFiltering(variants_per_gene_222,
                                                     self.family_XX_both_aff,
                                                     self.genes_X_linked_dominant,
                                                     None,
                                                     None, self.candidate_variants,
                                                      self.inhreport)
        inheritancefilter_222.inheritance_filter_genes()
        test_candidate_variants_222 = {'single_variants': {},
                                       'compound_hets': {}}

        self.assertEqual(inheritancefilter_222.candidate_variants,
                         test_candidate_variants_222)

        # mum aff, mum 0/0, dad 0/0 fail
        variants_per_gene_200 = create_test_variants_per_gene(self.variants_200,
                                                              self.family_XX_mum_aff)
        inheritancefilter_200 = InheritanceFiltering(variants_per_gene_200,
                                                     self.family_XX_mum_aff,
                                                     self.genes_X_linked_dominant,
                                                     None,
                                                     None, self.candidate_variants,
                                                      self.inhreport)
        inheritancefilter_200.inheritance_filter_genes()
        test_candidate_variants_200 = {'single_variants': {},
                                       'compound_hets': {}}

        self.assertEqual(inheritancefilter_200.candidate_variants,
                         test_candidate_variants_200)
        # mum aff, mum 0/1, dad 0/0 fail
        variants_per_gene_210 = create_test_variants_per_gene(self.variants_210,
                                                              self.family_XX_mum_aff)
        inheritancefilter_210 = InheritanceFiltering(variants_per_gene_210,
                                                     self.family_XX_mum_aff,
                                                     self.genes_X_linked_dominant,
                                                     None,
                                                     None, self.candidate_variants,
                                                      self.inhreport)
        inheritancefilter_210.inheritance_filter_genes()
        test_candidate_variants_210 = {'single_variants': {},
                                       'compound_hets': {}}

        self.assertEqual(inheritancefilter_210.candidate_variants,
                         test_candidate_variants_210)
        # mum aff, mum 1/1, dad 0/0 fail
        variants_per_gene_220 = create_test_variants_per_gene(self.variants_220,
                                                              self.family_XX_mum_aff)
        inheritancefilter_220 = InheritanceFiltering(variants_per_gene_220,
                                                     self.family_XX_mum_aff,
                                                     self.genes_X_linked_dominant,
                                                     None,
                                                     None, self.candidate_variants,
                                                      self.inhreport)
        inheritancefilter_220.inheritance_filter_genes()
        test_candidate_variants_220 = {'single_variants': {},
                                       'compound_hets': {}}

        self.assertEqual(inheritancefilter_220.candidate_variants,
                         test_candidate_variants_220)
        # mum aff, mum 0/0, dad 1/1 fail
        variants_per_gene_202 = create_test_variants_per_gene(self.variants_202,
                                                              self.family_XX_mum_aff)
        inheritancefilter_202 = InheritanceFiltering(variants_per_gene_202,
                                                     self.family_XX_mum_aff,
                                                     self.genes_X_linked_dominant,
                                                     None,
                                                     None, self.candidate_variants,
                                                      self.inhreport)
        inheritancefilter_202.inheritance_filter_genes()
        test_candidate_variants_202 = {'single_variants': {},
                                       'compound_hets': {}}

        self.assertEqual(inheritancefilter_202.candidate_variants,
                         test_candidate_variants_202)
        # mum aff, mum 0/1, dad 1/1 fail
        variants_per_gene_212 = create_test_variants_per_gene(self.variants_212,
                                                              self.family_XX_mum_aff)
        inheritancefilter_212 = InheritanceFiltering(variants_per_gene_212,
                                                     self.family_XX_mum_aff,
                                                     self.genes_X_linked_dominant,
                                                     None,
                                                     None, self.candidate_variants,
                                                      self.inhreport)
        inheritancefilter_212.inheritance_filter_genes()
        test_candidate_variants_212 = {'single_variants': {},
                                       'compound_hets': {}}

        self.assertEqual(inheritancefilter_212.candidate_variants,
                         test_candidate_variants_212)
        # mum aff, mum 1/1, dad 1/1 fail
        variants_per_gene_222 = create_test_variants_per_gene(self.variants_222,
                                                              self.family_XX_mum_aff)
        inheritancefilter_222 = InheritanceFiltering(variants_per_gene_222,
                                                     self.family_XX_mum_aff,
                                                     self.genes_X_linked_dominant,
                                                     None,
                                                     None, self.candidate_variants,
                                                      self.inhreport)
        inheritancefilter_222.inheritance_filter_genes()
        test_candidate_variants_222 = {'single_variants': {},
                                       'compound_hets': {}}

        self.assertEqual(inheritancefilter_222.candidate_variants,
                         test_candidate_variants_222)

        # dad aff, mum 0/0, dad 0/0 fail
        variants_per_gene_200 = create_test_variants_per_gene(self.variants_200,
                                                              self.family_XX_dad_aff)
        inheritancefilter_200 = InheritanceFiltering(variants_per_gene_200,
                                                     self.family_XX_dad_aff,
                                                     self.genes_X_linked_dominant,
                                                     None,
                                                     None, self.candidate_variants,
                                                      self.inhreport)
        inheritancefilter_200.inheritance_filter_genes()
        test_candidate_variants_200 = {'single_variants': {},
                                       'compound_hets': {}}

        self.assertEqual(inheritancefilter_200.candidate_variants,
                         test_candidate_variants_200)
        # dad aff, mum 0/1, dad 0/0 fail
        variants_per_gene_210 = create_test_variants_per_gene(self.variants_210,
                                                              self.family_XX_dad_aff)
        inheritancefilter_210 = InheritanceFiltering(variants_per_gene_210,
                                                     self.family_XX_dad_aff,
                                                     self.genes_X_linked_dominant,
                                                     None,
                                                     None, self.candidate_variants,
                                                      self.inhreport)
        inheritancefilter_210.inheritance_filter_genes()
        test_candidate_variants_210 = {'single_variants': {},
                                       'compound_hets': {}}

        self.assertEqual(inheritancefilter_210.candidate_variants,
                         test_candidate_variants_210)
        # dad aff, mum 1/1, dad 0/0 fail
        variants_per_gene_220 = create_test_variants_per_gene(self.variants_220,
                                                              self.family_XX_dad_aff)
        inheritancefilter_220 = InheritanceFiltering(variants_per_gene_220,
                                                     self.family_XX_dad_aff,
                                                     self.genes_X_linked_dominant,
                                                     None,
                                                     None, self.candidate_variants,
                                                      self.inhreport)
        inheritancefilter_220.inheritance_filter_genes()
        test_candidate_variants_220 = {'single_variants': {},
                                       'compound_hets': {}}

        self.assertEqual(inheritancefilter_220.candidate_variants,
                         test_candidate_variants_220)
        # dad aff, mum 0/0, dad 1/1 fail
        variants_per_gene_202 = create_test_variants_per_gene(self.variants_202,
                                                              self.family_XX_dad_aff)
        inheritancefilter_202 = InheritanceFiltering(variants_per_gene_202,
                                                     self.family_XX_dad_aff,
                                                     self.genes_X_linked_dominant,
                                                     None,
                                                     None, self.candidate_variants,
                                                      self.inhreport)
        inheritancefilter_202.inheritance_filter_genes()
        test_candidate_variants_202 = {'single_variants': {},
                                       'compound_hets': {}}

        self.assertEqual(inheritancefilter_202.candidate_variants,
                         test_candidate_variants_202)
        # dad aff, mum 0/1, dad 1/1 fail
        variants_per_gene_212 = create_test_variants_per_gene(self.variants_212,
                                                              self.family_XX_dad_aff)
        inheritancefilter_212 = InheritanceFiltering(variants_per_gene_212,
                                                     self.family_XX_dad_aff,
                                                     self.genes_X_linked_dominant,
                                                     None,
                                                     None, self.candidate_variants,
                                                      self.inhreport)
        inheritancefilter_212.inheritance_filter_genes()
        test_candidate_variants_212 = {'single_variants': {},
                                       'compound_hets': {}}

        self.assertEqual(inheritancefilter_212.candidate_variants,
                         test_candidate_variants_212)
        # dad aff, mum 1/1, dad 1/1 fail
        variants_per_gene_222 = create_test_variants_per_gene(self.variants_222,
                                                              self.family_XX_dad_aff)
        inheritancefilter_222 = InheritanceFiltering(variants_per_gene_222,
                                                     self.family_XX_dad_aff,
                                                     self.genes_X_linked_dominant,
                                                     None,
                                                     None, self.candidate_variants,
                                                      self.inhreport)
        inheritancefilter_222.inheritance_filter_genes()
        test_candidate_variants_222 = {'single_variants': {},
                                       'compound_hets': {}}

        self.assertEqual(inheritancefilter_222.candidate_variants,
                         test_candidate_variants_222)

        # both unaff, mum 0/0, dad 0/0 fail
        variants_per_gene_200 = create_test_variants_per_gene(self.variants_200,
                                                              self.family_XX_both_unaff)
        inheritancefilter_200 = InheritanceFiltering(variants_per_gene_200,
                                                     self.family_XX_both_unaff,
                                                     self.genes_X_linked_dominant,
                                                     None,
                                                     None, self.candidate_variants,
                                                      self.inhreport)
        inheritancefilter_200.inheritance_filter_genes()
        test_candidate_variants_200 = {'single_variants': {},
                                       'compound_hets': {}}

        self.assertEqual(inheritancefilter_200.candidate_variants,
                         test_candidate_variants_200)
        # both unaff, mum 0/1, dad 0/0 fail
        variants_per_gene_210 = create_test_variants_per_gene(self.variants_210,
                                                              self.family_XX_both_unaff)
        inheritancefilter_210 = InheritanceFiltering(variants_per_gene_210,
                                                     self.family_XX_both_unaff,
                                                     self.genes_X_linked_dominant,
                                                     None,
                                                     None, self.candidate_variants,
                                                      self.inhreport)
        inheritancefilter_210.inheritance_filter_genes()
        test_candidate_variants_210 = {'single_variants': {},
                                       'compound_hets': {}}

        self.assertEqual(inheritancefilter_210.candidate_variants,
                         test_candidate_variants_210)
        # both unaff, mum 1/1, dad 0/0 fail
        variants_per_gene_220 = create_test_variants_per_gene(self.variants_220,
                                                              self.family_XX_both_unaff)
        inheritancefilter_220 = InheritanceFiltering(variants_per_gene_220,
                                                     self.family_XX_both_unaff,
                                                     self.genes_X_linked_dominant,
                                                     None,
                                                     None, self.candidate_variants,
                                                      self.inhreport)
        inheritancefilter_220.inheritance_filter_genes()
        test_candidate_variants_220 = {'single_variants': {},
                                       'compound_hets': {}}

        self.assertEqual(inheritancefilter_220.candidate_variants,
                         test_candidate_variants_220)
        # both unaff, mum 0/0, dad 1/1 fail
        variants_per_gene_202 = create_test_variants_per_gene(self.variants_202,
                                                              self.family_XX_both_unaff)
        inheritancefilter_202 = InheritanceFiltering(variants_per_gene_202,
                                                     self.family_XX_both_unaff,
                                                     self.genes_X_linked_dominant,
                                                     None,
                                                     None, self.candidate_variants,
                                                      self.inhreport)
        inheritancefilter_202.inheritance_filter_genes()
        test_candidate_variants_202 = {'single_variants': {},
                                       'compound_hets': {}}

        self.assertEqual(inheritancefilter_202.candidate_variants,
                         test_candidate_variants_202)
        # both unaff, mum 0/1, dad 1/1 fail
        variants_per_gene_212 = create_test_variants_per_gene(self.variants_212,
                                                              self.family_XX_both_unaff)
        inheritancefilter_212 = InheritanceFiltering(variants_per_gene_212,
                                                     self.family_XX_both_unaff,
                                                     self.genes_X_linked_dominant,
                                                     None,
                                                     None, self.candidate_variants,
                                                      self.inhreport)
        inheritancefilter_212.inheritance_filter_genes()
        test_candidate_variants_212 = {'single_variants': {},
                                       'compound_hets': {}}

        self.assertEqual(inheritancefilter_212.candidate_variants,
                         test_candidate_variants_212)
        # both unaff, mum 1/1, dad 1/1 fail
        variants_per_gene_222 = create_test_variants_per_gene(self.variants_222,
                                                              self.family_XX_both_unaff)
        inheritancefilter_222 = InheritanceFiltering(variants_per_gene_222,
                                                     self.family_XX_both_unaff,
                                                     self.genes_X_linked_dominant,
                                                     None,
                                                     None, self.candidate_variants,
                                                      self.inhreport)
        inheritancefilter_222.inheritance_filter_genes()
        test_candidate_variants_222 = {'single_variants': {},
                                       'compound_hets': {}}

        self.assertEqual(inheritancefilter_222.candidate_variants,
                         test_candidate_variants_222)

    def test_gn_X_linked_over_dominant_gt_homozygous_parents_filter(self):
        # both aff, mum 0/0, dad 0/0 fail
        variants_per_gene_200 = create_test_variants_per_gene(self.variants_200,
                                                              self.family_XX_both_aff)
        inheritancefilter_200 = InheritanceFiltering(variants_per_gene_200,
                                                     self.family_XX_both_aff,
                                                     self.genes_X_linked_over_dominance,
                                                     None,
                                                     None, self.candidate_variants,
                                                      self.inhreport)
        inheritancefilter_200.inheritance_filter_genes()
        test_candidate_variants_200 = {'single_variants': {},
                                       'compound_hets': {}}

        self.assertEqual(inheritancefilter_200.candidate_variants,
                         test_candidate_variants_200)
        # both aff, mum 0/1, dad 0/0 fail
        variants_per_gene_210 = create_test_variants_per_gene(self.variants_210,
                                                              self.family_XX_both_aff)
        inheritancefilter_210 = InheritanceFiltering(variants_per_gene_210,
                                                     self.family_XX_both_aff,
                                                     self.genes_X_linked_over_dominance,
                                                     None,
                                                     None, self.candidate_variants,
                                                      self.inhreport)
        inheritancefilter_210.inheritance_filter_genes()
        test_candidate_variants_210 = {'single_variants': {},
                                       'compound_hets': {}}

        self.assertEqual(inheritancefilter_210.candidate_variants,
                         test_candidate_variants_210)
        # both aff, mum 1/1, dad 0/0 fail
        variants_per_gene_220 = create_test_variants_per_gene(self.variants_220,
                                                              self.family_XX_both_aff)
        inheritancefilter_220 = InheritanceFiltering(variants_per_gene_220,
                                                     self.family_XX_both_aff,
                                                     self.genes_X_linked_over_dominance,
                                                     None,
                                                     None, self.candidate_variants,
                                                      self.inhreport)
        inheritancefilter_220.inheritance_filter_genes()
        test_candidate_variants_220 = {'single_variants': {},
                                       'compound_hets': {}}

        self.assertEqual(inheritancefilter_220.candidate_variants,
                         test_candidate_variants_220)
        # both aff, mum 0/0, dad 1/1 fail
        variants_per_gene_202 = create_test_variants_per_gene(self.variants_202,
                                                              self.family_XX_both_aff)
        inheritancefilter_202 = InheritanceFiltering(variants_per_gene_202,
                                                     self.family_XX_both_aff,
                                                     self.genes_X_linked_over_dominance,
                                                     None,
                                                     None, self.candidate_variants,
                                                      self.inhreport)
        inheritancefilter_202.inheritance_filter_genes()
        test_candidate_variants_202 = {'single_variants': {},
                                       'compound_hets': {}}

        self.assertEqual(inheritancefilter_202.candidate_variants,
                         test_candidate_variants_202)
        # both aff, mum 0/1, dad 1/1 fail
        variants_per_gene_212 = create_test_variants_per_gene(self.variants_212,
                                                              self.family_XX_both_aff)
        inheritancefilter_212 = InheritanceFiltering(variants_per_gene_212,
                                                     self.family_XX_both_aff,
                                                     self.genes_X_linked_over_dominance,
                                                     None,
                                                     None, self.candidate_variants,
                                                      self.inhreport)
        inheritancefilter_212.inheritance_filter_genes()
        test_candidate_variants_212 = {'single_variants': {},
                                       'compound_hets': {}}

        self.assertEqual(inheritancefilter_212.candidate_variants,
                         test_candidate_variants_212)
        # both aff, mum 1/1, dad 1/1 fail
        variants_per_gene_222 = create_test_variants_per_gene(self.variants_222,
                                                              self.family_XX_both_aff)
        inheritancefilter_222 = InheritanceFiltering(variants_per_gene_222,
                                                     self.family_XX_both_aff,
                                                     self.genes_X_linked_over_dominance,
                                                     None,
                                                     None, self.candidate_variants,
                                                      self.inhreport)
        inheritancefilter_222.inheritance_filter_genes()
        test_candidate_variants_222 = {'single_variants': {},
                                       'compound_hets': {}}

        self.assertEqual(inheritancefilter_222.candidate_variants,
                         test_candidate_variants_222)

        # mum aff, mum 0/0, dad 0/0 fail
        variants_per_gene_200 = create_test_variants_per_gene(self.variants_200,
                                                              self.family_XX_mum_aff)
        inheritancefilter_200 = InheritanceFiltering(variants_per_gene_200,
                                                     self.family_XX_mum_aff,
                                                     self.genes_X_linked_over_dominance,
                                                     None,
                                                     None, self.candidate_variants,
                                                      self.inhreport)
        inheritancefilter_200.inheritance_filter_genes()
        test_candidate_variants_200 = {'single_variants': {},
                                       'compound_hets': {}}

        self.assertEqual(inheritancefilter_200.candidate_variants,
                         test_candidate_variants_200)
        # mum aff, mum 0/1, dad 0/0 fail
        variants_per_gene_210 = create_test_variants_per_gene(self.variants_210,
                                                              self.family_XX_mum_aff)
        inheritancefilter_210 = InheritanceFiltering(variants_per_gene_210,
                                                     self.family_XX_mum_aff,
                                                     self.genes_X_linked_over_dominance,
                                                     None,
                                                     None, self.candidate_variants,
                                                      self.inhreport)
        inheritancefilter_210.inheritance_filter_genes()
        test_candidate_variants_210 = {'single_variants': {},
                                       'compound_hets': {}}

        self.assertEqual(inheritancefilter_210.candidate_variants,
                         test_candidate_variants_210)
        # mum aff, mum 1/1, dad 0/0 fail
        variants_per_gene_220 = create_test_variants_per_gene(self.variants_220,
                                                              self.family_XX_mum_aff)
        inheritancefilter_220 = InheritanceFiltering(variants_per_gene_220,
                                                     self.family_XX_mum_aff,
                                                     self.genes_X_linked_over_dominance,
                                                     None,
                                                     None, self.candidate_variants,
                                                      self.inhreport)
        inheritancefilter_220.inheritance_filter_genes()
        test_candidate_variants_220 = {'single_variants': {},
                                       'compound_hets': {}}

        self.assertEqual(inheritancefilter_220.candidate_variants,
                         test_candidate_variants_220)
        # mum aff, mum 0/0, dad 1/1 fail
        variants_per_gene_202 = create_test_variants_per_gene(self.variants_202,
                                                              self.family_XX_mum_aff)
        inheritancefilter_202 = InheritanceFiltering(variants_per_gene_202,
                                                     self.family_XX_mum_aff,
                                                     self.genes_X_linked_over_dominance,
                                                     None,
                                                     None, self.candidate_variants,
                                                      self.inhreport)
        inheritancefilter_202.inheritance_filter_genes()
        test_candidate_variants_202 = {'single_variants': {},
                                       'compound_hets': {}}

        self.assertEqual(inheritancefilter_202.candidate_variants,
                         test_candidate_variants_202)
        # mum aff, mum 0/1, dad 1/1 fail
        variants_per_gene_212 = create_test_variants_per_gene(self.variants_212,
                                                              self.family_XX_mum_aff)
        inheritancefilter_212 = InheritanceFiltering(variants_per_gene_212,
                                                     self.family_XX_mum_aff,
                                                     self.genes_X_linked_over_dominance,
                                                     None,
                                                     None, self.candidate_variants,
                                                      self.inhreport)
        inheritancefilter_212.inheritance_filter_genes()
        test_candidate_variants_212 = {'single_variants': {},
                                       'compound_hets': {}}

        self.assertEqual(inheritancefilter_212.candidate_variants,
                         test_candidate_variants_212)
        # mum aff, mum 1/1, dad 1/1 fail
        variants_per_gene_222 = create_test_variants_per_gene(self.variants_222,
                                                              self.family_XX_mum_aff)
        inheritancefilter_222 = InheritanceFiltering(variants_per_gene_222,
                                                     self.family_XX_mum_aff,
                                                     self.genes_X_linked_over_dominance,
                                                     None,
                                                     None, self.candidate_variants,
                                                      self.inhreport)
        inheritancefilter_222.inheritance_filter_genes()
        test_candidate_variants_222 = {'single_variants': {},
                                       'compound_hets': {}}

        self.assertEqual(inheritancefilter_222.candidate_variants,
                         test_candidate_variants_222)

        # dad aff, mum 0/0, dad 0/0 fail
        variants_per_gene_200 = create_test_variants_per_gene(self.variants_200,
                                                              self.family_XX_dad_aff)
        inheritancefilter_200 = InheritanceFiltering(variants_per_gene_200,
                                                     self.family_XX_dad_aff,
                                                     self.genes_X_linked_over_dominance,
                                                     None,
                                                     None, self.candidate_variants,
                                                      self.inhreport)
        inheritancefilter_200.inheritance_filter_genes()
        test_candidate_variants_200 = {'single_variants': {},
                                       'compound_hets': {}}

        self.assertEqual(inheritancefilter_200.candidate_variants,
                         test_candidate_variants_200)
        # dad aff, mum 0/1, dad 0/0 fail
        variants_per_gene_210 = create_test_variants_per_gene(self.variants_210,
                                                              self.family_XX_dad_aff)
        inheritancefilter_210 = InheritanceFiltering(variants_per_gene_210,
                                                     self.family_XX_dad_aff,
                                                     self.genes_X_linked_over_dominance,
                                                     None,
                                                     None, self.candidate_variants,
                                                      self.inhreport)
        inheritancefilter_210.inheritance_filter_genes()
        test_candidate_variants_210 = {'single_variants': {},
                                       'compound_hets': {}}

        self.assertEqual(inheritancefilter_210.candidate_variants,
                         test_candidate_variants_210)
        # dad aff, mum 1/1, dad 0/0 fail
        variants_per_gene_220 = create_test_variants_per_gene(self.variants_220,
                                                              self.family_XX_dad_aff)
        inheritancefilter_220 = InheritanceFiltering(variants_per_gene_220,
                                                     self.family_XX_dad_aff,
                                                     self.genes_X_linked_over_dominance,
                                                     None,
                                                     None, self.candidate_variants,
                                                      self.inhreport)
        inheritancefilter_220.inheritance_filter_genes()
        test_candidate_variants_220 = {'single_variants': {},
                                       'compound_hets': {}}

        self.assertEqual(inheritancefilter_220.candidate_variants,
                         test_candidate_variants_220)
        # dad aff, mum 0/0, dad 1/1 fail
        variants_per_gene_202 = create_test_variants_per_gene(self.variants_202,
                                                              self.family_XX_dad_aff)
        inheritancefilter_202 = InheritanceFiltering(variants_per_gene_202,
                                                     self.family_XX_dad_aff,
                                                     self.genes_X_linked_over_dominance,
                                                     None,
                                                     None, self.candidate_variants,
                                                      self.inhreport)
        inheritancefilter_202.inheritance_filter_genes()
        test_candidate_variants_202 = {'single_variants': {},
                                       'compound_hets': {}}

        self.assertEqual(inheritancefilter_202.candidate_variants,
                         test_candidate_variants_202)
        # dad aff, mum 0/1, dad 1/1 fail
        variants_per_gene_212 = create_test_variants_per_gene(self.variants_212,
                                                              self.family_XX_dad_aff)
        inheritancefilter_212 = InheritanceFiltering(variants_per_gene_212,
                                                     self.family_XX_dad_aff,
                                                     self.genes_X_linked_over_dominance,
                                                     None,
                                                     None, self.candidate_variants,
                                                      self.inhreport)
        inheritancefilter_212.inheritance_filter_genes()
        test_candidate_variants_212 = {'single_variants': {},
                                       'compound_hets': {}}

        self.assertEqual(inheritancefilter_212.candidate_variants,
                         test_candidate_variants_212)
        # dad aff, mum 1/1, dad 1/1 fail
        variants_per_gene_222 = create_test_variants_per_gene(self.variants_222,
                                                              self.family_XX_dad_aff)
        inheritancefilter_222 = InheritanceFiltering(variants_per_gene_222,
                                                     self.family_XX_dad_aff,
                                                     self.genes_X_linked_over_dominance,
                                                     None,
                                                     None, self.candidate_variants,
                                                      self.inhreport)
        inheritancefilter_222.inheritance_filter_genes()
        test_candidate_variants_222 = {'single_variants': {},
                                       'compound_hets': {}}

        self.assertEqual(inheritancefilter_222.candidate_variants,
                         test_candidate_variants_222)

        # both unaff, mum 0/0, dad 0/0 fail
        variants_per_gene_200 = create_test_variants_per_gene(self.variants_200,
                                                              self.family_XX_both_unaff)
        inheritancefilter_200 = InheritanceFiltering(variants_per_gene_200,
                                                     self.family_XX_both_unaff,
                                                     self.genes_X_linked_over_dominance,
                                                     None,
                                                     None, self.candidate_variants,
                                                      self.inhreport)
        inheritancefilter_200.inheritance_filter_genes()
        test_candidate_variants_200 = {'single_variants': {},
                                       'compound_hets': {}}

        self.assertEqual(inheritancefilter_200.candidate_variants,
                         test_candidate_variants_200)
        # both unaff, mum 0/1, dad 0/0 fail
        variants_per_gene_210 = create_test_variants_per_gene(self.variants_210,
                                                              self.family_XX_both_unaff)
        inheritancefilter_210 = InheritanceFiltering(variants_per_gene_210,
                                                     self.family_XX_both_unaff,
                                                     self.genes_X_linked_over_dominance,
                                                     None,
                                                     None, self.candidate_variants,
                                                      self.inhreport)
        inheritancefilter_210.inheritance_filter_genes()
        test_candidate_variants_210 = {'single_variants': {},
                                       'compound_hets': {}}

        self.assertEqual(inheritancefilter_210.candidate_variants,
                         test_candidate_variants_210)
        # both unaff, mum 1/1, dad 0/0 fail
        variants_per_gene_220 = create_test_variants_per_gene(self.variants_220,
                                                              self.family_XX_both_unaff)
        inheritancefilter_220 = InheritanceFiltering(variants_per_gene_220,
                                                     self.family_XX_both_unaff,
                                                     self.genes_X_linked_over_dominance,
                                                     None,
                                                     None, self.candidate_variants,
                                                      self.inhreport)
        inheritancefilter_220.inheritance_filter_genes()
        test_candidate_variants_220 = {'single_variants': {},
                                       'compound_hets': {}}

        self.assertEqual(inheritancefilter_220.candidate_variants,
                         test_candidate_variants_220)
        # both unaff, mum 0/0, dad 1/1 fail
        variants_per_gene_202 = create_test_variants_per_gene(self.variants_202,
                                                              self.family_XX_both_unaff)
        inheritancefilter_202 = InheritanceFiltering(variants_per_gene_202,
                                                     self.family_XX_both_unaff,
                                                     self.genes_X_linked_over_dominance,
                                                     None,
                                                     None, self.candidate_variants,
                                                      self.inhreport)
        inheritancefilter_202.inheritance_filter_genes()
        test_candidate_variants_202 = {'single_variants': {},
                                       'compound_hets': {}}

        self.assertEqual(inheritancefilter_202.candidate_variants,
                         test_candidate_variants_202)
        # both unaff, mum 0/1, dad 1/1 fail
        variants_per_gene_212 = create_test_variants_per_gene(self.variants_212,
                                                              self.family_XX_both_unaff)
        inheritancefilter_212 = InheritanceFiltering(variants_per_gene_212,
                                                     self.family_XX_both_unaff,
                                                     self.genes_X_linked_over_dominance,
                                                     None,
                                                     None, self.candidate_variants,
                                                      self.inhreport)
        inheritancefilter_212.inheritance_filter_genes()
        test_candidate_variants_212 = {'single_variants': {},
                                       'compound_hets': {}}

        self.assertEqual(inheritancefilter_212.candidate_variants,
                         test_candidate_variants_212)
        # both unaff, mum 1/1, dad 1/1 fail
        variants_per_gene_222 = create_test_variants_per_gene(self.variants_222,
                                                              self.family_XX_both_unaff)
        inheritancefilter_222 = InheritanceFiltering(variants_per_gene_222,
                                                     self.family_XX_both_unaff,
                                                     self.genes_X_linked_over_dominance,
                                                     None,
                                                     None, self.candidate_variants,
                                                      self.inhreport)
        inheritancefilter_222.inheritance_filter_genes()
        test_candidate_variants_222 = {'single_variants': {},
                                       'compound_hets': {}}

        self.assertEqual(inheritancefilter_222.candidate_variants,
                         test_candidate_variants_222)

    def test_gn_hemizygous_gt_hemizygous_parents_filter(self):
        # both aff, mum 0/0, dad 0/0 pass
        variants_per_gene_200 = create_test_variants_per_gene(self.variants_200,
                                                              self.family_XY_both_aff)
        inheritancefilter_200 = InheritanceFiltering(variants_per_gene_200,
                                                     self.family_XY_both_aff,
                                                     self.genes_hemizygous,
                                                     None,
                                                     None, self.candidate_variants,
                                                      self.inhreport)
        inheritancefilter_200.inheritance_filter_genes()
        test_candidate_variants_200 = {'single_variants': {
            'X_1097183_A_GG': {
                'mode': {'hemizygous'},
                'variant':
                    variants_per_gene_200[
                        '1234'][
                        'X_1097183_A_GG'][
                        'child'],
                'hgncid': '1234'}}, 'compound_hets': {}}

        self.assertEqual(inheritancefilter_200.candidate_variants,
                         test_candidate_variants_200)
        # both aff, mum 0/1, dad 0/0 pass
        variants_per_gene_210 = create_test_variants_per_gene(self.variants_210,
                                                              self.family_XY_both_aff)
        inheritancefilter_210 = InheritanceFiltering(variants_per_gene_200,
                                                     self.family_XY_both_aff,
                                                     self.genes_hemizygous,
                                                     None,
                                                     None, self.candidate_variants,
                                                      self.inhreport)
        inheritancefilter_210.inheritance_filter_genes()
        test_candidate_variants_210 = {'single_variants': {
            'X_1097183_A_GG': {
                'mode': {'hemizygous'},
                'variant':
                    variants_per_gene_210[
                        '1234'][
                        'X_1097183_A_GG'][
                        'child'],
                'hgncid': '1234'}}, 'compound_hets': {}}

        self.assertEqual(inheritancefilter_210.candidate_variants,
                         test_candidate_variants_210)
        # both aff, mum 1/1, dad 0/0 pass
        variants_per_gene_220 = create_test_variants_per_gene(self.variants_220,
                                                              self.family_XY_both_aff)
        inheritancefilter_220 = InheritanceFiltering(variants_per_gene_220,
                                                     self.family_XY_both_aff,
                                                     self.genes_hemizygous,
                                                     None,
                                                     None, self.candidate_variants,
                                                      self.inhreport)
        inheritancefilter_220.inheritance_filter_genes()
        test_candidate_variants_220 = {'single_variants': {
            'X_1097183_A_GG': {
                'mode': {'hemizygous'},
                'variant':
                    variants_per_gene_220[
                        '1234'][
                        'X_1097183_A_GG'][
                        'child'],
                'hgncid': '1234'}}, 'compound_hets': {}}

        self.assertEqual(inheritancefilter_220.candidate_variants,
                         test_candidate_variants_220)
        # both aff, mum 0/0, dad 1/1 pass
        variants_per_gene_202 = create_test_variants_per_gene(self.variants_202,
                                                              self.family_XY_both_aff)
        inheritancefilter_202 = InheritanceFiltering(variants_per_gene_202,
                                                     self.family_XY_both_aff,
                                                     self.genes_hemizygous,
                                                     None,
                                                     None, self.candidate_variants,
                                                      self.inhreport)
        inheritancefilter_202.inheritance_filter_genes()
        test_candidate_variants_202 = {'single_variants': {
            'X_1097183_A_GG': {
                'mode': {'hemizygous'},
                'variant':
                    variants_per_gene_202[
                        '1234'][
                        'X_1097183_A_GG'][
                        'child'],
                'hgncid': '1234'}}, 'compound_hets': {}}

        self.assertEqual(inheritancefilter_202.candidate_variants,
                         test_candidate_variants_202)
        # both aff, mum 0/1, dad 1/1 pass
        variants_per_gene_212 = create_test_variants_per_gene(self.variants_212,
                                                              self.family_XY_both_aff)
        inheritancefilter_212 = InheritanceFiltering(variants_per_gene_212,
                                                     self.family_XY_both_aff,
                                                     self.genes_hemizygous,
                                                     None,
                                                     None, self.candidate_variants,
                                                      self.inhreport)
        inheritancefilter_212.inheritance_filter_genes()
        test_candidate_variants_212 = {'single_variants': {
            'X_1097183_A_GG': {
                'mode': {'hemizygous'},
                'variant':
                    variants_per_gene_212[
                        '1234'][
                        'X_1097183_A_GG'][
                        'child'],
                'hgncid': '1234'}}, 'compound_hets': {}}

        self.assertEqual(inheritancefilter_212.candidate_variants,
                         test_candidate_variants_212)
        # both aff, mum 1/1, dad 1/1 pass
        variants_per_gene_222 = create_test_variants_per_gene(self.variants_222,
                                                              self.family_XY_both_aff)
        inheritancefilter_222 = InheritanceFiltering(variants_per_gene_222,
                                                     self.family_XY_both_aff,
                                                     self.genes_hemizygous,
                                                     None,
                                                     None, self.candidate_variants,
                                                      self.inhreport)
        inheritancefilter_222.inheritance_filter_genes()
        test_candidate_variants_222 = {'single_variants': {
            'X_1097183_A_GG': {
                'mode': {'hemizygous'},
                'variant':
                    variants_per_gene_222[
                        '1234'][
                        'X_1097183_A_GG'][
                        'child'],
                'hgncid': '1234'}}, 'compound_hets': {}}

        self.assertEqual(inheritancefilter_222.candidate_variants,
                         test_candidate_variants_222)

        # mum aff, mum 0/0, dad 0/0 pass
        variants_per_gene_200 = create_test_variants_per_gene(self.variants_200,
                                                              self.family_XY_mum_aff)
        inheritancefilter_200 = InheritanceFiltering(variants_per_gene_200,
                                                     self.family_XY_mum_aff,
                                                     self.genes_hemizygous,
                                                     None,
                                                     None, self.candidate_variants,
                                                      self.inhreport)
        inheritancefilter_200.inheritance_filter_genes()
        test_candidate_variants_200 = {'single_variants': {
            'X_1097183_A_GG': {
                'mode': {'hemizygous'},
                'variant':
                    variants_per_gene_200[
                        '1234'][
                        'X_1097183_A_GG'][
                        'child'],
                'hgncid': '1234'}}, 'compound_hets': {}}

        self.assertEqual(inheritancefilter_200.candidate_variants,
                         test_candidate_variants_200)
        # mum aff, mum 0/1, dad 0/0 pass
        variants_per_gene_210 = create_test_variants_per_gene(self.variants_210,
                                                              self.family_XY_mum_aff)
        inheritancefilter_210 = InheritanceFiltering(variants_per_gene_210,
                                                     self.family_XY_mum_aff,
                                                     self.genes_hemizygous,
                                                     None,
                                                     None, self.candidate_variants,
                                                      self.inhreport)
        inheritancefilter_210.inheritance_filter_genes()
        test_candidate_variants_210 = {'single_variants': {
            'X_1097183_A_GG': {
                'mode': {'hemizygous'},
                'variant':
                    variants_per_gene_210[
                        '1234'][
                        'X_1097183_A_GG'][
                        'child'],
                'hgncid': '1234'}}, 'compound_hets': {}}

        self.assertEqual(inheritancefilter_210.candidate_variants,
                         test_candidate_variants_210)
        # mum aff, mum 1/1, dad 0/0 pass
        variants_per_gene_220 = create_test_variants_per_gene(self.variants_220,
                                                              self.family_XY_mum_aff)
        inheritancefilter_220 = InheritanceFiltering(variants_per_gene_220,
                                                     self.family_XY_mum_aff,
                                                     self.genes_hemizygous,
                                                     None,
                                                     None, self.candidate_variants,
                                                      self.inhreport)
        inheritancefilter_220.inheritance_filter_genes()
        test_candidate_variants_220 = {'single_variants': {
            'X_1097183_A_GG': {
                'mode': {'hemizygous'},
                'variant':
                    variants_per_gene_220[
                        '1234'][
                        'X_1097183_A_GG'][
                        'child'],
                'hgncid': '1234'}}, 'compound_hets': {}}

        self.assertEqual(inheritancefilter_220.candidate_variants,
                         test_candidate_variants_220)
        # mum aff, mum 0/0, dad 1/1 pass
        variants_per_gene_202 = create_test_variants_per_gene(self.variants_202,
                                                              self.family_XY_mum_aff)
        inheritancefilter_202 = InheritanceFiltering(variants_per_gene_202,
                                                     self.family_XY_mum_aff,
                                                     self.genes_hemizygous,
                                                     None,
                                                     None, self.candidate_variants,
                                                      self.inhreport)
        inheritancefilter_202.inheritance_filter_genes()
        test_candidate_variants_202 = {'single_variants': {
            'X_1097183_A_GG': {
                'mode': {'hemizygous'},
                'variant':
                    variants_per_gene_202[
                        '1234'][
                        'X_1097183_A_GG'][
                        'child'],
                'hgncid': '1234'}}, 'compound_hets': {}}

        self.assertEqual(inheritancefilter_202.candidate_variants,
                         test_candidate_variants_202)
        # mum aff, mum 0/1, dad 1/1 pass
        variants_per_gene_212 = create_test_variants_per_gene(self.variants_212,
                                                              self.family_XY_mum_aff)
        inheritancefilter_212 = InheritanceFiltering(variants_per_gene_212,
                                                     self.family_XY_mum_aff,
                                                     self.genes_hemizygous,
                                                     None,
                                                     None, self.candidate_variants,
                                                      self.inhreport)
        inheritancefilter_212.inheritance_filter_genes()
        test_candidate_variants_212 = {'single_variants': {
            'X_1097183_A_GG': {
                'mode': {'hemizygous'},
                'variant':
                    variants_per_gene_212[
                        '1234'][
                        'X_1097183_A_GG'][
                        'child'],
                'hgncid': '1234'}}, 'compound_hets': {}}

        self.assertEqual(inheritancefilter_212.candidate_variants,
                         test_candidate_variants_212)
        # mum aff, mum 1/1, dad 1/1 pass
        variants_per_gene_222 = create_test_variants_per_gene(self.variants_222,
                                                              self.family_XY_mum_aff)
        inheritancefilter_222 = InheritanceFiltering(variants_per_gene_222,
                                                     self.family_XY_mum_aff,
                                                     self.genes_hemizygous,
                                                     None,
                                                     None, self.candidate_variants,
                                                      self.inhreport)
        inheritancefilter_222.inheritance_filter_genes()
        test_candidate_variants_222 = {'single_variants': {
            'X_1097183_A_GG': {
                'mode': {'hemizygous'},
                'variant':
                    variants_per_gene_222[
                        '1234'][
                        'X_1097183_A_GG'][
                        'child'],
                'hgncid': '1234'}}, 'compound_hets': {}}

        self.assertEqual(inheritancefilter_222.candidate_variants,
                         test_candidate_variants_222)

        # dad aff, mum 0/0, dad 0/0 pass
        variants_per_gene_200 = create_test_variants_per_gene(self.variants_200,
                                                              self.family_XY_dad_aff)
        inheritancefilter_200 = InheritanceFiltering(variants_per_gene_200,
                                                     self.family_XY_dad_aff,
                                                     self.genes_hemizygous,
                                                     None,
                                                     None, self.candidate_variants,
                                                      self.inhreport)
        inheritancefilter_200.inheritance_filter_genes()
        test_candidate_variants_200 = {'single_variants': {
            'X_1097183_A_GG': {
                'mode': {'hemizygous'},
                'variant':
                    variants_per_gene_200[
                        '1234'][
                        'X_1097183_A_GG'][
                        'child'],
                'hgncid': '1234'}}, 'compound_hets': {}}

        self.assertEqual(inheritancefilter_200.candidate_variants,
                         test_candidate_variants_200)
        # dad aff, mum 0/1, dad 0/0 pass
        variants_per_gene_210 = create_test_variants_per_gene(self.variants_210,
                                                              self.family_XY_dad_aff)
        inheritancefilter_210 = InheritanceFiltering(variants_per_gene_210,
                                                     self.family_XY_dad_aff,
                                                     self.genes_hemizygous,
                                                     None,
                                                     None, self.candidate_variants,
                                                      self.inhreport)
        inheritancefilter_210.inheritance_filter_genes()
        test_candidate_variants_210 = {'single_variants': {
            'X_1097183_A_GG': {
                'mode': {'hemizygous'},
                'variant':
                    variants_per_gene_210[
                        '1234'][
                        'X_1097183_A_GG'][
                        'child'],
                'hgncid': '1234'}}, 'compound_hets': {}}

        self.assertEqual(inheritancefilter_210.candidate_variants,
                         test_candidate_variants_210)
        # dad aff, mum 1/1, dad 0/0 fail
        variants_per_gene_220 = create_test_variants_per_gene(self.variants_220,
                                                              self.family_XY_dad_aff)
        inheritancefilter_220 = InheritanceFiltering(variants_per_gene_220,
                                                     self.family_XY_dad_aff,
                                                     self.genes_hemizygous,
                                                     None,
                                                     None, self.candidate_variants,
                                                      self.inhreport)
        inheritancefilter_220.inheritance_filter_genes()
        test_candidate_variants_220 = {'single_variants': {},
                                       'compound_hets': {}}

        self.assertEqual(inheritancefilter_220.candidate_variants,
                         test_candidate_variants_220)
        # dad aff, mum 0/0, dad 1/1 pass
        variants_per_gene_202 = create_test_variants_per_gene(self.variants_202,
                                                              self.family_XY_dad_aff)
        inheritancefilter_202 = InheritanceFiltering(variants_per_gene_202,
                                                     self.family_XY_dad_aff,
                                                     self.genes_hemizygous,
                                                     None,
                                                     None, self.candidate_variants,
                                                      self.inhreport)
        inheritancefilter_202.inheritance_filter_genes()
        test_candidate_variants_202 = {'single_variants': {
            'X_1097183_A_GG': {
                'mode': {'hemizygous'},
                'variant':
                    variants_per_gene_202[
                        '1234'][
                        'X_1097183_A_GG'][
                        'child'],
                'hgncid': '1234'}}, 'compound_hets': {}}

        self.assertEqual(inheritancefilter_202.candidate_variants,
                         test_candidate_variants_202)
        # dad aff, mum 0/1, dad 1/1 pass
        variants_per_gene_212 = create_test_variants_per_gene(self.variants_212,
                                                              self.family_XY_dad_aff)
        inheritancefilter_212 = InheritanceFiltering(variants_per_gene_212,
                                                     self.family_XY_dad_aff,
                                                     self.genes_hemizygous,
                                                     None,
                                                     None, self.candidate_variants,
                                                      self.inhreport)
        inheritancefilter_212.inheritance_filter_genes()
        test_candidate_variants_212 = {'single_variants': {
            'X_1097183_A_GG': {
                'mode': {'hemizygous'},
                'variant':
                    variants_per_gene_212[
                        '1234'][
                        'X_1097183_A_GG'][
                        'child'],
                'hgncid': '1234'}}, 'compound_hets': {}}

        self.assertEqual(inheritancefilter_212.candidate_variants,
                         test_candidate_variants_212)
        # dad aff, mum 1/1, dad 1/1 fail
        variants_per_gene_222 = create_test_variants_per_gene(self.variants_222,
                                                              self.family_XY_dad_aff)
        inheritancefilter_222 = InheritanceFiltering(variants_per_gene_222,
                                                     self.family_XY_dad_aff,
                                                     self.genes_hemizygous,
                                                     None,
                                                     None, self.candidate_variants,
                                                      self.inhreport)
        inheritancefilter_222.inheritance_filter_genes()
        test_candidate_variants_222 = {'single_variants': {},
                                       'compound_hets': {}}

        self.assertEqual(inheritancefilter_222.candidate_variants,
                         test_candidate_variants_222)

        # both unaff, mum 0/0, dad 0/0 pass
        variants_per_gene_200 = create_test_variants_per_gene(self.variants_200,
                                                              self.family_XY_both_unaff)
        inheritancefilter_200 = InheritanceFiltering(variants_per_gene_200,
                                                     self.family_XY_both_unaff,
                                                     self.genes_hemizygous,
                                                     None,
                                                     None, self.candidate_variants,
                                                      self.inhreport)
        inheritancefilter_200.inheritance_filter_genes()
        test_candidate_variants_200 = {'single_variants': {
            'X_1097183_A_GG': {
                'mode': {'hemizygous'},
                'variant':
                    variants_per_gene_200[
                        '1234'][
                        'X_1097183_A_GG'][
                        'child'],
                'hgncid': '1234'}}, 'compound_hets': {}}

        self.assertEqual(inheritancefilter_200.candidate_variants,
                         test_candidate_variants_200)
        # both unaff, mum 0/1, dad 0/0 pass
        variants_per_gene_210 = create_test_variants_per_gene(self.variants_210,
                                                              self.family_XY_both_unaff)
        inheritancefilter_210 = InheritanceFiltering(variants_per_gene_210,
                                                     self.family_XY_both_unaff,
                                                     self.genes_hemizygous,
                                                     None,
                                                     None, self.candidate_variants,
                                                      self.inhreport)
        inheritancefilter_210.inheritance_filter_genes()
        test_candidate_variants_210 = {'single_variants': {
            'X_1097183_A_GG': {
                'mode': {'hemizygous'},
                'variant':
                    variants_per_gene_210[
                        '1234'][
                        'X_1097183_A_GG'][
                        'child'],
                'hgncid': '1234'}}, 'compound_hets': {}}

        self.assertEqual(inheritancefilter_210.candidate_variants,
                         test_candidate_variants_210)
        # both unaff, mum 1/1, dad 0/0 fail
        variants_per_gene_220 = create_test_variants_per_gene(self.variants_220,
                                                              self.family_XY_both_unaff)
        inheritancefilter_220 = InheritanceFiltering(variants_per_gene_220,
                                                     self.family_XY_both_unaff,
                                                     self.genes_hemizygous,
                                                     None,
                                                     None, self.candidate_variants,
                                                      self.inhreport)
        inheritancefilter_220.inheritance_filter_genes()
        test_candidate_variants_220 = {'single_variants': {},
                                       'compound_hets': {}}

        self.assertEqual(inheritancefilter_220.candidate_variants,
                         test_candidate_variants_220)
        # both unaff, mum 0/0, dad 1/1 pass
        variants_per_gene_202 = create_test_variants_per_gene(self.variants_202,
                                                              self.family_XY_both_unaff)
        inheritancefilter_202 = InheritanceFiltering(variants_per_gene_202,
                                                     self.family_XY_both_unaff,
                                                     self.genes_hemizygous,
                                                     None,
                                                     None, self.candidate_variants,
                                                      self.inhreport)
        inheritancefilter_202.inheritance_filter_genes()
        test_candidate_variants_202 = {'single_variants': {
            'X_1097183_A_GG': {
                'mode': {'hemizygous'},
                'variant':
                    variants_per_gene_202[
                        '1234'][
                        'X_1097183_A_GG'][
                        'child'],
                'hgncid': '1234'}}, 'compound_hets': {}}

        self.assertEqual(inheritancefilter_202.candidate_variants,
                         test_candidate_variants_202)
        # both unaff, mum 0/1, dad 1/1 pass
        variants_per_gene_212 = create_test_variants_per_gene(self.variants_212,
                                                              self.family_XY_both_unaff)
        inheritancefilter_212 = InheritanceFiltering(variants_per_gene_212,
                                                     self.family_XY_both_unaff,
                                                     self.genes_hemizygous,
                                                     None,
                                                     None, self.candidate_variants,
                                                      self.inhreport)
        inheritancefilter_212.inheritance_filter_genes()
        test_candidate_variants_212 = {'single_variants': {
            'X_1097183_A_GG': {
                'mode': {'hemizygous'},
                'variant':
                    variants_per_gene_200[
                        '1234'][
                        'X_1097183_A_GG'][
                        'child'],
                'hgncid': '1234'}}, 'compound_hets': {}}

        self.assertEqual(inheritancefilter_212.candidate_variants,
                         test_candidate_variants_212)
        # both unaff, mum 1/1, dad 1/1 fail
        variants_per_gene_222 = create_test_variants_per_gene(self.variants_222,
                                                              self.family_XY_both_unaff)
        inheritancefilter_222 = InheritanceFiltering(variants_per_gene_222,
                                                     self.family_XY_both_unaff,
                                                     self.genes_hemizygous,
                                                     None,
                                                     None, self.candidate_variants,
                                                      self.inhreport)
        inheritancefilter_222.inheritance_filter_genes()
        test_candidate_variants_222 = {'single_variants': {},
                                       'compound_hets': {}}

        self.assertEqual(inheritancefilter_222.candidate_variants,
                         test_candidate_variants_222)

    def test_gn_X_linked_dominant_gt_hemizygous_parents_filter(self):
        # both aff, mum 0/0, dad 0/0 pass
        variants_per_gene_200 = create_test_variants_per_gene(self.variants_200,
                                                              self.family_XY_both_aff)
        inheritancefilter_200 = InheritanceFiltering(variants_per_gene_200,
                                                     self.family_XY_both_aff,
                                                     self.genes_X_linked_dominant,
                                                     None,
                                                     None, self.candidate_variants,
                                                      self.inhreport)
        inheritancefilter_200.inheritance_filter_genes()
        test_candidate_variants_200 = {'single_variants': {
            'X_1097183_A_GG': {
                'mode': {'X-linked dominant'},
                'variant':
                    variants_per_gene_200[
                        '1234'][
                        'X_1097183_A_GG'][
                        'child'],
                'hgncid': '1234'}}, 'compound_hets': {}}

        self.assertEqual(inheritancefilter_200.candidate_variants,
                         test_candidate_variants_200)
        # both aff, mum 0/1, dad 0/0 pass
        variants_per_gene_210 = create_test_variants_per_gene(self.variants_210,
                                                              self.family_XY_both_aff)
        inheritancefilter_210 = InheritanceFiltering(variants_per_gene_210,
                                                     self.family_XY_both_aff,
                                                     self.genes_X_linked_dominant,
                                                     None,
                                                     None, self.candidate_variants,
                                                      self.inhreport)
        inheritancefilter_210.inheritance_filter_genes()
        test_candidate_variants_210 = {'single_variants': {
            'X_1097183_A_GG': {
                'mode': {'X-linked dominant'},
                'variant':
                    variants_per_gene_210[
                        '1234'][
                        'X_1097183_A_GG'][
                        'child'],
                'hgncid': '1234'}}, 'compound_hets': {}}

        self.assertEqual(inheritancefilter_210.candidate_variants,
                         test_candidate_variants_210)
        # both aff, mum 1/1, dad 0/0 pass
        variants_per_gene_220 = create_test_variants_per_gene(self.variants_220,
                                                              self.family_XY_both_aff)
        inheritancefilter_220 = InheritanceFiltering(variants_per_gene_220,
                                                     self.family_XY_both_aff,
                                                     self.genes_X_linked_dominant,
                                                     None,
                                                     None, self.candidate_variants,
                                                      self.inhreport)
        inheritancefilter_220.inheritance_filter_genes()
        test_candidate_variants_220 = {'single_variants': {
            'X_1097183_A_GG': {
                'mode': {'X-linked dominant'},
                'variant':
                    variants_per_gene_220[
                        '1234'][
                        'X_1097183_A_GG'][
                        'child'],
                'hgncid': '1234'}}, 'compound_hets': {}}

        self.assertEqual(inheritancefilter_220.candidate_variants,
                         test_candidate_variants_220)
        # both aff, mum 0/0, dad 1/1 pass
        variants_per_gene_202 = create_test_variants_per_gene(self.variants_202,
                                                              self.family_XY_both_aff)
        inheritancefilter_202 = InheritanceFiltering(variants_per_gene_202,
                                                     self.family_XY_both_aff,
                                                     self.genes_X_linked_dominant,
                                                     None,
                                                     None, self.candidate_variants,
                                                      self.inhreport)
        inheritancefilter_202.inheritance_filter_genes()
        test_candidate_variants_202 = {'single_variants': {
            'X_1097183_A_GG': {
                'mode': {'X-linked dominant'},
                'variant':
                    variants_per_gene_202[
                        '1234'][
                        'X_1097183_A_GG'][
                        'child'],
                'hgncid': '1234'}}, 'compound_hets': {}}

        self.assertEqual(inheritancefilter_202.candidate_variants,
                         test_candidate_variants_202)
        # both aff, mum 0/1, dad 1/1 pass
        variants_per_gene_212 = create_test_variants_per_gene(self.variants_212,
                                                              self.family_XY_both_aff)
        inheritancefilter_212 = InheritanceFiltering(variants_per_gene_212,
                                                     self.family_XY_both_aff,
                                                     self.genes_X_linked_dominant,
                                                     None,
                                                     None, self.candidate_variants,
                                                      self.inhreport)
        inheritancefilter_212.inheritance_filter_genes()
        test_candidate_variants_212 = {'single_variants': {
            'X_1097183_A_GG': {
                'mode': {'X-linked dominant'},
                'variant':
                    variants_per_gene_212[
                        '1234'][
                        'X_1097183_A_GG'][
                        'child'],
                'hgncid': '1234'}}, 'compound_hets': {}}

        self.assertEqual(inheritancefilter_212.candidate_variants,
                         test_candidate_variants_212)
        # both aff, mum 1/1, dad 1/1 pass
        variants_per_gene_222 = create_test_variants_per_gene(self.variants_222,
                                                              self.family_XY_both_aff)
        inheritancefilter_222 = InheritanceFiltering(variants_per_gene_222,
                                                     self.family_XY_both_aff,
                                                     self.genes_X_linked_dominant,
                                                     None,
                                                     None, self.candidate_variants,
                                                      self.inhreport)
        inheritancefilter_222.inheritance_filter_genes()
        test_candidate_variants_222 = {'single_variants': {
            'X_1097183_A_GG': {
                'mode': {'X-linked dominant'},
                'variant':
                    variants_per_gene_222[
                        '1234'][
                        'X_1097183_A_GG'][
                        'child'],
                'hgncid': '1234'}}, 'compound_hets': {}}

        self.assertEqual(inheritancefilter_222.candidate_variants,
                         test_candidate_variants_222)

        # mum aff, mum 0/0, dad 0/0 pass
        variants_per_gene_200 = create_test_variants_per_gene(self.variants_200,
                                                              self.family_XY_mum_aff)
        inheritancefilter_200 = InheritanceFiltering(variants_per_gene_200,
                                                     self.family_XY_mum_aff,
                                                     self.genes_X_linked_dominant,
                                                     None,
                                                     None, self.candidate_variants,
                                                      self.inhreport)
        inheritancefilter_200.inheritance_filter_genes()
        test_candidate_variants_200 = {'single_variants': {
            'X_1097183_A_GG': {
                'mode': {'X-linked dominant'},
                'variant':
                    variants_per_gene_200[
                        '1234'][
                        'X_1097183_A_GG'][
                        'child'],
                'hgncid': '1234'}}, 'compound_hets': {}}

        self.assertEqual(inheritancefilter_200.candidate_variants,
                         test_candidate_variants_200)
        # mum aff, mum 0/1, dad 0/0 pass
        variants_per_gene_210 = create_test_variants_per_gene(self.variants_210,
                                                              self.family_XY_mum_aff)
        inheritancefilter_210 = InheritanceFiltering(variants_per_gene_210,
                                                     self.family_XY_mum_aff,
                                                     self.genes_X_linked_dominant,
                                                     None,
                                                     None, self.candidate_variants,
                                                      self.inhreport)
        inheritancefilter_210.inheritance_filter_genes()
        test_candidate_variants_210 = {'single_variants': {
            'X_1097183_A_GG': {
                'mode': {'X-linked dominant'},
                'variant':
                    variants_per_gene_210[
                        '1234'][
                        'X_1097183_A_GG'][
                        'child'],
                'hgncid': '1234'}}, 'compound_hets': {}}

        self.assertEqual(inheritancefilter_210.candidate_variants,
                         test_candidate_variants_210)
        # mum aff, mum 1/1, dad 0/0 pass
        variants_per_gene_220 = create_test_variants_per_gene(self.variants_220,
                                                              self.family_XY_mum_aff)
        inheritancefilter_220 = InheritanceFiltering(variants_per_gene_220,
                                                     self.family_XY_mum_aff,
                                                     self.genes_X_linked_dominant,
                                                     None,
                                                     None, self.candidate_variants,
                                                      self.inhreport)
        inheritancefilter_220.inheritance_filter_genes()
        test_candidate_variants_220 = {'single_variants': {
            'X_1097183_A_GG': {
                'mode': {'X-linked dominant'},
                'variant':
                    variants_per_gene_220[
                        '1234'][
                        'X_1097183_A_GG'][
                        'child'],
                'hgncid': '1234'}}, 'compound_hets': {}}

        self.assertEqual(inheritancefilter_220.candidate_variants,
                         test_candidate_variants_220)
        # mum aff, mum 0/0, dad 1/1 pass
        variants_per_gene_202 = create_test_variants_per_gene(self.variants_202,
                                                              self.family_XY_mum_aff)
        inheritancefilter_202 = InheritanceFiltering(variants_per_gene_202,
                                                     self.family_XY_mum_aff,
                                                     self.genes_X_linked_dominant,
                                                     None,
                                                     None, self.candidate_variants,
                                                      self.inhreport)
        inheritancefilter_202.inheritance_filter_genes()
        test_candidate_variants_202 = {'single_variants': {
            'X_1097183_A_GG': {
                'mode': {'X-linked dominant'},
                'variant':
                    variants_per_gene_202[
                        '1234'][
                        'X_1097183_A_GG'][
                        'child'],
                'hgncid': '1234'}}, 'compound_hets': {}}

        self.assertEqual(inheritancefilter_202.candidate_variants,
                         test_candidate_variants_202)
        # mum aff, mum 0/1, dad 1/1 pass
        variants_per_gene_212 = create_test_variants_per_gene(self.variants_212,
                                                              self.family_XY_mum_aff)
        inheritancefilter_212 = InheritanceFiltering(variants_per_gene_212,
                                                     self.family_XY_mum_aff,
                                                     self.genes_X_linked_dominant,
                                                     None,
                                                     None, self.candidate_variants,
                                                      self.inhreport)
        inheritancefilter_212.inheritance_filter_genes()
        test_candidate_variants_212 = {'single_variants': {
            'X_1097183_A_GG': {
                'mode': {'X-linked dominant'},
                'variant':
                    variants_per_gene_212[
                        '1234'][
                        'X_1097183_A_GG'][
                        'child'],
                'hgncid': '1234'}}, 'compound_hets': {}}

        self.assertEqual(inheritancefilter_212.candidate_variants,
                         test_candidate_variants_212)
        # mum aff, mum 1/1, dad 1/1 pass
        variants_per_gene_222 = create_test_variants_per_gene(self.variants_222,
                                                              self.family_XY_mum_aff)
        inheritancefilter_222 = InheritanceFiltering(variants_per_gene_222,
                                                     self.family_XY_mum_aff,
                                                     self.genes_X_linked_dominant,
                                                     None,
                                                     None, self.candidate_variants,
                                                      self.inhreport)
        inheritancefilter_222.inheritance_filter_genes()
        test_candidate_variants_222 = {'single_variants': {
            'X_1097183_A_GG': {
                'mode': {'X-linked dominant'},
                'variant':
                    variants_per_gene_222[
                        '1234'][
                        'X_1097183_A_GG'][
                        'child'],
                'hgncid': '1234'}}, 'compound_hets': {}}

        self.assertEqual(inheritancefilter_222.candidate_variants,
                         test_candidate_variants_222)

        # dad aff, mum 0/0, dad 0/0 pass
        variants_per_gene_200 = create_test_variants_per_gene(self.variants_200,
                                                              self.family_XY_dad_aff)
        inheritancefilter_200 = InheritanceFiltering(variants_per_gene_200,
                                                     self.family_XY_dad_aff,
                                                     self.genes_X_linked_dominant,
                                                     None,
                                                     None, self.candidate_variants,
                                                      self.inhreport)
        inheritancefilter_200.inheritance_filter_genes()
        test_candidate_variants_200 = {'single_variants': {
            'X_1097183_A_GG': {
                'mode': {'X-linked dominant'},
                'variant':
                    variants_per_gene_200[
                        '1234'][
                        'X_1097183_A_GG'][
                        'child'],
                'hgncid': '1234'}}, 'compound_hets': {}}

        self.assertEqual(inheritancefilter_200.candidate_variants,
                         test_candidate_variants_200)
        # dad aff, mum 0/1, dad 0/0 pass
        variants_per_gene_210 = create_test_variants_per_gene(self.variants_210,
                                                              self.family_XY_dad_aff)
        inheritancefilter_210 = InheritanceFiltering(variants_per_gene_210,
                                                     self.family_XY_dad_aff,
                                                     self.genes_X_linked_dominant,
                                                     None,
                                                     None, self.candidate_variants,
                                                      self.inhreport)
        inheritancefilter_210.inheritance_filter_genes()
        test_candidate_variants_210 = {'single_variants': {
            'X_1097183_A_GG': {
                'mode': {'X-linked dominant'},
                'variant':
                    variants_per_gene_210[
                        '1234'][
                        'X_1097183_A_GG'][
                        'child'],
                'hgncid': '1234'}}, 'compound_hets': {}}

        self.assertEqual(inheritancefilter_210.candidate_variants,
                         test_candidate_variants_210)
        # dad aff, mum 1/1, dad 0/0 fail
        variants_per_gene_220 = create_test_variants_per_gene(self.variants_220,
                                                              self.family_XY_dad_aff)
        inheritancefilter_220 = InheritanceFiltering(variants_per_gene_220,
                                                     self.family_XY_dad_aff,
                                                     self.genes_X_linked_dominant,
                                                     None,
                                                     None, self.candidate_variants,
                                                      self.inhreport)
        inheritancefilter_220.inheritance_filter_genes()
        test_candidate_variants_220 = {'single_variants': {},
                                       'compound_hets': {}}

        self.assertEqual(inheritancefilter_220.candidate_variants,
                         test_candidate_variants_220)
        # dad aff, mum 0/0, dad 1/1 pass
        variants_per_gene_202 = create_test_variants_per_gene(self.variants_202,
                                                              self.family_XY_dad_aff)
        inheritancefilter_202 = InheritanceFiltering(variants_per_gene_202,
                                                     self.family_XY_dad_aff,
                                                     self.genes_X_linked_dominant,
                                                     None,
                                                     None, self.candidate_variants,
                                                      self.inhreport)
        inheritancefilter_202.inheritance_filter_genes()
        test_candidate_variants_202 = {'single_variants': {
            'X_1097183_A_GG': {
                'mode': {'X-linked dominant'},
                'variant':
                    variants_per_gene_202[
                        '1234'][
                        'X_1097183_A_GG'][
                        'child'],
                'hgncid': '1234'}}, 'compound_hets': {}}

        self.assertEqual(inheritancefilter_202.candidate_variants,
                         test_candidate_variants_202)
        # dad aff, mum 0/1, dad 1/1 pass
        variants_per_gene_212 = create_test_variants_per_gene(self.variants_212,
                                                              self.family_XY_dad_aff)
        inheritancefilter_212 = InheritanceFiltering(variants_per_gene_212,
                                                     self.family_XY_dad_aff,
                                                     self.genes_X_linked_dominant,
                                                     None,
                                                     None, self.candidate_variants,
                                                      self.inhreport)
        inheritancefilter_212.inheritance_filter_genes()
        test_candidate_variants_212 = {'single_variants': {
            'X_1097183_A_GG': {
                'mode': {'X-linked dominant'},
                'variant':
                    variants_per_gene_212[
                        '1234'][
                        'X_1097183_A_GG'][
                        'child'],
                'hgncid': '1234'}}, 'compound_hets': {}}

        self.assertEqual(inheritancefilter_212.candidate_variants,
                         test_candidate_variants_212)
        # dad aff, mum 1/1, dad 1/1 fail
        variants_per_gene_222 = create_test_variants_per_gene(self.variants_222,
                                                              self.family_XY_dad_aff)
        inheritancefilter_222 = InheritanceFiltering(variants_per_gene_222,
                                                     self.family_XY_dad_aff,
                                                     self.genes_X_linked_dominant,
                                                     None,
                                                     None, self.candidate_variants,
                                                      self.inhreport)
        inheritancefilter_222.inheritance_filter_genes()
        test_candidate_variants_222 = {'single_variants': {},
                                       'compound_hets': {}}

        self.assertEqual(inheritancefilter_222.candidate_variants,
                         test_candidate_variants_222)

        # both unaff, mum 0/0, dad 0/0 pass
        variants_per_gene_200 = create_test_variants_per_gene(self.variants_200,
                                                              self.family_XY_both_unaff)
        inheritancefilter_200 = InheritanceFiltering(variants_per_gene_200,
                                                     self.family_XY_both_unaff,
                                                     self.genes_X_linked_dominant,
                                                     None,
                                                     None, self.candidate_variants,
                                                      self.inhreport)
        inheritancefilter_200.inheritance_filter_genes()
        test_candidate_variants_200 = {'single_variants': {
            'X_1097183_A_GG': {
                'mode': {'X-linked dominant'},
                'variant':
                    variants_per_gene_200[
                        '1234'][
                        'X_1097183_A_GG'][
                        'child'],
                'hgncid': '1234'}}, 'compound_hets': {}}

        self.assertEqual(inheritancefilter_200.candidate_variants,
                         test_candidate_variants_200)
        # both unaff, mum 0/1, dad 0/0 pass
        variants_per_gene_210 = create_test_variants_per_gene(self.variants_210,
                                                              self.family_XY_both_unaff)
        inheritancefilter_210 = InheritanceFiltering(variants_per_gene_210,
                                                     self.family_XY_both_unaff,
                                                     self.genes_X_linked_dominant,
                                                     None,
                                                     None, self.candidate_variants,
                                                      self.inhreport)
        inheritancefilter_210.inheritance_filter_genes()
        test_candidate_variants_210 = {'single_variants': {
            'X_1097183_A_GG': {
                'mode': {'X-linked dominant'},
                'variant':
                    variants_per_gene_210[
                        '1234'][
                        'X_1097183_A_GG'][
                        'child'],
                'hgncid': '1234'}}, 'compound_hets': {}}

        self.assertEqual(inheritancefilter_210.candidate_variants,
                         test_candidate_variants_210)
        # both unaff, mum 1/1, dad 0/0 fail
        variants_per_gene_220 = create_test_variants_per_gene(self.variants_220,
                                                              self.family_XY_both_unaff)
        inheritancefilter_220 = InheritanceFiltering(variants_per_gene_220,
                                                     self.family_XY_both_unaff,
                                                     self.genes_X_linked_dominant,
                                                     None,
                                                     None, self.candidate_variants,
                                                      self.inhreport)
        inheritancefilter_220.inheritance_filter_genes()
        test_candidate_variants_220 = {'single_variants': {},
                                       'compound_hets': {}}

        self.assertEqual(inheritancefilter_220.candidate_variants,
                         test_candidate_variants_220)
        # both unaff, mum 0/0, dad 1/1 pass
        variants_per_gene_202 = create_test_variants_per_gene(self.variants_202,
                                                              self.family_XY_both_unaff)
        inheritancefilter_202 = InheritanceFiltering(variants_per_gene_202,
                                                     self.family_XY_both_unaff,
                                                     self.genes_X_linked_dominant,
                                                     None,
                                                     None, self.candidate_variants,
                                                      self.inhreport)
        inheritancefilter_202.inheritance_filter_genes()
        test_candidate_variants_202 = {'single_variants': {
            'X_1097183_A_GG': {
                'mode': {'X-linked dominant'},
                'variant':
                    variants_per_gene_202[
                        '1234'][
                        'X_1097183_A_GG'][
                        'child'],
                'hgncid': '1234'}}, 'compound_hets': {}}

        self.assertEqual(inheritancefilter_202.candidate_variants,
                         test_candidate_variants_202)
        # both unaff, mum 0/1, dad 1/1 pass
        variants_per_gene_212 = create_test_variants_per_gene(self.variants_212,
                                                              self.family_XY_both_unaff)
        inheritancefilter_212 = InheritanceFiltering(variants_per_gene_212,
                                                     self.family_XY_both_unaff,
                                                     self.genes_X_linked_dominant,
                                                     None,
                                                     None, self.candidate_variants,
                                                      self.inhreport)
        inheritancefilter_212.inheritance_filter_genes()
        test_candidate_variants_212 = {'single_variants': {
            'X_1097183_A_GG': {
                'mode': {'X-linked dominant'},
                'variant':
                    variants_per_gene_212[
                        '1234'][
                        'X_1097183_A_GG'][
                        'child'],
                'hgncid': '1234'}}, 'compound_hets': {}}

        self.assertEqual(inheritancefilter_212.candidate_variants,
                         test_candidate_variants_212)
        # both unaff, mum 1/1, dad 1/1 fail
        variants_per_gene_222 = create_test_variants_per_gene(self.variants_222,
                                                              self.family_XY_both_unaff)
        inheritancefilter_222 = InheritanceFiltering(variants_per_gene_222,
                                                     self.family_XY_both_unaff,
                                                     self.genes_X_linked_dominant,
                                                     None,
                                                     None, self.candidate_variants,
                                                      self.inhreport)
        inheritancefilter_222.inheritance_filter_genes()
        test_candidate_variants_222 = {'single_variants': {},
                                       'compound_hets': {}}

        self.assertEqual(inheritancefilter_222.candidate_variants,
                         test_candidate_variants_222)

    def test_gn_X_linked_over_dominant_gt_hemizygous_parents_filter(self):
        # both aff, mum 0/0, dad 0/0 fail
        variants_per_gene_200 = create_test_variants_per_gene(self.variants_200,
                                                              self.family_XY_both_aff)
        inheritancefilter_200 = InheritanceFiltering(variants_per_gene_200,
                                                     self.family_XY_both_aff,
                                                     self.genes_X_linked_over_dominance,
                                                     None,
                                                     None, self.candidate_variants,
                                                      self.inhreport)
        inheritancefilter_200.inheritance_filter_genes()
        test_candidate_variants_200 = {'single_variants': {},
                                       'compound_hets': {}}

        self.assertEqual(inheritancefilter_200.candidate_variants,
                         test_candidate_variants_200)
        # both aff, mum 0/1, dad 0/0 fail
        variants_per_gene_210 = create_test_variants_per_gene(self.variants_210,
                                                              self.family_XY_both_aff)
        inheritancefilter_210 = InheritanceFiltering(variants_per_gene_210,
                                                     self.family_XY_both_aff,
                                                     self.genes_X_linked_over_dominance,
                                                     None,
                                                     None, self.candidate_variants,
                                                      self.inhreport)
        inheritancefilter_210.inheritance_filter_genes()
        test_candidate_variants_210 = {'single_variants': {},
                                       'compound_hets': {}}

        self.assertEqual(inheritancefilter_210.candidate_variants,
                         test_candidate_variants_210)
        # both aff, mum 1/1, dad 0/0 fail
        variants_per_gene_220 = create_test_variants_per_gene(self.variants_220,
                                                              self.family_XY_both_aff)
        inheritancefilter_220 = InheritanceFiltering(variants_per_gene_220,
                                                     self.family_XY_both_aff,
                                                     self.genes_X_linked_over_dominance,
                                                     None,
                                                     None, self.candidate_variants,
                                                      self.inhreport)
        inheritancefilter_220.inheritance_filter_genes()
        test_candidate_variants_220 = {'single_variants': {},
                                       'compound_hets': {}}

        self.assertEqual(inheritancefilter_220.candidate_variants,
                         test_candidate_variants_220)
        # both aff, mum 0/0, dad 1/1 fail
        variants_per_gene_202 = create_test_variants_per_gene(self.variants_202,
                                                              self.family_XY_both_aff)
        inheritancefilter_202 = InheritanceFiltering(variants_per_gene_202,
                                                     self.family_XY_both_aff,
                                                     self.genes_X_linked_over_dominance,
                                                     None,
                                                     None, self.candidate_variants,
                                                      self.inhreport)
        inheritancefilter_202.inheritance_filter_genes()
        test_candidate_variants_202 = {'single_variants': {},
                                       'compound_hets': {}}

        self.assertEqual(inheritancefilter_202.candidate_variants,
                         test_candidate_variants_202)
        # both aff, mum 0/1, dad 1/1 fail
        variants_per_gene_212 = create_test_variants_per_gene(self.variants_212,
                                                              self.family_XY_both_aff)
        inheritancefilter_212 = InheritanceFiltering(variants_per_gene_212,
                                                     self.family_XY_both_aff,
                                                     self.genes_X_linked_over_dominance,
                                                     None,
                                                     None, self.candidate_variants,
                                                      self.inhreport)
        inheritancefilter_212.inheritance_filter_genes()
        test_candidate_variants_212 = {'single_variants': {},
                                       'compound_hets': {}}

        self.assertEqual(inheritancefilter_212.candidate_variants,
                         test_candidate_variants_212)
        # both aff, mum 1/1, dad 1/1 fail
        variants_per_gene_222 = create_test_variants_per_gene(self.variants_222,
                                                              self.family_XY_both_aff)
        inheritancefilter_222 = InheritanceFiltering(variants_per_gene_222,
                                                     self.family_XY_both_aff,
                                                     self.genes_X_linked_over_dominance,
                                                     None,
                                                     None, self.candidate_variants,
                                                      self.inhreport)
        inheritancefilter_222.inheritance_filter_genes()
        test_candidate_variants_222 = {'single_variants': {},
                                       'compound_hets': {}}

        self.assertEqual(inheritancefilter_222.candidate_variants,
                         test_candidate_variants_222)

        # mum aff, mum 0/0, dad 0/0 fail
        variants_per_gene_200 = create_test_variants_per_gene(self.variants_200,
                                                              self.family_XY_mum_aff)
        inheritancefilter_200 = InheritanceFiltering(variants_per_gene_200,
                                                     self.family_XY_mum_aff,
                                                     self.genes_X_linked_over_dominance,
                                                     None,
                                                     None, self.candidate_variants,
                                                      self.inhreport)
        inheritancefilter_200.inheritance_filter_genes()
        test_candidate_variants_200 = {'single_variants': {},
                                       'compound_hets': {}}

        self.assertEqual(inheritancefilter_200.candidate_variants,
                         test_candidate_variants_200)
        # mum aff, mum 0/1, dad 0/0 fail
        variants_per_gene_210 = create_test_variants_per_gene(self.variants_210,
                                                              self.family_XY_mum_aff)
        inheritancefilter_210 = InheritanceFiltering(variants_per_gene_210,
                                                     self.family_XY_mum_aff,
                                                     self.genes_X_linked_over_dominance,
                                                     None,
                                                     None, self.candidate_variants,
                                                      self.inhreport)
        inheritancefilter_210.inheritance_filter_genes()
        test_candidate_variants_210 = {'single_variants': {},
                                       'compound_hets': {}}

        self.assertEqual(inheritancefilter_210.candidate_variants,
                         test_candidate_variants_210)
        # mum aff, mum 1/1, dad 0/0 fail
        variants_per_gene_220 = create_test_variants_per_gene(self.variants_220,
                                                              self.family_XY_mum_aff)
        inheritancefilter_220 = InheritanceFiltering(variants_per_gene_220,
                                                     self.family_XY_mum_aff,
                                                     self.genes_X_linked_over_dominance,
                                                     None,
                                                     None, self.candidate_variants,
                                                      self.inhreport)
        inheritancefilter_220.inheritance_filter_genes()
        test_candidate_variants_220 = {'single_variants': {},
                                       'compound_hets': {}}

        self.assertEqual(inheritancefilter_220.candidate_variants,
                         test_candidate_variants_220)
        # mum aff, mum 0/0, dad 1/1 fail
        variants_per_gene_202 = create_test_variants_per_gene(self.variants_202,
                                                              self.family_XY_mum_aff)
        inheritancefilter_202 = InheritanceFiltering(variants_per_gene_202,
                                                     self.family_XY_mum_aff,
                                                     self.genes_X_linked_over_dominance,
                                                     None,
                                                     None, self.candidate_variants,
                                                      self.inhreport)
        inheritancefilter_202.inheritance_filter_genes()
        test_candidate_variants_202 = {'single_variants': {},
                                       'compound_hets': {}}

        self.assertEqual(inheritancefilter_202.candidate_variants,
                         test_candidate_variants_202)
        # mum aff, mum 0/1, dad 1/1 fail
        variants_per_gene_212 = create_test_variants_per_gene(self.variants_212,
                                                              self.family_XY_mum_aff)
        inheritancefilter_212 = InheritanceFiltering(variants_per_gene_212,
                                                     self.family_XY_mum_aff,
                                                     self.genes_X_linked_over_dominance,
                                                     None,
                                                     None, self.candidate_variants,
                                                      self.inhreport)
        inheritancefilter_212.inheritance_filter_genes()
        test_candidate_variants_212 = {'single_variants': {},
                                       'compound_hets': {}}

        self.assertEqual(inheritancefilter_212.candidate_variants,
                         test_candidate_variants_212)
        # mum aff, mum 1/1, dad 1/1 fail
        variants_per_gene_222 = create_test_variants_per_gene(self.variants_222,
                                                              self.family_XY_mum_aff)
        inheritancefilter_222 = InheritanceFiltering(variants_per_gene_222,
                                                     self.family_XY_mum_aff,
                                                     self.genes_X_linked_over_dominance,
                                                     None,
                                                     None, self.candidate_variants,
                                                      self.inhreport)
        inheritancefilter_222.inheritance_filter_genes()
        test_candidate_variants_222 = {'single_variants': {},
                                       'compound_hets': {}}

        self.assertEqual(inheritancefilter_222.candidate_variants,
                         test_candidate_variants_222)

        # dad aff, mum 0/0, dad 0/0 fail
        variants_per_gene_200 = create_test_variants_per_gene(self.variants_200,
                                                              self.family_XY_dad_aff)
        inheritancefilter_200 = InheritanceFiltering(variants_per_gene_200,
                                                     self.family_XY_dad_aff,
                                                     self.genes_X_linked_over_dominance,
                                                     None,
                                                     None, self.candidate_variants,
                                                      self.inhreport)
        inheritancefilter_200.inheritance_filter_genes()
        test_candidate_variants_200 = {'single_variants': {},
                                       'compound_hets': {}}

        self.assertEqual(inheritancefilter_200.candidate_variants,
                         test_candidate_variants_200)
        # dad aff, mum 0/1, dad 0/0 fail
        variants_per_gene_210 = create_test_variants_per_gene(self.variants_210,
                                                              self.family_XY_dad_aff)
        inheritancefilter_210 = InheritanceFiltering(variants_per_gene_210,
                                                     self.family_XY_dad_aff,
                                                     self.genes_X_linked_over_dominance,
                                                     None,
                                                     None, self.candidate_variants,
                                                      self.inhreport)
        inheritancefilter_210.inheritance_filter_genes()
        test_candidate_variants_210 = {'single_variants': {},
                                       'compound_hets': {}}

        self.assertEqual(inheritancefilter_210.candidate_variants,
                         test_candidate_variants_210)
        # dad aff, mum 1/1, dad 0/0 fail
        variants_per_gene_220 = create_test_variants_per_gene(self.variants_220,
                                                              self.family_XY_dad_aff)
        inheritancefilter_220 = InheritanceFiltering(variants_per_gene_220,
                                                     self.family_XY_dad_aff,
                                                     self.genes_X_linked_over_dominance,
                                                     None,
                                                     None, self.candidate_variants,
                                                      self.inhreport)
        inheritancefilter_220.inheritance_filter_genes()
        test_candidate_variants_220 = {'single_variants': {},
                                       'compound_hets': {}}

        self.assertEqual(inheritancefilter_220.candidate_variants,
                         test_candidate_variants_220)
        # dad aff, mum 0/0, dad 1/1 fail
        variants_per_gene_202 = create_test_variants_per_gene(self.variants_202,
                                                              self.family_XY_dad_aff)
        inheritancefilter_202 = InheritanceFiltering(variants_per_gene_202,
                                                     self.family_XY_dad_aff,
                                                     self.genes_X_linked_over_dominance,
                                                     None,
                                                     None, self.candidate_variants,
                                                      self.inhreport)
        inheritancefilter_202.inheritance_filter_genes()
        test_candidate_variants_202 = {'single_variants': {},
                                       'compound_hets': {}}

        self.assertEqual(inheritancefilter_202.candidate_variants,
                         test_candidate_variants_202)
        # dad aff, mum 0/1, dad 1/1 fail
        variants_per_gene_212 = create_test_variants_per_gene(self.variants_212,
                                                              self.family_XY_dad_aff)
        inheritancefilter_212 = InheritanceFiltering(variants_per_gene_212,
                                                     self.family_XY_dad_aff,
                                                     self.genes_X_linked_over_dominance,
                                                     None,
                                                     None, self.candidate_variants,
                                                      self.inhreport)
        inheritancefilter_212.inheritance_filter_genes()
        test_candidate_variants_212 = {'single_variants': {},
                                       'compound_hets': {}}

        self.assertEqual(inheritancefilter_212.candidate_variants,
                         test_candidate_variants_212)
        # dad aff, mum 1/1, dad 1/1 fail
        variants_per_gene_222 = create_test_variants_per_gene(self.variants_222,
                                                              self.family_XY_dad_aff)
        inheritancefilter_222 = InheritanceFiltering(variants_per_gene_222,
                                                     self.family_XY_dad_aff,
                                                     self.genes_X_linked_over_dominance,
                                                     None,
                                                     None, self.candidate_variants,
                                                      self.inhreport)
        inheritancefilter_222.inheritance_filter_genes()
        test_candidate_variants_222 = {'single_variants': {},
                                       'compound_hets': {}}

        self.assertEqual(inheritancefilter_222.candidate_variants,
                         test_candidate_variants_222)

        # both unaff, mum 0/0, dad 0/0 fail
        variants_per_gene_200 = create_test_variants_per_gene(self.variants_200,
                                                              self.family_XY_both_unaff)
        inheritancefilter_200 = InheritanceFiltering(variants_per_gene_200,
                                                     self.family_XY_both_unaff,
                                                     self.genes_X_linked_over_dominance,
                                                     None,
                                                     None, self.candidate_variants,
                                                      self.inhreport)
        inheritancefilter_200.inheritance_filter_genes()
        test_candidate_variants_200 = {'single_variants': {},
                                       'compound_hets': {}}

        self.assertEqual(inheritancefilter_200.candidate_variants,
                         test_candidate_variants_200)
        # both unaff, mum 0/1, dad 0/0 fail
        variants_per_gene_210 = create_test_variants_per_gene(self.variants_210,
                                                              self.family_XY_both_unaff)
        inheritancefilter_210 = InheritanceFiltering(variants_per_gene_210,
                                                     self.family_XY_both_unaff,
                                                     self.genes_X_linked_over_dominance,
                                                     None,
                                                     None, self.candidate_variants,
                                                      self.inhreport)
        inheritancefilter_210.inheritance_filter_genes()
        test_candidate_variants_210 = {'single_variants': {},
                                       'compound_hets': {}}

        self.assertEqual(inheritancefilter_210.candidate_variants,
                         test_candidate_variants_210)
        # both unaff, mum 1/1, dad 0/0 fail
        variants_per_gene_220 = create_test_variants_per_gene(self.variants_220,
                                                              self.family_XY_both_unaff)
        inheritancefilter_220 = InheritanceFiltering(variants_per_gene_220,
                                                     self.family_XY_both_unaff,
                                                     self.genes_X_linked_over_dominance,
                                                     None,
                                                     None, self.candidate_variants,
                                                      self.inhreport)
        inheritancefilter_220.inheritance_filter_genes()
        test_candidate_variants_220 = {'single_variants': {},
                                       'compound_hets': {}}

        self.assertEqual(inheritancefilter_220.candidate_variants,
                         test_candidate_variants_220)
        # both unaff, mum 0/0, dad 1/1 fail
        variants_per_gene_202 = create_test_variants_per_gene(self.variants_202,
                                                              self.family_XY_both_unaff)
        inheritancefilter_202 = InheritanceFiltering(variants_per_gene_202,
                                                     self.family_XY_both_unaff,
                                                     self.genes_X_linked_over_dominance,
                                                     None,
                                                     None, self.candidate_variants,
                                                      self.inhreport)
        inheritancefilter_202.inheritance_filter_genes()
        test_candidate_variants_202 = {'single_variants': {},
                                       'compound_hets': {}}

        self.assertEqual(inheritancefilter_202.candidate_variants,
                         test_candidate_variants_202)
        # both unaff, mum 0/1, dad 1/1 fail
        variants_per_gene_212 = create_test_variants_per_gene(self.variants_212,
                                                              self.family_XY_both_unaff)
        inheritancefilter_212 = InheritanceFiltering(variants_per_gene_212,
                                                     self.family_XY_both_unaff,
                                                     self.genes_X_linked_over_dominance,
                                                     None,
                                                     None, self.candidate_variants,
                                                      self.inhreport)
        inheritancefilter_212.inheritance_filter_genes()
        test_candidate_variants_212 = {'single_variants': {},
                                       'compound_hets': {}}

        self.assertEqual(inheritancefilter_212.candidate_variants,
                         test_candidate_variants_212)
        # both unaff, mum 1/1, dad 1/1 fail
        variants_per_gene_222 = create_test_variants_per_gene(self.variants_222,
                                                              self.family_XY_both_unaff)
        inheritancefilter_222 = InheritanceFiltering(variants_per_gene_222,
                                                     self.family_XY_both_unaff,
                                                     self.genes_X_linked_over_dominance,
                                                     None,
                                                     None, self.candidate_variants,
                                                      self.inhreport)
        inheritancefilter_222.inheritance_filter_genes()
        test_candidate_variants_222 = {'single_variants': {},
                                       'compound_hets': {}}

        self.assertEqual(inheritancefilter_222.candidate_variants,
                         test_candidate_variants_222)

    def test_gn_hemizygous_gt_heterozygous_parents_filter(self):
        # both aff, mum 0/0, dad 0/0 pass
        variants_per_gene_100 = create_test_variants_per_gene(self.variants_100,
                                                              self.family_XX_both_aff)
        inheritancefilter_100 = InheritanceFiltering(variants_per_gene_100,
                                                     self.family_XX_both_aff,
                                                     self.genes_hemizygous,
                                                     None,
                                                     None, self.candidate_variants,
                                                      self.inhreport)
        inheritancefilter_100.inheritance_filter_genes()
        test_candidate_variants_100 = {'single_variants': {
            'X_1097183_A_GG': {
                'mode': {'hemizygous'},
                'variant':
                    variants_per_gene_100[
                        '1234'][
                        'X_1097183_A_GG'][
                        'child'],
                'hgncid': '1234'}}, 'compound_hets': {}}

        self.assertEqual(inheritancefilter_100.candidate_variants,
                         test_candidate_variants_100)
        # both aff, mum 0/1, dad 0/0 pass
        variants_per_gene_110 = create_test_variants_per_gene(self.variants_110,
                                                              self.family_XX_both_aff)
        inheritancefilter_110 = InheritanceFiltering(variants_per_gene_110,
                                                     self.family_XX_both_aff,
                                                     self.genes_hemizygous,
                                                     None,
                                                     None, self.candidate_variants,
                                                      self.inhreport)
        inheritancefilter_110.inheritance_filter_genes()
        test_candidate_variants_110 = {'single_variants': {
            'X_1097183_A_GG': {
                'mode': {'hemizygous'},
                'variant':
                    variants_per_gene_110[
                        '1234'][
                        'X_1097183_A_GG'][
                        'child'],
                'hgncid': '1234'}}, 'compound_hets': {}}

        self.assertEqual(inheritancefilter_110.candidate_variants,
                         test_candidate_variants_110)
        # both aff, mum 1/1, dad 0/0 pass
        variants_per_gene_120 = create_test_variants_per_gene(self.variants_120,
                                                              self.family_XX_both_aff)
        inheritancefilter_120 = InheritanceFiltering(variants_per_gene_120,
                                                     self.family_XX_both_aff,
                                                     self.genes_hemizygous,
                                                     None,
                                                     None, self.candidate_variants,
                                                      self.inhreport)
        inheritancefilter_120.inheritance_filter_genes()
        test_candidate_variants_120 = {'single_variants': {
            'X_1097183_A_GG': {
                'mode': {'hemizygous'},
                'variant':
                    variants_per_gene_120[
                        '1234'][
                        'X_1097183_A_GG'][
                        'child'],
                'hgncid': '1234'}}, 'compound_hets': {}}

        self.assertEqual(inheritancefilter_120.candidate_variants,
                         test_candidate_variants_120)
        # both aff, mum 0/0, dad 1/1 pass
        variants_per_gene_102 = create_test_variants_per_gene(self.variants_102,
                                                              self.family_XX_both_aff)
        inheritancefilter_102 = InheritanceFiltering(variants_per_gene_102,
                                                     self.family_XX_both_aff,
                                                     self.genes_hemizygous,
                                                     None,
                                                     None, self.candidate_variants,
                                                      self.inhreport)
        inheritancefilter_102.inheritance_filter_genes()
        test_candidate_variants_102 = {'single_variants': {
            'X_1097183_A_GG': {
                'mode': {'hemizygous'},
                'variant':
                    variants_per_gene_102[
                        '1234'][
                        'X_1097183_A_GG'][
                        'child'],
                'hgncid': '1234'}}, 'compound_hets': {}}

        self.assertEqual(inheritancefilter_102.candidate_variants,
                         test_candidate_variants_102)
        # both aff, mum 0/1, dad 1/1 pass
        variants_per_gene_112 = create_test_variants_per_gene(self.variants_112,
                                                              self.family_XX_both_aff)
        inheritancefilter_112 = InheritanceFiltering(variants_per_gene_112,
                                                     self.family_XX_both_aff,
                                                     self.genes_hemizygous,
                                                     None,
                                                     None, self.candidate_variants,
                                                      self.inhreport)
        inheritancefilter_112.inheritance_filter_genes()
        test_candidate_variants_112 = {'single_variants': {
            'X_1097183_A_GG': {
                'mode': {'hemizygous'},
                'variant':
                    variants_per_gene_112[
                        '1234'][
                        'X_1097183_A_GG'][
                        'child'],
                'hgncid': '1234'}}, 'compound_hets': {}}

        self.assertEqual(inheritancefilter_112.candidate_variants,
                         test_candidate_variants_112)
        # both aff, mum 1/1, dad 1/1 fail
        variants_per_gene_122 = create_test_variants_per_gene(self.variants_122,
                                                              self.family_XX_both_aff)
        inheritancefilter_122 = InheritanceFiltering(variants_per_gene_122,
                                                     self.family_XX_both_aff,
                                                     self.genes_hemizygous,
                                                     None,
                                                     None, self.candidate_variants,
                                                      self.inhreport)
        inheritancefilter_122.inheritance_filter_genes()
        test_candidate_variants_122 = {'single_variants': {},
                                       'compound_hets': {}}

        self.assertEqual(inheritancefilter_122.candidate_variants,
                         test_candidate_variants_122)

        # mum aff, mum 0/0, dad 0/0 pass
        variants_per_gene_100 = create_test_variants_per_gene(self.variants_100,
                                                              self.family_XX_mum_aff)
        inheritancefilter_100 = InheritanceFiltering(variants_per_gene_100,
                                                     self.family_XX_mum_aff,
                                                     self.genes_hemizygous,
                                                     None,
                                                     None, self.candidate_variants,
                                                      self.inhreport)
        inheritancefilter_100.inheritance_filter_genes()
        test_candidate_variants_100 = {'single_variants': {
            'X_1097183_A_GG': {
                'mode': {'hemizygous'},
                'variant':
                    variants_per_gene_100[
                        '1234'][
                        'X_1097183_A_GG'][
                        'child'],
                'hgncid': '1234'}}, 'compound_hets': {}}

        self.assertEqual(inheritancefilter_100.candidate_variants,
                         test_candidate_variants_100)
        # mum aff, mum 0/1, dad 0/0 pass
        variants_per_gene_110 = create_test_variants_per_gene(self.variants_110,
                                                              self.family_XX_mum_aff)
        inheritancefilter_110 = InheritanceFiltering(variants_per_gene_110,
                                                     self.family_XX_mum_aff,
                                                     self.genes_hemizygous,
                                                     None,
                                                     None, self.candidate_variants,
                                                      self.inhreport)
        inheritancefilter_110.inheritance_filter_genes()
        test_candidate_variants_110 = {'single_variants': {
            'X_1097183_A_GG': {
                'mode': {'hemizygous'},
                'variant':
                    variants_per_gene_110[
                        '1234'][
                        'X_1097183_A_GG'][
                        'child'],
                'hgncid': '1234'}}, 'compound_hets': {}}

        self.assertEqual(inheritancefilter_110.candidate_variants,
                         test_candidate_variants_110)
        # mum aff, mum 1/1, dad 0/0 pass
        variants_per_gene_120 = create_test_variants_per_gene(self.variants_120,
                                                              self.family_XX_mum_aff)
        inheritancefilter_120 = InheritanceFiltering(variants_per_gene_120,
                                                     self.family_XX_mum_aff,
                                                     self.genes_hemizygous,
                                                     None,
                                                     None, self.candidate_variants,
                                                      self.inhreport)
        inheritancefilter_120.inheritance_filter_genes()
        test_candidate_variants_120 = {'single_variants': {
            'X_1097183_A_GG': {
                'mode': {'hemizygous'},
                'variant':
                    variants_per_gene_120[
                        '1234'][
                        'X_1097183_A_GG'][
                        'child'],
                'hgncid': '1234'}}, 'compound_hets': {}}

        self.assertEqual(inheritancefilter_120.candidate_variants,
                         test_candidate_variants_120)
        # mum aff, mum 0/0, dad 1/1 fail
        variants_per_gene_102 = create_test_variants_per_gene(self.variants_102,
                                                              self.family_XX_mum_aff)
        inheritancefilter_102 = InheritanceFiltering(variants_per_gene_102,
                                                     self.family_XX_mum_aff,
                                                     self.genes_hemizygous,
                                                     None,
                                                     None, self.candidate_variants,
                                                      self.inhreport)
        inheritancefilter_102.inheritance_filter_genes()
        test_candidate_variants_102 = {'single_variants': {},
                                       'compound_hets': {}}

        self.assertEqual(inheritancefilter_102.candidate_variants,
                         test_candidate_variants_102)
        # mum aff, mum 0/1, dad 1/1 fail
        variants_per_gene_112 = create_test_variants_per_gene(self.variants_112,
                                                              self.family_XX_mum_aff)
        inheritancefilter_112 = InheritanceFiltering(variants_per_gene_112,
                                                     self.family_XX_mum_aff,
                                                     self.genes_hemizygous,
                                                     None,
                                                     None, self.candidate_variants,
                                                      self.inhreport)
        inheritancefilter_112.inheritance_filter_genes()
        test_candidate_variants_112 = {'single_variants': {},
                                       'compound_hets': {}}

        self.assertEqual(inheritancefilter_112.candidate_variants,
                         test_candidate_variants_112)
        # mum aff, mum 1/1, dad 1/1 fail
        variants_per_gene_122 = create_test_variants_per_gene(self.variants_122,
                                                              self.family_XX_mum_aff)
        inheritancefilter_122 = InheritanceFiltering(variants_per_gene_122,
                                                     self.family_XX_mum_aff,
                                                     self.genes_hemizygous,
                                                     None,
                                                     None, self.candidate_variants,
                                                      self.inhreport)
        inheritancefilter_122.inheritance_filter_genes()
        test_candidate_variants_122 = {'single_variants': {},
                                       'compound_hets': {}}

        self.assertEqual(inheritancefilter_122.candidate_variants,
                         test_candidate_variants_122)

        # dad aff, mum 0/0, dad 0/0 pass
        variants_per_gene_100 = create_test_variants_per_gene(self.variants_100,
                                                              self.family_XX_dad_aff)
        inheritancefilter_100 = InheritanceFiltering(variants_per_gene_100,
                                                     self.family_XX_dad_aff,
                                                     self.genes_hemizygous,
                                                     None,
                                                     None, self.candidate_variants,
                                                      self.inhreport)
        inheritancefilter_100.inheritance_filter_genes()
        test_candidate_variants_100 = {'single_variants': {
            'X_1097183_A_GG': {
                'mode': {'hemizygous'},
                'variant':
                    variants_per_gene_100[
                        '1234'][
                        'X_1097183_A_GG'][
                        'child'],
                'hgncid': '1234'}}, 'compound_hets': {}}

        self.assertEqual(inheritancefilter_100.candidate_variants,
                         test_candidate_variants_100)
        # dad aff, mum 0/1, dad 0/0 pass
        variants_per_gene_110 = create_test_variants_per_gene(self.variants_110,
                                                              self.family_XX_dad_aff)
        inheritancefilter_110 = InheritanceFiltering(variants_per_gene_110,
                                                     self.family_XX_dad_aff,
                                                     self.genes_hemizygous,
                                                     None,
                                                     None, self.candidate_variants,
                                                      self.inhreport)
        inheritancefilter_110.inheritance_filter_genes()
        test_candidate_variants_110 = {'single_variants': {
            'X_1097183_A_GG': {
                'mode': {'hemizygous'},
                'variant':
                    variants_per_gene_110[
                        '1234'][
                        'X_1097183_A_GG'][
                        'child'],
                'hgncid': '1234'}}, 'compound_hets': {}}

        self.assertEqual(inheritancefilter_110.candidate_variants,
                         test_candidate_variants_110)
        # dad aff, mum 1/1, dad 0/0 fail
        variants_per_gene_120 = create_test_variants_per_gene(self.variants_120,
                                                              self.family_XX_dad_aff)
        inheritancefilter_120 = InheritanceFiltering(variants_per_gene_120,
                                                     self.family_XX_dad_aff,
                                                     self.genes_hemizygous,
                                                     None,
                                                     None, self.candidate_variants,
                                                      self.inhreport)
        inheritancefilter_120.inheritance_filter_genes()
        test_candidate_variants_120 = {'single_variants': {},
                                       'compound_hets': {}}

        self.assertEqual(inheritancefilter_120.candidate_variants,
                         test_candidate_variants_120)
        # dad aff, mum 0/0, dad 1/1 pass
        variants_per_gene_102 = create_test_variants_per_gene(self.variants_102,
                                                              self.family_XX_dad_aff)
        inheritancefilter_102 = InheritanceFiltering(variants_per_gene_102,
                                                     self.family_XX_dad_aff,
                                                     self.genes_hemizygous,
                                                     None,
                                                     None, self.candidate_variants,
                                                      self.inhreport)
        inheritancefilter_102.inheritance_filter_genes()
        test_candidate_variants_102 = {'single_variants': {
            'X_1097183_A_GG': {
                'mode': {'hemizygous'},
                'variant':
                    variants_per_gene_102[
                        '1234'][
                        'X_1097183_A_GG'][
                        'child'],
                'hgncid': '1234'}}, 'compound_hets': {}}

        self.assertEqual(inheritancefilter_102.candidate_variants,
                         test_candidate_variants_102)
        # dad aff, mum 0/1, dad 1/1 pass
        variants_per_gene_112 = create_test_variants_per_gene(self.variants_112,
                                                              self.family_XX_dad_aff)
        inheritancefilter_112 = InheritanceFiltering(variants_per_gene_112,
                                                     self.family_XX_dad_aff,
                                                     self.genes_hemizygous,
                                                     None,
                                                     None, self.candidate_variants,
                                                      self.inhreport)
        inheritancefilter_112.inheritance_filter_genes()
        test_candidate_variants_112 = {'single_variants': {
            'X_1097183_A_GG': {
                'mode': {'hemizygous'},
                'variant':
                    variants_per_gene_112[
                        '1234'][
                        'X_1097183_A_GG'][
                        'child'],
                'hgncid': '1234'}}, 'compound_hets': {}}

        self.assertEqual(inheritancefilter_112.candidate_variants,
                         test_candidate_variants_112)
        # dad aff, mum 1/1, dad 1/1 fail
        variants_per_gene_122 = create_test_variants_per_gene(self.variants_122,
                                                              self.family_XX_dad_aff)
        inheritancefilter_122 = InheritanceFiltering(variants_per_gene_122,
                                                     self.family_XX_dad_aff,
                                                     self.genes_hemizygous,
                                                     None,
                                                     None, self.candidate_variants,
                                                      self.inhreport)
        inheritancefilter_122.inheritance_filter_genes()
        test_candidate_variants_122 = {'single_variants': {},
                                       'compound_hets': {}}

        self.assertEqual(inheritancefilter_122.candidate_variants,
                         test_candidate_variants_122)

        # both unaff, mum 0/0, dad 0/0 pass
        variants_per_gene_100 = create_test_variants_per_gene(self.variants_100,
                                                              self.family_XX_both_unaff)
        inheritancefilter_100 = InheritanceFiltering(variants_per_gene_100,
                                                     self.family_XX_both_unaff,
                                                     self.genes_hemizygous,
                                                     None,
                                                     None, self.candidate_variants,
                                                      self.inhreport)
        inheritancefilter_100.inheritance_filter_genes()
        test_candidate_variants_100 = {'single_variants': {
            'X_1097183_A_GG': {
                'mode': {'hemizygous'},
                'variant':
                    variants_per_gene_100[
                        '1234'][
                        'X_1097183_A_GG'][
                        'child'],
                'hgncid': '1234'}}, 'compound_hets': {}}

        self.assertEqual(inheritancefilter_100.candidate_variants,
                         test_candidate_variants_100)
        # both unaff, mum 0/1, dad 0/0 pass
        variants_per_gene_110 = create_test_variants_per_gene(self.variants_110,
                                                              self.family_XX_both_unaff)
        inheritancefilter_110 = InheritanceFiltering(variants_per_gene_110,
                                                     self.family_XX_both_unaff,
                                                     self.genes_hemizygous,
                                                     None,
                                                     None, self.candidate_variants,
                                                      self.inhreport)
        inheritancefilter_110.inheritance_filter_genes()
        test_candidate_variants_110 = {'single_variants': {
            'X_1097183_A_GG': {
                'mode': {'hemizygous'},
                'variant':
                    variants_per_gene_110[
                        '1234'][
                        'X_1097183_A_GG'][
                        'child'],
                'hgncid': '1234'}}, 'compound_hets': {}}

        self.assertEqual(inheritancefilter_110.candidate_variants,
                         test_candidate_variants_110)
        # both unaff, mum 1/1, dad 0/0 fail
        variants_per_gene_120 = create_test_variants_per_gene(self.variants_120,
                                                              self.family_XX_both_unaff)
        inheritancefilter_120 = InheritanceFiltering(variants_per_gene_120,
                                                     self.family_XX_both_unaff,
                                                     self.genes_hemizygous,
                                                     None,
                                                     None, self.candidate_variants,
                                                      self.inhreport)
        inheritancefilter_120.inheritance_filter_genes()
        test_candidate_variants_120 = {'single_variants': {},
                                       'compound_hets': {}}

        self.assertEqual(inheritancefilter_120.candidate_variants,
                         test_candidate_variants_120)
        # both unaff, mum 0/0, dad 1/1 fail
        variants_per_gene_102 = create_test_variants_per_gene(self.variants_102,
                                                              self.family_XX_both_unaff)
        inheritancefilter_102 = InheritanceFiltering(variants_per_gene_102,
                                                     self.family_XX_both_unaff,
                                                     self.genes_hemizygous,
                                                     None,
                                                     None, self.candidate_variants,
                                                      self.inhreport)
        inheritancefilter_102.inheritance_filter_genes()
        test_candidate_variants_102 = {'single_variants': {},
                                       'compound_hets': {}}

        self.assertEqual(inheritancefilter_102.candidate_variants,
                         test_candidate_variants_102)
        # both unaff, mum 0/1, dad 1/1 fail
        variants_per_gene_112 = create_test_variants_per_gene(self.variants_112,
                                                              self.family_XX_both_unaff)
        inheritancefilter_112 = InheritanceFiltering(variants_per_gene_112,
                                                     self.family_XX_both_unaff,
                                                     self.genes_hemizygous,
                                                     None,
                                                     None, self.candidate_variants,
                                                      self.inhreport)
        inheritancefilter_112.inheritance_filter_genes()
        test_candidate_variants_112 = {'single_variants': {},
                                       'compound_hets': {}}

        self.assertEqual(inheritancefilter_112.candidate_variants,
                         test_candidate_variants_112)
        # both unaff, mum 1/1, dad 1/1 fail
        variants_per_gene_122 = create_test_variants_per_gene(self.variants_122,
                                                              self.family_XX_both_unaff)
        inheritancefilter_122 = InheritanceFiltering(variants_per_gene_122,
                                                     self.family_XX_both_unaff,
                                                     self.genes_hemizygous,
                                                     None,
                                                     None, self.candidate_variants,
                                                      self.inhreport)
        inheritancefilter_122.inheritance_filter_genes()
        test_candidate_variants_122 = {'single_variants': {},
                                       'compound_hets': {}}

        self.assertEqual(inheritancefilter_122.candidate_variants,
                         test_candidate_variants_122)

    def test_gn_X_linked_dominant_gt_heterozygous_parents_filter(self):
        # both aff, mum 0/0, dad 0/0 pass
        variants_per_gene_100 = create_test_variants_per_gene(self.variants_100,
                                                              self.family_XX_both_aff)
        inheritancefilter_100 = InheritanceFiltering(variants_per_gene_100,
                                                     self.family_XX_both_aff,
                                                     self.genes_X_linked_dominant,
                                                     None,
                                                     None, self.candidate_variants,
                                                      self.inhreport)
        inheritancefilter_100.inheritance_filter_genes()
        test_candidate_variants_100 = {'single_variants': {
            'X_1097183_A_GG': {
                'mode': {'X-linked dominant'},
                'variant':
                    variants_per_gene_100[
                        '1234'][
                        'X_1097183_A_GG'][
                        'child'],
                'hgncid': '1234'}}, 'compound_hets': {}}

        self.assertEqual(inheritancefilter_100.candidate_variants,
                         test_candidate_variants_100)
        # both aff, mum 0/1, dad 0/0 pass
        variants_per_gene_110 = create_test_variants_per_gene(self.variants_110,
                                                              self.family_XX_both_aff)
        inheritancefilter_110 = InheritanceFiltering(variants_per_gene_110,
                                                     self.family_XX_both_aff,
                                                     self.genes_X_linked_dominant,
                                                     None,
                                                     None, self.candidate_variants,
                                                      self.inhreport)
        inheritancefilter_110.inheritance_filter_genes()
        test_candidate_variants_110 = {'single_variants': {
            'X_1097183_A_GG': {
                'mode': {'X-linked dominant'},
                'variant':
                    variants_per_gene_110[
                        '1234'][
                        'X_1097183_A_GG'][
                        'child'],
                'hgncid': '1234'}}, 'compound_hets': {}}

        self.assertEqual(inheritancefilter_110.candidate_variants,
                         test_candidate_variants_110)
        # both aff, mum 1/1, dad 0/0 pass
        variants_per_gene_120 = create_test_variants_per_gene(self.variants_120,
                                                              self.family_XX_both_aff)
        inheritancefilter_120 = InheritanceFiltering(variants_per_gene_120,
                                                     self.family_XX_both_aff,
                                                     self.genes_X_linked_dominant,
                                                     None,
                                                     None, self.candidate_variants,
                                                      self.inhreport)
        inheritancefilter_120.inheritance_filter_genes()
        test_candidate_variants_120 = {'single_variants': {
            'X_1097183_A_GG': {
                'mode': {'X-linked dominant'},
                'variant':
                    variants_per_gene_120[
                        '1234'][
                        'X_1097183_A_GG'][
                        'child'],
                'hgncid': '1234'}}, 'compound_hets': {}}

        self.assertEqual(inheritancefilter_120.candidate_variants,
                         test_candidate_variants_120)
        # both aff, mum 0/0, dad 1/1 pass
        variants_per_gene_102 = create_test_variants_per_gene(self.variants_102,
                                                              self.family_XX_both_aff)
        inheritancefilter_102 = InheritanceFiltering(variants_per_gene_102,
                                                     self.family_XX_both_aff,
                                                     self.genes_X_linked_dominant,
                                                     None,
                                                     None, self.candidate_variants,
                                                      self.inhreport)
        inheritancefilter_102.inheritance_filter_genes()
        test_candidate_variants_102 = {'single_variants': {
            'X_1097183_A_GG': {
                'mode': {'X-linked dominant'},
                'variant':
                    variants_per_gene_102[
                        '1234'][
                        'X_1097183_A_GG'][
                        'child'],
                'hgncid': '1234'}}, 'compound_hets': {}}

        self.assertEqual(inheritancefilter_102.candidate_variants,
                         test_candidate_variants_102)
        # both aff, mum 0/1, dad 1/1 pass
        variants_per_gene_112 = create_test_variants_per_gene(self.variants_112,
                                                              self.family_XX_both_aff)
        inheritancefilter_112 = InheritanceFiltering(variants_per_gene_112,
                                                     self.family_XX_both_aff,
                                                     self.genes_X_linked_dominant,
                                                     None,
                                                     None, self.candidate_variants,
                                                      self.inhreport)
        inheritancefilter_112.inheritance_filter_genes()
        test_candidate_variants_112 = {'single_variants': {
            'X_1097183_A_GG': {
                'mode': {'X-linked dominant'},
                'variant':
                    variants_per_gene_112[
                        '1234'][
                        'X_1097183_A_GG'][
                        'child'],
                'hgncid': '1234'}}, 'compound_hets': {}}

        self.assertEqual(inheritancefilter_112.candidate_variants,
                         test_candidate_variants_112)
        # both aff, mum 1/1, dad 1/1 fail
        variants_per_gene_122 = create_test_variants_per_gene(self.variants_122,
                                                              self.family_XX_both_aff)
        inheritancefilter_122 = InheritanceFiltering(variants_per_gene_122,
                                                     self.family_XX_both_aff,
                                                     self.genes_X_linked_dominant,
                                                     None,
                                                     None, self.candidate_variants,
                                                      self.inhreport)
        inheritancefilter_122.inheritance_filter_genes()
        test_candidate_variants_122 = {'single_variants': {},
                                       'compound_hets': {}}

        self.assertEqual(inheritancefilter_122.candidate_variants,
                         test_candidate_variants_122)

        # mum aff, mum 0/0, dad 0/0 pass
        variants_per_gene_100 = create_test_variants_per_gene(self.variants_100,
                                                              self.family_XX_mum_aff)
        inheritancefilter_100 = InheritanceFiltering(variants_per_gene_100,
                                                     self.family_XX_mum_aff,
                                                     self.genes_X_linked_dominant,
                                                     None,
                                                     None, self.candidate_variants,
                                                      self.inhreport)
        inheritancefilter_100.inheritance_filter_genes()
        test_candidate_variants_100 = {'single_variants': {
            'X_1097183_A_GG': {
                'mode': {'X-linked dominant'},
                'variant':
                    variants_per_gene_100[
                        '1234'][
                        'X_1097183_A_GG'][
                        'child'],
                'hgncid': '1234'}}, 'compound_hets': {}}

        self.assertEqual(inheritancefilter_100.candidate_variants,
                         test_candidate_variants_100)
        # mum aff, mum 0/1, dad 0/0 pass
        variants_per_gene_110 = create_test_variants_per_gene(self.variants_110,
                                                              self.family_XX_mum_aff)
        inheritancefilter_110 = InheritanceFiltering(variants_per_gene_110,
                                                     self.family_XX_mum_aff,
                                                     self.genes_X_linked_dominant,
                                                     None,
                                                     None, self.candidate_variants,
                                                      self.inhreport)
        inheritancefilter_110.inheritance_filter_genes()
        test_candidate_variants_110 = {'single_variants': {
            'X_1097183_A_GG': {
                'mode': {'X-linked dominant'},
                'variant':
                    variants_per_gene_110[
                        '1234'][
                        'X_1097183_A_GG'][
                        'child'],
                'hgncid': '1234'}}, 'compound_hets': {}}

        self.assertEqual(inheritancefilter_110.candidate_variants,
                         test_candidate_variants_110)
        # mum aff, mum 1/1, dad 0/0 pass
        variants_per_gene_120 = create_test_variants_per_gene(self.variants_120,
                                                              self.family_XX_mum_aff)
        inheritancefilter_120 = InheritanceFiltering(variants_per_gene_120,
                                                     self.family_XX_mum_aff,
                                                     self.genes_X_linked_dominant,
                                                     None,
                                                     None, self.candidate_variants,
                                                      self.inhreport)
        inheritancefilter_120.inheritance_filter_genes()
        test_candidate_variants_120 = {'single_variants': {
            'X_1097183_A_GG': {
                'mode': {'X-linked dominant'},
                'variant':
                    variants_per_gene_120[
                        '1234'][
                        'X_1097183_A_GG'][
                        'child'],
                'hgncid': '1234'}}, 'compound_hets': {}}

        self.assertEqual(inheritancefilter_120.candidate_variants,
                         test_candidate_variants_120)
        # mum aff, mum 0/0, dad 1/1 fail
        variants_per_gene_102 = create_test_variants_per_gene(self.variants_102,
                                                              self.family_XX_mum_aff)
        inheritancefilter_102 = InheritanceFiltering(variants_per_gene_102,
                                                     self.family_XX_mum_aff,
                                                     self.genes_X_linked_dominant,
                                                     None,
                                                     None, self.candidate_variants,
                                                      self.inhreport)
        inheritancefilter_102.inheritance_filter_genes()
        test_candidate_variants_102 = {'single_variants': {}, 'compound_hets': {}}

        self.assertEqual(inheritancefilter_102.candidate_variants,
                         test_candidate_variants_102)
        # mum aff, mum 0/1, dad 1/1 fail
        variants_per_gene_112 = create_test_variants_per_gene(self.variants_112,
                                                              self.family_XX_mum_aff)
        inheritancefilter_112 = InheritanceFiltering(variants_per_gene_112,
                                                     self.family_XX_mum_aff,
                                                     self.genes_X_linked_dominant,
                                                     None,
                                                     None, self.candidate_variants,
                                                      self.inhreport)
        inheritancefilter_112.inheritance_filter_genes()
        test_candidate_variants_112 = {'single_variants': {},
                                       'compound_hets': {}}

        self.assertEqual(inheritancefilter_112.candidate_variants,
                         test_candidate_variants_112)
        # mum aff, mum 1/1, dad 1/1 fail
        variants_per_gene_122 = create_test_variants_per_gene(self.variants_122,
                                                              self.family_XX_mum_aff)
        inheritancefilter_122 = InheritanceFiltering(variants_per_gene_122,
                                                     self.family_XX_mum_aff,
                                                     self.genes_X_linked_dominant,
                                                     None,
                                                     None, self.candidate_variants,
                                                      self.inhreport)
        inheritancefilter_122.inheritance_filter_genes()
        test_candidate_variants_122 = {'single_variants': {},
                                       'compound_hets': {}}

        self.assertEqual(inheritancefilter_122.candidate_variants,
                         test_candidate_variants_122)

        # dad aff, mum 0/0, dad 0/0 pass
        variants_per_gene_100 = create_test_variants_per_gene(self.variants_100,
                                                              self.family_XX_dad_aff)
        inheritancefilter_100 = InheritanceFiltering(variants_per_gene_100,
                                                     self.family_XX_dad_aff,
                                                     self.genes_X_linked_dominant,
                                                     None,
                                                     None, self.candidate_variants,
                                                      self.inhreport)
        inheritancefilter_100.inheritance_filter_genes()
        test_candidate_variants_100 = {'single_variants': {
            'X_1097183_A_GG': {
                'mode': {'X-linked dominant'},
                'variant':
                    variants_per_gene_100[
                        '1234'][
                        'X_1097183_A_GG'][
                        'child'],
                'hgncid': '1234'}}, 'compound_hets': {}}

        self.assertEqual(inheritancefilter_100.candidate_variants,
                         test_candidate_variants_100)
        # dad aff, mum 0/1, dad 0/0 pass
        variants_per_gene_110 = create_test_variants_per_gene(self.variants_110,
                                                              self.family_XX_dad_aff)
        inheritancefilter_110 = InheritanceFiltering(variants_per_gene_110,
                                                     self.family_XX_dad_aff,
                                                     self.genes_X_linked_dominant,
                                                     None,
                                                     None, self.candidate_variants,
                                                      self.inhreport)
        inheritancefilter_110.inheritance_filter_genes()
        test_candidate_variants_110 = {'single_variants': {
            'X_1097183_A_GG': {
                'mode': {'X-linked dominant'},
                'variant':
                    variants_per_gene_110[
                        '1234'][
                        'X_1097183_A_GG'][
                        'child'],
                'hgncid': '1234'}}, 'compound_hets': {}}

        self.assertEqual(inheritancefilter_110.candidate_variants,
                         test_candidate_variants_110)
        # dad aff, mum 1/1, dad 0/0 fail
        variants_per_gene_120 = create_test_variants_per_gene(self.variants_120,
                                                              self.family_XX_dad_aff)
        inheritancefilter_120 = InheritanceFiltering(variants_per_gene_120,
                                                     self.family_XX_dad_aff,
                                                     self.genes_X_linked_dominant,
                                                     None,
                                                     None, self.candidate_variants,
                                                      self.inhreport)
        inheritancefilter_120.inheritance_filter_genes()
        test_candidate_variants_120 = {'single_variants': {},
                                       'compound_hets': {}}

        self.assertEqual(inheritancefilter_120.candidate_variants,
                         test_candidate_variants_120)
        # dad aff, mum 0/0, dad 1/1 pass
        variants_per_gene_102 = create_test_variants_per_gene(self.variants_102,
                                                              self.family_XX_dad_aff)
        inheritancefilter_102 = InheritanceFiltering(variants_per_gene_102,
                                                     self.family_XX_dad_aff,
                                                     self.genes_X_linked_dominant,
                                                     None,
                                                     None, self.candidate_variants,
                                                      self.inhreport)
        inheritancefilter_102.inheritance_filter_genes()
        test_candidate_variants_102 = {'single_variants': {
            'X_1097183_A_GG': {
                'mode': {'X-linked dominant'},
                'variant':
                    variants_per_gene_102[
                        '1234'][
                        'X_1097183_A_GG'][
                        'child'],
                'hgncid': '1234'}}, 'compound_hets': {}}

        self.assertEqual(inheritancefilter_102.candidate_variants,
                         test_candidate_variants_102)
        # dad aff, mum 0/1, dad 1/1 pass
        variants_per_gene_112 = create_test_variants_per_gene(self.variants_112,
                                                              self.family_XX_dad_aff)
        inheritancefilter_112 = InheritanceFiltering(variants_per_gene_112,
                                                     self.family_XX_dad_aff,
                                                     self.genes_X_linked_dominant,
                                                     None,
                                                     None, self.candidate_variants,
                                                      self.inhreport)
        inheritancefilter_112.inheritance_filter_genes()
        test_candidate_variants_112 = {'single_variants': {
            'X_1097183_A_GG': {
                'mode': {'X-linked dominant'},
                'variant':
                    variants_per_gene_112[
                        '1234'][
                        'X_1097183_A_GG'][
                        'child'],
                'hgncid': '1234'}}, 'compound_hets': {}}

        self.assertEqual(inheritancefilter_112.candidate_variants,
                         test_candidate_variants_112)
        # dad aff, mum 1/1, dad 1/1 fail
        variants_per_gene_122 = create_test_variants_per_gene(self.variants_122,
                                                              self.family_XX_dad_aff)
        inheritancefilter_122 = InheritanceFiltering(variants_per_gene_122,
                                                     self.family_XX_dad_aff,
                                                     self.genes_X_linked_dominant,
                                                     None,
                                                     None, self.candidate_variants,
                                                      self.inhreport)
        inheritancefilter_122.inheritance_filter_genes()
        test_candidate_variants_122 = {'single_variants': {},
                                       'compound_hets': {}}

        self.assertEqual(inheritancefilter_122.candidate_variants,
                         test_candidate_variants_122)

        # both unaff, mum 0/0, dad 0/0 pass
        variants_per_gene_100 = create_test_variants_per_gene(self.variants_100,
                                                              self.family_XX_both_unaff)
        inheritancefilter_100 = InheritanceFiltering(variants_per_gene_100,
                                                     self.family_XX_both_unaff,
                                                     self.genes_X_linked_dominant,
                                                     None,
                                                     None, self.candidate_variants,
                                                      self.inhreport)
        inheritancefilter_100.inheritance_filter_genes()
        test_candidate_variants_100 = {'single_variants': {
            'X_1097183_A_GG': {
                'mode': {'X-linked dominant'},
                'variant':
                    variants_per_gene_100[
                        '1234'][
                        'X_1097183_A_GG'][
                        'child'],
                'hgncid': '1234'}}, 'compound_hets': {}}

        self.assertEqual(inheritancefilter_100.candidate_variants,
                         test_candidate_variants_100)
        # both unaff, mum 0/1, dad 0/0 pass
        variants_per_gene_110 = create_test_variants_per_gene(self.variants_110,
                                                              self.family_XX_both_unaff)
        inheritancefilter_110 = InheritanceFiltering(variants_per_gene_110,
                                                     self.family_XX_both_unaff,
                                                     self.genes_X_linked_dominant,
                                                     None,
                                                     None, self.candidate_variants,
                                                      self.inhreport)
        inheritancefilter_110.inheritance_filter_genes()
        test_candidate_variants_110 = {'single_variants': {
            'X_1097183_A_GG': {
                'mode': {'X-linked dominant'},
                'variant':
                    variants_per_gene_110[
                        '1234'][
                        'X_1097183_A_GG'][
                        'child'],
                'hgncid': '1234'}}, 'compound_hets': {}}

        self.assertEqual(inheritancefilter_110.candidate_variants,
                         test_candidate_variants_110)
        # both unaff, mum 1/1, dad 0/0 fail
        variants_per_gene_120 = create_test_variants_per_gene(self.variants_120,
                                                              self.family_XX_both_unaff)
        inheritancefilter_120 = InheritanceFiltering(variants_per_gene_120,
                                                     self.family_XX_both_unaff,
                                                     self.genes_X_linked_dominant,
                                                     None,
                                                     None, self.candidate_variants,
                                                      self.inhreport)
        inheritancefilter_120.inheritance_filter_genes()
        test_candidate_variants_120 = {'single_variants': {},
                                       'compound_hets': {}}

        self.assertEqual(inheritancefilter_120.candidate_variants,
                         test_candidate_variants_120)
        # both unaff, mum 0/0, dad 1/1 fail
        variants_per_gene_102 = create_test_variants_per_gene(self.variants_102,
                                                              self.family_XX_both_unaff)
        inheritancefilter_102 = InheritanceFiltering(variants_per_gene_102,
                                                     self.family_XX_both_unaff,
                                                     self.genes_X_linked_dominant,
                                                     None,
                                                     None, self.candidate_variants,
                                                      self.inhreport)
        inheritancefilter_102.inheritance_filter_genes()
        test_candidate_variants_102 = {'single_variants': {},
                                       'compound_hets': {}}

        self.assertEqual(inheritancefilter_102.candidate_variants,
                         test_candidate_variants_102)
        # both unaff, mum 0/1, dad 1/1 fail
        variants_per_gene_112 = create_test_variants_per_gene(self.variants_112,
                                                              self.family_XX_both_unaff)
        inheritancefilter_112 = InheritanceFiltering(variants_per_gene_112,
                                                     self.family_XX_both_unaff,
                                                     self.genes_X_linked_dominant,
                                                     None,
                                                     None, self.candidate_variants,
                                                      self.inhreport)
        inheritancefilter_112.inheritance_filter_genes()
        test_candidate_variants_112 = {'single_variants': {},
                                       'compound_hets': {}}

        self.assertEqual(inheritancefilter_112.candidate_variants,
                         test_candidate_variants_112)
        # both unaff, mum 1/1, dad 1/1 fail
        variants_per_gene_122 = create_test_variants_per_gene(self.variants_122,
                                                              self.family_XX_both_unaff)
        inheritancefilter_122 = InheritanceFiltering(variants_per_gene_122,
                                                     self.family_XX_both_unaff,
                                                     self.genes_X_linked_dominant,
                                                     None,
                                                     None, self.candidate_variants,
                                                      self.inhreport)
        inheritancefilter_122.inheritance_filter_genes()
        test_candidate_variants_122 = {'single_variants': {},
                                       'compound_hets': {}}

        self.assertEqual(inheritancefilter_122.candidate_variants,
                         test_candidate_variants_122)

    def test_gn_X_linked_over_dominant_gt_heterozygous_parents_filter(self):
        # both aff, mum 0/0, dad 0/0 pass
        variants_per_gene_100 = create_test_variants_per_gene(self.variants_100,
                                                              self.family_XX_both_aff)
        inheritancefilter_100 = InheritanceFiltering(variants_per_gene_100,
                                                     self.family_XX_both_aff,
                                                     self.genes_X_linked_over_dominance,
                                                     None,
                                                     None, self.candidate_variants,
                                                      self.inhreport)
        inheritancefilter_100.inheritance_filter_genes()
        test_candidate_variants_100 = {'single_variants': {
            'X_1097183_A_GG': {
                'mode': {'X-linked over-dominance'},
                'variant':
                    variants_per_gene_100[
                        '1234'][
                        'X_1097183_A_GG'][
                        'child'],
                'hgncid': '1234'}}, 'compound_hets': {}}

        self.assertEqual(inheritancefilter_100.candidate_variants,
                         test_candidate_variants_100)
        # both aff, mum 0/1, dad 0/0 fail
        variants_per_gene_110 = create_test_variants_per_gene(self.variants_110,
                                                              self.family_XX_both_aff)
        inheritancefilter_110 = InheritanceFiltering(variants_per_gene_110,
                                                     self.family_XX_both_aff,
                                                     self.genes_X_linked_over_dominance,
                                                     None,
                                                     None, self.candidate_variants,
                                                      self.inhreport)
        inheritancefilter_110.inheritance_filter_genes()
        test_candidate_variants_110 = {'single_variants': {},
                                       'compound_hets': {}}

        self.assertEqual(inheritancefilter_110.candidate_variants,
                         test_candidate_variants_110)
        # both aff, mum 1/1, dad 0/0 fail
        variants_per_gene_120 = create_test_variants_per_gene(self.variants_120,
                                                              self.family_XX_both_aff)
        inheritancefilter_120 = InheritanceFiltering(variants_per_gene_120,
                                                     self.family_XX_both_aff,
                                                     self.genes_X_linked_over_dominance,
                                                     None,
                                                     None, self.candidate_variants,
                                                      self.inhreport)
        inheritancefilter_120.inheritance_filter_genes()
        test_candidate_variants_120 = {'single_variants': {},
                                       'compound_hets': {}}

        self.assertEqual(inheritancefilter_120.candidate_variants,
                         test_candidate_variants_120)
        # both aff, mum 0/0, dad 1/1 fail
        variants_per_gene_102 = create_test_variants_per_gene(self.variants_102,
                                                              self.family_XX_both_aff)
        inheritancefilter_102 = InheritanceFiltering(variants_per_gene_102,
                                                     self.family_XX_both_aff,
                                                     self.genes_X_linked_over_dominance,
                                                     None,
                                                     None, self.candidate_variants,
                                                      self.inhreport)
        inheritancefilter_102.inheritance_filter_genes()
        test_candidate_variants_102 = {'single_variants': {},
                                       'compound_hets': {}}

        self.assertEqual(inheritancefilter_102.candidate_variants,
                         test_candidate_variants_102)
        # both aff, mum 0/1, dad 1/1 fail
        variants_per_gene_112 = create_test_variants_per_gene(self.variants_112,
                                                              self.family_XX_both_aff)
        inheritancefilter_112 = InheritanceFiltering(variants_per_gene_112,
                                                     self.family_XX_both_aff,
                                                     self.genes_X_linked_over_dominance,
                                                     None,
                                                     None, self.candidate_variants,
                                                      self.inhreport)
        inheritancefilter_112.inheritance_filter_genes()
        test_candidate_variants_112 = {'single_variants': {},
                                       'compound_hets': {}}

        self.assertEqual(inheritancefilter_112.candidate_variants,
                         test_candidate_variants_112)
        # both aff, mum 1/1, dad 1/1 fail
        variants_per_gene_122 = create_test_variants_per_gene(self.variants_122,
                                                              self.family_XX_both_aff)
        inheritancefilter_122 = InheritanceFiltering(variants_per_gene_122,
                                                     self.family_XX_both_aff,
                                                     self.genes_X_linked_over_dominance,
                                                     None,
                                                     None, self.candidate_variants,
                                                      self.inhreport)
        inheritancefilter_122.inheritance_filter_genes()
        test_candidate_variants_122 = {'single_variants': {},
                                       'compound_hets': {}}

        self.assertEqual(inheritancefilter_122.candidate_variants,
                         test_candidate_variants_122)

        # mum aff, mum 0/0, dad 0/0 pass
        variants_per_gene_100 = create_test_variants_per_gene(self.variants_100,
                                                              self.family_XX_mum_aff)
        inheritancefilter_100 = InheritanceFiltering(variants_per_gene_100,
                                                     self.family_XX_mum_aff,
                                                     self.genes_X_linked_over_dominance,
                                                     None,
                                                     None, self.candidate_variants,
                                                      self.inhreport)
        inheritancefilter_100.inheritance_filter_genes()
        test_candidate_variants_100 = {'single_variants': {
            'X_1097183_A_GG': {
                'mode': {'X-linked over-dominance'},
                'variant':
                    variants_per_gene_100[
                        '1234'][
                        'X_1097183_A_GG'][
                        'child'],
                'hgncid': '1234'}}, 'compound_hets': {}}

        self.assertEqual(inheritancefilter_100.candidate_variants,
                         test_candidate_variants_100)
        # mum aff, mum 0/1, dad 0/0 pass
        variants_per_gene_110 = create_test_variants_per_gene(self.variants_110,
                                                              self.family_XX_mum_aff)
        inheritancefilter_110 = InheritanceFiltering(variants_per_gene_110,
                                                     self.family_XX_mum_aff,
                                                     self.genes_X_linked_over_dominance,
                                                     None,
                                                     None, self.candidate_variants,
                                                      self.inhreport)
        inheritancefilter_110.inheritance_filter_genes()
        test_candidate_variants_110 = {'single_variants': {
            'X_1097183_A_GG': {
                'mode': {'X-linked over-dominance'},
                'variant':
                    variants_per_gene_110[
                        '1234'][
                        'X_1097183_A_GG'][
                        'child'],
                'hgncid': '1234'}}, 'compound_hets': {}}

        self.assertEqual(inheritancefilter_110.candidate_variants,
                         test_candidate_variants_110)
        # mum aff, mum 1/1, dad 0/0 fail
        variants_per_gene_120 = create_test_variants_per_gene(self.variants_120,
                                                              self.family_XX_mum_aff)
        inheritancefilter_120 = InheritanceFiltering(variants_per_gene_120,
                                                     self.family_XX_mum_aff,
                                                     self.genes_X_linked_over_dominance,
                                                     None,
                                                     None, self.candidate_variants,
                                                      self.inhreport)
        inheritancefilter_120.inheritance_filter_genes()
        test_candidate_variants_120 = {'single_variants': {},
                                       'compound_hets': {}}

        self.assertEqual(inheritancefilter_120.candidate_variants,
                         test_candidate_variants_120)
        # mum aff, mum 0/0, dad 1/1 fail
        variants_per_gene_102 = create_test_variants_per_gene(self.variants_102,
                                                              self.family_XX_mum_aff)
        inheritancefilter_102 = InheritanceFiltering(variants_per_gene_102,
                                                     self.family_XX_mum_aff,
                                                     self.genes_X_linked_over_dominance,
                                                     None,
                                                     None, self.candidate_variants,
                                                      self.inhreport)
        inheritancefilter_102.inheritance_filter_genes()
        test_candidate_variants_102 = {'single_variants': {},
                                       'compound_hets': {}}

        self.assertEqual(inheritancefilter_102.candidate_variants,
                         test_candidate_variants_102)
        # mum aff, mum 0/1, dad 1/1 pass
        variants_per_gene_112 = create_test_variants_per_gene(self.variants_112,
                                                              self.family_XX_mum_aff)
        inheritancefilter_112 = InheritanceFiltering(variants_per_gene_112,
                                                     self.family_XX_mum_aff,
                                                     self.genes_X_linked_over_dominance,
                                                     None,
                                                     None, self.candidate_variants,
                                                      self.inhreport)
        inheritancefilter_112.inheritance_filter_genes()
        test_candidate_variants_112 = {'single_variants': {
            'X_1097183_A_GG': {
                'mode': {'X-linked over-dominance'},
                'variant':
                    variants_per_gene_112[
                        '1234'][
                        'X_1097183_A_GG'][
                        'child'],
                'hgncid': '1234'}}, 'compound_hets': {}}

        self.assertEqual(inheritancefilter_112.candidate_variants,
                         test_candidate_variants_112)
        # mum aff, mum 1/1, dad 1/1 fail
        variants_per_gene_122 = create_test_variants_per_gene(self.variants_122,
                                                              self.family_XX_mum_aff)
        inheritancefilter_122 = InheritanceFiltering(variants_per_gene_122,
                                                     self.family_XX_mum_aff,
                                                     self.genes_X_linked_over_dominance,
                                                     None,
                                                     None, self.candidate_variants,
                                                      self.inhreport)
        inheritancefilter_122.inheritance_filter_genes()
        test_candidate_variants_122 = {'single_variants': {},
                                       'compound_hets': {}}

        self.assertEqual(inheritancefilter_122.candidate_variants,
                         test_candidate_variants_122)

        # dad aff, mum 0/0, dad 0/0 pass
        variants_per_gene_100 = create_test_variants_per_gene(self.variants_100,
                                                              self.family_XX_dad_aff)
        inheritancefilter_100 = InheritanceFiltering(variants_per_gene_100,
                                                     self.family_XX_dad_aff,
                                                     self.genes_X_linked_over_dominance,
                                                     None,
                                                     None, self.candidate_variants,
                                                      self.inhreport)
        inheritancefilter_100.inheritance_filter_genes()
        test_candidate_variants_100 = {'single_variants': {
            'X_1097183_A_GG': {
                'mode': {'X-linked over-dominance'},
                'variant':
                    variants_per_gene_100[
                        '1234'][
                        'X_1097183_A_GG'][
                        'child'],
                'hgncid': '1234'}}, 'compound_hets': {}}

        self.assertEqual(inheritancefilter_100.candidate_variants,
                         test_candidate_variants_100)
        # dad aff, mum 0/1, dad 0/0 fail
        variants_per_gene_110 = create_test_variants_per_gene(self.variants_110,
                                                              self.family_XX_dad_aff)
        inheritancefilter_110 = InheritanceFiltering(variants_per_gene_110,
                                                     self.family_XX_dad_aff,
                                                     self.genes_X_linked_over_dominance,
                                                     None,
                                                     None, self.candidate_variants,
                                                      self.inhreport)
        inheritancefilter_110.inheritance_filter_genes()
        test_candidate_variants_110 = {'single_variants': {},
                                       'compound_hets': {}}

        self.assertEqual(inheritancefilter_110.candidate_variants,
                         test_candidate_variants_110)
        # dad aff, mum 1/1, dad 0/0 fail
        variants_per_gene_120 = create_test_variants_per_gene(self.variants_120,
                                                              self.family_XX_dad_aff)
        inheritancefilter_120 = InheritanceFiltering(variants_per_gene_120,
                                                     self.family_XX_dad_aff,
                                                     self.genes_X_linked_over_dominance,
                                                     None,
                                                     None, self.candidate_variants,
                                                      self.inhreport)
        inheritancefilter_120.inheritance_filter_genes()
        test_candidate_variants_120 = {'single_variants': {},
                                       'compound_hets': {}}

        self.assertEqual(inheritancefilter_120.candidate_variants,
                         test_candidate_variants_120)
        # dad aff, mum 0/0, dad 1/1 fail
        variants_per_gene_102 = create_test_variants_per_gene(self.variants_102,
                                                              self.family_XX_dad_aff)
        inheritancefilter_102 = InheritanceFiltering(variants_per_gene_102,
                                                     self.family_XX_dad_aff,
                                                     self.genes_X_linked_over_dominance,
                                                     None,
                                                     None, self.candidate_variants,
                                                      self.inhreport)
        inheritancefilter_102.inheritance_filter_genes()
        test_candidate_variants_102 = {'single_variants': {},
                                       'compound_hets': {}}

        self.assertEqual(inheritancefilter_102.candidate_variants,
                         test_candidate_variants_102)
        # dad aff, mum 0/1, dad 1/1 fail
        variants_per_gene_112 = create_test_variants_per_gene(self.variants_112,
                                                              self.family_XX_dad_aff)
        inheritancefilter_112 = InheritanceFiltering(variants_per_gene_112,
                                                     self.family_XX_dad_aff,
                                                     self.genes_X_linked_over_dominance,
                                                     None,
                                                     None, self.candidate_variants,
                                                      self.inhreport)
        inheritancefilter_112.inheritance_filter_genes()
        test_candidate_variants_112 = {'single_variants': {},
                                       'compound_hets': {}}

        self.assertEqual(inheritancefilter_112.candidate_variants,
                         test_candidate_variants_112)
        # dad aff, mum 1/1, dad 1/1 fail
        variants_per_gene_122 = create_test_variants_per_gene(self.variants_122,
                                                              self.family_XX_dad_aff)
        inheritancefilter_122 = InheritanceFiltering(variants_per_gene_122,
                                                     self.family_XX_dad_aff,
                                                     self.genes_X_linked_over_dominance,
                                                     None,
                                                     None, self.candidate_variants,
                                                      self.inhreport)
        inheritancefilter_122.inheritance_filter_genes()
        test_candidate_variants_122 = {'single_variants': {},
                                       'compound_hets': {}}

        self.assertEqual(inheritancefilter_122.candidate_variants,
                         test_candidate_variants_122)

        # both unaff, mum 0/0, dad 0/0 pass
        variants_per_gene_100 = create_test_variants_per_gene(self.variants_100,
                                                              self.family_XX_both_unaff)
        inheritancefilter_100 = InheritanceFiltering(variants_per_gene_100,
                                                     self.family_XX_both_unaff,
                                                     self.genes_X_linked_over_dominance,
                                                     None,
                                                     None, self.candidate_variants,
                                                      self.inhreport)
        inheritancefilter_100.inheritance_filter_genes()
        test_candidate_variants_100 = {'single_variants': {
            'X_1097183_A_GG': {
                'mode': {'X-linked over-dominance'},
                'variant':
                    variants_per_gene_100[
                        '1234'][
                        'X_1097183_A_GG'][
                        'child'],
                'hgncid': '1234'}}, 'compound_hets': {}}

        self.assertEqual(inheritancefilter_100.candidate_variants,
                         test_candidate_variants_100)
        # both unaff, mum 0/1, dad 0/0 fail
        variants_per_gene_110 = create_test_variants_per_gene(self.variants_110,
                                                              self.family_XX_both_unaff)
        inheritancefilter_110 = InheritanceFiltering(variants_per_gene_110,
                                                     self.family_XX_both_unaff,
                                                     self.genes_X_linked_over_dominance,
                                                     None,
                                                     None, self.candidate_variants,
                                                      self.inhreport)
        inheritancefilter_110.inheritance_filter_genes()
        test_candidate_variants_110 = {'single_variants': {},
                                       'compound_hets': {}}

        self.assertEqual(inheritancefilter_110.candidate_variants,
                         test_candidate_variants_110)
        # both unaff, mum 1/1, dad 0/0 pass
        variants_per_gene_120 = create_test_variants_per_gene(self.variants_120,
                                                              self.family_XX_both_unaff)
        inheritancefilter_120 = InheritanceFiltering(variants_per_gene_120,
                                                     self.family_XX_both_unaff,
                                                     self.genes_X_linked_over_dominance,
                                                     None,
                                                     None, self.candidate_variants,
                                                      self.inhreport)
        inheritancefilter_120.inheritance_filter_genes()
        test_candidate_variants_120 = {'single_variants': {
            'X_1097183_A_GG': {
                'mode': {'X-linked over-dominance'},
                'variant':
                    variants_per_gene_120[
                        '1234'][
                        'X_1097183_A_GG'][
                        'child'],
                'hgncid': '1234'}}, 'compound_hets': {}}

        self.assertEqual(inheritancefilter_120.candidate_variants,
                         test_candidate_variants_120)
        # both unaff, mum 0/0, dad 1/1 pass
        variants_per_gene_102 = create_test_variants_per_gene(self.variants_102,
                                                              self.family_XX_both_unaff)
        inheritancefilter_102 = InheritanceFiltering(variants_per_gene_102,
                                                     self.family_XX_both_unaff,
                                                     self.genes_X_linked_over_dominance,
                                                     None,
                                                     None, self.candidate_variants,
                                                      self.inhreport)
        inheritancefilter_102.inheritance_filter_genes()
        test_candidate_variants_102 = {'single_variants': {
            'X_1097183_A_GG': {
                'mode': {'X-linked over-dominance'},
                'variant':
                    variants_per_gene_102[
                        '1234'][
                        'X_1097183_A_GG'][
                        'child'],
                'hgncid': '1234'}}, 'compound_hets': {}}

        self.assertEqual(inheritancefilter_102.candidate_variants,
                         test_candidate_variants_102)
        # both unaff, mum 0/1, dad 1/1 fail
        variants_per_gene_112 = create_test_variants_per_gene(self.variants_112,
                                                              self.family_XX_both_unaff)
        inheritancefilter_112 = InheritanceFiltering(variants_per_gene_112,
                                                     self.family_XX_both_unaff,
                                                     self.genes_X_linked_over_dominance,
                                                     None,
                                                     None, self.candidate_variants,
                                                      self.inhreport)
        inheritancefilter_112.inheritance_filter_genes()
        test_candidate_variants_112 = {'single_variants': {},
                                       'compound_hets': {}}

        self.assertEqual(inheritancefilter_112.candidate_variants,
                         test_candidate_variants_112)
        # both unaff, mum 1/1, dad 1/1 pass
        variants_per_gene_122 = create_test_variants_per_gene(self.variants_122,
                                                              self.family_XX_both_unaff)
        inheritancefilter_122 = InheritanceFiltering(variants_per_gene_122,
                                                     self.family_XX_both_unaff,
                                                     self.genes_X_linked_over_dominance,
                                                     None,
                                                     None, self.candidate_variants,
                                                      self.inhreport)
        inheritancefilter_122.inheritance_filter_genes()
        test_candidate_variants_122 = {'single_variants': {
            'X_1097183_A_GG': {
                'mode': {'X-linked over-dominance'},
                'variant':
                    variants_per_gene_122[
                        '1234'][
                        'X_1097183_A_GG'][
                        'child'],
                'hgncid': '1234'}}, 'compound_hets': {}}

        self.assertEqual(inheritancefilter_122.candidate_variants,
                         test_candidate_variants_122)

    def test_gn_mono_y_hem_gt_hemizygous_parents_filter(self):
        # both aff, mum 0/0, dad 0/0 pass
        variants_per_gene_200 = create_test_variants_per_gene(self.variants_200Y,
                                                              self.family_XY_both_aff)
        inheritancefilter_200 = InheritanceFiltering(variants_per_gene_200,
                                                     self.family_XY_both_aff,
                                                     self.genes_monoallelic_Y_hemizygous,
                                                     None,
                                                     None, self.candidate_variants,
                                                      self.inhreport)
        inheritancefilter_200.inheritance_filter_genes()
        test_candidate_variants_200 = {'single_variants': {
            'Y_1097183_A_GG': {
                'mode': {'monoallelic_Y_hemizygous'},
                'variant':
                    variants_per_gene_200[
                        '1234'][
                        'Y_1097183_A_GG'][
                        'child'],
                'hgncid': '1234'}}, 'compound_hets': {}}

        self.assertEqual(inheritancefilter_200.candidate_variants,
                         test_candidate_variants_200)
        # both aff, mum 0/1, dad 0/0 pass
        variants_per_gene_210 = create_test_variants_per_gene(self.variants_210Y,
                                                              self.family_XY_both_aff)
        inheritancefilter_210 = InheritanceFiltering(variants_per_gene_200,
                                                     self.family_XY_both_aff,
                                                     self.genes_monoallelic_Y_hemizygous,
                                                     None,
                                                     None, self.candidate_variants,
                                                      self.inhreport)
        inheritancefilter_210.inheritance_filter_genes()
        test_candidate_variants_210 = {'single_variants': {
            'Y_1097183_A_GG': {
                'mode': {'monoallelic_Y_hemizygous'},
                'variant':
                    variants_per_gene_200[
                        '1234'][
                        'Y_1097183_A_GG'][
                        'child'],
                'hgncid': '1234'}}, 'compound_hets': {}}

        self.assertEqual(inheritancefilter_210.candidate_variants,
                         test_candidate_variants_210)
        # both aff, mum 1/1, dad 0/0 pass
        variants_per_gene_220 = create_test_variants_per_gene(self.variants_220Y,
                                                              self.family_XY_both_aff)
        inheritancefilter_220 = InheritanceFiltering(variants_per_gene_220,
                                                     self.family_XY_both_aff,
                                                     self.genes_monoallelic_Y_hemizygous,
                                                     None,
                                                     None, self.candidate_variants,
                                                      self.inhreport)
        inheritancefilter_220.inheritance_filter_genes()
        test_candidate_variants_220 = {'single_variants': {
            'Y_1097183_A_GG': {
                'mode': {'monoallelic_Y_hemizygous'},
                'variant':
                    variants_per_gene_200[
                        '1234'][
                        'Y_1097183_A_GG'][
                        'child'],
                'hgncid': '1234'}}, 'compound_hets': {}}

        self.assertEqual(inheritancefilter_220.candidate_variants,
                         test_candidate_variants_220)
        # both aff, mum 0/0, dad 1/1 pass
        variants_per_gene_202 = create_test_variants_per_gene(self.variants_202Y,
                                                              self.family_XY_both_aff)
        inheritancefilter_202 = InheritanceFiltering(variants_per_gene_202,
                                                     self.family_XY_both_aff,
                                                     self.genes_monoallelic_Y_hemizygous,
                                                     None,
                                                     None, self.candidate_variants,
                                                      self.inhreport)
        inheritancefilter_202.inheritance_filter_genes()
        test_candidate_variants_202 = {'single_variants': {
            'Y_1097183_A_GG': {
                'mode': {'monoallelic_Y_hemizygous'},
                'variant':
                    variants_per_gene_200[
                        '1234'][
                        'Y_1097183_A_GG'][
                        'child'],
                'hgncid': '1234'}}, 'compound_hets': {}}

        self.assertEqual(inheritancefilter_202.candidate_variants,
                         test_candidate_variants_202)
        # both aff, mum 0/1, dad 1/1 pass
        variants_per_gene_212 = create_test_variants_per_gene(self.variants_212Y,
                                                              self.family_XY_both_aff)
        inheritancefilter_212 = InheritanceFiltering(variants_per_gene_212,
                                                     self.family_XY_both_aff,
                                                     self.genes_monoallelic_Y_hemizygous,
                                                     None,
                                                     None, self.candidate_variants,
                                                      self.inhreport)
        inheritancefilter_212.inheritance_filter_genes()
        test_candidate_variants_212 = {'single_variants': {
            'Y_1097183_A_GG': {
                'mode': {'monoallelic_Y_hemizygous'},
                'variant':
                    variants_per_gene_200[
                        '1234'][
                        'Y_1097183_A_GG'][
                        'child'],
                'hgncid': '1234'}}, 'compound_hets': {}}

        self.assertEqual(inheritancefilter_212.candidate_variants,
                         test_candidate_variants_212)
        # both aff, mum 1/1, dad 1/1 pass
        variants_per_gene_222 = create_test_variants_per_gene(self.variants_222Y,
                                                              self.family_XY_both_aff)
        inheritancefilter_222 = InheritanceFiltering(variants_per_gene_222,
                                                     self.family_XY_both_aff,
                                                     self.genes_monoallelic_Y_hemizygous,
                                                     None,
                                                     None, self.candidate_variants,
                                                      self.inhreport)
        inheritancefilter_222.inheritance_filter_genes()
        test_candidate_variants_222 = {'single_variants': {
            'Y_1097183_A_GG': {
                'mode': {'monoallelic_Y_hemizygous'},
                'variant':
                    variants_per_gene_200[
                        '1234'][
                        'Y_1097183_A_GG'][
                        'child'],
                'hgncid': '1234'}}, 'compound_hets': {}}

        self.assertEqual(inheritancefilter_222.candidate_variants,
                         test_candidate_variants_222)

        # mum aff, mum 0/0, dad 0/0 pass
        variants_per_gene_200 = create_test_variants_per_gene(self.variants_200Y,
                                                              self.family_XY_mum_aff)
        inheritancefilter_200 = InheritanceFiltering(variants_per_gene_200,
                                                     self.family_XY_mum_aff,
                                                     self.genes_monoallelic_Y_hemizygous,
                                                     None,
                                                     None, self.candidate_variants,
                                                      self.inhreport)
        inheritancefilter_200.inheritance_filter_genes()
        test_candidate_variants_200 = {'single_variants': {
            'Y_1097183_A_GG': {
                'mode': {'monoallelic_Y_hemizygous'},
                'variant':
                    variants_per_gene_200[
                        '1234'][
                        'Y_1097183_A_GG'][
                        'child'],
                'hgncid': '1234'}}, 'compound_hets': {}}

        self.assertEqual(inheritancefilter_200.candidate_variants,
                         test_candidate_variants_200)
        # mum aff, mum 0/1, dad 0/0 pass
        variants_per_gene_210 = create_test_variants_per_gene(self.variants_210Y,
                                                              self.family_XY_mum_aff)
        inheritancefilter_210 = InheritanceFiltering(variants_per_gene_210,
                                                     self.family_XY_mum_aff,
                                                     self.genes_monoallelic_Y_hemizygous,
                                                     None,
                                                     None, self.candidate_variants,
                                                      self.inhreport)
        inheritancefilter_210.inheritance_filter_genes()
        test_candidate_variants_210 = {'single_variants': {
            'Y_1097183_A_GG': {
                'mode': {'monoallelic_Y_hemizygous'},
                'variant':
                    variants_per_gene_200[
                        '1234'][
                        'Y_1097183_A_GG'][
                        'child'],
                'hgncid': '1234'}}, 'compound_hets': {}}

        self.assertEqual(inheritancefilter_210.candidate_variants,
                         test_candidate_variants_210)
        # mum aff, mum 1/1, dad 0/0 pass
        variants_per_gene_220 = create_test_variants_per_gene(self.variants_220Y,
                                                              self.family_XY_mum_aff)
        inheritancefilter_220 = InheritanceFiltering(variants_per_gene_220,
                                                     self.family_XY_mum_aff,
                                                     self.genes_monoallelic_Y_hemizygous,
                                                     None,
                                                     None, self.candidate_variants,
                                                      self.inhreport)
        inheritancefilter_220.inheritance_filter_genes()
        test_candidate_variants_220 = {'single_variants': {
            'Y_1097183_A_GG': {
                'mode': {'monoallelic_Y_hemizygous'},
                'variant':
                    variants_per_gene_200[
                        '1234'][
                        'Y_1097183_A_GG'][
                        'child'],
                'hgncid': '1234'}}, 'compound_hets': {}}

        self.assertEqual(inheritancefilter_220.candidate_variants,
                         test_candidate_variants_220)
        # mum aff, mum 0/0, dad 1/1 fail
        variants_per_gene_202 = create_test_variants_per_gene(self.variants_202Y,
                                                              self.family_XY_mum_aff)
        inheritancefilter_202 = InheritanceFiltering(variants_per_gene_202,
                                                     self.family_XY_mum_aff,
                                                     self.genes_monoallelic_Y_hemizygous,
                                                     None,
                                                     None, self.candidate_variants,
                                                      self.inhreport)
        inheritancefilter_202.inheritance_filter_genes()
        test_candidate_variants_202 = {'single_variants': {}, 'compound_hets': {}}

        self.assertEqual(inheritancefilter_202.candidate_variants,
                         test_candidate_variants_202)
        # mum aff, mum 0/1, dad 1/1 fail
        variants_per_gene_212 = create_test_variants_per_gene(self.variants_212Y,
                                                              self.family_XY_mum_aff)
        inheritancefilter_212 = InheritanceFiltering(variants_per_gene_212,
                                                     self.family_XY_mum_aff,
                                                     self.genes_monoallelic_Y_hemizygous,
                                                     None,
                                                     None, self.candidate_variants,
                                                      self.inhreport)
        inheritancefilter_212.inheritance_filter_genes()
        test_candidate_variants_212 = {'single_variants': {}, 'compound_hets': {}}

        self.assertEqual(inheritancefilter_212.candidate_variants,
                         test_candidate_variants_212)
        # mum aff, mum 1/1, dad 1/1 fail
        variants_per_gene_222 = create_test_variants_per_gene(self.variants_222Y,
                                                              self.family_XY_mum_aff)
        inheritancefilter_222 = InheritanceFiltering(variants_per_gene_222,
                                                     self.family_XY_mum_aff,
                                                     self.genes_monoallelic_Y_hemizygous,
                                                     None,
                                                     None, self.candidate_variants,
                                                      self.inhreport)
        inheritancefilter_222.inheritance_filter_genes()
        test_candidate_variants_222 = {'single_variants': {}, 'compound_hets': {}}

        self.assertEqual(inheritancefilter_222.candidate_variants,
                         test_candidate_variants_222)

        # dad aff, mum 0/0, dad 0/0 pass
        variants_per_gene_200 = create_test_variants_per_gene(self.variants_200Y,
                                                              self.family_XY_dad_aff)
        inheritancefilter_200 = InheritanceFiltering(variants_per_gene_200,
                                                     self.family_XY_dad_aff,
                                                     self.genes_monoallelic_Y_hemizygous,
                                                     None,
                                                     None, self.candidate_variants,
                                                      self.inhreport)
        inheritancefilter_200.inheritance_filter_genes()
        test_candidate_variants_200 = {'single_variants': {
            'Y_1097183_A_GG': {
                'mode': {'monoallelic_Y_hemizygous'},
                'variant':
                    variants_per_gene_200[
                        '1234'][
                        'Y_1097183_A_GG'][
                        'child'],
                'hgncid': '1234'}}, 'compound_hets': {}}

        self.assertEqual(inheritancefilter_200.candidate_variants,
                         test_candidate_variants_200)
        # dad aff, mum 0/1, dad 0/0 pass
        variants_per_gene_210 = create_test_variants_per_gene(self.variants_210Y,
                                                              self.family_XY_dad_aff)
        inheritancefilter_210 = InheritanceFiltering(variants_per_gene_210,
                                                     self.family_XY_dad_aff,
                                                     self.genes_monoallelic_Y_hemizygous,
                                                     None,
                                                     None, self.candidate_variants,
                                                      self.inhreport)
        inheritancefilter_210.inheritance_filter_genes()
        test_candidate_variants_210 = {'single_variants': {
            'Y_1097183_A_GG': {
                'mode': {'monoallelic_Y_hemizygous'},
                'variant':
                    variants_per_gene_200[
                        '1234'][
                        'Y_1097183_A_GG'][
                        'child'],
                'hgncid': '1234'}}, 'compound_hets': {}}

        self.assertEqual(inheritancefilter_210.candidate_variants,
                         test_candidate_variants_210)
        # dad aff, mum 1/1, dad 0/0 pass
        variants_per_gene_220 = create_test_variants_per_gene(self.variants_220Y,
                                                              self.family_XY_dad_aff)
        inheritancefilter_220 = InheritanceFiltering(variants_per_gene_220,
                                                     self.family_XY_dad_aff,
                                                     self.genes_monoallelic_Y_hemizygous,
                                                     None,
                                                     None, self.candidate_variants,
                                                      self.inhreport)
        inheritancefilter_220.inheritance_filter_genes()
        test_candidate_variants_220 = {'single_variants': {
            'Y_1097183_A_GG': {
                'mode': {'monoallelic_Y_hemizygous'},
                'variant':
                    variants_per_gene_200[
                        '1234'][
                        'Y_1097183_A_GG'][
                        'child'],
                'hgncid': '1234'}},
                                       'compound_hets': {}}

        self.assertEqual(inheritancefilter_220.candidate_variants,
                         test_candidate_variants_220)
        # dad aff, mum 0/0, dad 1/1 pass
        variants_per_gene_202 = create_test_variants_per_gene(self.variants_202Y,
                                                              self.family_XY_dad_aff)
        inheritancefilter_202 = InheritanceFiltering(variants_per_gene_202,
                                                     self.family_XY_dad_aff,
                                                     self.genes_monoallelic_Y_hemizygous,
                                                     None,
                                                     None, self.candidate_variants,
                                                      self.inhreport)
        inheritancefilter_202.inheritance_filter_genes()
        test_candidate_variants_202 = {'single_variants': {
            'Y_1097183_A_GG': {
                'mode': {'monoallelic_Y_hemizygous'},
                'variant':
                    variants_per_gene_200[
                        '1234'][
                        'Y_1097183_A_GG'][
                        'child'],
                'hgncid': '1234'}}, 'compound_hets': {}}

        self.assertEqual(inheritancefilter_202.candidate_variants,
                         test_candidate_variants_202)
        # dad aff, mum 0/1, dad 1/1 pass
        variants_per_gene_212 = create_test_variants_per_gene(self.variants_212Y,
                                                              self.family_XY_dad_aff)
        inheritancefilter_212 = InheritanceFiltering(variants_per_gene_212,
                                                     self.family_XY_dad_aff,
                                                     self.genes_monoallelic_Y_hemizygous,
                                                     None,
                                                     None, self.candidate_variants,
                                                      self.inhreport)
        inheritancefilter_212.inheritance_filter_genes()
        test_candidate_variants_212 = {'single_variants': {
            'Y_1097183_A_GG': {
                'mode': {'monoallelic_Y_hemizygous'},
                'variant':
                    variants_per_gene_200[
                        '1234'][
                        'Y_1097183_A_GG'][
                        'child'],
                'hgncid': '1234'}}, 'compound_hets': {}}

        self.assertEqual(inheritancefilter_212.candidate_variants,
                         test_candidate_variants_212)
        # dad aff, mum 1/1, dad 1/1 pass
        variants_per_gene_222 = create_test_variants_per_gene(self.variants_222Y,
                                                              self.family_XY_dad_aff)
        inheritancefilter_222 = InheritanceFiltering(variants_per_gene_222,
                                                     self.family_XY_dad_aff,
                                                     self.genes_monoallelic_Y_hemizygous,
                                                     None,
                                                     None, self.candidate_variants,
                                                      self.inhreport)
        inheritancefilter_222.inheritance_filter_genes()
        test_candidate_variants_222 = {'single_variants': {
            'Y_1097183_A_GG': {
                'mode': {'monoallelic_Y_hemizygous'},
                'variant':
                    variants_per_gene_200[
                        '1234'][
                        'Y_1097183_A_GG'][
                        'child'],
                'hgncid': '1234'}},
                                       'compound_hets': {}}

        self.assertEqual(inheritancefilter_222.candidate_variants,
                         test_candidate_variants_222)

        # both unaff, mum 0/0, dad 0/0 pass
        variants_per_gene_200 = create_test_variants_per_gene(self.variants_200Y,
                                                              self.family_XY_both_unaff)
        inheritancefilter_200 = InheritanceFiltering(variants_per_gene_200,
                                                     self.family_XY_both_unaff,
                                                     self.genes_monoallelic_Y_hemizygous,
                                                     None,
                                                     None, self.candidate_variants,
                                                      self.inhreport)
        inheritancefilter_200.inheritance_filter_genes()
        test_candidate_variants_200 = {'single_variants': {
            'Y_1097183_A_GG': {
                'mode': {'monoallelic_Y_hemizygous'},
                'variant':
                    variants_per_gene_200[
                        '1234'][
                        'Y_1097183_A_GG'][
                        'child'],
                'hgncid': '1234'}}, 'compound_hets': {}}

        self.assertEqual(inheritancefilter_200.candidate_variants,
                         test_candidate_variants_200)
        # both unaff, mum 0/1, dad 0/0 pass
        variants_per_gene_210 = create_test_variants_per_gene(self.variants_210Y,
                                                              self.family_XY_both_unaff)
        inheritancefilter_210 = InheritanceFiltering(variants_per_gene_210,
                                                     self.family_XY_both_unaff,
                                                     self.genes_monoallelic_Y_hemizygous,
                                                     None,
                                                     None, self.candidate_variants,
                                                      self.inhreport)
        inheritancefilter_210.inheritance_filter_genes()
        test_candidate_variants_210 = {'single_variants': {
            'Y_1097183_A_GG': {
                'mode': {'monoallelic_Y_hemizygous'},
                'variant':
                    variants_per_gene_200[
                        '1234'][
                        'Y_1097183_A_GG'][
                        'child'],
                'hgncid': '1234'}}, 'compound_hets': {}}

        self.assertEqual(inheritancefilter_210.candidate_variants,
                         test_candidate_variants_210)
        # both unaff, mum 1/1, dad 0/0 pass
        variants_per_gene_220 = create_test_variants_per_gene(self.variants_220Y,
                                                              self.family_XY_both_unaff)
        inheritancefilter_220 = InheritanceFiltering(variants_per_gene_220,
                                                     self.family_XY_both_unaff,
                                                     self.genes_monoallelic_Y_hemizygous,
                                                     None,
                                                     None, self.candidate_variants,
                                                      self.inhreport)
        inheritancefilter_220.inheritance_filter_genes()
        test_candidate_variants_220 = {'single_variants': {
            'Y_1097183_A_GG': {
                'mode': {'monoallelic_Y_hemizygous'},
                'variant':
                    variants_per_gene_200[
                        '1234'][
                        'Y_1097183_A_GG'][
                        'child'],
                'hgncid': '1234'}},
                                       'compound_hets': {}}

        self.assertEqual(inheritancefilter_220.candidate_variants,
                         test_candidate_variants_220)
        # both unaff, mum 0/0, dad 1/1 fail
        variants_per_gene_202 = create_test_variants_per_gene(self.variants_202Y,
                                                              self.family_XY_both_unaff)
        inheritancefilter_202 = InheritanceFiltering(variants_per_gene_202,
                                                     self.family_XY_both_unaff,
                                                     self.genes_monoallelic_Y_hemizygous,
                                                     None,
                                                     None, self.candidate_variants,
                                                      self.inhreport)
        inheritancefilter_202.inheritance_filter_genes()
        test_candidate_variants_202 = {'single_variants': {}, 'compound_hets': {}}

        self.assertEqual(inheritancefilter_202.candidate_variants,
                         test_candidate_variants_202)
        # both unaff, mum 0/1, dad 1/1 fail
        variants_per_gene_212 = create_test_variants_per_gene(self.variants_212Y,
                                                              self.family_XY_both_unaff)
        inheritancefilter_212 = InheritanceFiltering(variants_per_gene_212,
                                                     self.family_XY_both_unaff,
                                                     self.genes_monoallelic_Y_hemizygous,
                                                     None,
                                                     None, self.candidate_variants,
                                                      self.inhreport)
        inheritancefilter_212.inheritance_filter_genes()
        test_candidate_variants_212 = {'single_variants': {}, 'compound_hets': {}}

        self.assertEqual(inheritancefilter_212.candidate_variants,
                         test_candidate_variants_212)
        # both unaff, mum 1/1, dad 1/1 fail
        variants_per_gene_222 = create_test_variants_per_gene(self.variants_222Y,
                                                              self.family_XY_both_unaff)
        inheritancefilter_222 = InheritanceFiltering(variants_per_gene_222,
                                                     self.family_XY_both_unaff,
                                                     self.genes_monoallelic_Y_hemizygous,
                                                     None,
                                                     None, self.candidate_variants,
                                                      self.inhreport)
        inheritancefilter_222.inheritance_filter_genes()
        test_candidate_variants_222 = {'single_variants': {},
                                       'compound_hets': {}}

        self.assertEqual(inheritancefilter_222.candidate_variants,
                         test_candidate_variants_222)

if __name__ == '__main__':
    unittest.main()